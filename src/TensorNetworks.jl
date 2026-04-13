function build_SSH_MPO_OBC(sites; v=1.0, w=0.5, Δ=0.0, V=0.0)
    N = length(sites)
    os = OpSum()

    # Hopping terms
    for i in 1:2:N-1
        os += -v, "Cdag", i, "C", i+1; os += -v, "Cdag", i+1, "C", i
    end
    for i in 2:2:N-1
        os += -w, "Cdag", i, "C", i+1; os += -w, "Cdag", i+1, "C", i
    end

    # Staggered Potential
    if Δ != 0.0
        for i in 1:N; os += Δ * (-1)^i, "N", i; end
    end

    # Nearest-neighbor interaction
    if V != 0.0
        for i in 1:N-1; os += V, "N", i, "N", i+1; end
    end

    return MPO(os, sites)
end

function measure_observables(psi, N, La_width)
    center = Int(N/2)
    start_site = center - Int(La_width/2) + 1
    end_site   = center + Int(La_width/2)
    subsystem_inds = start_site:end_site

    # Entanglement Entropy
    orthogonalize!(psi, center) 
    row_inds = uniqueinds(psi[center], psi[center+1])
    U, S_vals, V_mat = svd(psi[center], row_inds)
    S_EE = 0.0
    for λ in diag(S_vals)
        p = λ^2
        if p > 1e-12; S_EE -= p * log(p); end
    end

    # Particle Number Fluctuations
    NM = correlation_matrix(psi, "N", "N"; sites=subsystem_inds)
    N_A_mean = 0.0
    N_A_sq_mean = 0.0
    for i in 1:La_width
        for j in 1:La_width; N_A_sq_mean += NM[i, j]; end
        N_A_mean += NM[i, i]
    end
    F = N_A_sq_mean - (N_A_mean)^2

    return S_EE, F
end

function product_state_Nf(N::Int, Nf::Int; pattern::Symbol=:alternating)
    @assert 0<= Nf <= N
    st = fill("0", N)

    if pattern == :alternating 
        order = vcat(collect(1:2:N), collect(2:2:N))
    else
        error("Unknown pattern = $pattern")
    end

    for k in 1:Nf
        st[order[k]] = "1"
    end

    return st
end

function ground_energy(H::MPO, sites; init_state, nsweeps=12)
    # 1. Use randomMPS to break the strict 1,0,1,0 pinning
    psi0 = randomMPS(sites, init_state)
    
    # 2. Use the Sweeps object to apply enough noise to allow 
    # the edge states to tunnel across the chain and symmetrize
    sweeps = Sweeps(nsweeps)
    maxdim!(sweeps, 50, 100, 200, 400, 800)
    cutoff!(sweeps, 1e-10)
    noise!(sweeps, 1E-4, 1E-5, 1E-6, 1E-7, 0.0) # Slightly higher initial noise
    
    energy, psi = dmrg(H, psi0, sweeps; outputlevel=1)
    return energy, psi
end

function create_CDW_states(N::Int)
    # 1. Ground state: perfectly alternating |1 0 1 0 1 0 ...>
    st0 = fill("0", N)
    for i in 1:2:N; st0[i] = "1"; end
    
    # 2. Ep State: Add particle to an even site as close to the center as possible
    stp = copy(st0)
    center_even = (N ÷ 2) % 2 == 0 ? (N ÷ 2) : (N ÷ 2) + 1
    stp[center_even] = "1"
    
    # 3. Em State: Remove particle from an odd site as close to the center as possible
    stm = copy(st0)
    center_odd = (N ÷ 2) % 2 != 0 ? (N ÷ 2) : (N ÷ 2) + 1
    stm[center_odd] = "0"
    
    return st0, stp, stm
end

function run_dmrg_cdw(H, sites, st_array)
    # randomMPS introduces immediate quantum fluctuations to break symmetry
    psi0 = randomMPS(sites, st_array) 
    
    # Aggressive noise schedule to unpin domain walls!
    sweeps = Sweeps(15)
    maxdim!(sweeps, 50, 100, 200, 400, 800)
    cutoff!(sweeps, 1E-10)
    noise!(sweeps, 1E-3, 1E-4, 1E-4, 1E-5, 1E-5, 1E-6, 1E-7, 1E-8, 0.0)
    
    energy, psi = dmrg(H, psi0, sweeps; outputlevel=0)
    return energy
end