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
        for i in 1:N-1
            os += V, "N", i, "N", i+1
            os += -V/2.0, "N", i
            os += -V/2.0, "N", i+1
        end
    end

    return MPO(os, sites)
end

function measure_observables(psi, N, V)
    La_width = N ÷ 2
    subsystem_inds = 1:La_width

    # 1. Entanglement Entropy
    orthogonalize!(psi, La_width) 
    row_inds = uniqueinds(psi[La_width], psi[La_width+1])
    U, S_vals, V_mat = svd(psi[La_width], row_inds)

    S_EE = 0.0
    for n in 1:dim(S_vals, 1)
        λ = S_vals[n, n]
        p = λ^2
        if p > 1e-12 # Avoid log(0) issues
            S_EE -= p * log(p)
        end
    end
    S_EE = real(S_EE) 

    # 2. Particle Number Fluctuations
    if iszero(V)
        CM = correlation_matrix(psi, "Cdag", "C", sites = subsystem_inds)
        F = real(tr(CM * (I - CM)))
    else
        NM = correlation_matrix(psi, "N", "N", sites = subsystem_inds)
        N_A = real(tr(NM))
        F = real(sum(NM) - N_A^2)
    end
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

function build_tebd_gates(sites; v_f=1.0, w_f=0.5, dt=0.05)
    N = length(sites)
    gates = ITensor[]
    
    # Forward sweep: dt/2
    for j in 1:(N - 1)
        s1 = sites[j]
        s2 = sites[j + 1]
        
        # Alternating hopping parameters for SSH
        hop = (j % 2 != 0) ? v_f : w_f
        
        # Local Hamiltonian bond (Non-interacting)
        hj = -hop * op("Cdag", s1) * op("C", s2) - hop * op("Cdag", s2) * op("C", s1)
        
        # Exponentiate the local Hamiltonian to create the Trotter gate.
        # ITensors automatically pairs the primed/unprimed indices.
        Gj = exp(-im * dt / 2 * hj)
        push!(gates, Gj)
    end
    
    # Reverse sweep: dt/2 to complete the 2nd order Trotter step
    append!(gates, reverse(gates))

    println("Constructed ", length(gates), " TEBD gates for quench evolution.")
    
    return gates
end

function simulate_quench(N_sites, T_max, dt; v_i=1.0, w_i=0.5, V_i=1.0, v_f=1.0, w_f=0.5)
    # 1. Initialize sites with particle number conservation
    sites = siteinds("Fermion", N_sites; conserve_qns=true)
    
    # 2. Prepare the Initial State (t < 0)
    println("Preparing interacting ground state (t < 0)...")
    H_init = build_SSH_MPO_OBC(sites; v=v_i, w=w_i, V=V_i)
    
    # Using half-filling based on product state function
    init_state = product_state_Nf(N_sites, N_sites ÷ 2) 
    E0, psi = ground_energy(H_init, sites; init_state=init_state, nsweeps=12)
    println("Initial Ground State Energy: ", E0)
    
    # 3. Build TEBD gates for post-quench Hamiltonian (t > 0)
    println("Constructing TEBD gates for quench...")
    gates = build_tebd_gates(sites; v_f=v_f, w_f=w_f, dt=dt)
    
    # 4. Time Evolution Loop
    times = 0.0:dt:T_max
    S_EE_vals = Float64[]
    F_vals = Float64[]
    La_width = N_sites ÷ 2 # Subsystem size for measurements
    
    println("Starting time evolution...")
    for t in times
        # Measure observables before applying the time step
        S_EE, F = measure_observables(psi, N_sites, V_i)
        push!(S_EE_vals, S_EE)
        push!(F_vals, F)
        
        # Apply the Trotter gates to evolve the state by dt
        # 'cutoff' and 'maxdim' are critical here to manage entanglement growth
        psi = apply(gates, psi; cutoff=1e-10, maxdim=1000)
        normalize!(psi) # Normalize after each full Trotter step
        println("Time: ", round(t, digits=3), " | S_EE: ", round(S_EE, digits=4), " | F: ", round(F, digits=4))
    end
    
    return times, S_EE_vals, F_vals
end

function build_SSH_MPO_OBC_disordered(sites; disorder_strength=0.5, v=1.0, w=0.5, Δ=0.0, V=0.0)
    N = length(sites)
    os = OpSum()

    # Hopping terms with disorder
    for i in 1:2:N-1
        v_disorder = v + disorder_strength * (rand() - 0.5) 
        os += -v_disorder, "Cdag", i, "C", i+1; os += -v_disorder, "Cdag", i+1, "C", i
    end
    for i in 2:2:N-1
        w_disorder = w + disorder_strength * (rand() - 0.5)
        os += -w_disorder, "Cdag", i, "C", i+1; os += -w_disorder, "Cdag", i+1, "C", i
    end

    # Staggered Potential
    if Δ != 0.0
        for i in 1:N; os += Δ * (-1)^i, "N", i; end
    end

    # Nearest-neighbor interaction
    if V != 0.0
        for i in 1:N-1
            os += V, "N", i, "N", i+1
            os += -V/2.0, "N", i
            os += -V/2.0, "N", i+1
        end
    end

    return MPO(os, sites)

end

function simulate_quench_initial_disorder(N_sites, T_max, dt; v_i=1.0, w_i=0.5, V_i=1.0, v_f=1.0, w_f=0.5, W, bond_dim=1000)
    # 1. Initialize sites with particle number conservation
    sites = siteinds("Fermion", N_sites; conserve_qns=true)
    
    # 2. Prepare the Initial State (t < 0)
    println("Preparing interacting ground state (t < 0)...")
    H_init = build_SSH_MPO_OBC_disordered(sites; disorder_strength=W, v=v_i, w=w_i, Δ=0.0, V=V_i)
    
    # Using half-filling based on product state function
    init_state = product_state_Nf(N_sites, N_sites ÷ 2) 
    E0, psi = ground_energy(H_init, sites; init_state=init_state, nsweeps=12)
    println("Initial Ground State Energy: ", E0)
    
    # 3. Build TEBD gates for post-quench Hamiltonian (t > 0)
    println("Constructing TEBD gates for quench...")
    gates = build_tebd_gates(sites; v_f=v_f, w_f=w_f, dt=dt)
    
    # 4. Time Evolution Loop
    times = 0.0:dt:T_max
    S_EE_vals = Float64[]
    F_vals = Float64[]
    La_width = N_sites ÷ 2 # Subsystem size for measurements
    
    println("Starting time evolution...")
    for t in times
        # Measure observables before applying the time step
        S_EE, F = measure_observables(psi, N_sites, V_i)
        push!(S_EE_vals, S_EE)
        push!(F_vals, F)
        
        # Apply the Trotter gates to evolve the state by dt
        # 'cutoff' and 'maxdim' are critical here to manage entanglement growth
        psi = apply(gates, psi; cutoff=1e-10, maxdim=bond_dim)
        normalize!(psi) # Normalize after each full Trotter step
        println("Time: ", round(t, digits=3), " | S_EE: ", round(S_EE, digits=4), " | F: ", round(F, digits=4))
    end
    
    return times, S_EE_vals, F_vals
end