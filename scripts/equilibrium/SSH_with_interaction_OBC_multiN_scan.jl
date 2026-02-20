using LinearAlgebra
using ITensors
using ITensorMPS
using Statistics
using ITensorGaussianMPS
using Printf

# =========================================================
# 1. SETUP FUNCTIONS (Standard)
# =========================================================

function buil_SSH_MPO_OBC(sites; v=1.0, w=0.5, Δ=0.0, V=0.0)
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

function ground_energy(H::MPO, sites; init_state, nsweeps=12, maxdim=[50, 100, 200, 400, 800], cutoff=1e-10, outputlevel=0)
    psi0 = MPS(sites, init_state)
    energy, psi = dmrg(H, psi0; nsweeps=nsweeps, maxdim=maxdim, noise=[1E-6, 1E-7, 0.0, 0.0, 0.0], cutoff=cutoff, outputlevel=outputlevel)
    return energy , psi
end 

# =========================================================
# 2. MAIN EXPERIMENT: "Charge Gap Scan"
# =========================================================

function scan_charge_gap_multiN()
    v_fixed = 1.0
    w_fixed = 1.0
    Δ_fixed = 0.0

    Ns = [50, 100, 200, 400]
    V_values = 0.0:0.2:4.0

    open("data/equilibrium/charge_gap_multiN_probing_phase_transition.csv", "w") do file 
        write(file, "N,V,ChargeGap\n")

        for N in Ns
            @printf("=== N = %d ===\n", N)

            sites = siteinds("Fermion", N; conserve_qns=true)

            Nf0 = N ÷ 2
            st0 = product_state_Nf(N, Nf0)
            stp = product_state_Nf(N, Nf0 + 1)
            stm = product_state_Nf(N, Nf0 - 1)

            for V in V_values
                H = buil_SSH_MPO_OBC(sites; v=v_fixed, w=w_fixed, Δ=Δ_fixed, V=V)

                E0 , _ = ground_energy(H, sites; init_state=st0)
                Ep , _ = ground_energy(H, sites; init_state=stp)
                Em , _ = ground_energy(H, sites; init_state=stm)

                charge_gap = (Ep + Em - 2.0 * E0)

                write(file, "$N, $V, $charge_gap\n")
                @printf("N=%d | V=%4.2f | ChargeGap=%.4e\n", N, V, charge_gap)
            end
        end
    end
end

# Run the scan
scan_charge_gap_multiN()