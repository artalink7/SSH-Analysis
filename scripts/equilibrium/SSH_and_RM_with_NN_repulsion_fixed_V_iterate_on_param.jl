using LinearAlgebra
using ITensors
using ITensorMPS
using Statistics
using ITensorGaussianMPS
using Printf

# =========================================================
# 1. SETUP FUNCTIONS (Standard)
# =========================================================

function build_SSH_MPO_PBC(sites; v=1.0, w=0.5, Δ=0.0, V=0.0)
    N = length(sites)
    os = OpSum()
    # Bulk Terms
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
    # Interactions
    if V != 0.0
        for i in 1:N-1; os += V, "N", i, "N", i+1; end
    end
    # Periodic Boundary Terms
    os += -w, "Cdag", N, "C", 1; os += -w, "Cdag", 1, "C", N
    if V != 0.0; os += V, "N", N, "N", 1; end
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

# =========================================================
# 2. MAIN EXPERIMENT: "Gap Sweep"
# =========================================================

function main_gap_sweep()
    # System Setup
    N = 600
    La_width = 300 
    sites = siteinds("Fermion", N; conserve_qns=true) 
    state_pattern = [isodd(n) ? "1" : "0" for n in 1:N]

    # --- KEY PARAMETER: FIXED INTERACTION ---
    V_fixed = 1.0  # We fix the interaction strength!

    # ==================================================
    # EXPERIMENT 1: SSH (Varying w to tune the Gap)
    # ==================================================
    println("\n=== EXPERIMENT 1: SSH (Varying w, V=$V_fixed) ===")
    
    w_values = 1.2:0.1:2.0

    open("data/equilibrium/data_SSH_gap_sweep_new.csv", "w") do file
        write(file, "w,S,F,Ratio_SF,GapNorm,RawGap\n")
        
        v_param = 1.0
        Delta_param = 0.0
        
        for w_param in w_values
            # 1. Update Bandwidth (It changes as we change w!)
            W_SSH = 2.0 * (abs(v_param) + abs(w_param))
            
            # 2. Build Hamiltonian
            H = build_SSH_MPO_PBC(sites, v=v_param, w=w_param, Δ=Delta_param, V=V_fixed)

            # 3. Ground State (E0)
            psi0_init = randomMPS(sites, state_pattern, linkdims=10)
            energy0, psi0 = dmrg(H, psi0_init; nsweeps=8, maxdim=[20, 60, 100, 200], noise=[1E-6, 1E-7, 0.0, 0.0], cutoff=1E-10, outputlevel=0)

            # 4. Excited State (E1)
            psi1_init = randomMPS(sites, state_pattern, linkdims=10)
            energy1, psi1 = dmrg(H, [psi0], psi1_init; nsweeps=8, maxdim=[20, 60, 100, 200], noise=[1E-6, 1E-7, 0.0, 0.0], cutoff=1E-10, weight=20.0, outputlevel=0)

            # 5. Measure Gap
            raw_gap = real(energy1 - energy0)
            #if raw_gap < 0; raw_gap = abs(raw_gap); end # Safety
            norm_gap = raw_gap / W_SSH

            # 6. Observables
            S, F = measure_observables(psi0, N, La_width)
            ratio = (F > 1e-10) ? S / F : 0.0
            
            write(file, "$w_param, $S, $F, $ratio, $norm_gap, $raw_gap\n")
            @printf("SSH | w=%-4.2f | GapNorm=%-7.4f | S=%-7.4f\n", w_param, norm_gap, S)
        end
    end 

    # ==================================================
    # EXPERIMENT 2: Rice-Mele (Varying Delta to tune the Gap)
    # ==================================================
    #println("\n=== EXPERIMENT 2: Rice-Mele (Varying Delta, V=$V_fixed) ===")
    
    # We scan Delta from 0.1 to 2.0. (Avoid exactly 0.0)
    Delta_values = 0.1:0.2:2.0

    open("data/equilibrium/data_RM_gap_sweep_new.csv", "w") do file
        write(file, "Delta,S,F,Ratio_SF,GapNorm,RawGap\n")
        
        v_param = 1.0
        w_param = 1.0 # Standard RM typically has v=w
        
        for Delta_param in Delta_values
            # 1. Update Bandwidth (Changes with Delta!)
            W_RM = 2.0 * sqrt((abs(v_param) + abs(w_param))^2 + Delta_param^2)
            
            # 2. Build Hamiltonian
            H = build_SSH_MPO_PBC(sites, v=v_param, w=w_param, Δ=Delta_param, V=V_fixed)

            # 3. Ground State
            psi0_init = randomMPS(sites, state_pattern, linkdims=10)
            energy0, psi0 = dmrg(H, psi0_init; nsweeps=8, maxdim=[20, 60, 100, 200], noise=[1E-6, 1E-7, 0.0, 0.0], cutoff=1E-10, outputlevel=0)

            # 4. Excited State
            psi1_init = randomMPS(sites, state_pattern, linkdims=10)
            energy1, psi1 = dmrg(H, [psi0], psi1_init; nsweeps=8, maxdim=[20, 60, 100, 200], noise=[1E-6, 1E-7, 0.0, 0.0], cutoff=1E-10, weight=20.0, outputlevel=0)

            # 5. Measure Gap
            raw_gap = real(energy1 - energy0)
            #if raw_gap < 0; raw_gap = abs(raw_gap); end
            norm_gap = raw_gap / W_RM

            # 6. Observables
            S, F = measure_observables(psi0, N, La_width)
            ratio = (F > 1e-10) ? S / F : 0.0
            
            write(file, "$Delta_param, $S, $F, $ratio, $norm_gap, $raw_gap\n")
            @printf("RM  | D=%-4.2f | GapNorm=%-7.4f | S=%-7.4f\n", Delta_param, norm_gap, S)
        end
    end 
    
    println("\nData collection complete.")
end

# Run it
main_gap_sweep()