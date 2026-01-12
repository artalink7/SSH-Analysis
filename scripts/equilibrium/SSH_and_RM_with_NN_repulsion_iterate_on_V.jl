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
    subsystem_inds = start_site:(center + Int(La_width/2))

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
# 2. MAIN EXPERIMENT: "Interaction V Scan"
# =========================================================

function main_interaction_scan()
    # System Setup
    N = 600
    La_width = 300 
    sites = siteinds("Fermion", N; conserve_qns=true) 
    state_pattern = [isodd(n) ? "1" : "0" for n in 1:N]

    # --- SCAN PARAMETER: INTERACTION V ---
    V_values = 0.1:0.1:2.0  # Incrementing V to open the gap

    # ==================================================
    # EXPERIMENT 1: SSH (Fixed w=1.0, Varying V)
    # ==================================================
    println("\n=== EXPERIMENT 1: SSH (Fixed w=1.0, Varying V) ===")
    
    v_fixed = 1.0
    w_fixed = 1.3 
    Delta_fixed = 0.0
    
    # Bandwidth is constant here because v, w are fixed
    W_SSH = 2.0 * (abs(v_fixed) + abs(w_fixed)) 

    open("data/equilibrium/data_SSH_interaction_scan.csv", "w") do file
        write(file, "V,S,F,Ratio_SF,GapNorm,RawGap\n")
        
        for V_int in V_values
            # 1. Build Hamiltonian with varying V
            H = build_SSH_MPO_PBC(sites, v=v_fixed, w=w_fixed, Δ=Delta_fixed, V=V_int)

            # 2. Ground State (E0)
            psi0_init = randomMPS(sites, state_pattern, linkdims=10)
            energy0, psi0 = dmrg(H, psi0_init; nsweeps=8, maxdim=[20, 60, 100, 200], noise=[1E-6, 1E-7, 0.0, 0.0], cutoff=1E-10, outputlevel=0)

            # 3. Excited State (E1)
            psi1_init = randomMPS(sites, state_pattern, linkdims=10)
            energy1, psi1 = dmrg(H, [psi0], psi1_init; nsweeps=8, maxdim=[20, 60, 100, 200], noise=[1E-6, 1E-7, 0.0, 0.0], cutoff=1E-10, weight=20.0, outputlevel=0)

            # 4. Measure Gap
            raw_gap = real(energy1 - energy0)
            # if raw_gap < 0; raw_gap = abs(raw_gap); end
            norm_gap = raw_gap / W_SSH

            # 5. Observables
            S, F = measure_observables(psi0, N, La_width)
            ratio = (F > 1e-10) ? S / F : 0.0
            
            write(file, "$V_int, $S, $F, $ratio, $norm_gap, $raw_gap\n")
            @printf("SSH | V=%-4.2f | GapNorm=%-7.4f | S=%-7.4f\n", V_int, norm_gap, S)
        end
    end 

    # ==================================================
    # EXPERIMENT 2: Rice-Mele (Fixed Delta=0.5, Varying V)
    # ==================================================
    println("\n=== EXPERIMENT 2: Rice-Mele (Fixed Delta=0.5, Varying V) ===")
    
    v_fixed = 1.0
    w_fixed = 1.0 
    Delta_fixed = 0.5 # Fixed Trivial Gap
    
    # Constant Bandwidth
    W_RM = 2.0 * sqrt((abs(v_fixed) + abs(w_fixed))^2 + Delta_fixed^2)

    open("data/equilibrium/data_RM_interaction_scan.csv", "w") do file
        write(file, "V,S,F,Ratio_SF,GapNorm,RawGap\n")
        
        for V_int in V_values
            # 1. Build Hamiltonian with varying V
            H = build_SSH_MPO_PBC(sites, v=v_fixed, w=w_fixed, Δ=Delta_fixed, V=V_int)

            # 2. Ground State
            psi0_init = randomMPS(sites, state_pattern, linkdims=10)
            energy0, psi0 = dmrg(H, psi0_init; nsweeps=8, maxdim=[20, 60, 100, 200], noise=[1E-6, 1E-7, 0.0, 0.0], cutoff=1E-10, outputlevel=0)

            # 3. Excited State
            psi1_init = randomMPS(sites, state_pattern, linkdims=10)
            energy1, psi1 = dmrg(H, [psi0], psi1_init; nsweeps=8, maxdim=[20, 60, 100, 200], noise=[1E-6, 1E-7, 0.0, 0.0], cutoff=1E-10, weight=20.0, outputlevel=0)

            # 4. Measure Gap
            raw_gap = real(energy1 - energy0)
            # if raw_gap < 0; raw_gap = abs(raw_gap); end
            norm_gap = raw_gap / W_RM

            # 5. Observables
            S, F = measure_observables(psi0, N, La_width)
            ratio = (F > 1e-10) ? S / F : 0.0
            
            write(file, "$V_int, $S, $F, $ratio, $norm_gap, $raw_gap\n")
            @printf("RM  | V=%-4.2f | GapNorm=%-7.4f | S=%-7.4f\n", V_int, norm_gap, S)
        end
    end
    
    println("\nInteraction Scan Complete.")
end

# Run it
main_interaction_scan()