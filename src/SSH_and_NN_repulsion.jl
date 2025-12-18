using LinearAlgebra
using ITensors
using ITensorMPS
using Statistics
using ITensorGaussianMPS
using Printf

# 1) Build single-particle SSH Hamiltonian 
function build_SSH_single_particle(Ncells; v= 1.0, w=0.5, Δ=0.0)
    N = 2 * Ncells  # total number of sites
    h = zeros(Float64, N, N)

    # --- hopping terms ---
    for i in 1:Ncells
        #intracell hopping
        a = 2*i - 1
        b = 2*i
        h[a, b] = -v
        h[b, a] = -v
    end
    for i in 1:(Ncells - 1)
        #intercell hopping
        b = 2*i
        a_next = 2*i + 1
        h[b, a_next] = -w
        h[a_next, b] = -w
    end

    # --- staggered onsite potential ---
    if Δ != 0.0
        for i in 1:N
            h[i, i] += Δ * (-1)^i
        end
    end

    return h
end

function build_SSH_MPO_PBC(sites; v=1.0, w=0.5, Δ=0.0, V=0.0)
    N = length(sites)
    os = OpSum()

    # --- 1. Bulk Terms (Same as before) ---
    for i in 1:2:N-1
        # Intracell (v): 1-2, 3-4, ...
        os += -v, "Cdag", i, "C", i+1
        os += -v, "Cdag", i+1, "C", i
    end
    for i in 2:2:N-1
        # Intercell (w): 2-3, 4-5, ...
        os += -w, "Cdag", i, "C", i+1
        os += -w, "Cdag", i+1, "C", i
    end

    # --- 2. Staggered Potential (Same as before) ---
    if Δ != 0.0
        for i in 1:N
            os += Δ * (-1)^i, "N", i
        end
    end

    # --- 3. Interactions (Bulk) ---
    if V != 0.0
        for i in 1:N-1
            os += V, "N", i, "N", i+1
        end
    end

    # =========================================================
    # --- 4. PERIODIC BOUNDARY TERMS (Wrapping N -> 1) ---
    # =========================================================
    
    # Hopping: Connects N (sublattice B) to 1 (sublattice A)
    # This closes the intercell gap, so we use 'w'
    os += -w, "Cdag", N, "C", 1
    os += -w, "Cdag", 1, "C", N

    # Interaction: Repulsion across the boundary
    if V != 0.0
        os += V, "N", N, "N", 1
    end

    # Convert to MPO 
    H = MPO(os, sites)
    return H
end

function measure_observables(psi, N, La_width)
    # Cut in the middle of the chain
    center = Int(N/2)
    # Define Subsystem A: from (center - La/2) to (center + La/2)
    start_site = center - Int(La_width/2) + 1
    end_site   = center + Int(La_width/2)
    subsystem_inds = start_site:end_site

    # --- 1. Entanglement Entropy (S) ---
    # We cut the bond at 'end_site'
    orthogonalize!(psi, end_site)
    # Perform SVD across the bond
    row_inds = uniqueinds(psi[end_site], psi[end_site+1])
    U, S_vals, V = svd(psi[end_site], row_inds)
    
    S_EE = 0.0
    for λ in diag(S_vals)
        p = λ^2
        if p > 1e-12
            S_EE -= p * log(p)
        end
    end

    # --- 2. Particle Number Fluctuations (F) ---
    # F = <N_A^2> - <N_A>^2
    # Expanding this: F = sum_{i,j in A} (<ni nj> - <ni><nj>)
    
    # We first measure the Correlation Matrix for the subsystem
    # This returns C[i,j] = <ni nj>
    # Note: For fermions, order matters, but for density operators ni and nj commute.
    NM = correlation_matrix(psi, "N", "N"; sites=subsystem_inds)
    
    N_A_mean = 0.0
    N_A_sq_mean = 0.0
    
    # Sum the elements of the correlation matrix
    # The matrix NM indices run from 1 to La_width
    for i in 1:La_width
        for j in 1:La_width
            N_A_sq_mean += NM[i, j]
        end
        # The diagonal elements NM[i,i] are <ni^2> = <ni> for fermions
        N_A_mean += NM[i, i]
    end

    F = N_A_sq_mean - (N_A_mean)^2

    return S_EE, F
end

function main_scan()
    # System Parameters
    N = 600
    La_width = 300 
    
    sites = siteinds("Fermion", N; conserve_qns=true) 

    # Scan parameters
    V_values = 0.0:0.2:4.0

    # ==========================================
    # EXPERIMENT 1: SSH
    # ==========================================
    println("\n=== EXPERIMENT 1: Topological SSH (Saving to data_SSH.csv) ===")
    
    open("data_SSH.csv", "w") do file
        write(file, "V,S,F\n")
        
        for V_int in V_values
            H = build_SSH_MPO_PBC(sites, v=1.0, w=1.3, Δ=0.0, V=V_int)

            state_pattern = [isodd(n) ? "1" : "0" for n in 1:N] 
            psi0 = randomMPS(sites, state_pattern, linkdims=10)

            # --- FIXED DMRG CALL ---
            energy, psi = dmrg(H, psi0; nsweeps=6, maxdim=[20, 60, 100, 200], cutoff=1E-10, outputlevel=0)

            S, F = measure_observables(psi, N, La_width)
            
            write(file, "$V_int, $S, $F\n")
            @printf("SSH | V=%-4.2f | S=%-8.5f | F=%-8.5f\n", V_int, S, F)
        end
    end 

    # ==========================================
    # EXPERIMENT 2: Rice-Mele
    # ==========================================
    println("\n=== EXPERIMENT 2: Trivial Rice-Mele (Saving to data_RM.csv) ===")
    
    open("data_RM.csv", "w") do file
        write(file, "V,S,F\n")
        
        for V_int in V_values
            H = build_SSH_MPO_PBC(sites, v=1.0, w=1.0, Δ=0.5, V=V_int)

            state_pattern = [isodd(n) ? "1" : "0" for n in 1:N] 
            psi0 = randomMPS(sites, state_pattern, linkdims=10)

            # --- FIXED DMRG CALL ---
            energy, psi = dmrg(H, psi0; nsweeps=6, maxdim=[20, 60, 100, 200], cutoff=1E-10, outputlevel=0)

            S, F = measure_observables(psi, N, La_width)
            
            write(file, "$V_int, $S, $F\n")
            @printf("RM  | V=%-4.2f | S=%-8.5f | F=%-8.5f\n", V_int, S, F)
        end
    end
    
    println("\nData collection complete.")
end

# Run it
main_scan()
