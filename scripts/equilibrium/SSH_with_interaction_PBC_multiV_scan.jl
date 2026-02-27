using LinearAlgebra
using ITensors
using ITensorMPS
using Statistics
using ITensorGaussianMPS
using Printf
# =========================================================
# 1. SETUP FUNCTIONS 
# =========================================================

function build_SSH_MPO_PBC(sites; v=1.0, w=1.0, Δ=0.0, V=0.0)
    N = length(sites)
    os = OpSum()

    # Bulk Hopping terms
    for i in 1:2:N-1
        os += -v, "Cdag", i, "C", i+1; os += -v, "Cdag", i+1, "C", i
    end
    for i in 2:2:N-1
        os += -w, "Cdag", i, "C", i+1; os += -w, "Cdag", i+1, "C", i
    end

    # Boundary Hopping (PBC connection N -> 1)
    bond_type = (N % 2 == 0) ? w : v
    # Fermionic PBC requires careful handling of the boundary sign depending on particle number,
    # but ITensor handles the Jordan-Wigner strings automatically for OpSums.
    os += bond_type, "Cdag", N, "C", 1
    os += bond_type, "Cdag", 1, "C", N

    # Staggered Potential (Keep 0 for criticality)
    if Δ != 0.0
        for i in 1:N; os += Δ * (-1)^i, "N", i; end
    end

    # Bulk Nearest-neighbor interaction
    if V != 0.0
        for i in 1:N-1; os += V, "N", i, "N", i+1; end
        # Boundary Interaction
        os += V, "N", N, "N", 1
    end

    return MPO(os, sites)
end

function product_state_Nf(N::Int, Nf::Int; pattern::Symbol=:alternating)
    @assert 0 <= Nf <= N
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

function measure_fluctuations_vs_ell(psi)
    N = length(psi)
    
    # 1. Get the full 2-point correlation matrix: <N_i N_j>
    NN_matrix = correlation_matrix(psi, "N", "N")
    
    # 2. Get the single particle density: <N_i>
    N_expected = expect(psi, "N")
    
    fluctuations = zeros(Float64, N-1)
    
    # 3. Calculate fluctuation for each subsystem size ell
    for ell in 1:N-1
        # Sum the sub-block of the correlation matrix and density vector
        N_val = sum(N_expected[1:ell])
        N_sq_val = sum(NN_matrix[1:ell, 1:ell])
        
        # Variance: <N_A^2> - <N_A>^2
        fluctuations[ell] = N_sq_val - (N_val)^2
    end
    
    return fluctuations
end

# =========================================================
# 2. MAIN EXPERIMENT: Subsystem Fluctuations Scan
# =========================================================

function scan_fluctuations_PBC()
    N = 100 # Match the system size L=100 from the paper
    v_fixed = 1.0
    w_fixed = 1.0
    Δ_fixed = 0.0
    
    # V values corresponding to Delta = -0.2 to 0.9 in steps of 0.1
    # V = 2 * v * Delta 
    V_values = -0.4:0.4:0.8
    
    # We only need one set of sites for the whole scan
    sites = siteinds("Fermion", N; conserve_qns=true)
    Nf0 = N ÷ 2
    st0 = product_state_Nf(N, Nf0)
    
    # Open CSV for writing
    open("fluctuations_PBC_vs_ell_full.csv", "w") do file 
        write(file, "V,ell,Fluctuation\n")
        
        for V in V_values
            @printf("=== Running V = %.2f ===\n", V)
            
            # Build Hamiltonian
            H = build_SSH_MPO_PBC(sites; v=v_fixed, w=w_fixed, Δ=Δ_fixed, V=V)
            
            # Initialize state. Using randomMPS with the product state array 
            # helps avoid getting stuck in local minima with PBCs.
            psi0 = randomMPS(sites, st0)
            
            # Define DMRG parameters (PBC requires slightly more aggressive sweeps and noise)
            sweeps = Sweeps(10)
            maxdim!(sweeps, 50, 100, 200, 400, 800)
            cutoff!(sweeps, 1E-9)
            noise!(sweeps, 1E-5, 1E-6, 1E-7, 1E-8, 0.0)
            
            # Run DMRG
            energy, psi = dmrg(H, psi0, sweeps; outputlevel=1)
            
            # Measure fluctuations
            fluctuations = measure_fluctuations_vs_ell(psi)
            
            # Write results to file using proper Julia string interpolation
            for ell in 1:N-1
                write(file, "$V,$ell,$(fluctuations[ell])\n")
            end
        end
    end
    println("Scan complete! Data saved to fluctuations_PBC_vs_ell.csv")
end

# Run the experiment
scan_fluctuations_PBC()