push!(LOAD_PATH, joinpath(@__DIR__,"..","..","src"))
using SSHAnalysis

# =========================================================
# 1. SETUP FUNCTIONS (Standard)
# =========================================================

function measure_fluctuations_vs_ell(psi)
    N = length(psi)
    Nc = N ÷ 2 # Total number of unit cells in the chain
    
    # 1. Get the 2-point correlation matrix and density
    NN_matrix = correlation_matrix(psi, "N", "N")
    N_expected = expect(psi, "N")
    
    fluctuations = Dict{Int, Float64}()
    
    # 2. Determine the starting subsystem size (in cells) based on chain parity
    # If Nc is even, start with 2 cells. If Nc is odd, start with 1 cell.
    k_start = (Nc % 2 == 0) ? 2 : 1
    
    # We stop at Nc - 2 to always leave at least one cell on both the left and right boundaries
    k_max = Nc - 2 
    
    # 3. Expand the subsystem symmetrically (adding 2 cells / 4 sites at a time)
    for k in k_start:2:k_max
        ell = k * 2 # Convert number of subsystem cells back to number of sites
        
        # Calculate exactly which cells make up the subsystem
        c_start = (Nc - k) ÷ 2 + 1
        c_end   = c_start + k - 1
        
        # Map the cell indices to the exact physical site indices
        # Cell c always occupies sites (2c - 1) and (2c)
        start_site = 2 * c_start - 1
        end_site   = 2 * c_end
        
        # Sum the sub-block of the correlation matrix and density vector
        N_val = sum(N_expected[start_site:end_site])
        N_sq_val = sum(NN_matrix[start_site:end_site, start_site:end_site])
        
        # Variance: <N_A^2> - <N_A>^2
        fluctuations[ell] = N_sq_val - (N_val)^2
    end
    println("Calculated fluctuations for subsystem sizes: ", sort(collect(keys(fluctuations))))
    
    return fluctuations
end
# =========================================================
# 2. MAIN EXPERIMENT: Scan fluctuations in the gapped phase of the interacting SSH model
# =========================================================

function scan_fluctuations_gapped_phase_interacting()
    N = 100
    v_fixed = 1.0
    w_fixed = 3.0
    Δ_fixed = 0.0
    V = 6.0

    sites = siteinds("Fermion", N; conserve_qns=true)
    Nf0 = N ÷ 2 
    st0 = product_state_Nf(N, Nf0)

    open("data/equilibrium/scan_fluctuations_gapped_SSH_topological_V6.csv", "w") do file 
        write(file, "ell,fluctuations\n")
        
        H = build_SSH_MPO_OBC(sites; v=v_fixed, w=w_fixed, Δ=Δ_fixed, V=V)
        E0 , psi = ground_energy(H, sites; init_state=st0)
        
        fluctuations = measure_fluctuations_vs_ell(psi)
        
        # Safely extract and sort the valid subsystem sizes
        ell_vals = sort(collect(keys(fluctuations)))
        
        for ell in ell_vals
            write(file, "$ell, $(fluctuations[ell])\n")
        end
    end
end

# Run the scan
scan_fluctuations_gapped_phase_interacting()