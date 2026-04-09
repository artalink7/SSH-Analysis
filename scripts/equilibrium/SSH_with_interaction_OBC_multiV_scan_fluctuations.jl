push!(LOAD_PATH, joinpath(@__DIR__,"..","..","src"))
using SSHAnalysis

# =========================================================
# 1. SETUP FUNCTIONS (Standard)
# =========================================================

function measure_fluctuations_vs_ell(psi)
    N = length(psi)
    
    # 1. Get the full 2-point correlation matrix: <N_i N_j>
    NN_matrix = correlation_matrix(psi, "N", "N")
    
    # 2. Get the single particle density: <N_i>
    N_expected = expect(psi, "N")

    for i in 1:N
        NN_matrix[i, i] = N_expected[i]
    end
    
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
# 2. MAIN EXPERIMENT: "Charge Gap Scan"
# =========================================================

function scan_fluctuations_gapped_phase()
    N = 100
    v_fixed = 1.0
    w_fixed = 1.0
    Δ_fixed = 0.0
    V_values = [4.0, 5.0, 6.0, 7.0]

    sites = siteinds("Fermion", N; conserve_qns=true)
    Nf0 = N ÷ 2
    st0 = product_state_Nf(N, Nf0)

    open("data/equilibrium/SSH_with_interaction_OBC_multiV_scan_fluctuations_v_equal_w.csv", "w") do file 
        write(file, "V,ell,fluctuations\n")

        for V in V_values
            @printf("=== V = %.2f ===\n", V)

            H = build_SSH_MPO_OBC(sites; v=v_fixed, w=w_fixed, Δ=Δ_fixed, V=V)
            E0 , psi = ground_energy(H, sites; init_state=st0)

            fluctuations = measure_fluctuations_vs_ell(psi)

            for ell in 1:N-1
                write(file, "$V, $ell, $(fluctuations[ell])\n")
                @printf("V=%.2f | ell=%d | fluctuations=%.4e\n", V, ell, fluctuations[ell])
            end
        end
    end
end

# Run the scan
scan_fluctuations_gapped_phase()