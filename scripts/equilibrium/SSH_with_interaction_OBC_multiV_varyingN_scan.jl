push!(LOAD_PATH, joinpath(@__DIR__,"..","..","src"))
using SSHAnalysis

function scan_charge_gap_vs_N_multiV()
    v_fixed = 1.0
    w_fixed = 1.0
    Δ_fixed = 0.0
     

    Ns = 310:20:610
    V_vals = [3.0, 4.0, 5.0, 6.0]

    open("data/equilibrium/SSH_gapped_criticalparam_charge_gap_vs_N_multiV.csv", "w") do file 
        write(file, "V,N,ChargeGap\n")
        for V in V_vals
            @printf("=== V = %.1f === \n", V)

            for N in Ns
                
                sites = siteinds("Fermion", N; conserve_qns=true)
                st0, stp, stm = create_CDW_states(N)
                
                H = build_SSH_MPO_OBC(sites; v=v_fixed, w=w_fixed, Δ=Δ_fixed, V=V)
                
                E0 = run_dmrg_cdw(H, sites, st0)
                Ep = run_dmrg_cdw(H, sites, stp)
                Em = run_dmrg_cdw(H, sites, stm)
                
                charge_gap = (Ep + Em - 2.0 * E0)

                write(file, "$V, $N, $charge_gap\n")
                @printf("V=%.1f | N=%d | ChargeGap=%.6e\n", V, N, charge_gap)
            end
        end
    end
end

scan_charge_gap_vs_N_multiV()