push!(LOAD_PATH, joinpath(@__DIR__,"..","..","src"))
using SSHAnalysis

# =========================================================
# 1. MAIN EXPERIMENT: "Charge Gap Scan"
# =========================================================

function scan_charge_gap_multiN()
    v_fixed = 1.0
    w_fixed = 1.0
    Δ_fixed = 0.0

    Ns = [300, 400, 500, 600]
    V_values = 0.0:0.2:4.0

    open("data/equilibrium/charge_gap_multiN2_probing_phase_transition.csv", "w") do file 
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