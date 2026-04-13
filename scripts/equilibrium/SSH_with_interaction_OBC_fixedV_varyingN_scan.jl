push!(LOAD_PATH, joinpath(@__DIR__,"..","..","src"))
using SSHAnalysis

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

function run_dmrg_cdw(H, sites, st_array)
    # randomMPS introduces immediate quantum fluctuations to break symmetry
    psi0 = randomMPS(sites, st_array) 
    
    # Aggressive noise schedule to unpin domain walls!
    sweeps = Sweeps(15)
    maxdim!(sweeps, 50, 100, 200, 400, 800)
    cutoff!(sweeps, 1E-10)
    noise!(sweeps, 1E-3, 1E-4, 1E-4, 1E-5, 1E-5, 1E-6, 1E-7, 1E-8, 0.0)
    
    energy, psi = dmrg(H, psi0, sweeps; outputlevel=0)
    return energy
end

function scan_charge_gap_vs_N_fixed()
    v_fixed = 1.0
    w_fixed = 1.0
    Δ_fixed = 0.0
    V_fixed = 4.0 

    Ns = 310:20:610

    open("data/equilibrium/SSH_gapped_charge_gap_vs_N.csv", "w") do file 
        write(file, "N,ChargeGap\n")

        for N in Ns 
            @printf("=== N = %d ===\n", N)

            sites = siteinds("Fermion", N; conserve_qns=true)
            st0, stp, stm = create_CDW_states(N)

            H = build_SSH_MPO_OBC(sites; v=v_fixed, w=w_fixed, Δ=Δ_fixed, V=V_fixed)
            
            E0 = run_dmrg_cdw(H, sites, st0)
            Ep = run_dmrg_cdw(H, sites, stp)
            Em = run_dmrg_cdw(H, sites, stm)

            charge_gap = (Ep + Em - 2.0 * E0)

            write(file, "$N, $charge_gap\n")
            @printf("N=%d | ChargeGap=%.6e\n", N, charge_gap)
        end
    end
end

scan_charge_gap_vs_N_fixed()