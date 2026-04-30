push!(LOAD_PATH, joinpath(@__DIR__, "..","..", "src"))
using SSHAnalysis

data_dir = joinpath(ProjectRoot(), "data", "quench")


# Define parameters

pre_quench  = (φ=0, w=0.7, v=1, V=0.0)
post_quench = (φ=0, w=1, v=0.7, V=0.0)
style = "discrete"
L_cells = 400
N_cells = 800

# Generate Data

open(data_dir * "/run_quench_gapped_to_gapped_tracking_heat_SSH_400cells_200sec.csv", "w") do file 
    write(file , "Time,Entropy,Fluctuations,Energy,EnergyDiff\n")
    time = 0:0.4:199.6
    E_0_plus = ExpectationValue_HamiltonianSub(pre_quench, post_quench, style, L_cells, N_cells, 0.0)
    for t in time
        S, V, E_exp = EnergyFluctuationsEntropy(pre_quench, post_quench, style, L_cells, N_cells, t)
        energy_diff = E_exp - E_0_plus
        write(file, "$t,$S,$V,$E_exp,$energy_diff\n")
        @printf("t= %.1f, Entropy = %.4f, Fluctuations = %.4f, Energy = %.4f, EnergyDiff = %.4f\n", t, S, V, E_exp, energy_diff)
    end
end
