push!(LOAD_PATH, joinpath(@__DIR__, "..","..", "src"))
using SSHAnalysis

# Define parameters
style = "continuos"
N_cells = 1
w = 0.400:0.01:2.000

# Generate Data
entropy = zeros(length(w))
@showprogress for (i, w_val) in enumerate(w)
    entropy_equilibrium = GenerateDataEntropyEquilibrium([0, w_val, 1, 0.0], style, N_cells)
    entropy[i] = entropy_equilibrium[1]
end

# Save Data
data_dir = joinpath(ProjectRoot(), "data", "equilibrium")
SaveWithTimeStamp(entropy, data_dir, "entropy_1_cell_thermodynamic_comparison_hop_coeff")

