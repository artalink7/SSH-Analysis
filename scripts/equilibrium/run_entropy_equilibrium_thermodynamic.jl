push!(LOAD_PATH, joinpath(@__DIR__, "..","..", "src"))
using SSHAnalysis

# Define parameters

pre_quench  = (φ=0, w=0.7, v=1, V=0.0)
style = "continuos"
N_cells = 20


# Generate Data
entropy = GenerateDataEntropyEquilibrium(pre_quench, style, N_cells)

# Save Data
data_dir = joinpath(ProjectRoot(), "data", "equilibrium")
SaveWithTimeStamp(entropy, data_dir, "entropy_10_cells_thermodynamic")

