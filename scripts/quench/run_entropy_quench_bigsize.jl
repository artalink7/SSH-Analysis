push!(LOAD_PATH, joinpath(@__DIR__, "..","..", "src"))
using SSHAnalysis

# Define parameters

pre_quench  = (φ=0, w=0.7, v=1, V=0.0)
post_quench = (φ=0, w=1, v=0.7, V=0.0)
style = "discrete"
L_cells = 100
N_cells = 200

# Generate Data
entropy = GenerateDataEntropy(pre_quench, post_quench, style, L_cells, N_cells)

# Save Data
data_dir = joinpath(ProjectRoot(), "data", "quench")
SaveWithTimeStamp(entropy, data_dir, "entropy_200_cells_noflux")

