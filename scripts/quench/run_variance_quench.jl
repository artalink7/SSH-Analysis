push!(LOAD_PATH, joinpath(@__DIR__, "..","..", "src"))
using SSHAnalysis

# Define parameters

pre_quench  = (φ=0, w=0.7, v=1, V=0.0)
post_quench = (φ=0, w=1, v=0.7, V=0.0)
style = "discrete"
L_cells = 20
N_cells = 40

# Generate Data
variance = GenerateDataVariance(pre_quench, post_quench, style, L_cells, N_cells)

# Save Data
data_dir = joinpath(ProjectRoot(), "data", "quench")
SaveWithTimeStamp(variance, data_dir, "variance_40_cells_noflux")

