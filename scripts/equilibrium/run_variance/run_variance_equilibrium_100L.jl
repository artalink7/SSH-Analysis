push!(LOAD_PATH, joinpath(@__DIR__,"..","..","..","src"))
using SSHAnalysis

# Define parameters

pre_quench  = (φ=0, w=3, v=1, V=0.0)
style = "discrete"
N_cells = 100

data_variance = zeros(Float64, N_cells)

# Generate Data
for n in 1:N_cells
    data_variance[n] = EntropyandVariance_sub(pre_quench, pre_quench, style, n, N_cells, 0.0)[2]
end

# Save Data
data_dir = joinpath(ProjectRoot(), "data", "equilibrium")
SaveWithTimeStamp(data_variance, data_dir, "variance_100_cells_try")

