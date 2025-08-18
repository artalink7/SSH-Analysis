push!(LOAD_PATH, joinpath(@__DIR__, "..","..", "src"))
using SSHAnalysis

# Define parameters

pre_quench  = (φ=0, w=0.7, v=1, V=0.0)
post_quench = (φ=0, w=1, v=0.7, V=0.0)
style = "continuos"
N_cells = 200

# Computing the modulus squared of the correlation, as the distance between cells increases
correlation_list = zeros(Float64, 51)

@showprogress for (i, r) in enumerate(0:50)
    ρ = RealSpaceDensityM(pre_quench, post_quench, r, style, N_cells, 0.0)
    correlation_list[i] = abs(ρ[1, 2])
end

# Save Data
data_dir = joinpath(ProjectRoot(), "data", "equilibrium")
SaveWithTimeStamp(correlation_list, data_dir, "correlation_200_cells_noflux")
