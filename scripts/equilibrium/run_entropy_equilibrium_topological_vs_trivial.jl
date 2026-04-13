push!(LOAD_PATH, joinpath(@__DIR__,"..","..","src"))
using SSHAnalysis

# Define parameters
w_1 =  0.2
w_2 = 1.5
φ=0    # phase of the flux
v=1    # intracell hopping
V=0.0  # staggering potential 
style = "discrete"
N_cells = 400
data_dir = joinpath(ProjectRoot(), "data", "equilibrium")

topological_params = [φ, w_2, v, V]
trivial_params = [φ, w_1, v, V]

# Generate data for topological and trivial cases
eigen_corr_topological = EigenvaluesDensity_sub(topological_params, [0, 1, 1, 0], style, div(N_cells, 2), N_cells, 0.0)
eigen_corr_trivial = EigenvaluesDensity_sub(trivial_params, [0, 1, 1, 0], style, div(N_cells, 2), N_cells, 0.0)



data_dir = joinpath(ProjectRoot(), "data", "equilibrium")
SaveWithTimeStamp(eigen_corr_topological, data_dir, "eigen_corr_topological")
SaveWithTimeStamp(eigen_corr_trivial, data_dir, "eigen_corr_trivial")




