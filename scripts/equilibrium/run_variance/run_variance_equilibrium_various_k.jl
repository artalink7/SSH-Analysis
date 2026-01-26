push!(LOAD_PATH, joinpath(@__DIR__,"..","..","...","src"))
using SSHAnalysis

# Define parameters
w_vals = vcat(collect(0.5:0.1:0.9), collect(0.95:0.01:0.99)) # various intercell hoppings
φ=0    # phase of the flux
v=1    # intracell hopping
V=0.0  # staggering potential 
style = "discrete"
N_cells = 200
data_dir = joinpath(ProjectRoot(), "data", "equilibrium")


@showprogress 1 "Computing entropies..." for  w in w_vals
    pre_quench = [φ, w, v, V]    
    # Generate Data
    variance = GenerateDataVarianceEquilibrium(pre_quench, style, N_cells)
    # Save Data
    SaveWithTimeStamp(variance, data_dir, "variance_200_cells_noflux_k$(round(w,digits=3))")
end
