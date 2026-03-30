push!(LOAD_PATH, joinpath(@__DIR__, "..","..", "src"))
using SSHAnalysis

# Saving path
data_dir = joinpath(ProjectRoot(), "data", "quench")


# Define parameters

pre_quench  = (φ=0, w=0.5, v=1, V=0.0)
post_quench = (φ=0, w=0.5, v=1, V=1.0)
style = "discrete"
L_cells = 100
N_cells = 200

# Generate Data
open(data_dir * "/run_quench_on_stag_RM_200cells.csv", "w") do file 
    write(file , "Time,Entropy,Variance\n")
    time = 0:0.4:199.6
    for t in time
        entropy = EntanglementEntropy(pre_quench, post_quench, style, L_cells, N_cells, t)
        variance = SingleVariance(pre_quench, post_quench, style, L_cells, N_cells, t)
        write(file, "$t,$entropy,$variance\n")
        @printf("t= %.1f, S_A= %.4f, ΔN_A = %.4f\n", t, entropy, variance)
    end
end

