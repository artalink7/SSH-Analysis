push!(LOAD_PATH, joinpath(@__DIR__, "..","..", "src"))
using SSHAnalysis

# Saving path
data_dir = joinpath(ProjectRoot(), "data", "quench")


# Define parameters

pre_quench  = (φ=0, w=1, v=1, V=4.0)
post_quench = (φ=0, w=1, v=1, V=0.0)
style = "discrete"
L_cells = 200
N_cells = 400

# Generate Data
open(data_dir * "/run_quench_on_stag_RM_gapped_to_gapless_400cells_1000sec.csv", "w") do file 
    write(file , "Time,Entropy,Variance\n")
    time = 0:1.0:1000.0
    for t in time
        entropy, variance = EntropyandVariance_sub(pre_quench, post_quench, style, L_cells, N_cells, t) 
        write(file, "$t,$entropy,$variance\n")
        @printf("t= %.1f, S_A= %.4f, ΔN_A = %.4f\n", t, entropy, variance)
    end
end

