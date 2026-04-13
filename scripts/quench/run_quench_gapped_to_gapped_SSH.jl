push!(LOAD_PATH, joinpath(@__DIR__, "..","..", "src"))
using SSHAnalysis

data_dir = joinpath(ProjectRoot(), "data", "quench")


# Define parameters

pre_quench  = (φ=0, w=0.7, v=1, V=0.0)
post_quench = (φ=0, w=1, v=0.7, V=0.0)
style = "discrete"
L_cells = 200
N_cells = 400

# Generate Data

open(data_dir * "/run_quench_gapped_to_gapped_SSH_400cells_200sec.csv", "w") do file 
    write(file , "Time,Entropy,Variance\n")
    time = 0:0.4:199.6
    for t in time
        S, V = EntropyandVariance_sub(pre_quench, post_quench, style, L_cells, N_cells, t)
        write(file, "$t,$S,$V\n")
        @printf("t= %.1f, S_A= %.4f, ΔN_A = %.4f\n", t, S, V)
    end
end
