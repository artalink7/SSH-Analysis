push!(LOAD_PATH, joinpath(@__DIR__, "..","..", "src"))
using SSHAnalysis

# Saving path
data_dir = joinpath(ProjectRoot(), "data", "quench")


# Define parameters

pre_quench  = (φ=0, w=1.0, v=1.0, V=-1.0)
post_quench = (φ=0, w=1.0, v=1.0, V=1.0)
style = "discrete"
L_cells = 200
N_cells = 400

# Generate Data
open(data_dir * "/run_quench_on_stag_change_sign_RM_400cells.csv", "w") do file 
    write(file , "Time,Entropy,Variance\n")
    time = 0:0.4:199.6
    for t in time
        S, V = EntropyandVariance_sub(pre_quench, post_quench, style, L_cells, N_cells, t)
        write(file, "$t,$S,$V\n")
        @printf("t= %.1f, S_A= %.4f, ΔN_A = %.4f\n", t, S, V)
    end
end

