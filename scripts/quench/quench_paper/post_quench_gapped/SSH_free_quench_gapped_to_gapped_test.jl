push!(LOAD_PATH, joinpath(@__DIR__, "..","..","..","..", "src"))
using SSHAnalysis

# Saving path
data_dir = joinpath(ProjectRoot(), "data", "quench", "quench_paper")


# Define parameters

style = "discrete"
L_cells = 200
N_cells = 400
 

pre_quench  = (φ=0, w=1.0, v=0.5, V=0.0)
post_quench = (φ=0, w=0.5, v=1.0, V=0.0)


# Generate Data
open(data_dir * "/run_SSH_free_quench_gapped_to_gapped_no_flux_L200cells_100sec.csv", "w") do file 
    write(file , "Time,Entropy,Variance\n")
    time = 0:0.05:100
    for t in time
        entropy, variance = EntropyandVariance_sub(pre_quench, post_quench, style, L_cells, N_cells, t) 
        write(file, "$t,$entropy,$variance\n")
        @printf("t= %.1f, S_A= %.4f, ΔN_A = %.4f\n", t, entropy, variance)
    end
end

