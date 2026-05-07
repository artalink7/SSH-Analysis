push!(LOAD_PATH, joinpath(@__DIR__, "..","..","..","..", "src"))
using SSHAnalysis

# Saving path
data_dir = joinpath(ProjectRoot(), "data", "quench", "quench_paper")


# Define parameters

style = "discrete"
L_cells = 100
N_cells = 200
flux_param = π / N_cells

pre_quench  = (φ=flux_param, w=1.0, v=1.0, V=0.5)
post_quench = (φ=0, w=1.0, v=1.0, V=0.5)


# Generate Data
open(data_dir * "/run_RM_free_quench_gapped_with_flux_to_gapped_no_flux_L100cells_200sec.csv", "w") do file 
    write(file , "Time,Entropy,Variance\n")
    time = 0:0.1:200
    for t in time
        entropy, variance = EntropyandVariance_sub(pre_quench, post_quench, style, L_cells, N_cells, t) 
        write(file, "$t,$entropy,$variance\n")
        @printf("t= %.1f, S_A= %.4f, ΔN_A = %.4f\n", t, entropy, variance)
    end
end

