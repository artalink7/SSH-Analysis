push!(LOAD_PATH, joinpath(@__DIR__,"..","..","src"))
using SSHAnalysis
using LinearAlgebra
using ProgressMeter
using Base.Threads
using DelimitedFiles

gap = vcat(
    [0.005, 0.008, 0.01, 0.02],
    10 .^ range(log10(0.02), log10(0.3), length=10)[2:end]
)
N_cells = 600
N_k_points = 20000  

Ratio_RM = zeros(Float64, length(gap))
Ratio_SSH = zeros(Float64, length(gap))
#p = Progress(length(gap))
Threads.@threads for i in eachindex(gap)
    r = gap[i]
    S_ssh, V_ssh = EntropyandVariance_sub([0, 1 - r/2, 1 , 0], [0, 0.5, 1, 0], "discrete", N_cells, N_k_points, 0)
    S_rm, V_rm = EntropyandVariance_sub([0, 1 , 1, r/2], [0, 1 , 1, r/2], "discrete", N_cells, N_k_points, 0)
    Ratio_RM[i] = S_rm / V_rm
    Ratio_SSH[i] = S_ssh / V_ssh
    #next!(p)
end 

data_dir = joinpath(ProjectRoot(), "data", "equilibrium")
SaveWithTimeStamp(Ratio_SSH, data_dir, "ratio_entropy_variance_vs_gap_ssh_600cells_newgap")
SaveWithTimeStamp(Ratio_RM, data_dir, "ratio_entropy_variance_vs_gap_rm_600cells_newgap")
println("Done!")
