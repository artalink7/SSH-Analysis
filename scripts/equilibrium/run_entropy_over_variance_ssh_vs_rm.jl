push!(LOAD_PATH, joinpath(@__DIR__,"..","..","src"))
using SSHAnalysis
using LinearAlgebra
using ProgressMeter
using Base.Threads
using DelimitedFiles

gap = 10 .^ range(-4, -1; length=20)
N_cells = 5000  
#Entropy_SSH = zeros(Float64, length(gap))
#Variance_SSH = zeros(Float64, length(gap))
Ratio_SSH = zeros(Float64, length(gap))
#Entropy_RM = zeros(Float64, length(gap))
#Variance_RM = zeros(Float64, length(gap))
#Ratio_RM = zeros(Float64, length(gap))
#p = Progress(length(gap))
Threads.@threads for i in eachindex(gap)
    r = gap[i]
    # Entropy_SSH[i] = EntanglementEntropy([0, 1 - r/2, 1, 0], [0, 1 - r/2, 1, 0], "discrete", N_cells, 2*N_cells, "1", 0)
    # Variance_SSH[i] = SingleVariance([0, 1 - r/2, 1, 0], [0, 1 - r/2, 1, 0], "discrete", N_cells, 2*N_cells, "1",0)
    # Ratio_SSH[i] = Entropy_SSH[i] / Variance_SSH[i]
    # Entropy_RM[i] = EntanglementEntropy([0, 1 , 1, r/2], [0, 1 , 1, r/2], "discrete", N_cells, 2*N_cells, "1", 0)
    # Variance_RM[i] = SingleVariance([0, 1 , 1, r/2], [0, 1 , 1, r/2], "discrete", N_cells, 2*N_cells, "1", 0)
    # Ratio_RM[i] = Entropy_RM[i] / Variance_RM[i]
    S, V = EntropyandVariance_sub([0, 1 - r/2, 1, 0], [0, 1- r/2, 1, 0], "discrete", N_cells, 2*N_cells, 0)
    Ratio_SSH[i] = S / V
    #next!(p)
end 

data_dir = joinpath(ProjectRoot(), "data", "equilibrium")
SaveWithTimeStamp(Ratio_SSH, data_dir, "ratio_entropy_variance_vs_gap_ssh_1e-4_half_filling_5000cells")
# SaveWithTimeStamp(Ratio_RM, data_dir, "ratio_entropy_variance_vs_gap_rm_1e-4_half_filling")
# SaveWithTimeStamp(Entropy_SSH, data_dir, "entropy_vs_gap_ssh_1e-4_half_filling")
# SaveWithTimeStamp(Variance_SSH, data_dir, "variance_vs_gap_ssh_1e-4_half_filling")
# SaveWithTimeStamp(Entropy_RM, data_dir, "entropy_vs_gap_rm_1e-4_half_filling")
# SaveWithTimeStamp(Variance_RM, data_dir, "variance_vs_gap_rm_1e-4_half_filling")