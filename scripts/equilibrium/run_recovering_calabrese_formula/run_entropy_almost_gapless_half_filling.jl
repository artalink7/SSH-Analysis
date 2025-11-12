push!(LOAD_PATH, joinpath(@__DIR__,"..","..","..","src"))
using SSHAnalysis
using LinearAlgebra
using StatsPlots    
using ProgressMeter
using Plots
using DelimitedFiles

N_cells = 200:40:800
Entropy = zeros(Float64, length(N_cells))
p = Progress(length(N_cells))
for (i, n) in enumerate(N_cells)
    Entropy[i] = EntanglementEntropy([0, 1, 0.99999, 0], [0, 0.9999, 1, 0], "discrete", n, 2n, "1", 0)
    next!(p)
end

# Save Data
data_dir = joinpath(ProjectRoot(), "data", "equilibrium")
SaveWithTimeStamp(Entropy, data_dir, "entropy_vs_Ncells_almostgapless_half_filling")


plot(N_cells, Entropy; xscale = :log10 ,label="Entropy", marker=:circle, markersize=3, linestyle=:solid)
xlabel!("N cells")
ylabel!("Entropy")
title!("Entropy vs N cells for gap ~ 0 - semilogx")
savefig("plot/entropy_vs_Ncells_almostgapless_half_filling.pdf")

