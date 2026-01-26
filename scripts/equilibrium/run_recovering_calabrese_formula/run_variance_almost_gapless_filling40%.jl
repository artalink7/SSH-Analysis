push!(LOAD_PATH, joinpath(@__DIR__,"..","..","..","src"))
using SSHAnalysis
using LinearAlgebra
using StatsPlots    
using ProgressMeter
using Plots
using DelimitedFiles

N_cells = 200:40:800
Variance = zeros(Float64, length(N_cells))
p = Progress(length(N_cells))
for (i, n) in enumerate(N_cells)
    Variance[i] = SingleVariance([0, 1, 0.99999, 0], [0, 0.9999, 1, 0], "discrete", n, 2n, "1", 0; filling_fraction=0.4)
    next!(p)
end

# Save Data
data_dir = joinpath(ProjectRoot(), "data", "equilibrium")
SaveWithTimeStamp(Variance, data_dir, "variance_vs_Ncells_almostgapless_filling40percent")


plot(N_cells, Variance; xscale = :log10 ,label="Variance", marker=:circle, markersize=3, linestyle=:solid)
xlabel!("N cells")
ylabel!("Variance")
title!("Variance vs N cells for gap ~ 0 - semilogx")
savefig("plot/variance_vs_Ncells_almostgapless_filling40percent.pdf")

