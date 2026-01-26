push!(LOAD_PATH, joinpath(@__DIR__,"..","..","..","src"))
using SSHAnalysis
using LinearAlgebra
using StatsPlots    
using ProgressMeter
using Plots
using DelimitedFiles

using DataFrames, GLM, Statistics

# Suppose you already have these:
N_cells = 200:40:800
Entropy = vec(readdlm("data/equilibrium/entropy_vs_Ncells_almostgapless_filling20percent.txt")) 

# 1. Create DataFrame with log(N_cells)
df = DataFrame(
    logN = log.(N_cells),     # natural log of N_cells
    Ent = Entropy
)

# 2. Fit a linear model: Ent = a * logN + b
model = lm(@formula(Ent ~ logN), df)

# 3. Show results
println(coef(model))          # [b, a]
println(r2(model))            # goodness of fit
display(coeftable(model))

# 4. Plot data + fitted line
scatter(df.logN, df.Ent; label="Data", xlabel="log(N_cells)", ylabel="Entropy")
plot!(df.logN, predict(model); label="Linear fit", lw=2)
savefig("plot/entropy_vs_logNcells_almostgapless_filling20percent_fit.pdf")
