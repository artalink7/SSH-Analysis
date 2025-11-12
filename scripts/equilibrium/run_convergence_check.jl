push!(LOAD_PATH, joinpath(@__DIR__,"..","..","src"))
using SSHAnalysis
using LinearAlgebra
using ProgressMeter
using Printf

# === PARAMETERS ===
r = 1e-8                      # very small mass → metallic limit
t = 0.0                       # equilibrium (pre-quench ground state)
filling_fraction = 1.0        # fully fill lower band
style = "discrete"            # use discrete-k integration

# SSH parameters
pre_quench = [0.0, 1 - r/2, 1.0, 0.0]
post_quench = pre_quench      # irrelevant at t=0

# number of cells (momentum resolution) and subsystem size
N_cells_list = 2*[100, 200, 400, 800, 1600]
L_over_N = 0.1               # fraction of total system for subsystem
println("Testing convergence toward π²/3 ≈ ", π^2/3, "\n")

for N_cells in N_cells_list
    L_cells = round(Int, L_over_N * N_cells)

    S = EntanglementEntropy(pre_quench, post_quench, style, L_cells, N_cells, "1", t)
    Var = SingleVariance(pre_quench, post_quench, style, L_cells, N_cells, "1", t)
    ratio = S / Var

    @printf "N = %4d, L = %4d → S/Var = %.6f (deviation = %.3e)\n" N_cells L_cells ratio abs(ratio - π^2/3)
end
