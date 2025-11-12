push!(LOAD_PATH, joinpath(@__DIR__, "..", "..", "src"))
using SSHAnalysis
using LinearAlgebra
using StatsPlots
using LsqFit
using ProgressMeter
using Base.Threads

gap = 1e-10
pre_quench_ssh = [0, 1 - gap/2, 1, 0]
post_quench_ssh = pre_quench_ssh
t = 0
L_cells = [200, 400, 600, 800, 1000, 1300, 1500, 2000, 2500, 3000]
Ratio = zeros(Float64, length(L_cells))


Threads.@threads for i in eachindex(L_cells)
    n = L_cells[i]
    Entropy = EntanglementEntropy(pre_quench_ssh, post_quench_ssh, "discrete", n, 2n, "1", t)
    Variance = SingleVariance(pre_quench_ssh, post_quench_ssh, "discrete", n, 2n, "1", t)
    Ratio[i] = Entropy / Variance
end
data_dir = joinpath(ProjectRoot(), "data", "equilibrium")
SaveWithTimeStamp(Ratio, data_dir, "ratio_entropy_variance_vs_gap_ssh_1e-10_for_fit")

# N_data::Vector, R_data::Vector already computed
x = 1.0 ./ log.(L_cells)
X = hcat(ones(length(x)), x, x.^2)   # design matrix [1  x  x^2]
β = X \ Ratio                       # linear least-squares solution -> [A, B, C]
resid = Ratio - X * β
dof = length(Ratio) - length(β)
σ2 = sum(resid.^2) / dof
covβ = σ2 * inv(X' * X)              # covariance matrix of β
stderr = sqrt.(diag(covβ))

A_est, B_est, C_est = β
println("A = $A_est ± $(stderr[1])")
println("B = $B_est ± $(stderr[2])")
println("C = $C_est ± $(stderr[3])")
println("Residual norm: ", norm(resid))
