push!(LOAD_PATH, joinpath(@__DIR__, "..", "..", "src"))
using SSHAnalysis
using LinearAlgebra
using StatsPlots
using LsqFit
using ProgressMeter

gap = 1e-10
pre_quench_ssh = [0, 1 , 1, gap/2]
post_quench_ssh = pre_quench_ssh
t = 0
L_cells = [200, 300, 400, 500, 600, 700, 800, 900, 1000]
Ratio = zeros(Float64, length(L_cells))
filling_fraction = 0.4

p = Progress(length(L_cells))

for (i, n) in enumerate(L_cells)
    Entropy = EntanglementEntropy(pre_quench_ssh, post_quench_ssh, "discrete", n, 2n, "1", t, filling_fraction = filling_fraction)
    Variance = SingleVariance(pre_quench_ssh, post_quench_ssh, "discrete", n, 2n, "1", t, filling_fraction= filling_fraction)
    Ratio[i] = Entropy / Variance
    next!(p)
end

# Fit R(N) = π²/3 + B / ln(N) + C / (ln(N))²
model(N, p) = (π)^2/3 .+ p[1] ./ log.(N) .+ p[2] ./ (log.(N)).^2
p0 = [1.0, 1.0]
fit = curve_fit(model, L_cells, Ratio, p0)

B, C = fit.param
println("Fit parameters: B = $(B), C = $(C)")
println("Residual norm: ", norm(fit.resid))
