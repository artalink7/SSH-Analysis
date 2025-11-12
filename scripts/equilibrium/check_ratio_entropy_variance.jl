push!(LOAD_PATH, joinpath(@__DIR__,"..","..","src"))
using SSHAnalysis
using LinearAlgebra
using ProgressMeter
using DelimitedFiles
using Printf

using Plots
π2over3 = π^2/3

# Inputs you already produced
gap = 10 .^ range(-8, -2; length=30)
Ratio_SSH = vec(readdlm("data/equilibrium/ratio_entropy_variance_vs_gap_ssh_filling40percent.txt"))
Ratio_RM = vec(readdlm("data/equilibrium/ratio_entropy_variance_vs_gap_rm_filling40percent.txt"))
ℓ = 100

# correlation lengths (exact for your parametrizations)
xi_RM = 4 ./(gap)               # since t1=t2=1, u=gap/2 -> xi=4/gap
xi_SSH = 4 ./(gap) .- 1         # since t1=1-gap/2, t2=1 -> xi = 4/gap - 1

# dimensionless scaling variable
x_RM = ℓ ./ xi_RM
x_SSH = ℓ ./ xi_SSH

# Plot S/Var vs ℓ/ξ
p1 = plot(x_RM, Ratio_RM, label="RM", xlabel="ℓ/ξ", ylabel="S/Var", legend=:topleft)
plot!(p1, x_SSH, Ratio_SSH, label="SSH", marker=:diamond)
hline!(p1, [π2over3], linestyle=:dash, label="π²/3")

# Plot deviation from π²/3 (percent)
dev_RM = 100 .* (Ratio_RM .- π2over3) ./ π2over3
dev_SSH = 100 .* (Ratio_SSH .- π2over3) ./ π2over3
p2 = plot(x_RM, dev_RM, label="RM % dev", xlabel="ℓ/ξ", ylabel="% deviation", legend=:topleft)
plot!(p2, x_SSH, dev_SSH, label="SSH % dev", marker=:diamond)
hline!(p2, [0.0], linestyle=:dash, color=:black)

# show both
plot(p1, p2, layout=(2,1), size=(800,700))

# Print summary table to console
println("gap\tℓ/ξ_RM\tℓ/ξ_SSH\tRatio_RM\tRatio_SSH\t%dev_RM\t%dev_SSH")
for i in eachindex(gap)
    @printf("%1.1e\t%0.4f\t%0.4f\t%0.6f\t%0.6f\t%0.3f\t%0.3f\n",
            gap[i], x_RM[i], x_SSH[i], Ratio_RM[i], Ratio_SSH[i], dev_RM[i], dev_SSH[i])
end

savefig(joinpath(@__DIR__, "check_ratio_entropy_variance_ssh_and_rm_percent_filling.pdf"))