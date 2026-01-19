push!(LOAD_PATH, joinpath(@__DIR__,"..","..","src"))
using SSHAnalysis

using LinearAlgebra
using DelimitedFiles
using Printf
using Base.Threads

# --- Option B driver: fix x0 = Δ * L_sub (small) and send L_sub -> ∞ ---
# Uses YOUR existing EntropyandVariance_sub(...) function.

push!(LOAD_PATH, joinpath(@__DIR__,"..","..","src"))
using SSHAnalysis
using LinearAlgebra
using DelimitedFiles
using Printf

const CFT_RATIO = π^2 / 3

# SSH with fixed t2 = 1 and variable t1 = 1 - δ
# Gap: Δ = 2δ
function SSH_params_fixed(δ)
    t1 = 1.0 - δ
    t2 = 1.0
    pre  = [0, t1, t2, 0]
    post = pre           # consistent with your t=0 "no quench" usage
    return pre, post
end

# ---------------- user knobs ----------------
x0 = 0.2                       # target scaling variable x = Δ L_sub (pick 0.1-0.3)
Lsubs = [100, 200, 400, 800]   # increase if feasible (e.g. add 1200, 1600)
Nk_factor = 4                  # N_k_points = Nk_factor * L_sub (use 4-8 for small gaps)
style = "discrete"
t = 0.0
# --------------------------------------------

# We will store rows: L_sub, delta, Delta, x, ratio, deviation
rows = Matrix{Float64}(undef, length(Lsubs), 6)

for (j, L_sub) in enumerate(Lsubs)
    Δ = x0 / L_sub          # enforce x = Δ L_sub = x0
    δ = Δ / 2               # because Δ = 2δ in this parametrization

    N_k_points = Nk_factor * L_sub

    pre, post = SSH_params_fixed(δ)
    S, V = EntropyandVariance_sub(pre, post, style, L_sub, N_k_points, t)

    ratio = S / V
    dev = ratio - CFT_RATIO

    rows[j, :] .= (L_sub, δ, Δ, x0, ratio, dev)

    @printf("L_sub=%4d | δ=%.4e | Δ=%.4e | x=%.3f | S/V=%.8f | dev=%.3e\n",
            L_sub, δ, Δ, x0, ratio, dev)
end

# Save
outdir = joinpath(ProjectRoot(), "data", "ssh_optionB_fixed_x")
mkpath(outdir)
outfile = joinpath(outdir, @sprintf("ssh_fixed_x0_%.3f.dat", x0))
writedlm(outfile, rows)

println("\nSaved: ", outfile)
println("Columns: L_sub, delta, Delta, x0, ratio, ratio_minus_pi2_over3")
