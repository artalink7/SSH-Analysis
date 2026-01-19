push!(LOAD_PATH, joinpath(@__DIR__,"..","..","src"))
using SSHAnalysis

using LinearAlgebra
using DelimitedFiles
using Printf
using Base.Threads

const CFT_RATIO = π^2 / 3

function SSH_params_fixed(δ)
    t1 = 1.0 - δ
    t2 = 1.0
    pre  = [0, t1, t2, 0]
    post = pre
    return pre, post
end

Lsubs  = [200, 400, 800]
deltas = [0.00025, 0.0005, 0.001, 0.002, 0.004]

Nk_factor = 2
style = "discrete"
t = 0.0

results = Dict{Int, Vector{Tuple{Float64, Float64}}}()

for L_sub in Lsubs
    println("\nRunning L_sub = $L_sub")
    N_k_points = Nk_factor * L_sub

    data = Vector{Tuple{Float64, Float64}}()

    for δ in deltas
        Δ = 2δ
        x = Δ * L_sub

        if x > 2.0
            continue
        end

        pre, post = SSH_params_fixed(δ)

        S, V = EntropyandVariance_sub(pre, post, style, L_sub, N_k_points, t)

        ratio = S / V
        deviation = ratio - CFT_RATIO

        push!(data, (x, deviation))

        @printf("  δ = %.4e | ΔL = %.3f | S/V = %.6f | dev = %.3e\n",
                δ, x, ratio, deviation)
    end

    results[L_sub] = data
end

outdir = "data/ssh_fixed_hopping_scaling"
mkpath(outdir)

for (L_sub, data) in results
    fname = joinpath(outdir, "ssh_scaling_fixed_t2_Lsub_$(L_sub).dat")
    writedlm(fname, data)
end
println("\nAll done!")