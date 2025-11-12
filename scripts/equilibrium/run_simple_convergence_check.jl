push!(LOAD_PATH, joinpath(@__DIR__,"..","..","src"))
using SSHAnalysis
using LinearAlgebra
using ProgressMeter
using Printf
using Statistics
using LinearAlgebra
using Base.Threads
using Plots
using DelimitedFiles
using FFTW
using Base.Threads

#-------------------
# Replace these with your exact data
# If you printed them while running the direct-eigen code, collect them here:
# e.g. ℓs = [5000, 10000, 15000, 20000]
#      ratios_exact = [3.261761619082826, 3.2634585495259536, 3.2643594435181207, 3.264962244337278]
ℓs = [5000, 10000, 15000, 20000]               # <- replace/extend if you have more
ratios_exact = [3.261761619082826, 3.2634585495259536, 3.2643594435181207, 3.264962244337278]  # <- your exact results

# -----------------------------
# Build design matrix for model: R(ℓ) = A + B / ln(ℓ) + C / ln(ℓ)^2
# Linear in [A, B, C]
xs = log.(ℓs)
X = [ones(length(xs)) 1.0 ./ xs 1.0 ./ (xs.^2)]
y = ratios_exact

# Least squares solve
coeffs = X \ y  # [A, B, C]
A_est, B_est, C_est = coeffs

# Predict and compute residuals
y_pred = X * coeffs
rms = sqrt(mean((y - y_pred).^2))

println("Fitted asymptote A = $A_est")
println("B = $B_est, C = $C_est")
println("RMS residual = $rms")
println("π^2 / 3 = $(π^2/3)")
println("Difference (A - π^2/3) = $(A_est - π^2/3)")

# -----------------------------
# Optional: plot data + fit (use many points for smooth curve)
ℓs_dense = collect(minimum(ℓs):100:maximum(ℓs))
xs_dense = log.(ℓs_dense)
X_dense = [ones(length(xs_dense)) 1.0 ./ xs_dense 1.0 ./ (xs_dense.^2)]
y_dense = X_dense * coeffs

plot(ℓs, ratios_exact, seriestype=:scatter, label="exact data", xlabel="ℓ", ylabel="S/Var(N)", legend=:bottomright)
plot!(ℓs_dense, y_dense, lw=2, label="fit: A + B/lnℓ + C/(lnℓ)^2")
hline!([π^2/3], linestyle=:dash, color=:red, label="π^2/3")
title!("Extrapolated asymptote A ≈ $(round(A_est; digits=6))")
savefig("S_over_Var_extrapolation.pdf")
println("Plot saved as S_over_Var_extrapolation.pdf")
