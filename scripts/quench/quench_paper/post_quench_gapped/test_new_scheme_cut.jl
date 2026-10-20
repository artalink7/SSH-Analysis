using ITensors
using LinearAlgebra
using ITensorMPS
push!(LOAD_PATH, joinpath(@__DIR__, "..","..","..","..", "src"))
using SSHAnalysis

# --- assumes single_cut_entropy, bulk_entanglement_entropy, and
#     product_state_Nf are already defined/included from your codebase ---

# ---------------------------------------------------------------------
# Spectrum-extracting variants of both methods, so we compare the full
# entanglement spectrum, not just the entropy (a spectrum match is a much
# stronger check -- two different sign/ordering bugs could in principle
# cancel in the entropy sum but won't survive an eigenvalue-by-eigenvalue
# comparison).
# ---------------------------------------------------------------------

function single_cut_spectrum(psi::MPS, cut::Int)
    orthogonalize!(psi, cut)
    row_inds = uniqueinds(psi[cut], psi[cut + 1])
    _, S_vals, _ = svd(psi[cut], row_inds)
    λs = [S_vals[n, n]^2 for n in 1:dim(S_vals, 1)]
    return sort(λs, rev=true)
end

function bulk_spectrum(psi::MPS, a::Int, b::Int)
    orthogonalize!(psi, a)
    ket = psi[a]
    bra = prime(dag(psi[a]), "Link")
    E = ket * bra
    for j in (a + 1):b
        ket = psi[j]
        bra = prime(dag(psi[j]), "Link")
        E = E * ket
        E = E * bra
    end
    row_inds = filter(i -> plev(i) == 0, inds(E))
    col_inds = filter(i -> plev(i) == 1, inds(E))
    D, _ = eigen(E, row_inds, col_inds)
    λs = [real(D[n, n]) for n in 1:dim(D, 1)]
    return sort(λs, rev=true)
end

# ---------------------------------------------------------------------
# Build a small, cheap test state. A random MPS at fixed particle number
# is enough -- this test is about the linear-algebra/index bookkeeping,
# not about the physics, so it doesn't need to be a real ground state.
# ---------------------------------------------------------------------

N_test = 20
sites = siteinds("Fermion", N_test; conserve_qns=true)
init_state = product_state_Nf(N_test, N_test ÷ 2)
psi_test = randomMPS(sites, init_state; linkdims=8)  # nontrivial entanglement

# ---------------------------------------------------------------------
# The actual check: two different code paths computing the SAME
# bipartition ([1,cut] vs [cut+1,N]) should agree to numerical precision.
# ---------------------------------------------------------------------

cut = 15

S_single = single_cut_entropy(psi_test, cut)
S_bulk   = bulk_entanglement_entropy(psi_test, 1, cut)

println("S (single-cut SVD)   = ", S_single)
println("S (bulk two-sided)   = ", S_bulk)
println("|difference|         = ", abs(S_single - S_bulk))

λ_single = single_cut_spectrum(psi_test, cut)
λ_bulk   = bulk_spectrum(psi_test, 1, cut)

n = min(length(λ_single), length(λ_bulk))
println("\nlength(spectrum): single=$(length(λ_single))  bulk=$(length(λ_bulk))")
println("max |Δλ| over first $n matched eigenvalues: ",
        maximum(abs.(λ_single[1:n] .- λ_bulk[1:n])))