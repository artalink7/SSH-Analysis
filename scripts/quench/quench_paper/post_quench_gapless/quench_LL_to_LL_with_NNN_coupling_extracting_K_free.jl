push!(LOAD_PATH, joinpath(@__DIR__, "..","..","..","..", "src"))
using SSHAnalysis
using ITensors
using LinearAlgebra
using DataFrames
using CSV
using ITensors.NDTensors: Strided
using Dates

n_cores = Threads.nthreads()
if n_cores == 1
    @warn "Julia started with only 1 thread. Relaunch with `julia -t 12`!"
end

ITensors.disable_threaded_blocksparse()
BLAS.set_num_threads(n_cores)
Strided.set_num_threads(n_cores)

println("Julia threads: $n_cores, BLAS threads: $(BLAS.get_num_threads()), " *
        "Strided threads: $(Strided.get_num_threads())")

using LinearAlgebra

"""
Exact free-fermion bipartite fluctuations for a CENTERED subsystem [a,b],
symmetric about N_sites÷2, matching the same centering convention used
in extract_bipartite_fluctuations_centered_new (a = div(N-l,2)+1, b = a+l-1).
"""
function exact_free_fermion_F_centered(N_sites::Int, t0::Float64, t2::Float64,
                                        l_min::Int, l_max::Int)
    # Build single-particle Hamiltonian matrix (OBC)
    h = zeros(Float64, N_sites, N_sites)
    for i in 1:N_sites-1
        h[i, i+1] = -t0
        h[i+1, i] = -t0
    end
    if t2 != 0.0
        for i in 1:N_sites-2
            h[i, i+2] = -t2
            h[i+2, i] = -t2
        end
    end

    evals, evecs = eigen(Symmetric(h))
    N_fill = N_sites ÷ 2
    occ_orbitals = evecs[:, 1:N_fill]
    G = occ_orbitals * occ_orbitals'    # correlation matrix G_ij = <c_i^dag c_j>

    l_vals = Int[]
    F_vals = Float64[]
    a_vals = Int[]
    b_vals = Int[]
    for l in l_min:l_max
        a = div(N_sites - l, 2) + 1
        b = a + l - 1
        Gsub = G[a:b, a:b]
        F = tr(Gsub) - sum(abs2, Gsub)
        push!(l_vals, l)
        push!(F_vals, F)
        push!(a_vals, a)
        push!(b_vals, b)
    end
    return l_vals, F_vals, a_vals, b_vals
end

# --- Parameters ---
N_sites = 1000
l_min = 4
l_max = 500
t0 = 1.0
t2 = 0.3

# --- Paths ---
data_dir = joinpath(ProjectRoot(), "data", "fluctuations", "quench_paper")
log_dir = joinpath(ProjectRoot(), "logs", "fluctuations", "quench_paper")
mkpath(data_dir)
mkpath(log_dir)

log_file = joinpath(log_dir, "bipartite_fluctuations_N$(N_sites)_l$(l_min)-$(l_max)_V00_ED_t203_CENTERED.log")

l_pre, F_pre, a_pre, b_pre = exact_free_fermion_F_centered(N_sites, t0, t2, l_min, l_max)

# quick sanity print, mirrors the per-l log style you've used elsewhere
open(log_file, "w") do io
    println(io, "="^60)
    println(io, "Exact free-fermion CENTERED subsystem check: $(now())")
    println(io, "N_sites=$N_sites, t0=$t0, t2=$t2, l=$l_min:$l_max")
    println(io, "="^60)
    for i in eachindex(l_pre)
        println(io, "  -> l=$(l_pre[i]) (sites $(a_pre[i]) to $(b_pre[i])): F=$(round(F_pre[i], digits=6))")
    end
end

# --- Save Data ---
df_pre = DataFrame(length_sub=l_pre, Variance=F_pre, a=a_pre, b=b_pre)
filename = "bipartite_fluctuations_post_quench_N$(N_sites)_l$(l_min)-$(l_max)_V00_ED_t203_CENTERED.csv"
full_save_path = joinpath(data_dir, filename)

CSV.write(full_save_path, df_pre)
println("\n>>> Data successfully saved to $full_save_path")
println(">>> Log saved to $log_file")