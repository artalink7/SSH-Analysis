push!(LOAD_PATH, joinpath(@__DIR__, "..","..","..","..", "src"))
using SSHAnalysis
using ITensors
using LinearAlgebra

# ---------------------------------------------------------------------
# Threading setup — do this BEFORE building any sites/tensors.
#
# You have 24 logical cores essentially idle (load avg ~1.2). This won't
# change *when* the entanglement wall hits, but it changes how much wall
# clock each χ costs you, which is what lets you push the wall further
# out in practice.
#
#   * Launch this script with:   julia -t 8 run_quench_paper_patched.jl
#     (or `julia -t auto`, or set JULIA_NUM_THREADS=8 in the environment)
#
#   * enable_threaded_blocksparse is specifically for QN-conserving
#     tensors (which you already have via conserve_qns=true) — it lets
#     ITensor parallelize contractions across the QN blocks.
#
#   * IMPORTANT: block-sparse threading parallelizes ACROSS your QN blocks
#     using Julia threads. It is coarse-grained parallelism, one block per
#     thread. If BLAS/Strided.jl then ALSO try to multithread the (small)
#     linear algebra work inside each individual block, you get nested
#     parallelism fighting over the same cores (8 Julia threads x 3 BLAS
#     threads = 24 threads contending, with only 24 physical cores and no
#     room left for anything else, plus constant thread-handoff overhead).
#     ITensors' own advice — confirmed by the warnings you saw — is to set
#     BOTH BLAS and Strided.jl to single-threaded and let the Julia threads
#     do all the parallel work. Do NOT split 24 cores proportionally here;
#     that was wrong in an earlier version of this script.
# ---------------------------------------------------------------------
using ITensors.NDTensors: Strided

n_cores = Threads.nthreads()
if n_cores == 1
    @warn "Julia started with only 1 thread. Relaunch with `julia -t 12`!"
end

# 1. Turn OFF block-sparse threading globally to banish the [13] crash
ITensors.disable_threaded_blocksparse()

# 2. Turn ON BLAS threading so dense matrix math uses all 12 cores!
BLAS.set_num_threads(n_cores)

# 3. Turn ON Strided threading to help with tensor permutations
Strided.set_num_threads(n_cores)

println("Julia threads: $n_cores, BLAS threads: $(BLAS.get_num_threads()), " *
        "Strided threads: $(Strided.get_num_threads())")


data_dir = joinpath(ProjectRoot(), "data", "quench", "quench_paper")


using LinearAlgebra

"""
Single-particle hopping matrix (real, symmetric) for NN + NNN hopping only.
Valid whenever V=0 (interaction term absent).
"""
function build_h_matrix(N_sites::Int, t0::Float64, t2::Float64)
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
    return h
end

"""
Ground-state correlation matrix G0_ij = <c_i^dagger c_j> for the N_fill
lowest-energy orbitals of h (real symmetric single-particle Hamiltonian).
"""
function ground_state_correlation_matrix(h::Matrix{Float64}, N_fill::Int)
    evals, evecs = eigen(Symmetric(h))   # ascending eigenvalues -> lowest N_fill = Fermi sea
    occ = evecs[:, 1:N_fill]
    return occ * occ'                    # real symmetric
end

"""
Exact free-fermion time evolution of the correlation matrix under h_post,
starting from G0. Returns G(t) as a complex Hermitian matrix.
"""
function evolve_correlation_matrix(G0::Matrix{Float64}, h_post::Matrix{Float64}, t::Float64)
    evals, V = eigen(Symmetric(h_post))
    M0 = V' * G0 * V
    phase = [exp(im * (evals[m] - evals[n]) * t) for m in eachindex(evals), n in eachindex(evals)]
    Mt = M0 .* phase
    return V * Mt * V'
end

"""
Bipartite charge fluctuation F for subsystem sites a:b, from a
(possibly complex Hermitian) correlation matrix G.
"""
function F_from_correlation_matrix(G::AbstractMatrix, a::Int, b::Int)
    Gsub = G[a:b, a:b]
    return real(tr(Gsub) - sum(abs2, Gsub))
end

"""
Exact free-fermion quench: F(t) for subsystem [a,b], over a time grid,
for a V=0 quench where only (t0,t2) change between pre- and post-quench H.
"""
function exact_free_fermion_quench_F(N_sites::Int, N_fill::Int,
                                      t0_pre::Float64, t2_pre::Float64,
                                      t0_post::Float64, t2_post::Float64,
                                      a::Int, b::Int, times::AbstractVector{Float64};
                                      logfile::Union{String,Nothing}=nothing)
    log_io = logfile === nothing ? nothing : open(logfile, "w")
    function log_msg(msg)
        println(msg)
        if log_io !== nothing; println(log_io, msg); flush(log_io); end
    end

    log_msg("="^60)
    log_msg("Exact free-fermion quench validation")
    log_msg("N_sites=$N_sites, N_fill=$N_fill, subsystem=[$a,$b]")
    log_msg("Pre:  t0=$t0_pre, t2=$t2_pre")
    log_msg("Post: t0=$t0_post, t2=$t2_post")
    log_msg("="^60)

    h_pre  = build_h_matrix(N_sites, t0_pre, t2_pre)
    h_post = build_h_matrix(N_sites, t0_post, t2_post)
    G0 = ground_state_correlation_matrix(h_pre, N_fill)

    # Sanity check: t=0 must reproduce the equilibrium F exactly (self-consistency of the code)
    F0_direct = F_from_correlation_matrix(complex.(G0), a, b)
    Gt0 = evolve_correlation_matrix(G0, h_post, 0.0)
    F0_evolved = F_from_correlation_matrix(Gt0, a, b)
    log_msg("Sanity check -- F(0) direct: $F0_direct, F(0) via evolve at t=0: $F0_evolved, diff = $(abs(F0_direct-F0_evolved))")

    F_vals = Float64[]
    for t in times
        Gt = evolve_correlation_matrix(G0, h_post, t)
        push!(F_vals, F_from_correlation_matrix(Gt, a, b))
    end

    if log_io !== nothing; close(log_io); end
    return collect(times), F_vals
end

using CSV, DataFrames

N_sites = 400
N_fill  = N_sites ÷ 2
a, b    = 161, 240          # same centered subsystem as your TEBD run
T_max   = 7.0
dt_grid = 0.05
times   = 0.0:dt_grid:T_max

# V=0 throughout; quench t2 instead, e.g. 0.3 -> 0.1 (pick whatever you like,
# just make sure both sides stay away from t2=t0/2 where the dispersion
# develops extra Fermi points / band structure changes qualitatively)
t_vals, F_exact = exact_free_fermion_quench_F(
    N_sites, N_fill,
    1.0, 0.3,    # pre:  t0=1.0, t2=0.3
    1.0, 0.1,    # post: t0=1.0, t2=0.1
    a, b, times;
    logfile = "exact_ff_quench_check.log"
)

df = DataFrame(Time = t_vals, F_exact = F_exact)
filename = "quench_LL_to_LL_NNN_ED_no_V.csv"
full_save_path = joinpath(data_dir, filename)
CSV.write(full_save_path, df)
println(">>> Saved exact free-fermion F(t) for comparison against TEBD. to $(full_save_path)")