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

# Order matters: enable_threaded_blocksparse checks BLAS/Strided thread
# counts AT CALL TIME and warns if they're not 1 yet. Set them first.
BLAS.set_num_threads(1)
Strided.disable_threads()   # equivalent to Strided.set_num_threads(1)
ITensors.enable_threaded_blocksparse(true)

n_julia_threads = Threads.nthreads()
if n_julia_threads == 1
    @warn "Julia started with only 1 thread. Relaunch with `julia -t 8` " *
          "(or similar) to actually use enable_threaded_blocksparse."
end
println("Julia threads: $n_julia_threads, BLAS threads: $(BLAS.get_num_threads()), " *
        "Strided threads: $(Strided.get_num_threads())")

# Parameters
N_sites = 200
T_max = 10.0
dt = 0.05

# Paths
data_dir = joinpath(ProjectRoot(), "data", "quench", "quench_paper")
log_dir = joinpath(ProjectRoot(), "logs", "quench", "quench_paper")
ckpt_dir = joinpath(ProjectRoot(), "checkpoints", "quench", "quench_paper")
mkpath(log_dir)
mkpath(ckpt_dir)
log_file = joinpath(log_dir, "quench_from_LL_to_LL_decrease_V_DV02_$(T_max)_centered_subsystem_l70.log")
checkpoint_file = joinpath(ckpt_dir, "quench_from_LL_to_LL_decrease_V_DV02_$(T_max)_centered_subsystem_l70.h5")

# Bond dim: bumped up from 1700. Ramp this rather than jumping straight to
# something huge — go up, watch memory in htop/free, confirm it's stable,
# then push further if you still have headroom and still need it.
# With conserve_qns=true, actual memory use is well below the naive dense
# estimate (N * chi^2 * d), since only the physically populated QN blocks
# are stored — so don't be scared off by a large dense estimate alone,
# just verify empirically at your actual filling/parameters.
bond_dimension = 2000
cutoff = 1e-8

pre_quench = SSHParams(v=1.0, w=1.0, Δ=0, V=0.7)
post_quench = SSHParams(v=1.0, w=1.0, Δ=0, V=0.5)

# --- Fresh run with checkpointing enabled ---
# checkpoint_every=20 steps * dt=0.05 => a checkpoint roughly every 1.0
# time unit. Adjust to trade "protection against losing progress" against
# checkpoint I/O overhead (writing chi~2000+ MPS tensors isn't free).
times, variances = simulate_quench_only_fluctuations(
    N_sites, T_max, dt, pre_quench, post_quench;
    bond_dim=bond_dimension, cutoff=cutoff, logfile=log_file,
    checkpoint_file=checkpoint_file, checkpoint_every=20, subsystem_a = 65, subsystem_b = 134
)

# --- Example: resuming a run that stopped or that you want to extend at
#     higher bond_dim, instead of the call above ---
#
# times, entropies, variances = simulate_quench(
#     N_sites, T_max, dt, pre_quench, post_quench;
#     bond_dim=3000, cutoff=cutoff, logfile=log_file,
#     checkpoint_file=checkpoint_file, checkpoint_every=20,
#     resume_from=checkpoint_file,
# )

df = DataFrame(Time=times, Variance=variances)
filename = "run_LL_quench_from_LL_to_LL_decrease_V_DV02_200sites_time$(T_max)_bonddim$(bond_dimension)_cutoff$(cutoff)_centered_subsystem_l70.csv"
full_save_path = joinpath(data_dir, filename)
CSV.write(full_save_path, df)
println("Quench simulation completed and data saved to $full_save_path")
