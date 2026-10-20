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

# Parameters
N_sites = 1200
T_max = 7.0
dt = 0.05

# Paths
data_dir = joinpath(ProjectRoot(), "data", "quench", "quench_paper")
log_dir = joinpath(ProjectRoot(), "logs", "quench", "quench_paper")
ckpt_dir = joinpath(ProjectRoot(), "checkpoints", "quench", "quench_paper")
mkpath(log_dir)
mkpath(ckpt_dir)
log_file = joinpath(log_dir, "quench_from_LL_to_LL_NNN_coupling_bd1600_$(N_sites)_decrease_V_DV02_$(T_max)_edge_subsystem_l80.log")
checkpoint_file = joinpath(ckpt_dir, "quench_from_LL_to_LL_NNN_coupling_bd1600_$(N_sites)_decrease_V_DV02_$(T_max)_edge_subsystem_l80.h5")

# Bond dim: bumped up from 1700. Ramp this rather than jumping straight to
# something huge — go up, watch memory in htop/free, confirm it's stable,
# then push further if you still have headroom and still need it.
# With conserve_qns=true, actual memory use is well below the naive dense
# estimate (N * chi^2 * d), since only the physically populated QN blocks
# are stored — so don't be scared off by a large dense estimate alone,
# just verify empirically at your actual filling/parameters.
bond_dimension = 1600
cutoff = 1e-10

pre_quench =  ExtendedTBParams(t0=1.0, t2=0.3, V=0.7)
post_quench = ExtendedTBParams(t0=1.0, t2=0.3, V=0.5)

# --- Fresh run with checkpointing enabled ---
# checkpoint_every=20 steps * dt=0.05 => a checkpoint roughly every 1.0
# time unit. Adjust to trade "protection against losing progress" against
# checkpoint I/O overhead (writing chi~2000+ MPS tensors isn't free).
times, variances = simulate_quench_only_fluctuations_extended_tebd_new(
    N_sites, T_max, dt, pre_quench, post_quench;
    bond_dim=bond_dimension, cutoff=cutoff, logfile=log_file,
    checkpoint_file=checkpoint_file, checkpoint_every=20, subsystem_a = 561, subsystem_b = 640, 
    resume_from=checkpoint_file
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
filename = "quench_from_LL_to_LL_NNN_coupling_decrease_V_DV02_sites$(N_sites)_time$(T_max)_bonddim$(bond_dimension)_cutoff$(cutoff)_central_subsystem_l80.csv"
full_save_path = joinpath(data_dir, filename)
CSV.write(full_save_path, df)
println("Quench simulation completed and data saved to $full_save_path")
