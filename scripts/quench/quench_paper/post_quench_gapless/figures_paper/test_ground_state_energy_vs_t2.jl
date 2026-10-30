push!(LOAD_PATH, joinpath(@__DIR__, "..","..","..","..","..", "src"))
using SSHAnalysis
using ITensors
using LinearAlgebra
using Dates

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

# ---------------------------------------------------------------------
# System Parameters
# ---------------------------------------------------------------------
N_sites = 600
N_particles = 300  # Half-filling
t0 = 1.0
V  = 0.5
bond_dimension = 1800 # Passed implicitly or set in your SSHAnalysis environment

# ---------------------------------------------------------------------
# Paths & Directories
# ---------------------------------------------------------------------
# Simplified directory structure for a static ground state sweep
data_dir = joinpath(ProjectRoot(), "data", "ground_state_sweep")
log_dir  = joinpath(ProjectRoot(), "logs", "ground_state_sweep")
gs_dir   = joinpath(ProjectRoot(), "ground_state", "sweep")

mkpath(data_dir)
mkpath(log_dir)
mkpath(gs_dir)

# ---------------------------------------------------------------------
# t2 Sweep Execution
# ---------------------------------------------------------------------
t2_vals = 0.0:0.02:0.2
E0_array = Float64[]

println("Starting t2 sweep at half-filling (N=$N_sites, V=$V)...")

for t2 in t2_vals
    println("--> Computing ground state for t2 = $t2")
    
    # Update parameters for the current step
    params = ExtendedTBParams(t0=t0, t2=t2, V=V)    
    
    # Create unique filenames for this specific t2 run
    gs_file  = joinpath(gs_dir,  "gs_N$(N_sites)_t2$(t2)_V$(V).h5")
    log_file = joinpath(log_dir, "gs_N$(N_sites)_t2$(t2)_V$(V).log")
    
    # Compute ground state (fixed 'pre_quench' typo here)
    E, psi, quality = compute_and_save_ground_state(
        N_sites, 
        params, 
        gs_file;
        Nf = N_particles,
        logfile = log_file
    )
    
    push!(E0_array, E)
end

# ---------------------------------------------------------------------
# Data Export
# ---------------------------------------------------------------------
# Store the arrays in a DataFrame (fixed the scalar 'E' typo here)
df = DataFrame(t2 = t2_vals, GroundEnergy = E0_array)

filename = "energy_vs_t2_N$(N_sites)_V$(V).csv"
full_save_path = joinpath(data_dir, filename)

CSV.write(full_save_path, df)
println("Sweep completed successfully! Data saved to: $full_save_path")