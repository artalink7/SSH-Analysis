push!(LOAD_PATH, joinpath(@__DIR__, "..","..","..","..", "src"))
using SSHAnalysis
using ITensors
using LinearAlgebra
using DataFrames
using CSV
using ITensors.NDTensors: Strided
using Dates  # For timestamping the log

# ---------------------------------------------------------------------
# 1. Fluctuations Measurement Function with Integrated Logging
# ---------------------------------------------------------------------

# ---------------------------------------------------------------------
# 2. Main Script Execution
# ---------------------------------------------------------------------
n_cores = Threads.nthreads()
if n_cores == 1
    @warn "Julia started with only 1 thread. Relaunch with `julia -t 12`!"
end

ITensors.disable_threaded_blocksparse()
BLAS.set_num_threads(n_cores)
Strided.set_num_threads(n_cores)

println("Julia threads: $n_cores, BLAS threads: $(BLAS.get_num_threads()), " *
        "Strided threads: $(Strided.get_num_threads())")

# --- Parameters ---
N_sites = 400
l_min = 4
l_max = 200

# --- Paths ---
# Separated from quench data for better organization
data_dir = joinpath(ProjectRoot(), "data", "fluctuations", "quench_paper")
log_dir = joinpath(ProjectRoot(), "logs", "fluctuations", "quench_paper")
mkpath(data_dir)
mkpath(log_dir)

log_file = joinpath(log_dir, "bipartite_fluctuations_N$(N_sites)_l$(l_min)-$(l_max)_V07.log")

# --- Physics Parameters ---
pre_quench = ExtendedTBParams(t0=1.0, t2=0.3, V=0.7)

# --- Run Calculation ---
l_pre, F_pre = extract_bipartite_fluctuations_edge(
    N_sites, 
    pre_quench; 
    l_min=l_min, 
    l_max=l_max, 
    logfile=log_file
)

# --- Save Data ---
df_pre = DataFrame(length_sub = l_pre, Variance=F_pre)
filename = "bipartite_fluctuations_post_quench_N$(N_sites)_l$(l_min)-$(l_max)_V07.csv"
full_save_path = joinpath(data_dir, filename)

CSV.write(full_save_path, df_pre)
println("\n>>> Data successfully saved to $full_save_path")
println(">>> Execution log saved to $log_file")