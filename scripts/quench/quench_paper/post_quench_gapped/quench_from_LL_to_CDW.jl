push!(LOAD_PATH, joinpath(@__DIR__, "..","..","..","..", "src"))
using SSHAnalysis

# Saving path
data_dir = joinpath(ProjectRoot(), "data", "quench", "quench_paper")

# log file
log_dir = joinpath(ProjectRoot(), "logs", "quench", "quench_paper")
mkpath(log_dir)  # Create the log directory if it doesn't exist
log_file = joinpath(log_dir, "quench_from_LL_to_CDW.log")

# Define parameters

N_sites = 200 # Number of cells is 100, so number of sites is 200
T_max = 10.0
dt = 0.05
bond_dimension = 1700
cutoff = 1e-8

pre_quench = SSHParams(v = 1.0, w= 1.0, Δ = 0, V = 1.0)
post_quench = SSHParams(v = 1.0, w= 1.0, Δ = 0, V = 4.0)

# Run the quench simulation
times, entropies, variances = simulate_quench(N_sites, T_max, dt, pre_quench, post_quench; 
                                                bond_dim=bond_dimension, cutoff=cutoff, logfile=log_file)

# Save the data to a CSV file
df = DataFrame(Time=times, Entropy=entropies, Variance=variances)
filename = "run_LL_quench_from_LL_to_CDW_200sites_10sec_bonddim$(bond_dimension)_cutoff$(cutoff).csv"
full_save_path = joinpath(data_dir, filename)
CSV.write(full_save_path, df)

println("Quench simulation completed and data saved to $full_save_path")

