push!(LOAD_PATH, joinpath(@__DIR__, "..","..","..","..", "src"))
using SSHAnalysis

# Saving path
data_dir = joinpath(ProjectRoot(), "data", "quench", "quench_paper")

# log file
log_dir = joinpath(ProjectRoot(), "logs", "quench", "quench_paper")
mkpath(log_dir)  # Create the log directory if it doesn't exist
log_file = joinpath(log_dir, "quench_from_Luttinger_liquid_to_Luttinger_liquid.log")

# Define parameters

N_sites = 200 # Number of cells is 100, so number of sites is 200
T_max = 10.0
dt = 0.05
bond_dimension = 3000

pre_quench = SSHParams(v = 1.0, w= 1.0, Δ = 0, V = 0.5)
post_quench = SSHParams(v = 1.0, w= 1.0, Δ = 0, V = 1.0)

# Run the quench simulation
times, entropies, variances = simulate_quench(N_sites, T_max, dt, pre_quench, post_quench; 
                                                bond_dim=bond_dimension, logfile=log_file)

# Save the data to a CSV file
df = DataFrame(Time=times, Entropy=entropies, Variance=variances)
filename = "run_LL_quench_from_Luttinger_liquid_to_Luttinger_liquid_200sites_10sec_BondDim3000.csv"
CSV.write(data_dir * "/$filename", df)

println("Quench simulation completed and data saved to $data_dir/$filename")

