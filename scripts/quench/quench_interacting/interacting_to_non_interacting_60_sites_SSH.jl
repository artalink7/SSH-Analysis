push!(LOAD_PATH, joinpath(@__DIR__, "..","..","..", "src"))
using SSHAnalysis

# Saving path
data_dir = joinpath(ProjectRoot(), "data", "quench")


# Define parameters

N_sites = 80
T_max = 10.0
dt = 0.05
v_i = 1.0
w_i = 0.5
V_i = 1.0
v_f = 1.0
w_f = 0.5

# Run the quench simulation
times, entropies, variances = simulate_quench(N_sites, T_max, dt; v_i=v_i, w_i=w_i, V_i=V_i, v_f=v_f, w_f=w_f)

# Save the data to a CSV file
df = DataFrame(Time=times, Entropy=entropies, Variance=variances)
filename = "run_quench_interacting_to_non_interacting_SSH_80sites_10sec.csv"
CSV.write(data_dir * "/$filename", df)

println("Quench simulation completed and data saved to $data_dir/$filename")

