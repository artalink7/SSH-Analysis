push!(LOAD_PATH, joinpath(@__DIR__, "..","..","..","..", "src"))
using SSHAnalysis

# Saving path
data_dir = joinpath(ProjectRoot(), "data", "quench", "quench_paper")


# Define parameters

N_sites = 200 # Number of cells is 100, so number of sites is 200
T_max = 10.0
dt = 0.05
v_i = 1.0
w_i = 1.0
V_i = 0.0
v_f = 1.0
w_f = 1.0
disorder_strength = 0.5
num_runs = 10
bond_dimension = 1500
println("Starting $num_runs runs of quench simulation with disorder strength W = $disorder_strength...")

for run in 1:num_runs
    println("\n=== Starting Run $run of $num_runs ===")


times, entropies, variances = simulate_quench_initial_disorder(
    N_sites, T_max, dt; 
    v_i=v_i, w_i=w_i, V_i=V_i, v_f=v_f, w_f=w_f, W=disorder_strength, bond_dim=bond_dimension
)

df = DataFrame(Time=times, Entropy=entropies, Variance=variances)
filename = "run_TB_quench_disordered_to_non_interacting_gapless_200sites_10sec_BondDim$(bond_dimension)_run$(run).csv"
CSV.write(joinpath(data_dir, filename), df)
println("Run $run completed and data saved to $data_dir/$filename")
end

println("\nAll $num_runs runs completed.")

