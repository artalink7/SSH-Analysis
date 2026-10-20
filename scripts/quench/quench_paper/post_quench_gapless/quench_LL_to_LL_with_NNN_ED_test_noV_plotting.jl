push!(LOAD_PATH, joinpath(@__DIR__, "..","..","..","..", "src"))
using SSHAnalysis
using ITensors
using LinearAlgebra

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df_exact = pd.read_csv("data/quench/quench_paper/quench_LL_to_LL_NNN_ED_no_V.csv")
df_tebd  = pd.read_csv("data/quench/quench_paper/quench_from_LL_to_LL_NNN_coupling_decrease_V_DV02_sites400_time7.0_bonddim1600_cutoff1.0e-10_central_subsystem_l80_testing_noV.csv")   # Time, Variance columns

# Align on time grid (should match exactly if dt is the same)
merged = pd.merge(df_exact, df_tebd, on="Time", suffixes=("_exact", "_tebd"))
merged["abs_diff"] = np.abs(merged["F_exact"] - merged["Variance"])

print(merged[["Time", "F_exact", "Variance", "abs_diff"]].to_string(index=False))
print(f"\nMax abs diff: {merged['abs_diff'].max():.2e}")
print(f"Mean abs diff: {merged['abs_diff'].mean():.2e}")

plt.figure(figsize=(8,5))
plt.semilogy(merged["Time"], merged["abs_diff"], marker='o', markersize=3)
plt.xlabel("Time t")
plt.ylabel("|F_exact - F_TEBD|")
plt.title("TEBD vs. exact free-fermion quench: absolute error")
plt.grid(True, linestyle='-', alpha=0.6)
plt.tight_layout()
plt.show()