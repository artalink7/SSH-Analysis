
# SSH-Analysis

A Julia-based framework to simulate and analyze the Su-Schrieffer–Heeger (SSH) and Rice-Mele models, with post-processing and visualization in Python. This project is designed for entanglement entropy and variance studies in equilibrium and quench dynamics.

---

## 📁 Repository Layout

```
SSH-Analysis/
│
├── data/                  # Output data (entropy, variance, etc.)
│   ├── equilibrium/
│   └── quench/
│
├── plot/                  # Python notebooks for plotting
│   └── analysis.ipynb
│
├── scripts/               # Julia scripts for running simulations
│   ├── equilibrium/
│   └── quench/
│
├── src/                   # Julia source files
│   ├── SSHAnalysis.jl     # Main module
│   ├── DensityMatrix.jl
│   ├── Entropy.jl
│   ├── Variance.jl
│   ├── BlochVectors.jl
│   └── Utils.jl
│
└── README.md              # This file
```

---

## Getting Started

### 1. Clone the repository

```bash
git clone https://github.com/artalink7/SSH-Analysis.git
cd SSH-Analysis
```

### 2. Run Julia Scripts

All scripts assume relative paths based on the project root. You can run them using:

```bash
julia scripts/equilibrium/run_entropy_equilibrium.jl
```

Make sure your data directories exist. If not, they will be created automatically.

### 3. Visualize Results

Use the Python notebook in `plot/`:

```bash
cd plot
jupyter notebook analysis.ipynb
```
---

## ⚙️ Dependencies

### Julia
- Julia ≥ 1.6
- Python
- Libraries used inside src and plot
