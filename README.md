# ACE-TransferInts.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://Frost-group.github.io/ACE-TransferInts.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://Frost-group.github.io/ACE-TransferInts.jl/dev/)
[![Build Status](https://github.com/Frost-group/ACE-TransferInts.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Frost-group/ACE-TransferInts.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/Frost-group/ACE-TransferInts.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/Frost-group/ACE-TransferInts.jl)

A Julia research package for predicting **intermolecular transfer integrals** in organic molecular dimers using **Atomic Cluster Expansion (ACE)** machine-learning potentials.

Transfer integrals (electronic coupling matrix elements) determine charge- and exciton-transport rates in organic semiconductors. This project replaces expensive quantum-chemistry calculations with fast ACE surrogate models trained on dimer geometry libraries.

---

## Repository layout

```
ACE-TransferInts.jl/
│
├── src/                          # Julia library — importable functions
│   ├── ACE-TransferInts.jl       # Module entry point (ACETransferInts)
│   ├── dimer.jl                  # rot_x/y/z, centroid, construct_dimer
│   ├── fitting.jl                # extract_energies, rmse, mae, sign_accuracy
│   └── io.jl                     # read_J, write_xyz, read_lc_csv
│
├── experiments/                  # One folder per molecule
│   ├── ethylene/
│   │   ├── datagen/              # Scripts that generate training structures
│   │   ├── fitting/              # ACE fit scripts and hyperparameter sweeps
│   │   └── results/              # CSV learning curves, plots
│   ├── thiophene/
│   │   ├── datagen/
│   │   ├── fitting/
│   │   └── results/
│   └── naphthalene/
│       ├── datagen/
│       ├── fitting/
│       ├── results/
│       ├── heavy_atom/           # All-heavy-atom representation
│       │   ├── datagen/
│       │   ├── fitting/
│       │   └── results/
│       └── bead/                 # Coarse-grained bead representation
│           ├── datagen/
│           ├── fitting/
│           └── results/
│
├── notebooks/                    # Exploratory and publication notebooks
│   ├── learning_curves.jl        # Interactive learning-curve exploration
│   ├── learning_curves_paper.jl  # Final paper figures (multi-molecule)
│   └── ACE_Package.jl            # Environment setup helper
│
├── data/
│   ├── reference/                # Si_dataset.xyz and other reference data
│   ├── raw/                      # Original extxyz from quantum chemistry
│   │   ├── ethylene/
│   │   ├── thiophene/
│   │   └── naphthalene/
│   └── processed/                # Train/test splits, modified structures
│
├── archive/
│   └── Notebook_2024-2025/       # Chronological research notebooks (read-only)
│
├── test/
│   └── runtests.jl
├── docs/
│   ├── make.jl
│   └── src/index.md
│
├── Project.toml
├── Manifest.toml
└── LICENSE
```

---

## Background

**Atomic Cluster Expansion (ACE)** builds a systematically improvable, many-body descriptor for atomic environments from rotationally invariant polynomials of interatomic distances and angles. Fitting an ACE model to reference transfer-integral data gives a fast surrogate that:

- is **invariant** to translation, rotation, and permutation of equivalent atoms;
- scales with **body order** (ν) and **polynomial degree** (td) so accuracy is controllable;
- requires only atomic positions at inference time.

---

## Molecules studied

| Molecule | Atom types in model | Notes |
|---|---|---|
| Ethylene (C₂H₄) | C + pseudo-Si | Smallest π-conjugated test case |
| Thiophene (C₄H₄S) | C, Si or C, Si, S, O | Sulphur-containing heterocycle |
| Naphthalene (C₁₀H₈) | Heavy atoms or beads | Fused-ring aromatic; two representations |

H atoms are replaced by heavier pseudo-Si atoms so the ACE radial basis covers the relevant interatomic distances cleanly. A coarse-grained **bead** model of naphthalene is also studied as a more compact representation.

---

## Library (`src/`)

```julia
using ACETransferInts

# Dimer geometry
R = rot_y(π/3)                          # 3×3 rotation matrix
c = centroid(atoms.X)                   # centre of geometry
construct_dimer("mol.extxyz", "dimer.extxyz", params, π/3, J)

# Metrics
rmse(predicted, reference)
mae(predicted, reference)
sign_accuracy(predicted, reference)     # fraction with correct sign

# I/O
params, J, logJ = read_J("out.dat")
write_xyz("mol.xyz", atoms)
n, rmse_v, mae_v, sign_acc = read_lc_csv("LC_Ethylene.csv")
```

---

## ACE basis parameters

| Parameter | Meaning | Typical range |
|---|---|---|
| `ν` (order) | Body order (2 = pair, 3 = 3-body, …) | 2 – 5 |
| `td` (total degree) | Maximum polynomial degree | 5 – 40 |
| `rcut` | Cutoff radius (Å) | 7 – 15 |
| `p`, `q` | Agnesi transform shape | tuned per molecule |

Fitting scripts and hyperparameter sweep results live in `experiments/<molecule>/fitting/`.

---

## Learning curves

Learning-curve data are stored as row-keyed CSV files:

```
n_train,  10, 50, 100, 200, 500, 1000
rmse,     ...
mae,      ...
sign_acc, ...
```

Final paper plots (RMSE, MAE, and sign accuracy on log–log axes for all molecules) are in `notebooks/learning_curves_paper.jl` and the per-molecule PNGs in `experiments/*/results/`.

---

## Installation

```julia
import Pkg
Pkg.develop(path="path/to/ACE-TransferInts.jl")
```

The main external dependency is [ACEpotentials.jl](https://github.com/ACEsuit/ACEpotentials.jl) (requires the ACEsuit registry). All other dependencies are in `Project.toml`.

---

## Archive

`archive/Notebook_2024-2025/` preserves the original chronological research notebooks (July 2024 – April 2026) exactly as written. These are the historical record; new work should go in `experiments/` or `notebooks/`.

---

## License

MIT — see [LICENSE](LICENSE).  
Copyright © Jarvist Moore Frost and contributors.
