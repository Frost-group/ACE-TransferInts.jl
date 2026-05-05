# ACE-TransferInts.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://Frost-group.github.io/ACE-TransferInts.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://Frost-group.github.io/ACE-TransferInts.jl/dev/)
[![Build Status](https://github.com/Frost-group/ACE-TransferInts.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Frost-group/ACE-TransferInts.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/Frost-group/ACE-TransferInts.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/Frost-group/ACE-TransferInts.jl)

A Julia research package for predicting **intermolecular transfer integrals** in organic molecular dimers using **Atomic Cluster Expansion (ACE)** machine-learning potentials.

Transfer integrals (electronic coupling matrix elements) are key quantities for charge- and exciton-transport in organic semiconductors. Direct quantum-chemistry calculations are expensive; this project trains ACE-based models on a library of dimer geometries so that transfer integrals can be predicted from atomic positions alone.

---

## Background

The **Atomic Cluster Expansion** (ACE) is a systematically improvable, many-body descriptor for atomic environments. It expands the potential-energy surface (or, here, an electronic property) in a complete basis of body-ordered, rotationally invariant polynomials of interatomic distances and angles.

By fitting an ACE model to reference transfer-integral data (obtained, e.g., from DFT or tight-binding calculations), one obtains a fast surrogate that:

- is **invariant** to translation, rotation, and permutation of equivalent atoms;
- is **systematically improvable** by raising the polynomial degree or body order;
- requires only atomic positions as input at inference time.

---

## Molecules Studied

| Molecule | Atoms in model | Notes |
|---|---|---|
| Ethylene (C₂H₄) | C, (pseudo-)Si | Smallest π-conjugated testcase |
| Thiophene (C₄H₄S) | C, Si or C, Si, S, O | Sulphur-containing heterocycle |
| Naphthalene (C₁₀H₈) | Heavy atoms or bead model | Fused-ring aromatic |

Hydrogen atoms are often replaced by pseudo-silicon atoms (heavier ghost atoms) so that the ACE radial basis covers the relevant interatomic distances without ambiguity. A coarse-grained **bead** representation of naphthalene is also explored.

---

## Repository Layout

```
ACE-TransferInts.jl/
├── src/
│   ├── ACE-TransferInts.jl      # Package module entry point
│   ├── ProcessFiles.jl          # Convert JuLIP Atoms → XYZ files
│   └── process_dimer.jl         # Centroid / centre-of-mass utilities
│
├── Notebook/                    # Research notebooks (chronological)
│   ├── 2024_07/ … 2026_04/      # Julia scripts and Pluto notebooks
│   └── …
│
├── test/
│   └── runtests.jl              # Test suite (placeholder)
├── docs/
│   ├── make.jl
│   └── src/index.md
│
├── Project.toml                 # Julia package dependencies
├── Manifest.toml                # Locked dependency versions
└── LICENSE                      # MIT
```

The main research workflow lives in `Notebook/`, organised by date. Key subdirectories:

| Path | Content |
|---|---|
| `Notebook/2024_07–08/` | Initial dimer construction, extxyz reading, spherical harmonics visualisation |
| `Notebook/2024_09–11/` | First ACE fits, radial basis analysis, hyperparameter scans |
| `Notebook/2025_03–04/` | Large-scale fitting, model comparisons, cutoff / degree sweeps |
| `Notebook/2025_05–06/` | Transfer integral prediction and validation |
| `Notebook/2026_04/` | Learning curves (RMSE, MAE, sign accuracy) for the paper |

---

## Key Workflows

### 1. Dimer Construction

Molecular dimers are built programmatically by applying rotation matrices and translations to a monomer geometry. Rotation functions `rot_x`, `rot_y`, `rot_z` generate all required orientations. The centroid (centre of mass) of each fragment is computed with the utility in `src/process_dimer.jl`.

### 2. ACE Basis Setup

The ACE1 basis is configured with:

| Parameter | Meaning | Typical range |
|---|---|---|
| `ν` (order) | Body order (2 = pair, 3 = 3-body, …) | 2 – 5 |
| `td` (total degree) | Maximum polynomial degree | 5 – 40 |
| `rcut` | Cutoff radius (Å) | 7 – 15 |
| `p`, `q` | Agnesi transform shape | tuned per molecule |

The radial distribution function (RDF) of the training set is analysed first to inform the choice of cutoff and transform parameters.

### 3. Model Fitting

ACE models are fitted by linear regression against reference transfer-integral values stored in extended-XYZ files. Performance is measured with:

- **RMSE** and **MAE** (in μeV)
- **Sign accuracy** — fraction of dimers whose coupling sign is correctly predicted (important for coherent transport)

### 4. Learning Curves

Learning-curve scripts in `Notebook/2026_04/` train models on subsets of increasing size and record RMSE, MAE, and sign accuracy versus number of training structures. Results are saved as CSV files and plotted on a log–log scale to reveal the convergence rate.

Example CSV format:

```
n_train, 10, 50, 100, 200, 500, 1000
rmse,    ...
mae,     ...
sign_acc,...
```

---

## Installation

```julia
# From the Julia REPL
import Pkg
Pkg.develop(path="path/to/ACE-TransferInts.jl")
```

Dependencies are listed in `Project.toml`. The main external requirement is [ACEpotentials.jl](https://github.com/ACEsuit/ACEpotentials.jl) together with JuLIP and ASE (accessed via a Python bridge).

---

## Dependencies

| Package | Role |
|---|---|
| `ACEpotentials` | ACE basis construction and fitting |
| `JuLIP` | Atomic structure representation |
| `ASE` | Python ASE bridge for reading/writing structures |
| `Rotations`, `StaticArrays` | Geometry manipulation |
| `MultivariateStats` | Dimensionality reduction utilities |
| `Plots`, `LaTeXStrings` | Visualisation |
| `CSV`, `DataFrames`, `XLSX` | Data I/O |

---

## License

MIT — see [LICENSE](LICENSE).  
Copyright © Jarvist Moore Frost and contributors.
