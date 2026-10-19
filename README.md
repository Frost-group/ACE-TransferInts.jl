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

| Molecule | Atom types in model | rcut (Å) | RMSE at convergence | Sign accuracy |
|---|---|---|---|---|
| Ethylene (C₂H₄) | Si, C | 7.5 | ~2.85 μeV/mol | ~86% |
| Thiophene (C₄H₄S) | Si, C, S, O | 9.0 | ~1.95 μeV/mol | ~86% |
| Naphthalene (C₁₀H₈) | Si, C (heavy atom) | 13.6 | ~12 meV/atom | ~97% |
| Naphthalene (bead) | Si, C (4-bead CG) | 13.6 | ~18 meV/bead | ~97% |

H atoms are replaced by heavier pseudo-Si atoms so the ACE radial basis covers the relevant interatomic distances cleanly. A coarse-grained **bead** model of naphthalene is also studied as a more compact representation. All models use ν=2, td=10 with 8,000 randomly sampled training structures; convergence is reached at ~1,500–2,000 structures.

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

| Parameter | Meaning | Small molecules | Large molecules |
|---|---|---|---|
| `ν` (order) | Body order (2 = pair, 3 = 3-body, …) | 2 | 2 |
| `td` (total degree) | Maximum polynomial degree | 10 | 10 |
| `rcut` | Cutoff radius (Å) | 7.5 – 9.0 | 13.6 |
| `r₀` | Reference radius (Å) | 6.625 | 12.0 |
| `p` | Agnesi transform (numerator) | 19 | 19 |
| `q` | Agnesi transform (denominator) | 45 | 43 |

The canonical values above were fixed through systematic sweeps in April 2025 (rcut, p, q, td each varied independently) and validated across all three molecules. Fitting scripts and per-molecule sweep results live in `experiments/<molecule>/fitting/`.

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

## Uncertainty quantification

Models are fitted with a **Bayesian Linear Regression (BLR) committee** via `ACEfit.BLR(committee_size=500)`. Each committee member is drawn from the posterior over the linear coefficients, and the predictive uncertainty for a configuration is the committee standard deviation

```
σ(x) = std{ J⁽ᵐ⁾(x) },   m = 1 … M      (M = 500)
```

evaluated with `ACE1.co_energy`. These σ values appear as error bars on the per-configuration parity plots in the results figures.

### Basis stability under high molecular overlap

The element-dependent radial envelope `f_env = (y − 1)²(y − y_cut)²` drives the radial basis smoothly to zero as `r_ij → 0`, so predictions stay bounded even when monomers strongly overlap. Driving a **cofacial dimer inward past the sampled separation band** keeps the committee-mean J finite and smooth, while σ grows by several orders of magnitude — demonstrating both basis stability and honest extrapolation flagging. The demonstration sweep (ethylene; J ± σ vs intermolecular distance) is in `archive/Notebook_2024-2025/2025_06/05/basis_test/basis_stability_sweep.jl`; a parallel heteroatom (thiophene) sweep script is under `2025_06/10/basis_test/`.

---

## QM backend

Training labels are effective transfer integrals computed with [xTB](https://github.com/grimme-lab/xtb) via the `--dipro` flag (GFN1-xTB level). Raw DIPRO output gives J and the overlap integral S; the effective coupling is:

```
J_eff = (J × 1928 − S × E_avg) / (1 − S²)
```

where E_avg is the average site energy (meV). The Python script `script.py` in each `datagen/` directory batch-processes dimer XYZ files and attaches `energy_eff` to each extxyz frame.

---

## Key findings

- **Convergence is fast.** All systems reach diminishing returns by ~1,500–2,000 training structures; 5,000 structures yield only marginal improvement.
- **Sign accuracy is the practical target.** Models achieve 80–97% correct sign prediction depending on system, which is sufficient to rank dimer orientations for transport calculations.
- **Including real heteroatoms helps.** For thiophene, using Si, C, S, O (4-element model) reduces RMSE by ~35% compared to a 2-element Si, C proxy at large training set sizes.
- **Larger molecules need larger cutoffs but are easier to learn.** Naphthalene's sign accuracy (~97%) exceeds ethylene's (~86%) because its π–π stacking geometry is more regular and predictable.
- **The bead model is fragile at small N.** The 4-bead coarse-grained naphthalene representation needs >900 training structures to stabilise; below that, RMSE diverges.
- **Predictive uncertainty is available.** A 500-member BLR committee gives per-configuration error bars; the committee spread grows sharply outside the training domain, flagging extrapolation.
- **The ACE basis is stable at high overlap.** Pushing a cofacial dimer inward past the sampled separation range keeps J finite and smooth while σ inflates — no divergence, by construction of the radial envelope.

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
