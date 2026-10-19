# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

ACE-TransferInts.jl is a Julia package for predicting intermolecular transfer integrals (electronic coupling matrix elements) in organic molecular dimers using Atomic Cluster Expansion (ACE) machine-learning potentials. Transfer integrals govern charge/exciton transport in organic semiconductors.

## Common Commands

```julia
# Start Julia with the project environment
julia --project=.

# Run tests
julia --project=. -e "using Pkg; Pkg.test()"

# Run a single test (from REPL)
julia --project=.
# then: include("test/runtests.jl")

# Instantiate dependencies (first time or after Project.toml changes)
julia --project=. -e "using Pkg; Pkg.instantiate()"

# Build documentation
julia --project=docs docs/make.jl
```

For experiment scripts, navigate to the relevant subdirectory and activate that project if it has its own `Project.toml`.

## Architecture

The package is a thin library layer (`src/`) over the ACEpotentials.jl/JuLIP.jl ecosystem. Heavy scientific work lives in `experiments/` and `archive/`.

### `src/` — Library

- **ACETransferInts.jl** — module entry point (`module ACETransferInts`), re-exports everything
- **dimer.jl** — `rot_x/y/z(θ)` rotation matrices; `centroid(atoms)`; `construct_dimer(monomer, params)` where `params = [θx, θy, θz, dx, dy, dz]` (angles + displacement)
- **fitting.jl** — `extract_energies(dataset)`, `rmse`, `mae`, `sign_accuracy` (fraction of structures where predicted/true sign agree)
- **io.jl** — `read_J(file)` for transfer-integral output files (columns: 6 params, J, logJ); `write_xyz(atoms)`; `read_lc_csv(file)` for row-keyed learning-curve CSVs

### `experiments/` — Active Research

Naphthalene has `datagen/`, `fitting/`, `bead/` and `heavy_atom/`; ethylene and thiophene currently hold only
`results/` (learning-curve CSVs) — their data generation and fits live in `archive/`.

Molecules: **ethylene** (C2H4, simplest), **thiophene** (C4H4S, heterocycle), **naphthalene** (C10H8, fused-ring; includes heavy-atom and two-bead variants). Hydrogens are dropped; the two molecules' sites are labelled as different pseudo-elements (C vs Si) so ACE can tell them apart.

### `data/`

- `reference/` — Si_dataset.xyz (the only tracked data file; `data/` is git-ignored)

### `archive/`

Read-only chronological research notebooks (`archive/Notebook_2024-2025/YYYY_MM/DD/`). Do not modify these; treat them as historical record.

### `notebooks/`

Interactive Julia scripts for exploratory analysis and paper figures.

## Key ACE Hyperparameters

ACE models are defined by:
- `ν` — body order (2=pair, 3=three-body, …; typical: 2–5)
- `td` (total degree) — polynomial truncation (typical: 5–40)
- `rcut` — cutoff radius in Å (typical: 7–15)
- `p`, `q` — Agnesi transform parameters (tuned per molecule)

## Learning Curve CSV Format

Row-keyed format used by `read_lc_csv`:
```
n_train, 10, 50, 100, 200, 500, 1000
rmse,    ...
mae,     ...
sign_acc,...
```

## Dependencies

Core: **ACEpotentials**, **JuLIP**, **ASE** (Python interop via PyCall/Conda). Data: **CSV**, **DataFrames**. Math: **Rotations**, **StaticArrays**. Viz: **Plots**.

The `.gitignore` excludes `*.extxyz`, `*.xyz`, `data/`, `*.png`, `*.dat`, and `pixi.lock/toml` — these are large or generated artifacts.
