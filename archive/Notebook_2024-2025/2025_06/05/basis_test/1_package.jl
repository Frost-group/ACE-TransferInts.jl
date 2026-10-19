using Pkg
cd("C:/Users/Keerati/Documents/GitHub/ACE-TransferInts.jl/archive/Notebook_2024-2025/2025_06/05/basis_test/")
#cd("/Users/kk419/Documents/GitHub/ACE-TransferInts.jl/Notebook/2025_04/11/Datagen")
Pkg.activate(".")
Base.active_project()
Pkg.instantiate()
Pkg.add(["JuLIP", "StaticArrays", "LinearAlgebra", "Plots","LaTeXStrings", "Rotations", "Random", "Printf"])

using JuLIP
using StaticArrays
using LinearAlgebra
using Plots
using LaTeXStrings
using Rotations  # for Euler rotation matrices
using Random
using Printf