using Pkg
cd("C:/Users/Keerati/Documents/GitHub/ACE-TransferInts.jl/Notebook/2025_06/03/Datagen") # Window
#cd("/Users/kk419/Documents/GitHub/ACE-TransferInts.jl/Notebook/2025_04/11/Datagen") # MacOS
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