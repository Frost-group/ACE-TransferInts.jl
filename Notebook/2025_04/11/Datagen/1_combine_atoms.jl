using Pkg
cd("C:/Users/Keerati/Documents/GitHub/ACE-TransferInts.jl/Notebook/2025_04/11/Datagen")
#cd("/Users/kk419/Documents/GitHub/ACE-TransferInts.jl/Notebook/2025_04/11/Datagen")
Pkg.activate(".")
Base.active_project()
Pkg.instantiate()
Pkg.add(["JuLIP", "StaticArrays", "LinearAlgebra", "Plots","LaTeXStrings"])

using JuLIP
using StaticArrays
using LinearAlgebra
using Plots
using LaTeXStrings

file1 = "all_CSi_modified625.extxyz"
file2 = "all_CSi_modified650.extxyz"
file3 = "all_CSi_modified675.extxyz"
file4 = "all_CSi_modified700.extxyz"
# List of .extxyz files to combine
files = [file1, file2, file3, file4]
all_atoms = Vector{Atoms{Float64}}([])

for file in files
    # Read all frames in this file
    atoms_list = read_extxyz(file)
    append!(all_atoms, atoms_list)
end

write_extxyz("all_CSi_modified_all.extxyz", all_atoms)