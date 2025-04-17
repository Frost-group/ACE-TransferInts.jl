using Pkg
cd("C:/Users/Keerati/Documents/GitHub/ACE-TransferInts.jl/Notebook/2025_04/12/Datagen")
Pkg.activate(".")
Base.active_project()
Pkg.instantiate()
Pkg.add(["JuLIP", "StaticArrays", "LinearAlgebra", "Plots","LaTeXStrings"])

using JuLIP
using StaticArrays
using LinearAlgebra
using Plots
using LaTeXStrings

using JuLIP, StaticArrays

# --- Function to map molecules to centroid beads ---
function coarse_grain_centroid(frames)
    new_frames = Vector{Atoms{Float64}}([])
    for atoms in frames
        # Group 1: atoms 1–9
        pos_group1 = atoms.X[1:6]
        centroid1 = sum(pos_group1) / length(pos_group1)
        
        # Group 2: atoms 10–18
        pos_group2 = atoms.X[7:12]
        centroid2 = sum(pos_group2) / length(pos_group2)
        
        # Coarse-grained species (e.g., use Carbon (C) for both beads)
        species = [6, 14]  # atomic number 6 = Carbon (or choose whatever you want)
        
        # Positions of the two beads
        positions = [centroid1, centroid2]
        
        # Create new coarse-grained Atoms object
        cg_atoms = JuLIP.Atoms(species, positions)
        
        # Copy over any metadata like :theta, :phi
        set_cell!(cg_atoms, 20.0 * I(3)) 
        cg_atoms["energy"] = atoms["energy_eff"]*198.2
        push!(new_frames, cg_atoms)
    end
    return new_frames
end

# --- Load your file ---
frames = read_extxyz("output.extxyz")

# --- Apply coarse-graining ---
cg_frames = coarse_grain_centroid(frames)

# --- Save the new coarse-grained frames ---
write_extxyz("output_coarse_grained.extxyz", cg_frames)
