using Pkg
cd("C:/Users/Keerati/Documents/GitHub/ACE-TransferInts.jl/Notebook/2025_04/11/Datagen")
Pkg.activate(".")
Base.active_project()
Pkg.instantiate()
Pkg.add(["JuLIP", "StaticArrays", "LinearAlgebra", "Plots","LaTeXStrings"])

using JuLIP
using StaticArrays
using LinearAlgebra
using Plots
using LaTeXStrings

function change_Si_to_C(frames)
    new_frames = Vector{Atoms{Float64}}([])
    for atoms in frames
        # Find indices of all Si atoms (Z == 14)
        C_indices = findall(z -> z == 6, atoms.Z)
        
        # Check if there are at least 4 Si atoms
        if length(C_indices) >= 4
            # Change the 3rd and 4th Si atoms to Carbon (Z = 6)
            atoms.Z[C_indices[3]] = 14
            atoms.Z[C_indices[4]] = 14
        else
            @warn "Frame has fewer than 4 Si atoms. Skipping change."
        end

        # Keep properties like :theta, :phi, etc.
        push!(new_frames, atoms)
    end
    return new_frames
end

# --- Load the file ---
frames = read_extxyz("all_C_dimers700.extxyz")

# --- Perform the Si -> C substitution ---
frames_modified = change_Si_to_C(frames)



# --- Save the modified frames --Cwrite_extxyz("all_bead_modified.extxyz", frames_modified)
write_extxyz("all_C_modified700.extxyz", frames_modified)