using Pkg
cd("C:/Users/Keerati/Documents/GitHub/ACE-TransferInts.jl/Notebook/2025_04/13/Datagen")
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
        si_indices = findall(z -> z == 14, atoms.Z)
        
        # Check if there are at least 4 Si atoms
        if length(si_indices) >= 4
            # Change the 3rd and 4th Si atoms to Carbon (Z = 6)
            atoms.Z[si_indices[3]] = 6
            atoms.Z[si_indices[4]] = 6
        else
            @warn "Frame has fewer than 4 Si atoms. Skipping change."
        end

        # Keep properties like :theta, :phi, etc.
        push!(new_frames, atoms)
    end
    return new_frames
end

# --- Load the file ---
frames = read_extxyz("all_bead.extxyz")

# --- Perform the Si -> C substitution ---
frames_modified = change_Si_to_C(frames)


# --- Save the modified frames ---
write_extxyz("all_bead_modified.extxyz", frames_modified)


function change_atom(frames, A::Int64, B::Int64)
    new_frames = Vector{Atoms{Float64}}([])
    for atoms in frames
        # Find indices of all Si atoms (Z == 14)
        indices = findall(z -> z == A, atoms.Z)
        
        # Check if there are at least 4 Si atoms
        if length(indices) >= 20
            # Change the 3rd and 4th Si atoms to Carbon (Z = 6)
            for i in range(1,10)
                atoms.Z[indices[i]] = B
            end
        else
            @warn "Frame has fewer than 4 Si atoms. Skipping change."
        end

        # Keep properties like :theta, :phi, etc.
        push!(new_frames, atoms)
    end
    return new_frames
end

# --- Load the file ---
frames = read_extxyz("all_C_dimers.extxyz")

# --- Perform the Si -> C substitution ---
frames_modified = change_atom(frames, 6, 14)

# --- Save the modified frames ---
write_extxyz("all_C_modified.extxyz", frames_modified)


