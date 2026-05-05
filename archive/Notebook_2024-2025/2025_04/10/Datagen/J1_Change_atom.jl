using Pkg
cd("C:/Users/Keerati/Documents/GitHub/ACE-TransferInts.jl/Notebook/2025_04/10/Datagen")
Pkg.activate(".")
Base.active_project()
Pkg.instantiate()
Pkg.add(["JuLIP", "StaticArrays", "LinearAlgebra", "Plots","LaTeXStrings"])

using JuLIP
using StaticArrays
using LinearAlgebra
using Plots
using LaTeXStrings

function change_atom(frames, A::Int64, B::Int64, N::Int64)
    new_frames = Vector{Atoms{Float64}}([])
    for atoms in frames
        # Find indices of all Si atoms (Z == 14)
        indices = findall(z -> z == A, atoms.Z)
        
        # Check if there are at least 4 Si atoms
        if length(indices) >= N+1
            # Change the 3rd and 4th Si atoms to Carbon (Z = 6)
            for i in range(1, N)
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
frames_modified = change_atom(frames, 6, 14, 4)
frames_modified2 = change_atom(frames_modified, 16, 8, 1)

# --- Save the modified frames --Cwrite_extxyz("all_bead_modified.extxyz", frames_modified)
write_extxyz("all_C_modified.extxyz", frames_modified2)

frames = read_extxyz("all_C_dimers_removedS.extxyz")
frames_modified = change_atom(frames, 6, 14, 4)

# --- Save the modified frames --Cwrite_extxyz("all_bead_modified.extxyz", frames_modified)
write_extxyz("all_CSi_modified.extxyz", frames_modified)

function change_C_to_Si(frames)
    new_frames = Vector{Atoms{Float64}}([])
    for atoms in frames
        # Find indices of all Si atoms (Z == 14)
        si_indices = findall(z -> z == 6, atoms.Z)
        
        # Check if there are at least 4 Si atoms
        if length(si_indices) >= 4
            # Change the 3rd and 4th Si atoms to Carbon (Z = 6)
            atoms.Z[si_indices[1]] = 14
            atoms.Z[si_indices[2]] = 14
            atoms.Z[si_indices[3]] = 14
            atoms.Z[si_indices[4]] = 14
        else
            @warn "Frame has fewer than 4 Si atoms. Skipping change."
        end

        # Keep properties like :theta, :phi, etc.
        push!(new_frames, atoms)
    end
    return new_frames
end

function change_S_to_O(frames)
    new_frames = Vector{Atoms{Float64}}([])
    for atoms in frames
        # Find indices of all Si atoms (Z == 14)
        si_indices = findall(z -> z == 16, atoms.Z)
        
        # Check if there are at least 4 Si atoms
        if length(si_indices) >= 1
            # Change the 3rd and 4th Si atoms to Carbon (Z = 6)
            atoms.Z[si_indices[1]] = 14
        else
            @warn "Frame has fewer than 4 Si atoms. Skipping change."
        end

        # Keep properties like :theta, :phi, etc.
        push!(new_frames, atoms)
    end
    return new_frames
end

function change_O_to_S(frames)
    new_frames = Vector{Atoms{Float64}}([])
    for atoms in frames
        # Find indices of all Si atoms (Z == 14)
        si_indices = findall(z -> z == 16, atoms.Z)
        
        # Check if there are at least 4 Si atoms
        if length(si_indices) >= 1
            # Change the 3rd and 4th Si atoms to Carbon (Z = 6)
            atoms.Z[si_indices[1]] = 6
        else
            @warn "Frame has fewer than 4 Si atoms. Skipping change."
        end

        # Keep properties like :theta, :phi, etc.
        push!(new_frames, atoms)
    end
    return new_frames
end


frames = read_extxyz("all_C_dimers.extxyz")
frames_modified = change_C_to_Si(frames)
frames_modified2 = change_S_to_O(frames_modified)
frames_modified3 = change_O_to_S(frames_modified2)
write_extxyz("all_CSi.extxyz", frames_modified3)