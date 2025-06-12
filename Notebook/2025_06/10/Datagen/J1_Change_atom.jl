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
        if length(C_indices) >= 8
            # Change the 3rd and 4th Si atoms to Carbon (Z = 6)
            atoms.Z[C_indices[5]] = 14
            atoms.Z[C_indices[6]] = 14
            atoms.Z[C_indices[7]] = 14
            atoms.Z[C_indices[8]] = 14
        else
            @warn "Frame has fewer than 4 Si atoms. Skipping change."
        end

        # Keep properties like :theta, :phi, etc.
        push!(new_frames, atoms)
    end
    return new_frames
end

function change_Si_to_C_2(frames)
    new_frames = Vector{Atoms{Float64}}([])
    for atoms in frames
        # Find indices of all Si atoms (Z == 14)
        C_indices = findall(z -> z == 6, atoms.Z)
        
        # Check if there are at least 4 Si atoms
        if length(C_indices) >= 8
            # Change the 3rd and 4th Si atoms to Carbon (Z = 6)
            atoms.Z[C_indices[5]] = 14
            atoms.Z[C_indices[6]] = 14
            atoms.Z[C_indices[7]] = 14
            atoms.Z[C_indices[8]] = 14
        else
            @warn "Frame has fewer than 4 Si atoms. Skipping change."
        end

        S_indices = findall(z -> z == 16, atoms.Z)
        if length(C_indices) >= 2
            # Change the 3rd and 4th Si atoms to Carbon (Z = 6)
            atoms.Z[S_indices[2]] = 8
        else
            @warn "Frame has fewer than 4 Si atoms. Skipping change."
        end

        # Keep properties like :theta, :phi, etc.
        push!(new_frames, atoms)
    end
    return new_frames
end


Nr, Ntheta, Nphi, Nx, Ny, Nz = 20, 20, 20, 1, 1, 1 

fname = @sprintf "all_C_dimers_r%d_t%d_p%d_rx%d_ry%d_rz%d_2.extxyz" Nr Ntheta Nphi Nx Ny Nz
# --- Load the file ---
frames = read_extxyz(fname)

# --- Perform the Si -> C substitution ---
frames_modified = change_Si_to_C(frames)

fname2 = @sprintf "all_CSi_dimers_r%d_t%d_p%d_rx%d_ry%d_rz%d_2.extxyz" Nr Ntheta Nphi Nx Ny Nz
# --- Save the modified frames --Cwrite_extxyz("all_bead_modified.extxyz", frames_modified)
write_extxyz(fname2, frames_modified)

######################################
seed = 1000
N = 8000
# --- Load the file ---
frames = read_extxyz("all_C_dimers_rand"*string(seed)*"_Num"*string(N)*".extxyz")

# --- Perform the Si -> C substitution ---
frames_modified = change_Si_to_C_2(frames)


# --- Save the modified frames --Cwrite_extxyz("all_bead_modified.extxyz", frames_modified)
write_extxyz("all_CSi_SO_modified_rand"*string(seed)*"_Num"*string(N)*".extxyz", frames_modified)

#######################################

seed = 1000
N = 8000
# --- Load the file ---
frames = read_extxyz("all_C_dimers_removedS_rand"*string(seed)*"_Num"*string(N)*".extxyz")

# --- Perform the Si -> C substitution ---
frames_modified = change_Si_to_C(frames)


# --- Save the modified frames --Cwrite_extxyz("all_bead_modified.extxyz", frames_modified)
write_extxyz("all_CSi_modified_rand"*string(seed)*"_Num"*string(N)*".extxyz", frames_modified)

#######################################

# --- Load the file ---
frames = read_extxyz("output"*string(num)*".extxyz")

# --- Perform the Si -> C substitution ---
frames_modified = change_Si_to_C(frames)


# --- Save the modified frames --Cwrite_extxyz("all_bead_modified.extxyz", frames_modified)
write_extxyz("all_CH_modified"*string(num)*".extxyz", frames_modified)


frames = read_extxyz("all_CH_modified"*string(num)*".extxyz")

# --- Perform the Si -> C substitution ---
frames_modified = change_H_to_C_Si(frames)


# --- Save the modified frames --Cwrite_extxyz("all_bead_modified.extxyz", frames_modified)
write_extxyz("all_CH_modified"*string(num)*"_all.extxyz", frames_modified)