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

function change_H_to_C_Si(frames)
    new_frames = Vector{Atoms{Float64}}([])
    for atoms in frames
        # Find indices of all Si atoms (Z == 14)
        C_indices = findall(z -> z == 1, atoms.Z)
        
        # Check if there are at least 4 Si atoms
        atoms.Z[C_indices[1]] = 6
        atoms.Z[C_indices[2]] = 6
        atoms.Z[C_indices[3]] = 6
        atoms.Z[C_indices[4]] = 6
        
        atoms.Z[C_indices[5]] = 14
        atoms.Z[C_indices[6]] = 14
        atoms.Z[C_indices[7]] = 14
        atoms.Z[C_indices[8]] = 14

        # Keep properties like :theta, :phi, etc.
        push!(new_frames, atoms)
    end
    return new_frames
end

num = 625
rot = 60
# --- Load the file ---
frames = read_extxyz("all_C_dimers"*string(num)*"_rot"*string(rot)*".extxyz")

# --- Perform the Si -> C substitution ---
frames_modified = change_Si_to_C(frames)


# --- Save the modified frames --Cwrite_extxyz("all_bead_modified.extxyz", frames_modified)
write_extxyz("all_CSi_modified"*string(num)*"_rot"*string(rot)*".extxyz", frames_modified)


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