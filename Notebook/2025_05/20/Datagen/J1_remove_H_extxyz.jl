cd("C:/Users/Keerati/Documents/GitHub/ACE-TransferInts.jl/Notebook/2025_05/20/Datagen")
# --- Function to remove H atoms from all frames ---
function remove_hydrogen(frames)
    new_frames = Vector{Atoms{Float64}}([])
    for atoms in frames
        # Identify non-Hydrogen atoms (atomic number Z != 1)
        keep = findall(z -> z != 1, atoms.Z)
        
        # Filter species and positions
        new_Z = atoms.Z[keep]
        new_X = atoms.X[keep]
        
        # Create new Atoms object
        new_atoms = JuLIP.Atoms(new_Z, new_X)
        
        set_cell!(new_atoms, 20.0 * I(3)) 
        new_atoms["energy"] = atoms["energy_eff"]
        push!(new_frames, new_atoms)
    end
    return new_frames
end
num = 625
rot = 30
# --- Load your file ---
frames = read_extxyz("output"*string(num)*"_rot"*string(rot)*".extxyz")

# --- Remove Hydrogens ---
frames_noH = remove_hydrogen(frames)

frames_noH[1]["energy"]

# --- Save the new file ---
write_extxyz("all_C_dimers"*string(num)*"_rot"*string(rot)*".extxyz", frames_noH)

