#cd("C:/Users/Keerati/Documents/GitHub/ACE-TransferInts.jl/Notebook/2025_05/20/Datagen")
cd("/Users/kk419/Documents/GitHub/ACE-TransferInts.jl/Notebook/2025_05/20/Datagen")

# --- Load molecule ---
mol = read_extxyz("ethylene.extxyz")[1]


# --- Geometric center function ---
function centroid(positions)
    n = length(positions)
    return sum(positions) / n
end

# --- Center first molecule at origin ---
center1 = centroid(mol.X)
mol1 = JuLIP.Atoms(mol.Z, [pos - center1 for pos in mol.X])

# --- Generate a unit vector from spherical angles θ, φ ---
function spherical_to_cartesian(θ, φ)
    x = sin(θ) * cos(φ)
    y = sin(θ) * sin(φ)
    z = cos(θ)
    return SVector(z, y, x)
end
# Define rotation matrix for 60 degrees about z-axis
rot = 60
θ_rot = deg2rad(rot)
Rz = [cos(θ_rot) -sin(θ_rot) 0.0;
      sin(θ_rot)  cos(θ_rot) 0.0;
      0.0         0.0        1.0]


# --- Prepare θ and φ grids ---
θ_list = range(0.01, stop=π, length=17);
φ_list = range(0.01, stop=2π, length=33);

# --- Create and collect all dimers ---
all_dimers = Vector{Atoms{Float64}}([])
#JuLIP.Atoms[]
scale = 625
for θ in θ_list, φ in φ_list
    # Shift direction
    direction = spherical_to_cartesian(θ, φ)
    shift_vector = scale/100 * direction  # 6 Å along rotated unit vector
    # Create second molecule
    mol2 = JuLIP.Atoms(mol1.Z, [Rz*pos + shift_vector for pos in mol1.X])
    
    # Combine both molecules
    all_species = vcat(mol1.Z, mol2.Z)
    all_positions = vcat(mol1.X, mol2.X)
    dimer = JuLIP.Atoms(all_species, all_positions)
    set_cell!(dimer, 20.0 * I(3))  # Large enough box

    # Optionally add properties
    dimer[:theta] = θ
    dimer[:phi] = φ

    # Collect
    push!(all_dimers, dimer)
end

write_extxyz("dimer_ethylene"*string(scale)*"_rot"*string(rot)*".extxyz", all_dimers)


frames = read_extxyz("dimer_ethylene"*string(scale)*"_rot"*string(rot)*".extxyz")
function process_xyz_all(filename, frames)
    open(filename, "w") do io
        for atoms in frames
            println(io, length(atoms))  # Number of atoms

            # Write theta and phi into the comment line
            theta = get_data(atoms, "theta")
            phi = get_data(atoms, "phi")
            println(io, "theta=$(theta) phi=$(phi)")
            #println(io, "")  # Blank comment line

            # Atom lines
            for (z, x) in zip(atoms.Z, atoms.X)
                element = chemical_symbol(z)  # Atomic number to symbol
                println(io, "$(element) $(x[1]) $(x[2]) $(x[3])")
            end
        end
    end
end

# --- Save all 36 frames into a single file ---
process_xyz_all("dimer_ethylene"*string(scale)*"_rot"*string(rot)*".xyz", frames)

frames[1].X == frames[2].X
#get_data(frames[1], "theta")