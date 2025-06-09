using JuLIP
using StaticArrays
using LinearAlgebra
using Rotations
using Random

# === Set random seed for reproducibility ===
seed = 1000
Random.seed!(seed)

# === Load molecule and center ===
mol = read_extxyz("ethylene.extxyz")[1]

function centroid(positions)
    return sum(positions) / length(positions)
end

center1 = centroid(mol.X)
mol1 = JuLIP.Atoms(mol.Z, [pos - center1 for pos in mol.X])

# === Random vector of magnitude in [rmin, rmax] ===
function random_vector_in_shell(rmin, rmax)
    direction = normalize(randn(SVector{3, Float64}))  # random unit vector
    r = rand() * (rmax - rmin) + rmin
    return r * direction
end

# === Apply rotation to molecule around its own centroid ===
function rotate_atoms(atoms::Atoms{T}, θx, θy, θz) where T
    R = RotXYZ(θx, θy, θz)  # Rotation object from Rotations.jl
    c = centroid(atoms.X)
    rotated_positions = [SVector{3, T}(R * (pos - c)) + c for pos in atoms.X]
    return JuLIP.Atoms(atoms.Z, rotated_positions)
end

# === Generate N random dimers ===
N = 50000  # <--- Set desired number of dimers
rmin = 6.25
rmax = 7.00
all_dimers = Vector{Atoms{Float64}}([])

for i in 1:N
    # Random translation vector
    shift_vec = random_vector_in_shell(rmin, rmax)
    x, y, z = shift_vec
    r = norm(shift_vec)
    θ = acos(z / r)
    φ = atan(y, x)
    # Random Euler angles in radians
    θx = 2π #* rand()
    θy = 2π #* rand()
    θz = π/6 #+ rand() * (5π/6 - π/6)

    # Rotate and translate mol2
    mol2 = rotate_atoms(mol1, θx, θy, θz)
    mol2_shifted = JuLIP.Atoms(mol2.Z, [pos + shift_vec for pos in mol2.X])

    # Combine molecules
    species = vcat(mol1.Z, mol2_shifted.Z)
    positions = vcat(mol1.X, mol2_shifted.X)
    dimer = JuLIP.Atoms(species, positions)

    # Set box and metadata
    set_cell!(dimer, 20.0 * I(3))
    dimer[:r] = r
    dimer[:theta] = θ
    dimer[:phi] = φ
    dimer[:θx] = θx
    dimer[:θy] = θy
    dimer[:θz] = θz
    dimer[:shift_vec] = shift_vec  # store full vecto
    #dimer[:r] = norm(shift_vec)

    push!(all_dimers, dimer)
end
string(seed)
# === Write to extxyz ===
write_extxyz("dimer_ethylene_rand"*string(seed)*"_Num"*string(N)*".extxyz", all_dimers)

frames = read_extxyz("dimer_ethylene_rand"*string(seed)*"_Num"*string(N)*".extxyz")
# === Optional: write to .xyz with metadata ===
function process_xyz_all(filename, frames)
    open(filename, "w") do io
        for atoms in frames
            println(io, length(atoms))
            
            # Extract rotation and shift metadata
            θx = get_data(atoms, "θx")
            θy = get_data(atoms, "θy")
            θz = get_data(atoms, "θz")
            shift_vec = get_data(atoms, "shift_vec")  # Expecting SVector
            
            # Write metadata to comment line
            println(io, "theta_x=$(θx) theta_y=$(θy) theta_z=$(θz) shift_x=$(shift_vec[1]) shift_y=$(shift_vec[2]) shift_z=$(shift_vec[3])")
            
            # Atom lines
            for (z, x) in zip(atoms.Z, atoms.X)
                println(io, "$(chemical_symbol(z)) $(x[1]) $(x[2]) $(x[3])")
            end
        end
    end
end


process_xyz_all("dimer_ethylene_rand"*string(seed)*"_Num"*string(N)*".xyz", frames)
