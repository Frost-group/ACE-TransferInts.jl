using JuLIP
using StaticArrays
using LinearAlgebra
using Rotations
using Random

# === Set random seed for reproducibility ===
seed = 1000
Random.seed!(seed)

# === Load and center molecule ===
mol = read_extxyz("ethylene.extxyz")[1]

function centroid(positions)
    sum(positions) / length(positions)
end


center1 = centroid(mol.X)
mol1 = JuLIP.Atoms(mol.Z, [pos - center1 for pos in mol.X])

# === Apply rotation to molecule around its own centroid ===
function rotate_atoms(atoms::Atoms{T}, θx, θy, θz) where T
    R = RotXYZ(θx, θy, θz)  # Rotation object from Rotations.jl
    c = centroid(atoms.X)
    rotated_positions = [SVector{3, T}(R * (pos - c)) + c for pos in atoms.X]
    return JuLIP.Atoms(atoms.Z, rotated_positions)
end

# === Generate grid samples ===
function linspace_range(minv, maxv, N)
    N == 1 ? [(minv + maxv)/2] : range(minv, maxv; length=N)
end

function generate_dataset_grid(Nr, Ntheta, Nphi, Nx, Ny, Nz)
    r_vals = linspace_range(6.25, 7.00, Nr)
    θ_vals = linspace_range(0.1, π, Ntheta)
    φ_vals = linspace_range(0.1, 2π, Nphi)
    θx_vals = linspace_range(0.1, π, Nx)
    θy_vals = linspace_range(0.1, π, Ny)
    θz_vals = linspace_range(π/6, 5*π/6, Nz)

    all_dimers = Vector{Atoms{Float64}}([])

    for r in r_vals, θ in θ_vals, φ in φ_vals, θx in θx_vals, θy in θy_vals, θz in θz_vals
        # Spherical to Cartesian conversion for shift vector
        x = r * sin(θ) * cos(φ)
        y = r * sin(θ) * sin(φ)
        z = r * cos(θ)
        shift_vec = SVector(x, y, z)

        # Rotate and shift mol2
        mol2 = rotate_atoms(mol1, θx, θy, θz)
        mol2_shifted = JuLIP.Atoms(mol2.Z, [pos + shift_vec for pos in mol2.X])

        # Combine
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

        push!(all_dimers, dimer)
    end

    return all_dimers
end

# === Main execution ===
Nr, Ntheta, Nphi, Nx, Ny, Nz = 20, 20, 20, 1, 1, 1  # <-- Customize here

all_dimers = generate_dataset_grid(Nr, Ntheta, Nphi, Nx, Ny, Nz)
# === Write to file ===
fname = @sprintf "dimer_grid_r%d_t%d_p%d_rx%d_ry%d_rz%d_2.extxyz" Nr Ntheta Nphi Nx Ny Nz
write_extxyz(fname, all_dimers)

frames = read_extxyz(fname)
# === Optional: write to .xyz with metadata ===
function process_xyz_all(filename, frames)
    open(filename, "w") do io
        for atoms in frames
            println(io, length(atoms))
            
            # Extract rotation and shift metadata
            r = get_data(atoms, "r")
            θ = get_data(atoms, "theta")
            ϕ = get_data(atoms, "phi")
            θx = get_data(atoms, "θx")
            θy = get_data(atoms, "θy")
            θz = get_data(atoms, "θz")
            

            # Write metadata to comment line
            println(io, "r=$(r) theta=$(θ) phi=$(ϕ) theta_x=$(θx) theta_y=$(θy) theta_z=$(θz)")
            
            # Atom lines
            for (z, x) in zip(atoms.Z, atoms.X)
                println(io, "$(chemical_symbol(z)) $(x[1]) $(x[2]) $(x[3])")
            end
        end
    end
end

fname1 = @sprintf "dimer_grid_r%d_t%d_p%d_rx%d_ry%d_rz%d_2.xyz" Nr Ntheta Nphi Nx Ny Nz
process_xyz_all(fname1 , frames)
