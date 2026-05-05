"""
Rotation matrices and dimer construction utilities.
"""

function rot_x(θ)
    return [[1.0, 0.0, 0.0] [0.0, cos(θ), sin(θ)] [0.0, -sin(θ), cos(θ)]]
end

function rot_y(θ)
    return [[cos(θ), 0.0, -sin(θ)] [0.0, 1.0, 0.0] [sin(θ), 0.0, cos(θ)]]
end

function rot_z(θ)
    return [[cos(θ), sin(θ), 0.0] [-sin(θ), cos(θ), 0.0] [0.0, 0.0, 1.0]]
end

"""
    centroid(positions) -> SVector

Arithmetic mean position of a collection of atomic coordinates.
"""
function centroid(positions)
    return sum(positions) / length(positions)
end

"""
    construct_dimer(monomer_file, save_file, params, rot_angle_y, J_values)

Build a set of dimer configurations by rotating and translating a monomer.
`params` rows are [θ_x, θ_y, θ_z, dx, dy, dz]; only the translation columns
[4:6] and `rot_angle_y` are used here (extend as needed).
"""
function construct_dimer(monomer_file, save_file, params, rot_angle_y, J_values)
    mol1   = read_extxyz(monomer_file)[1]
    natoms = size(mol1.M, 1)
    ndata  = size(params, 1)
    dimers = Vector{Atoms{Float64}}()
    for jj in 1:ndata
        mol2 = read_extxyz(monomer_file)[1]
        for ii in 1:natoms
            mol2.X[ii] = rot_y(rot_angle_y) * mol1.X[ii] .+ params[jj, 4:6]
        end
        set_data!(mol2, "energy", J_values[jj])
        push!(dimers, mol2)
    end
    write_extxyz(save_file, dimers)
end
