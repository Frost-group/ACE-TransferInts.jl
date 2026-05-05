"""
I/O utilities: reading transfer-integral output files and writing XYZ.
"""

"""
    read_J(filename) -> (params, J, logJ)

Read a whitespace-delimited transfer-integral file with columns
[param1..param6, J, logJ] as produced by the ORCA/TB workflow.
"""
function read_J(filename)
    data  = CSV.File(filename, delim=" ") |> Tables.matrix
    params = data[:, 1:6]
    J      = data[:, 7]
    logJ   = data[:, 8]
    return params, J, logJ
end

"""
    write_xyz(filename, atoms)

Write a JuLIP Atoms object to plain XYZ format.
"""
function write_xyz(filename, atoms)
    open(filename, "w") do io
        println(io, length(atoms))
        println(io, "")
        for (z, x) in zip(atoms.Z, atoms.X)
            elem = chemical_symbol(z)
            println(io, "$(elem) $(x[1]) $(x[2]) $(x[3])")
        end
    end
end

"""
    read_lc_csv(filepath) -> (n_train, rmse, mae, sign_acc)

Read a learning-curve CSV with row-keyed format:
  n_train, v1, v2, ...
  rmse,    v1, v2, ...
  mae,     v1, v2, ...
  sign_acc,v1, v2, ...   (optional)
"""
function read_lc_csv(filepath)
    data = Dict{String,Vector{Float64}}()
    open(filepath, "r") do f
        for line in eachline(f)
            parts = split(strip(line), ",")
            key   = strip(parts[1])
            vals  = parse.(Float64, strip.(parts[2:end]))
            data[key] = vals
        end
    end
    n        = Int.(data["n_train"])
    rmse_v   = data["rmse"]
    mae_v    = data["mae"]
    sign_acc = haskey(data, "sign_acc") ? data["sign_acc"] : fill(NaN, length(n))
    return n, rmse_v, mae_v, sign_acc
end
