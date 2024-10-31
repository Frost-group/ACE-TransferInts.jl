begin
    using Pkg
    Pkg.activate(".")
    Pkg.add("CSV")
    Pkg.add("DataFrames")
    Pkg.add("DelimitedFiles")
    Pkg.add("LaTeXStrings")
    Pkg.add("MultivariateStats")
    Pkg.add("Plots")
    Pkg.add("Suppressor")
    Pkg.add("ExtXYZ")
    Pkg.add("JuLIP")
    Pkg.add("Plots")
    Pkg.add("LaTeXStrings")
    
    using CSV
    using DataFrames
    using DelimitedFiles
    using JuLIP
    using Plots
    using LaTeXStrings
    using Plots
    using LaTeXStrings
end

function extract_dimer(oldfile, n)
    olddimer = read_extxyz(oldfile)
    newdimer = Vector{Atoms{Float64}}([])

    dimer = olddimer[n]
    atoms = Atoms(X = dimer.X, Z = dimer.Z)
    push!(newdimer, atoms)

    write_extxyz("dimer_"*string(n)*".xyz", newdimer)
end

extract_dimer("all_CH_dimers.extxyz", 1)