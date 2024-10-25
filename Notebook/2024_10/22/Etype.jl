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

Z_bead = [14, 14, 8, 8] #[Si O Si O]

Z_bead_Si = [6, 6, 6, 6]

Z_1234 = [1, 2, 3, 4]

Z_Csame = [6, 6, 14, 14, 6, 6, 6, 6, 6, 6, 6, 6, 14, 14, 6, 6, 6, 6, 6, 6]

Z_Cdiff = [6, 6, 14, 14, 6, 6, 14, 14, 14, 14, 6, 6, 14, 14, 6, 6, 14, 14, 14, 14]

Z_Cdiff2 = [6, 6, 14, 14, 6, 6, 8, 8, 8, 8, 6, 6, 14, 14, 6, 6, 8, 8, 8, 8]

Z_Cdiff3 = [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14]

Z_Cdiff4 = [6, 6, 8, 8, 6, 6, 6, 6, 6, 6, 14, 14, 16, 16, 14, 14, 14, 14, 14, 14]

function element_change(oldfile, newfile, Zatomic)
    olddimer = read_extxyz(oldfile)
    newdimer = Vector{Atoms{Float64}}([])	

    for nn in 1:length(olddimer)
        dimer = olddimer[nn]
        dimer.Z = Zatomic

        push!(newdimer, dimer)
    end
    write_extxyz(newfile, newdimer)
end

element_change("train_bead.extxyz", "train_bead.extxyz", Z_bead_Si)
element_change("test_bead.extxyz", "test_bead.extxyz", Z_bead_Si)

element_change("train_bead.extxyz", "train_bead_SiO2.extxyz", Z_bead)
element_change("test_bead.extxyz", "test_bead_SiO2.extxyz", Z_bead)

element_change("train_bead.extxyz", "train_bead_1234.extxyz", Z_1234)
element_change("test_bead.extxyz", "test_bead_1234.extxyz", Z_1234)

element_change("train_carbon.extxyz", "train_carbon_CSiC.extxyz", Z_Csame)
element_change("test_carbon.extxyz", "test_carbon_CSiC.extxyz", Z_Csame)

element_change("train_carbon.extxyz", "train_carbon_CSiSi.extxyz", Z_Cdiff)
element_change("test_carbon.extxyz", "test_carbon_CSiSi.extxyz", Z_Cdiff)

element_change("train_carbon.extxyz", "train_carbon_CSiO.extxyz", Z_Cdiff2)
element_change("test_carbon.extxyz", "test_carbon_CSiO.extxyz", Z_Cdiff2)

element_change("train_carbon.extxyz", "train_carbon_CSi3.extxyz", Z_Cdiff3)
element_change("test_carbon.extxyz", "test_carbon_CSi3.extxyz", Z_Cdiff3)

element_change("train_carbon.extxyz", "train_carbon_COC_SiSSi.extxyz", Z_Cdiff4)
element_change("test_carbon.extxyz", "test_carbon_COC_SiSSi.extxyz", Z_Cdiff4)