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


using CSV
using DataFrames
using DelimitedFiles
using JuLIP
using ACEpotentials
using Plots
using LaTeXStrings

function read_J(namefile)
	data_J = (CSV.File(namefile, delim=" ")|> Tables.matrix)
	para, J, logJ = data_J[:, 1:6], data_J[:, 7], data_J[:, 8]
	return para, J, logJ
end

para, J, logJ = read_J("out.dat");

function rot_x(θ)
    return [[1.0, 0.0, 0.0] [0.0, cos(θ), sin(θ)] [0.0, -sin(θ), cos(θ)]]
end
function rot_y(θ)
    return [[cos(θ), 0.0, -sin(θ)] [0.0, 1.0, 0.0] [sin(θ), 0.0, cos(θ)]]
end
function rot_z(θ)
    return [[cos(θ), sin(θ), 0.0] [-sin(θ), cos(θ), 0.0] [0.0, -sin(θ), 1.0]]
end

function construct_dimer(namefile, savefile, para, rot_ang_y, J)
	mol1 = read_extxyz(namefile)[1] 	#Read Molecular Structure from file.extxyz
	ndata = size(para, 1) 				#Collect number of data from TranInt out.dat
	natoms = size(mol1.M)[1] 			#Collect number of atoms in molecules
	Set_Mol2 = Vector{Atoms{Float64}}([])						#Set of 2nd Molecules
	for jj in 1:ndata
		mol2 = read_extxyz(namefile)[1]
		#println("Dimer: "*string(jj))
		for ii in 1:natoms
			mol2.X[ii] = rot_y(rot_ang_y) * mol1.X[ii] + para[jj, 4:6] 		#Rotation + Translation
			#mol2.X[ii] = rot_z(para[jj, 2])*rot_y(para[jj, 1]) * mol1.X[ii] + para[jj, 4:6] 		#Rotation + Translation
		#println(mol2.X[ii])
			set_data!(mol2, "energy", J[jj])
		end
		push!(Set_Mol2, mol2) 	
		#println(Set_Mol2)
		#println("_________________________________________________________")
	end
	#return Set_Mol2
	write_extxyz(savefile, Set_Mol2)
end

construct_dimer("naphthalene.extxyz", "dimer_naphthalene.extxyz", para, 2*pi/3, J);

function construct_all_dimer(namefile1, namefile2, savefile)
	mol1 = read_extxyz(namefile1)[1] 			#Read Molecular Structure from 													file.extxyz
	Set_Mol2 =  read_extxyz(namefile2)	
	Set_dimer = Vector{Atoms{Float64}}([])		#Set of dimer Molecules
	
	for nn in 1:length(Set_Mol2)
		mol2 = Set_Mol2[nn]
		#println("Dimer: "*string(jj))
		dimer = Atoms(X = vcat(mol2.X, mol1.X),
			M = vcat(mol2.M, mol1.M),
			Z = vcat(mol2.Z, mol1.Z),
			cell = mol1.cell*100,
			pbc = fill(false, 3)
		)
		
		set_data!(dimer, "energy", 
			get_data(mol2, "energy")
		)
		
		push!(Set_dimer, dimer) 	
		#println(Set_Mol2)
		#println("_________________________________________________________")
	end
	#return Set_Mol2
	write_extxyz(savefile, Set_dimer)
end

construct_all_dimer("naphthalene.extxyz", 
	"dimer_naphthalene.extxyz", 
	"Large_dimers.extxyz"
)

function construct_all_dimer(namefile1, namefile2, savefile)
	mol1 = read_extxyz(namefile1)[1] 			#Read Molecular Structure from 													file.extxyz
	Set_Mol2 =  read_extxyz(namefile2)	
	Set_dimer = Vector{Atoms{Float64}}([])		#Set of dimer Molecules
	
	for nn in 1:length(Set_Mol2)
		mol2 = Set_Mol2[nn]
		#println("Dimer: "*string(jj))
		dimer = Atoms(X = vcat(mol2.X, mol1.X),
			M = vcat(mol2.M, mol1.M),
			Z = vcat(mol2.Z, mol1.Z),
			cell = mol1.cell*100,
			pbc = fill(false, 3)
		)
		
		set_data!(dimer, "energy", 
			get_data(mol2, "energy")
		)
		
		push!(Set_dimer, dimer) 	
		#println(Set_Mol2)
		#println("_________________________________________________________")
	end
	#return Set_Mol2
	write_extxyz(savefile, Set_dimer)
end


construct_all_dimer("naphthalene.extxyz", 
	"dimer_naphthalene.extxyz", 
	"All_CH_dimers.extxyz"
)

function train_test(all_structure_file, num_train, train_file, test_file)
	s_dim =  read_extxyz(all_structure_file)
	
	num_test = length(s_dim) - num_train

	trian_in = ones((1))

	for nn in 1:num_train
		ind = rand(1:length(s_dim))
		if !(ind in trian_in)
			append!(trian_in, ind)
		end
	end
	println(length(trian_in))
	train_index = Int.(trian_in)
	test_index = setdiff(collect(1:length(s_dim)), train_index)
	
	train_dat = Vector{Atoms{Float64}}([])	
	test_dat = Vector{Atoms{Float64}}([])

	for ii in train_index
		confi = s_dim[ii]
		push!(train_dat, confi)
	end
	write_extxyz(train_file, train_dat)
	for jj in test_index
		confi = s_dim[jj]
		push!(test_dat, confi)
	end
	write_extxyz(test_file, test_dat)
	
end
train_test("all_CH_dimers.extxyz", 900, "train_all.extxyz", "test_all.extxyz")

