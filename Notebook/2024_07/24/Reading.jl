### A Pluto.jl notebook ###
# v0.19.45

using Markdown
using InteractiveUtils

# ╔═╡ 4c6fc2f8-44f8-11ef-30b4-6118ceefb650
begin
	using CSV
	using DataFrames
	using DelimitedFiles
end

# ╔═╡ ab81d644-20d3-4617-b6cd-fb92849f334b
begin
	using JuLIP
	using ACEpotentials
	#using ExtXYZ
end

# ╔═╡ ae2a999b-a9e0-4c32-b3f9-c45cb38d56cf
using StaticArrays, LinearAlgebra

# ╔═╡ e55a98a6-924c-4bf8-a79c-c12e46532b0d
begin
	using Plots
	using LaTeXStrings
end

# ╔═╡ e7919cb0-09e4-4a5f-a0ce-4c5088ceb813
#df_J = DataFrame(CSV.File("out.dat", delim=" "))

# ╔═╡ 738ee49a-d3b4-4ac7-9076-d96d900404f9
# ╠═╡ disabled = true
#=╠═╡
#c10 = DataFrame(CSV.File("napthalene.tex", delim=" "))
  ╠═╡ =#

# ╔═╡ 4b2668af-48ac-49ba-8dc8-f14ea08006c8
function read_J(namefile)
	data_J = (CSV.File(namefile, delim=" ")|> Tables.matrix)
	para, J, logJ = data_J[:, 1:6], data_J[:, 7], data_J[:, 8]
	return para, J, logJ
end

# ╔═╡ 2436ae5d-f680-4322-8059-c1761688d4a5
begin
	#Rotation about x, y, z-axis
	function rot_x(θ)
		return [[1.0, 0.0, 0.0] [0.0, cos(θ), sin(θ)] [0.0, -sin(θ), cos(θ)]]
	end
	function rot_y(θ)
		return [[cos(θ), 0.0, -sin(θ)] [0.0, 1.0, 0.0] [sin(θ), 0.0, cos(θ)]]
	end
	function rot_z(θ)
		return [[cos(θ), sin(θ), 0.0] [-sin(θ), cos(θ), 0.0] [0.0, -sin(θ), 1.0]]
	end
end

# ╔═╡ 15c07f11-a3af-41b4-abeb-ddd110ce5c12
function construct_dimer(namefile, savefile, para, rot_ang_y)
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
		end
		push!(Set_Mol2, mol2) 	
		#println(Set_Mol2)
		#println("_________________________________________________________")
	end
	#return Set_Mol2
	write_extxyz(savefile, Set_Mol2)
end

# ╔═╡ 6fec2a43-14fa-42be-bfb0-9530de667a83
para, J, logJ = read_J("out.dat")

# ╔═╡ 61187086-ab9c-4aeb-8dc6-04227e862270
#Find element
function find_element(Molecule, Z)
	return findall(x->x==Z, Molecule.Z)
end

# ╔═╡ 9f36b8e9-2b50-4819-902d-8c14fc67bfed
# ╠═╡ disabled = true
#=╠═╡
function collect_allC_ring(num_C, shared_C)
	ring_C = zeros((2,6))
	for jj in 1:2
		for ii in 1:6
			if
		
		end
	end
	return
end
  ╠═╡ =#

# ╔═╡ 83ece228-eae5-4805-beac-218c459cbbc8
#C10 = read_frames("napthalene.xyz")

# ╔═╡ 6201194b-bafd-4fee-b353-85e6b880bc4f
#write_frames("naphthalene.extxyz", C10, append=true)

# ╔═╡ c39b69ee-51d9-4df0-941d-f1d477c8a067


# ╔═╡ 1f3f1281-6539-4bc7-8874-a952b5089767
function Centroid_bead(Molecule, shared_C, num_C)
	cent_Vec = [] 					#collect all bead centroids
	
	ring_C = zeros(Int64,(2,6))		#collect Carbon atoms (Shared + Unshared)
	ring_C[1, 1:6] = num_C[1:6]
	ring_C[2, 1:2] = shared_C
	ring_C[2, 3:6] = num_C[7:10]
	
	for jj in 1:2
		vec = zeros((3))
		bead_mass = 0
		for ii in 1:6
			vec += Molecule.X[ring_C[jj, ii]]*Molecule.M[ring_C[jj, ii]]
			bead_mass += Molecule.M[ring_C[jj, ii]]
		end
		#println(vec)
		#println(bead_mass)
		#println(vec/bead_mass)
		#println(vec/78.11)
		push!(cent_Vec, vec/bead_mass) #collect centroid positions with the CM of the beads
	end
	return cent_Vec
end

# ╔═╡ 8b6804f2-0d25-4cf8-b4da-2f3e638399a6
const JVec = SVector{3}

# ╔═╡ 919a4e3f-a121-4504-aa96-01967b967c6b


# ╔═╡ 2d935524-f013-4ad5-b6d8-de6f394fc713


# ╔═╡ 32ac6e32-50a4-4d43-9464-60f2ecb60c43


# ╔═╡ b1dfd94d-94c9-450a-bb38-36e098a7162a


# ╔═╡ 591c8375-9601-41aa-a540-b6b323b05ccc


# ╔═╡ 8892ca90-4b2e-404c-810c-73397a88b3b1
begin
	mol1 = read_extxyz("naphthalene.extxyz")[1]
	Set_Mol2 =  read_extxyz("dimer_naphthalene.extxyz")
	nn = 400
	mol2 = Set_Mol2[nn]
end

# ╔═╡ 817d7fc2-d928-4bf4-9059-dfd09e12a2ca
begin
	num_C = find_element(mol1, 6)
	shared_C = [3, 4];
end;

# ╔═╡ 0015d7e0-3fe7-4605-bdc9-2d47e32533d1
begin
	ring_C = zeros(Int64,(2,6))
	ring_C[1, 1:6] = num_C[1:6]
	ring_C[2, 1:2] = shared_C
	ring_C[2, 3:6] = num_C[7:10]
end;

# ╔═╡ b0d0efb5-9570-4d80-8fc7-d1849bbfa1f2
ring_C

# ╔═╡ ca265187-c7f9-42d6-bf73-360011477509
x1 = Centroid_bead(mol1, shared_C, num_C)

# ╔═╡ 9d3d33f4-d1c5-4360-8f53-29d900ca5b5e
norm(x1[1,:]-x1[2,:])

# ╔═╡ 587ab7d6-588a-4375-8f69-1e00333c0efc
x2 = Centroid_bead(mol2, shared_C, num_C)

# ╔═╡ 1c674758-faaa-4017-9f8c-833253ae6778
norm(x2[1,:]-x2[2,:])

# ╔═╡ a238c2e8-e75b-4484-8277-ccd04bf40a50


# ╔═╡ 088ad824-ed48-45d8-8bea-849b42b5cdf7
mol2.X

# ╔═╡ b932c8f9-d29a-4dc8-8bd9-9a4335f2bb54
begin
	J_mat = zeros((49,Int(size(J,1)/49)))
	logJ_mat = zeros((49,Int(size(J,1)/49)))
	θ = zeros((Int(size(J,1)/49)))
	ϕ = zeros((49))
end

# ╔═╡ 33feb269-3620-43f1-8d31-66e8b9472481
begin
	for bb in 1:Int(size(J,1)/49)
		J_mat[:, bb] = J[49*(bb-1)+1: 49*bb]	
		logJ_mat[:, bb] = logJ[49*(bb-1)+1: 49*bb]
		θ[bb] = para[49*bb,1]
	end
	for aa in 1:49
		ϕ[aa] = para[aa,2]
	end
end

# ╔═╡ 1f58b44f-bae8-420a-857d-2a1a5fbea8bf
heatmap(ϕ[1:49],
	    θ[1:97], J_mat[:,:],
	    c=cgrad(["#064635","#519259", "#96BB7C", "#F0BB62", "#FAD586","#F4EEA9"]),
	    xlabel= L"$\phi$",
		ylabel=L"$\theta$",
		#aspectratio=:equal,
		size=(700, 700),
	    title="Transfer Integrals" )

# ╔═╡ a9877a72-8601-4c24-9a26-10d51e3b96b6
heatmap(ϕ[1:49],
	    θ[1:97], logJ_mat[:,:],
	    c=cgrad(["#064635","#519259", "#96BB7C", "#F0BB62", "#FAD586","#F4EEA9"]),
	    xlabel="feature coord. (n x d)",
		ylabel="feature coord. (n x d)",
		#aspectratio=:equal,
		size=(700, 700),
	    title="PbTe_FC2 (Traning Data = " *string(199) *")" )

# ╔═╡ 9e546e9c-a521-4b88-b8ad-7ed0ecd49966
begin
	θ_grid = [θ for _ in ϕ]  # Duplicate θ for each φ
	ϕ_grid = [ϕ for _ in θ]'  # Duplicate φ for each θ and transpose
end

# ╔═╡ 42212782-5b70-4364-b2bb-a75ee65e00c5
plot(θ,ϕ,
	J_mat,
	st=:surface,
	xlabel= L"$\phi$",
	ylabel=L"$\theta$",
	zlabel=L"J",
	c=cgrad(["#064635","#519259", "#96BB7C", "#F0BB62", "#FAD586","#F4EEA9"]))

# ╔═╡ ad03089d-4b66-4f93-814e-b379e16da537
plot(θ,ϕ,logJ_mat,st=:surface,
	xlabel= L"$\phi$",
	ylabel=L"$\theta$",
	zlabel=L"log(J)",
	c=cgrad(["#064635","#519259", "#96BB7C", "#F0BB62", "#FAD586","#F4EEA9"]))

# ╔═╡ Cell order:
# ╠═4c6fc2f8-44f8-11ef-30b4-6118ceefb650
# ╠═e7919cb0-09e4-4a5f-a0ce-4c5088ceb813
# ╠═738ee49a-d3b4-4ac7-9076-d96d900404f9
# ╠═4b2668af-48ac-49ba-8dc8-f14ea08006c8
# ╠═2436ae5d-f680-4322-8059-c1761688d4a5
# ╠═15c07f11-a3af-41b4-abeb-ddd110ce5c12
# ╠═6fec2a43-14fa-42be-bfb0-9530de667a83
# ╠═ab81d644-20d3-4617-b6cd-fb92849f334b
# ╠═61187086-ab9c-4aeb-8dc6-04227e862270
# ╠═817d7fc2-d928-4bf4-9059-dfd09e12a2ca
# ╠═9f36b8e9-2b50-4819-902d-8c14fc67bfed
# ╠═0015d7e0-3fe7-4605-bdc9-2d47e32533d1
# ╠═b0d0efb5-9570-4d80-8fc7-d1849bbfa1f2
# ╠═83ece228-eae5-4805-beac-218c459cbbc8
# ╠═6201194b-bafd-4fee-b353-85e6b880bc4f
# ╠═c39b69ee-51d9-4df0-941d-f1d477c8a067
# ╠═1f3f1281-6539-4bc7-8874-a952b5089767
# ╠═ca265187-c7f9-42d6-bf73-360011477509
# ╠═587ab7d6-588a-4375-8f69-1e00333c0efc
# ╠═ae2a999b-a9e0-4c32-b3f9-c45cb38d56cf
# ╠═8b6804f2-0d25-4cf8-b4da-2f3e638399a6
# ╠═1c674758-faaa-4017-9f8c-833253ae6778
# ╠═9d3d33f4-d1c5-4360-8f53-29d900ca5b5e
# ╠═919a4e3f-a121-4504-aa96-01967b967c6b
# ╠═2d935524-f013-4ad5-b6d8-de6f394fc713
# ╠═32ac6e32-50a4-4d43-9464-60f2ecb60c43
# ╠═b1dfd94d-94c9-450a-bb38-36e098a7162a
# ╠═591c8375-9601-41aa-a540-b6b323b05ccc
# ╠═8892ca90-4b2e-404c-810c-73397a88b3b1
# ╠═a238c2e8-e75b-4484-8277-ccd04bf40a50
# ╠═088ad824-ed48-45d8-8bea-849b42b5cdf7
# ╠═b932c8f9-d29a-4dc8-8bd9-9a4335f2bb54
# ╠═33feb269-3620-43f1-8d31-66e8b9472481
# ╠═e55a98a6-924c-4bf8-a79c-c12e46532b0d
# ╠═1f58b44f-bae8-420a-857d-2a1a5fbea8bf
# ╠═a9877a72-8601-4c24-9a26-10d51e3b96b6
# ╠═9e546e9c-a521-4b88-b8ad-7ed0ecd49966
# ╠═42212782-5b70-4364-b2bb-a75ee65e00c5
# ╠═ad03089d-4b66-4f93-814e-b379e16da537
