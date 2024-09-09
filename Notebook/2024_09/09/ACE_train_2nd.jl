### A Pluto.jl notebook ###
# v0.19.45

using Markdown
using InteractiveUtils

# ╔═╡ 239d539c-6eb6-11ef-343d-bb46d309ea33
begin
	# add and load general packages used in this notebook.
	using Pkg
	Pkg.activate(".")
	Pkg.add("LaTeXStrings")
	Pkg.add("MultivariateStats")
	Pkg.add("Plots")
	Pkg.add("Suppressor")
	using LaTeXStrings, MultivariateStats, Plots, Printf, Statistics, Suppressor
end;

# ╔═╡ 426292c8-f314-46d9-895f-d18a41c108fc
begin
	Pkg.activate(".")
		Pkg.Registry.add("General")  # only needed when installing Julia for the first time
		Pkg.Registry.add(RegistrySpec(url="https://github.com/ACEsuit/ACEregistry"))
		Pkg.add("ACEpotentials")
		using ACEpotentials
end

# ╔═╡ 78beae24-2848-43b7-8c6e-8068f01c769f
bead_tiny_dataset = read_extxyz("Naphthalene/2nd_train_bead.extxyz")

# ╔═╡ c57a246a-fbdd-44dc-965d-a71de8a841ac
bead_dataset = read_extxyz("Naphthalene/2nd_test_bead.extxyz")

# ╔═╡ 6abe9b78-2d15-471f-a64a-3aa3529b3d54
begin
	r_cut = 5.0
	rdf_tiny = ACEpotentials.get_rdf(bead_tiny_dataset, r_cut; rescale = true)
	plt_rdf_1 = histogram(rdf_tiny[(:Si, :Si)], bins=150, label = "rdf",
		                      title="Bead_tiny_dataset", titlefontsize=10,
		                      xlabel = L"r[\AA]", ylabel = "RDF", yticks = [],
		                      xlims=(2.4,2.5), ylims=(0.0, 5), size=(400,200), left_margin = 2Plots.mm)
		#vline!(rnn(:Ge)*[1.0, 1.633, 1.915, 2.3, 2.5], label = "r1, r2, ...", lw=5, color = "black")
	plot!(plt_rdf_1, size=(600,400))
end

# ╔═╡ 512b58a3-10f6-4409-89e1-e721789eaf71
rdf_tiny

# ╔═╡ d8622929-c1da-4b26-b3db-c066a92b19c5
begin
	model1 = acemodel(elements = [:Si,],
	                 order = 4,
	                 totaldegree = 3,
	                 rcut = 2.5,
	                 Eref = [:Si => 0.0])
	@show length(model1.basis);
end

# ╔═╡ 7f101e48-e78b-4d50-86a1-66c44ba3f95a
solver = ACEfit.BLR(committee_size=200, factorization=:svd)

# ╔═╡ b42a2fc9-8c8c-4c15-9465-e423bfd494eb
data_keys1 = (energy_key = "energy", force_key = "force")

# ╔═╡ 4dcd8478-abdb-40ce-9dfd-f7b0619b5b53
acefit!(model1, bead_tiny_dataset; solver=solver, data_keys1...);

# ╔═╡ 21ac6f21-2cc2-4d6e-99b9-a1017bb36eb9
begin
	@info("Training Errors")
	ACEpotentials.linear_errors(bead_tiny_dataset, model1; data_keys1...);
	
	@info("Test Error")
	ACEpotentials.linear_errors(bead_dataset, model1; data_keys1...);
end

# ╔═╡ Cell order:
# ╠═239d539c-6eb6-11ef-343d-bb46d309ea33
# ╠═426292c8-f314-46d9-895f-d18a41c108fc
# ╠═78beae24-2848-43b7-8c6e-8068f01c769f
# ╠═c57a246a-fbdd-44dc-965d-a71de8a841ac
# ╠═512b58a3-10f6-4409-89e1-e721789eaf71
# ╠═6abe9b78-2d15-471f-a64a-3aa3529b3d54
# ╠═d8622929-c1da-4b26-b3db-c066a92b19c5
# ╠═7f101e48-e78b-4d50-86a1-66c44ba3f95a
# ╠═b42a2fc9-8c8c-4c15-9465-e423bfd494eb
# ╠═4dcd8478-abdb-40ce-9dfd-f7b0619b5b53
# ╠═21ac6f21-2cc2-4d6e-99b9-a1017bb36eb9
