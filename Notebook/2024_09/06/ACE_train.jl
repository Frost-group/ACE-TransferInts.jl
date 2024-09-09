### A Pluto.jl notebook ###
# v0.19.45

using Markdown
using InteractiveUtils

# ╔═╡ e212233e-6c4c-11ef-115e-476408bfabc1
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

# ╔═╡ 018914e2-b683-4ce1-b38d-c0a3f1a418b0
begin
	Pkg.activate(".")
		Pkg.Registry.add("General")  # only needed when installing Julia for the first time
		Pkg.Registry.add(RegistrySpec(url="https://github.com/ACEsuit/ACEregistry"))
		Pkg.add("ACEpotentials")
		using ACEpotentials
end

# ╔═╡ 953fd3da-8994-4c78-baa9-f1588ad3e5bd
Pkg.status()

# ╔═╡ eea2c9b4-94bd-4947-a4c3-54023c90e6d9
bead_tiny_dataset = read_extxyz("train_bead.extxyz")

# ╔═╡ aa54b142-89cc-48e9-bc0a-746fffc4809c
Dataset_energies 

# ╔═╡ fb701de3-ddd4-4d43-8aac-cdb58e7c3a32
begin
	r_cut = 20.0
	rdf_tiny = ACEpotentials.get_rdf(bead_tiny_dataset, r_cut; rescale = true)
	plt_rdf_1 = histogram(rdf_tiny[(:Ge, :Ge)], bins=150, label = "rdf",
		                      title="Si_tiny_dataset", titlefontsize=10,
		                      xlabel = L"r[\AA]", ylabel = "RDF", yticks = [],
		                      xlims=(1.5,40), size=(400,200), left_margin = 2Plots.mm)
		vline!(rnn(:Ge)*[1.0, 1.633, 1.915, 2.3, 2.5], label = "r1, r2, ...", lw=5, color = "black")
end

# ╔═╡ 899acd5f-ebde-4e34-804d-3d805274e10e
rdf_tiny

# ╔═╡ d54bfe21-dc89-41e6-8fa2-f1f1138ef006
bead_dataset = read_extxyz("test_bead.extxyz")

# ╔═╡ 74b64aa9-80f6-4493-869f-b1bddb9c808b
begin
	model = acemodel(elements = [:Si,],
	                 order = 5,
	                 totaldegree = 10,
	                 rcut = 5.5,
	                 Eref = [:Si => 0.0])
	@show length(model.basis);
end

# ╔═╡ b497ddf4-9b01-43f7-9b71-c32a42b1c09e
model

# ╔═╡ 94ee83f4-262f-4476-b7a7-340ee806e874
solver = ACEfit.QR(lambda=1e-1)

# ╔═╡ 88c6460b-05f3-4be6-9929-2d503bdc68d2
data_keys = (energy_key = "energy")

# ╔═╡ 76290b18-7c23-43d0-a567-e1e7b47fde60
acefit!(model, bead_tiny_dataset; solver=solver, data_keys...);

# ╔═╡ Cell order:
# ╠═e212233e-6c4c-11ef-115e-476408bfabc1
# ╠═018914e2-b683-4ce1-b38d-c0a3f1a418b0
# ╠═953fd3da-8994-4c78-baa9-f1588ad3e5bd
# ╠═eea2c9b4-94bd-4947-a4c3-54023c90e6d9
# ╠═899acd5f-ebde-4e34-804d-3d805274e10e
# ╠═aa54b142-89cc-48e9-bc0a-746fffc4809c
# ╠═fb701de3-ddd4-4d43-8aac-cdb58e7c3a32
# ╠═d54bfe21-dc89-41e6-8fa2-f1f1138ef006
# ╠═74b64aa9-80f6-4493-869f-b1bddb9c808b
# ╠═b497ddf4-9b01-43f7-9b71-c32a42b1c09e
# ╠═94ee83f4-262f-4476-b7a7-340ee806e874
# ╠═88c6460b-05f3-4be6-9929-2d503bdc68d2
# ╠═76290b18-7c23-43d0-a567-e1e7b47fde60
