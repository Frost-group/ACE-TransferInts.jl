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

# ╔═╡ d54bfe21-dc89-41e6-8fa2-f1f1138ef006
bead_dataset = read_extxyz("test_bead.extxyz")

# ╔═╡ 74b64aa9-80f6-4493-869f-b1bddb9c808b
begin
	model = acemodel(elements = [:Si,],
		                 order = 3,
		                 totaldegree = 4,
		                 rcut = 13.5,
		                 Eref = [:Si => 0.0125]
						 )
		@show length(model.basis);
end

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
# ╠═d54bfe21-dc89-41e6-8fa2-f1f1138ef006
# ╠═74b64aa9-80f6-4493-869f-b1bddb9c808b
# ╠═94ee83f4-262f-4476-b7a7-340ee806e874
# ╠═88c6460b-05f3-4be6-9929-2d503bdc68d2
# ╠═76290b18-7c23-43d0-a567-e1e7b47fde60
