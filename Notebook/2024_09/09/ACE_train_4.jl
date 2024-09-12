### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ 5a004da0-703d-11ef-210b-3ff7e0bcff4b
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

# ╔═╡ 1ec31829-6c47-4b70-a1ed-fe01af7dd29b
begin
	Pkg.activate(".")
	Pkg.Registry.add("General")  # only needed when installing Julia for the first time
	Pkg.Registry.add(RegistrySpec(url="https://github.com/ACEsuit/ACEregistry"))
	Pkg.add("ACEpotentials")
end

# ╔═╡ 149b3b50-00c8-4982-9328-2b7e69b09791
begin
	pkg"activate pathtoproject"
	Pkg.instantiate()
end

# ╔═╡ ed8a8bd6-c78a-4387-b285-7967f829f8c8
using ACEpotentials

# ╔═╡ 4c09ee04-7f74-4eb7-aea0-18f2222e43c9
Pkg.instantiate()

# ╔═╡ Cell order:
# ╠═4c09ee04-7f74-4eb7-aea0-18f2222e43c9
# ╠═5a004da0-703d-11ef-210b-3ff7e0bcff4b
# ╠═1ec31829-6c47-4b70-a1ed-fe01af7dd29b
# ╠═149b3b50-00c8-4982-9328-2b7e69b09791
# ╠═ed8a8bd6-c78a-4387-b285-7967f829f8c8
