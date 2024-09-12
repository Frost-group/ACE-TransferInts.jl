### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ 5e1c2860-7046-11ef-3249-13209bbf0420
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

# ╔═╡ c74523e5-3efa-4719-b92e-8a3c6f1d3c68
Pkg.instantiate()

# ╔═╡ Cell order:
# ╠═c74523e5-3efa-4719-b92e-8a3c6f1d3c68
# ╠═5e1c2860-7046-11ef-3249-13209bbf0420
