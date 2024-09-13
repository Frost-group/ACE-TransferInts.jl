### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ 1f79d0ea-710d-11ef-1d60-375671b1b654
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

# ╔═╡ 229beb96-de55-4fff-84ee-7848945c9528
begin
	Pkg.activate(".")
		Pkg.Registry.add("General")  # only needed when installing Julia for the first time
		Pkg.Registry.add(RegistrySpec(url="https://github.com/ACEsuit/ACEregistry"))
		Pkg.add("ACEpotentials")
		using ACEpotentials
end

# ╔═╡ 7cf292f2-8872-4abc-9754-589e9a71d42c
function get_Pr(p, q, r0, rcut)
    basis = ACE1x.ace_basis(; elements = [:Si,],
	order = 3,
	totaldegree = 9,
	r0 = r0,
	rcut = rcut, 
    pair_transform = (:agnesi, p, q), 
    pair_envelope = (:r, 1)
	)
    Pr = basis.BB[1].J[1,1]
    env = Pr.envelope
    rr = range(0.01, rcut, length=400)
    vals = ACE1.evaluate.(Ref(Pr), rr)
    vals2 = [ [vals[j][i] / ACE1.evaluate(Pr.envelope, rr[j]) for j = 1:length(vals) ] 
              for i = 1:length(vals[1]) ]
    return rr, vals2 
end

# ╔═╡ f21bfb7d-5cbc-4848-8746-e30edc599f4d
begin
	r0 = 10.0; rcut = 17.0 
	rr, Pr = get_Pr(1, 4, r0, rcut)
	plt = plot()
	for n = 1:length(Pr) 
	    plot!(plt, rr, Pr[n], lw=2, label = "P_$n")
	end
	vline!([r0,], lw=2, ls = :dash, c= :black, label = "r0")
	plt
	#savefig("Rnl.png")
end

# ╔═╡ a03cd246-cb6a-4acc-9667-ebc972d03fa2


# ╔═╡ Cell order:
# ╠═1f79d0ea-710d-11ef-1d60-375671b1b654
# ╠═229beb96-de55-4fff-84ee-7848945c9528
# ╠═7cf292f2-8872-4abc-9754-589e9a71d42c
# ╠═f21bfb7d-5cbc-4848-8746-e30edc599f4d
# ╠═a03cd246-cb6a-4acc-9667-ebc972d03fa2
