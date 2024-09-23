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
		#Pkg.Registry.add("General")  # only needed when installing Julia for the first time
		Pkg.Registry.add(RegistrySpec(url="https://github.com/ACEsuit/ACEregistry"))
		Pkg.add("ACEpotentials")
		using ACEpotentials
end

# ╔═╡ 7cf292f2-8872-4abc-9754-589e9a71d42c
function get_Pr(p, q, r0, rcut)
    basis = ACE1x.ace_basis(; elements = [:Si,],
	order = 3,
	totaldegree = 14,
	r0 = r0,
	rcut = rcut, 
    pair_transform = (:agnesi, p, q),
    pair_envelope = (:x, 2, 2)
	)
    Pr = basis.BB[1].J[1,1]
    env = Pr.envelope
    rr = range(0.01, rcut, length=400)
    vals = ACE1.evaluate.(Ref(Pr), rr)
    vals2 = [ [vals[j][i] / ACE1.evaluate(Pr.envelope, rr[j]) for j = 1:length(vals) ] 
              for i = 1:length(vals[1]) ]
    return rr, vals2 
end

# ╔═╡ a03cd246-cb6a-4acc-9667-ebc972d03fa2
begin
	bead_tiny_dataset = read_extxyz("train_bead.extxyz")
	bead_dataset = read_extxyz("test_bead.extxyz")
end

# ╔═╡ 0bb80fdb-7de9-404a-843d-633f2c9f52d5
begin
	r_cut = 20.0
	rdf_tiny = ACEpotentials.get_rdf(bead_tiny_dataset, r_cut; rescale = true)
	plt_rdf_1 = stephist(rdf_tiny[(:Si, :Si)], 
		bins=150, label = "rdf",
		title="Bead_tiny_dataset", titlefontsize=10,
		#xlabel = L"r[\AA]",
		ylabel = "RDF",
		yticks = [0,50,100], xlims=(1.5,15), ylims=(0.0, 100),
		size=(400,200), left_margin = 2Plots.mm, fill=true, fillalpha=0.5, color= "#F0BB62")
	#vline!(rnn(:Si)*[1.0, 4.65, 5.55], label = "r1, r2, ...", lw=2, color = "black")
	
	rdf = ACEpotentials.get_rdf(bead_dataset, r_cut; rescale = true);
	plt_rdf_2 = stephist(rdf[(:Si, :Si)],
		bins=150, label = "rdf",
		title="Bead_dataset", titlefontsize=10,
		#xlabel = L"r[\AA]",
		ylabel = "RDF",
		yticks = [0,50,100], xlims=(1.5,15), ylims=(0.0, 100),
		size=(400,200), left_margin = 2Plots.mm, fill=true, fillalpha=0.5, color= "#F0BB62")
	#vline!(rnn(:Si)*[1.0, 4.65, 5.55], label = "r1, r2, ...", lw=2, color = "black")
	r0 = 12.0; rcut = 14.5
	rr, Pr = get_Pr(11, 30, r0, rcut)
	plt = plot()
	for n = 1:length(Pr) 
	    plot!(plt, rr, Pr[n], lw=2, label = "P_$n", 
			xlabel = L"r[\AA]", ylabel = "RBF", legend=false, xlims=(1.5,15))
	end
	vline!([r0,], lw=2, ls = :dash, c= :black, label = "r0")
	plt
	
	plot!(plt_rdf_1, plt_rdf_2, plt, layout=(3,1), size=(900,600))
	#savefig("old_RDF_RBF_bead.png") 
end

# ╔═╡ cd4c712a-d07b-4705-8818-508a7ee7bbca
begin
	carbon_tiny_dataset = read_extxyz("train_d_Car.extxyz")
	carbon_dataset = read_extxyz("test_d_Car.extxyz")
end

# ╔═╡ 255b1dec-9293-49d9-9e57-3dc7678ab60d
begin
	r_cut2 = 20.0
	rdf_tiny2 = ACEpotentials.get_rdf(carbon_tiny_dataset, r_cut2; rescale = true)
	plt_rdf_3 = stephist(rdf_tiny2[(:C, :C)], 
		bins=150, label = "rdf",
		title="Carbon_tiny_dataset", titlefontsize=10,
		#xlabel = L"r[\AA]",
		ylabel = "RDF",
		yticks = [0,50,100], xlims=(1.5,17), ylims=(0.0, 200),
		size=(400,200), left_margin = 2Plots.mm, fill=true, fillalpha=0.5, color= "#F0BB62")
	#vline!(rnn(:Si)*[1.0, 4.65, 5.55], label = "r1, r2, ...", lw=2, color = "black")
	
	rdf2 = ACEpotentials.get_rdf(carbon_dataset, r_cut2; rescale = true);
	plt_rdf_4 = stephist(rdf2[(:C, :C)],
		bins=150, label = "rdf",
		title="Carbon_dataset", titlefontsize=10,
		#xlabel = L"r[\AA]",
		ylabel = "RDF",
		yticks = [0,50,100], xlims=(1.5,17), ylims=(0.0, 200),
		size=(400,200), left_margin = 2Plots.mm, fill=true, fillalpha=0.5, color= "#F0BB62")
	#vline!(rnn(:Si)*[1.0, 4.65, 5.55], label = "r1, r2, ...", lw=2, color = "black")
	r01 = 10.0; rcut1 = 16.5
	rr1, Pr1 = get_Pr(11, 30, r01, rcut1)
	plt2 = plot()
	for n = 1:length(Pr)
	    plot!(plt2, rr1, Pr1[n], lw=2, label = "P_$n",
			xlabel = L"r[\AA]", ylabel = "RBF", legend=false, xlims=(1.5,17))
	end
	vline!([r01,], lw=2, ls = :dash, c= :black, label = "r0")
	plt2
	
	plot!(plt_rdf_3, plt_rdf_4, plt2, layout=(3,1), size=(900,600))
	savefig("old_RDF_RBF_Carbon.png") 
end

# ╔═╡ Cell order:
# ╠═1f79d0ea-710d-11ef-1d60-375671b1b654
# ╠═229beb96-de55-4fff-84ee-7848945c9528
# ╠═7cf292f2-8872-4abc-9754-589e9a71d42c
# ╠═a03cd246-cb6a-4acc-9667-ebc972d03fa2
# ╠═0bb80fdb-7de9-404a-843d-633f2c9f52d5
# ╠═cd4c712a-d07b-4705-8818-508a7ee7bbca
# ╠═255b1dec-9293-49d9-9e57-3dc7678ab60d
