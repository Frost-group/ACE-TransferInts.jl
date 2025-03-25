using Pkg
Pkg.activate(".")
Base.active_project()
Pkg.instantiate()
begin# add and load general packages used in this notebook.
	Pkg.add(["LaTeXStrings", "MultivariateStats", "Plots", "Suppressor"
         ])
	using LaTeXStrings, MultivariateStats, Plots, Printf, Statistics, Suppressor
end;
begin
	#Pkg.Registry.add("General")  # only needed when installing Julia for the first time
	Pkg.Registry.add(RegistrySpec(url="https://github.com/ACEsuit/ACEregistry"))
	Pkg.add(PackageSpec(name="ACEpotentials", version="0.6.7"))
end

using LaTeXStrings, MultivariateStats, Plots, Printf,
	Statistics, Suppressor, ACEpotentials

begin
    bead_tiny_dataset = read_extxyz("train_carbon_CSi3.extxyz")
    bead_dataset = read_extxyz("test_carbon_CSi3.extxyz")
end

function get_Pr(p, q, r0, rcut)
    basis = ACE1x.ace_basis(; elements = [:Si,:C,],
	order = 2,
	totaldegree = 10,
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

begin
	bead_tiny_dataset = read_extxyz("train_bead_SiO2.extxyz")
	bead_dataset = read_extxyz("test_bead_SiO2.extxyz")
end

begin
	r_cut = 20.0
	rdf_tiny = ACEpotentials.get_rdf(bead_tiny_dataset, r_cut; rescale = true)
	plt_rdf_1 = stephist(rdf_tiny[(:Si, :C)], 
		bins=150, label = "rdf",
		title="Bead_tiny_dataset", titlefontsize=10,
		xlabel = L"r[\AA]", ylabel = "RDF",
		yticks = [0,300,600,900,1200], xlims=(5,18), ylims=(0.0, 1200),
		size=(400,200), left_margin = 2Plots.mm, fill=true, fillalpha=0.5, color= "#F0BB62")
	#vline!(rnn(:Si)*[1.0, 4.65, 5.55], label = "r1, r2, ...", lw=2, color = "black")
	
	rdf = ACEpotentials.get_rdf(bead_dataset, r_cut; rescale = true);
	plt_rdf_2 = stephist(rdf[(:Si, :C)],
		bins=150, label = "rdf",
		title="Bead_dataset", titlefontsize=10,
		xlabel = L"r[\AA]", ylabel = "RDF",
		yticks = [0,900,1800,2700,3600], xlims=(5,18), ylims=(0.0, 3600),
		size=(400,200), left_margin = 2Plots.mm, fill=true, fillalpha=0.5, color= "#F0BB62")
	#vline!(rnn(:Si)*[1.0, 4.65, 5.55], label = "r1, r2, ...", lw=2, color = "black")
	
	r0 = 12.0; rcut = 17.5
	rr, Pr = get_Pr(11, 15, r0, rcut)
	plt1 = plot()
	for n = 1:length(Pr) 
	    plot!(plt1, rr, Pr[n], lw=2, label = "P_$n", 
			xlabel = L"r[\AA]", ylabel = "RBF", legend=false, xlims=(5,18))
	end
	vline!([r0,], lw=2, ls = :dash, c= :black, label = "r0")

	r0 = 12.0; rcut = 17.5
	rr, Pr = get_Pr(11, 30, r0, rcut)
	plt2 = plot()
	for n = 1:length(Pr) 
	    plot!(plt2, rr, Pr[n], lw=2, label = "P_$n", 
			xlabel = L"r[\AA]", ylabel = "RBF", legend=false, xlims=(5,18))
	end
	vline!([r0,], lw=2, ls = :dash, c= :black, label = "r0")
	
    r0 = 12.0; rcut = 14.5
	rr, Pr = get_Pr(11, 30, r0, rcut)
	plt3 = plot()
	for n = 1:length(Pr) 
	    plot!(plt3, rr, Pr[n], lw=2, label = "P_$n", 
			xlabel = L"r[\AA]", ylabel = "RBF", legend=false, xlims=(5,18))
	end
	vline!([r0,], lw=2, ls = :dash, c= :black, label = "r0")

	plot!(plt_rdf_1, plt_rdf_2, plt1, plt2, plt3, layout=(5,1), size=(1500,900))
	savefig("_RDF_RBF_carbon.png") 
	#savefig("RDF_bead.png") 
end


function get_Pr(p, q, r0, rcut)
    basis = ACE1x.ace_basis(; elements = [:Si,:O,],
	order = 2,
	totaldegree = 10,
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


begin
	r_cut = 20.0
	rdf_tiny = ACEpotentials.get_rdf(bead_tiny_dataset, r_cut; rescale = true)
	plt_rdf_1 = stephist(rdf_tiny[(:Si, :O)], 
		bins=150, label = "rdf",
		title="Bead_tiny_dataset", titlefontsize=10,
		xlabel = L"r[\AA]", ylabel = "RDF",
		yticks = [0,25,50,75,100], xlims=(8,16), ylims=(0.0, 100),
		size=(400,200), left_margin = 2Plots.mm, fill=true, fillalpha=0.5, color= "#F0BB62")
	#vline!(rnn(:Si)*[1.0, 4.65, 5.55], label = "r1, r2, ...", lw=2, color = "black")
	
	rdf = ACEpotentials.get_rdf(bead_dataset, r_cut; rescale = true);
	plt_rdf_2 = stephist(rdf[(:Si, :O)],
		bins=150, label = "rdf",
		title="Bead_dataset", titlefontsize=10,
		xlabel = L"r[\AA]", ylabel = "RDF",
		yticks = [0,100,200,300,400], xlims=(8,16), ylims=(0.0, 400),
		size=(400,200), left_margin = 2Plots.mm, fill=true, fillalpha=0.5, color= "#F0BB62")
	#vline!(rnn(:Si)*[1.0, 4.65, 5.55], label = "r1, r2, ...", lw=2, color = "black")
	
	r0 = 12.0; rcut = 14.5
	rr, Pr = get_Pr(11, 30, r0, rcut)
	plt = plot()
	for n = 1:length(Pr) 
	    plot!(plt, rr, Pr[n], lw=2, label = "P_$n", 
			xlabel = L"r[\AA]", ylabel = "RBF", legend=false, xlims=(8,16))
	end
	vline!([r0,], lw=2, ls = :dash, c= :black, label = "r0")

	r0 = 12.0; rcut = 13.5
	rr, Pr = get_Pr(16, 35, r0, rcut)
	plt1 = plot()
	for n = 1:length(Pr) 
	    plot!(plt1, rr, Pr[n], lw=2, label = "P_$n", 
			xlabel = L"r[\AA]", ylabel = "RBF", legend=false, xlims=(8,16))
	end
	vline!([r0,], lw=2, ls = :dash, c= :black, label = "r0")
	
	r0 = 12.0; rcut = 13.3
	rr, Pr = get_Pr(26, 45, r0, rcut)
	plt2 = plot()
	for n = 1:length(Pr) 
	    plot!(plt2, rr, Pr[n], lw=2, label = "P_$n", 
			xlabel = L"r[\AA]", ylabel = "RBF", legend=false, xlims=(8,16))
	end
	vline!([r0,], lw=2, ls = :dash, c= :black, label = "r0")
	
	plot!(plt_rdf_1, plt_rdf_2, plt, plt1, plt2, layout=(5,1), size=(1500,900))
	savefig("_RDF_RBF_bead1.png") 
	#savefig("RDF_bead.png") 
end