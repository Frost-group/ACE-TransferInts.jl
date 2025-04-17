using Pkg
#cd("/Users/kk419/Documents/GitHub/ACE-TransferInts.jl/Notebook/2025_04/03")
cd("C:/Users/Keerati/Documents/GitHub/ACE-TransferInts.jl/Notebook/2025_04/13/Fitting/")

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

println(pwd())
 
begin
	# --- Load all frames ---
	frames = read_extxyz("all_CSi_modified.extxyz")

	# Total number of frames
	n_total = length(frames)

	# Number of frames to select for training
	n_train = 200

	# Uniformly spaced indices
	train_indices = round.(Int, range(1, stop=n_total, length=n_train))

	# Create training and testing sets
	train_data = frames[train_indices]
	test_data = [frames[i] for i in 1:n_total if i ∉ train_indices]
end


begin
	p = 1
	q = 4
	model1 = acemodel(elements = [:Si,:C,],
	        order = 2,
	        totaldegree = 10, 
			r0 = 12.0,
	        rcut = 13.6,
			#transform = (:agnesi, p, q),
			pair_transform = (:agnesi, p, q),
			pair_envelope = (:x, 2, 2),
			#envelope = (:x, 2, 2),
	        Eref = [:Si => 0.0, :C => 0.0]
	)
	@show length(model1.basis);
end

begin
	p = 19
	q = 43
	model2 = acemodel(elements = [:Si,:C,],
	        order = 2,
	        totaldegree = 10, #{2,3,4}-body 
			r0 = 12.0,
	        rcut = 13.6,
			#transform = (:agnesi, p, q),
			pair_transform = (:agnesi, p, q),
			pair_envelope = (:r, 2, 2),
			#envelope = (:x, 2, 2),
	        Eref = [:Si => 0.0, :C => 0.0]
	)
	@show length(model2.basis);
end


function get_Pr(p, q, r0, rcut, td)
    basis = ACE1x.ace_basis(;elements = [:Si,:C,],
	order = 2,
	totaldegree = td, #{2,3,4}-body 
	r0 = r0,
	rcut = rcut,
	#transform = (:agnesi, p, q),
	pair_transform = (:agnesi, p, q),
	pair_envelope = (:r, 2, 2),
	#envelope = (:x, 2, 2),
	)
    Pr = basis.BB[1].J[1,1]
    env = Pr.envelope
    rr = range(0.01, rcut, length=400)
    vals = ACE1.evaluate.(Ref(Pr), rr)
    vals2 = [ [vals[j][i] / ACE1.evaluate(Pr.envelope, rr[j]) for j = 1:length(vals) ] 
              for i = 1:length(vals[1]) ]
    return rr, vals2 
end

function get_Pr2(p, q, r0, rcut, td)
	basis = ACE1x.ace_basis(;elements = [:Si,:C,],
	order = 2,
	totaldegree = td, #{2,3,4}-body 
	r0 = r0,
	rcut = rcut,
	#transform = (:agnesi, p, q),
	pair_transform = (:agnesi, p, q),
	pair_envelope = (:x, 2, 2),
	#envelope = (:x, 2, 2),
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
	gr(dpi=1200)
	endg = 18.
	r_cut = 20.0; a = 0.2
	rdf_tiny = ACEpotentials.get_rdf(train_data, r_cut; rescale = true)
	plt_rdf_1 = plot()
	stephist!(rdf_tiny[(:Si,:Si)], 
		bins=50, label = "Intra-RDF",
		title="d-Si_trainset", titlefontsize=10,
		yticks = [0,25,50,75,100], xlims=(1,endg), ylims=(0.0, 100),
		size=(450,150), left_margin = 2Plots.mm, fill=true, fillalpha=0.5, color= "#F0BB62")
	#vline!(rnn(:Si)*Si_scale, label = L"r_{1}, r_{2}, ..., r_{5}", lw=1, color = "black", ls = :dash)
	
	b = 5
	rdf = ACEpotentials.get_rdf(train_data, r_cut; rescale = true);
	stephist!(rdf[(:Si,:C)],
		bins=50, label = "Inter-RDF",
		title="Training Set", titlefontsize=10,
		#xlabel = L"r[\AA]",
		ylabel = "RDF",
		yticks = [0*b, 25*b, 50*b, 75*b, 100*b], xlims=(1,endg), ylims=(0.0, 100*b),
		size=(300,200), left_margin = 2Plots.mm, fill=true, fillalpha=0.5, color= "#C64756")
	vline!(12.0*[1.0,], label = "Centre of Mass Separation", lw=1,
		 color = "black", ls = :dash, legend=(0.3, 0.9),size=(300,200))
	
	td = 10
	r0 = 12.0; rcut = 15.0
	rr, Pr = get_Pr(1, 4, r0, rcut, td)

	plt = plot()
	for n = 1:length(Pr)
	    plot!(plt, rr, Pr[n], lw=2, label = "", 
		 	#label = "P_$n", 
			xlabel = L"r[\AA]", ylabel = "RBF", legend=false, xlims=(0,endg),
			color= "#F0BB62")
	end
	vline!(r0*[1.0,], label = "r0", lw=1,
		color = "black", ls = :dash, legend=(0.3, 0.9))
	vline!(rcut*[1.0,], label = "rcut", lw=2,
		color = "#F0BB62", ls = :dash, legend=(0.3, 0.9),
		yticks = [0,25,50,75,100], xlims=(1,endg))	 

	r0 = 12.0; rcut = 14.5
	rr2, Pr2 = get_Pr2(11, 30, r0, rcut, td)	 
	plt1 = plot()
	for n = 1:length(Pr)
		plot!(plt1, rr2, Pr2[n], lw=2, label = "", 
			#label = "P_$n", 
			xlabel = L"r[\AA]", ylabel = "RBF", legend=false, xlims=(0,endg),
			color= "#C64756")	
	end
	vline!(r0*[1.0,], label = "r0", lw=1,
		 color = "black", ls = :dash, legend=(0.3, 0.9))
	vline!(rcut*[1.0,], label = "rcut", lw=2,
		 color = "#C64756", ls = :dash, legend=(0.3, 0.9),
		 yticks = [0,25,50,75,100], xlims=(1,endg), size=(450,150))

	#vline!([r0*s2,], lw=2, ls = :dash, c= "#F0BB62", label = L"r_{1}", legend=true)
	plot!(plt_rdf_1, 
	#plt_rdf_2,
	#plt, 
	plt1, layout=(2,1), size=(875,400))
	#savefig("_RDF_RBF_Si.png") 
	savefig("RDF_NH.png") 
end



#Bead

od = [2, 3, 4, 5]

RMSE_bead_od = [ 12.629,  10.732, 10.543, 10.295].*4 ./20
MAE_bead_od = [ 8.967,  7.435 , 7.237, 7.074].*4 ./20

model_s_od = [338, 1138, 2106, 2858 ]

RMSE_bead_td = [ 40.660,  12.629, 10.755, 10.793, 11.195, 14.471].*4 ./20
MAE_bead_td = [ 34.639,  8.967, 7.308, 7.407, 7.438, 8.304 ].*4 ./20

model_s_td = [65, 338, 987, 2178, 6850, 15688]

#NH

order = [2, 3, 4, 5]
RMSE_NH_od = [ 1.775,   1.632, 1.598, 1.579]
MAE_NH_od = [ 1.236,  1.074 , 1.033, 1.051]


td = [5, 10, 15, 20, 30, 40]

RMSE_NH_td = [ 3.883,  1.775, 1.663, 1.787, 2.339, 2.019]
MAE_NH_td = [ 2.895,  1.236, 1.186, 1.186, 1.526, 1.262]

begin
gr(dpi=1200)
plot(layout=(1,1))
plot!(size=(850,700))
plot!(od,
		[RMSE_bead_od MAE_bead_od RMSE_NH_od MAE_NH_od], 
		color = ["#F0BB62" "#6E9A50" "#F0BB62" "#6E9A50"],
		linestyle = [:solid :solid :dash :dash],
		labels = ["RMSE"*L"_{\mathbf{bead}}" "MAE"*L"_{\mathbf{bead}}" "RMSE"*L"_{\mathbf{NH}}" "MAE"*L"_{\mathbf{NH}}"],
		label="Line with Markers",
		marker=[:circle :utriangle],
		linewidth=3.5,
		xlabel = "Interaction Order "*L"(\nu)",
		ylims=(0.8, 2.8), ylabel = "RMSE & MAE meV/atom",
		left_margin = 4Plots.mm,
		title = "Fine-tuning with Interaction Order "*L"(\nu)" , 
		titlefontsize=20,
		legendfontsize=12,
		xguidefontsize=14,
		yguidefontsize=14,
)
savefig("od_error.png") 
end

begin
	gr(dpi=1200)
	plot(layout=(1,1))
	plot!(size=(850,700));
	plot!(td,
			[RMSE_bead_td MAE_bead_td RMSE_NH_td MAE_NH_td], 
			color = ["#F0BB62" "#6E9A50" "#F0BB62" "#6E9A50"],
			linestyle = [:solid :solid :dash :dash],
			labels = ["RMSE"*L"_{\mathbf{bead}}" "MAE"*L"_{\mathbf{bead}}" "RMSE"*L"_{\mathbf{NH}}" "MAE"*L"_{\mathbf{NH}}"],
			label="Line with Markers",
			marker=[:circle :utriangle],
			linewidth=3.5,
			xlabel = "Total Degrees of Polynomials" ,
			ylims=(0., 9), ylabel = "RMSE & MAE meV/atom",
			left_margin = 4Plots.mm,
			title="Fine-tuning with Total Degrees",
			titlefontsize=20,
			legendfontsize=12,
			xguidefontsize=14,
			yguidefontsize=14,
	)
savefig("td_error.png") 
end


begin
	gr(dpi=1200)
	plot(layout=(1,1));
	plot!(size=(850,700));
	plot!(od , model_s_od, 
			color = ["#F0BB62" ],
			linestyle = [:solid],
			labels = "",
			label="Line with Markers",
			marker=[:utriangle],
			linewidth=3.5,
			xlabel = "Interacting-or|der "*L"(\nu)" ,
			ylims=(100, 3000), ylabel = "Basis size", yaxis=:right,
			left_margin = 4Plots.mm,
			title="Basis size with Interacting-order "*L"(\nu)",
			titlefontsize=20,
			legendfontsize=12,
			xguidefontsize=14,
			yguidefontsize=14,
	)
savefig("od_basis.png") 
end

begin
	gr(dpi=1200)
	plot(layout=(1,1));
	plot!(size=(850,700));
	plot!(td, model_s_td, 
			color = ["#6E9A50"],
			linestyle = [:solid],
			labels = "",
			label="Line with Markers",
			marker=[:utriangle],
			linewidth=3.5,
			xlabel = "Total Degrees of Polynomials",
			#ylims=(0.8, 2.8),
			ylabel = "Basis size", yaxis=:right, yscale=:log10,
			left_margin = 4Plots.mm,
			title="Basis size with Total Degrees",
			titlefontsize=20,
			legendfontsize=12,
			xguidefontsize=14,
			yguidefontsize=14,
	)
savefig("td_basis.png") 
end
