using Pkg
#cd("/Users/kk419/Documents/GitHub/ACE-TransferInts.jl/Notebook/2025_04/03")
cd("C:/Users/Keerati/Documents/GitHub/ACE-TransferInts.jl/Notebook/2025_04/10/Fitting/")

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
	frames = read_extxyz("all_CS_SiO_modified.extxyz")
	frames = read_extxyz("all_CSi.extxyz")
	frames = read_extxyz("all_CSi_modified.extxyz")
	# Total number of frames
	n_total = length(frames)

	# Number of frames to select for training
	n_train = 900

	# Uniformly spaced indices
	train_indices = round.(Int, range(1, stop=n_total, length=n_train))

	# Create training and testing sets
	train_data = frames[train_indices]
	test_data = [frames[i] for i in 1:n_total if i ∉ train_indices]
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
ACEpotentials.get_rdf(train_data, r_cut; rescale = true)

begin
	gr(dpi=1200)
	endg = 8.
	r_cut = 7.6; a = 0.2
	rdf_tiny = ACEpotentials.get_rdf(train_data, r_cut; rescale = true)
	plt_rdf_1 = plot()
	stephist!(rdf_tiny[(:C,:C)], 
		bins=150, label = "Intra-RDF" ,
		title="d-Si_trainset", titlefontsize=10,
		yticks = [0,25,50,75,100], xlims=(1,endg), ylims=(0.0, 100),
		size=(450,150), left_margin = 2Plots.mm, fill=true, fillalpha=0.5, color= "#F0BB62")
	#vline!(rnn(:Si)*Si_scale, label = L"r_{1}, r_{2}, ..., r_{5}", lw=1, color = "black", ls = :dash)
	
	b = 8
	stephist!(rdf_tiny[(:Si,:C)],
		bins=50, label = "Inter-RDF",
		title="Training Set", titlefontsize=10,
		xlabel = L"r[\AA]",
		ylabel = "RDF",
		yticks = [0*b, 25*b, 50*b, 75*b, 100*b], xlims=(1,endg), ylims=(0.0, 100*b),
		size=(300,200), left_margin = 4Plots.mm, fill=true, fillalpha=0.5, color= "#C64756")
	vline!(6.0*[1.0,], label = "Centre of Mass Separation", lw=1,
		color = "black", ls = :dash, top_margin = 4Plots.mm,
		legend=(0.3, 0.9),size=(300,200))
	
	td = 10
	r0 = 6.0; rcut = 15.0
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

	r0 = 6.0; rcut = 7.6
	rr2, Pr2 = get_Pr2(19, 43, r0, rcut, td)	 
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
		 color = "#C64756", ls = :dash, legend=(0.3, 0.9),  bottom_margin = 4Plots.mm,
		 yticks = [0,25,50,75,100], xlims=(1,endg), size=(450,150)
		 )

	#vline!([r0*s2,], lw=2, ls = :dash, c= "#F0BB62", label = L"r_{1}", legend=true)
	plot!(plt_rdf_1, 
	#plt_rdf_2,
	#plt, 
	plt1, layout=(2,1), size=(875,425))
	#savefig("_RDF_RBF_Si.png") 
	savefig("RDF_NH.png")
end



#Bead 1

od = [2, 3, 4, 5]
td = [5, 10, 15, 20, 30, 40]
model_s_od = [338, 1138, 2106, 2858]
model_s_od_2 = [2484, 14180, 38616, 67544]

model_s_td = [65, 338, 987, 2178, 6850, 15688]
model_s_td_2 = [426, 2484, 7518, 16848, 53740, 123824]

RMSE_bead_od = [418.176,  418.176, 10.543, 10.295].*2 ./10
MAE_bead_od = [254.612,  254.612 , 7.237, 7.074].*2 ./10

RMSE_bead_td = [ 40.660,  13.636, 10.755, 10.793, 11.195, 14.471].*2 ./10
MAE_bead_td = [ 34.639,  11.061, 7.308, 7.407, 7.438, 8.304 ].*2 ./10


#NH 1

RMSE_NH_od_1 = [13.636, 13.616, 13.743, 13.762].*4 ./5
MAE_NH_od_1 = [11.061, 11.024, 11.137, 11.132].*4 ./5

RMSE_NH_td_1 = [ 41.113, 13.636, 13.733, 13.953, 13.133, 13.067].*4 ./5
MAE_NH_td_1 = [ 29.729, 11.061, 11.095, 11.100, 10.095, 10.482].*4 ./5

#NH 2

RMSE_NH_od_2 = [ 11.274,  10.758, 10.902, 10.936]
MAE_NH_od_2 = [ 8.864, 8.487, 8.630, 8.652]

RMSE_NH_td_2 = [ 13.630,  11.274, 11.078, 10.871, 10.545, 10.446]
MAE_NH_td_2 = [ 10.508,  8.864, 8.681, 8.527, 8.300, 8.170]

begin
gr(dpi=1200)
plot(layout=(1,1))
plot!(size=(850,700))
plot!(od,
		[RMSE_NH_od_1 MAE_NH_od_1 RMSE_NH_od_2 MAE_NH_od_2], 
		color = ["#F0BB62" "#F0BB62" "#6E9A50" "#6E9A50" ],
			linestyle = [:solid :dashdot :solid  :dashdot],
			labels = ["RMSE"*L"_{\mathbf{HA-C}}" "MAE"*L"_{\mathbf{HA-C}}" "RMSE"*L"_{\mathbf{HA-CS}}" "MAE"*L"_{\mathbf{HA-CS}}"],
			label="Line with Markers",
			marker=[:circle :utriangle :circle :utriangle],
			linewidth=[3.5 3.5 3.5 3.5],
		xlabel = "Interaction Order "*L"(\nu)",
		ylims=(8, 12), ylabel = "RMSE & MAE meV/atom",
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
			[RMSE_NH_td_1 MAE_NH_td_1 RMSE_NH_td_2 MAE_NH_td_2], 
			color = ["#F0BB62" "#F0BB62" "#6E9A50" "#6E9A50" ],
			linestyle = [:solid :dashdot :solid  :dashdot],
			labels = ["RMSE"*L"_{\mathbf{HA-C}}" "MAE"*L"_{\mathbf{HA-C}}" "RMSE"*L"_{\mathbf{HA-CS}}" "MAE"*L"_{\mathbf{HA-CS}}"],
			label="Line with Markers",
			marker=[:circle :utriangle :circle :utriangle],
			linewidth=[3.5 3.5 3.5 3.5],
			xlabel = "Total Degrees of Polynomials" ,
			xlims=(0., 45), ylims=(7., 35), ylabel = "RMSE & MAE meV/atom",
			left_margin = 4Plots.mm,
			title="Fine-tuning with Total Degrees",
			titlefontsize=20,
			legendfontsize=12,
			xguidefontsize=14,
			yguidefontsize=14,
	)
savefig("td_error.png") 
end

