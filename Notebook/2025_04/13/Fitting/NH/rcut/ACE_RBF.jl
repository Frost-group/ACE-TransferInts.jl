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
	# --- Load all frames ---
	frames = read_extxyz("all_bead_modified.extxyz")

	# Total number of frames
	n_total = length(frames)

	# Number of frames to select for training
	n_train = 900

	# Uniformly spaced indices
	train_indices = round.(Int, range(1, stop=n_total, length=n_train))

	# Create training and testing sets
	train_data1 = frames[train_indices]
	bead_dataset = [frames[i] for i in 1:n_total if i ∉ train_indices]

end


function get_Pr(p, q, r0, rcut, td)
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
	b = 2
	endg = 18.
	r_cut = 20.0; a = 0.2
	rdf_tiny = ACEpotentials.get_rdf(train_data1, r_cut; rescale = true)
	plt_rdf_1 = plot()
	stephist!(rdf_tiny[(:Si,:C)], 
		bins=50, label = "Intra-RDF",
		title="Bead Training Set", titlefontsize=10,
		ylabel = "RDF",
		yticks = [0*b, 25*b, 50*b, 75*b, 100*b], xlims=(1,endg), ylims=(0.0, 100*b),
		size=(450,150), left_margin = 4Plots.mm, fill=true, fillalpha=0.5, color= "#C64756")
		vline!(12.0*[1.0,], label = "Centre of Mass Separation", lw=1,
		color = "black", ls = :dash, legend=(0.3, 0.9),size=(300,200),
		titlefontsize=14)
		
	b = 5
	rdf = ACEpotentials.get_rdf(train_data, r_cut; rescale = true);
	plt_rdf_2 = plot()
	stephist!(rdf[(:Si,:C)],
		bins=50, label = "Inter-RDF",
		title="Heavy-atom Training Set", titlefontsize=10,
		#xlabel = L"r[\AA]",
		ylabel = "RDF",
		yticks = [0*b, 25*b, 50*b, 75*b, 100*b], xlims=(1,endg), ylims=(0.0, 100*b),
		size=(300,200), left_margin = 6Plots.mm, fill=true, fillalpha=0.5, color= "#C64756")
	vline!(12.0*[1.0,], label = "Centre of Mass Separation", lw=1,
		 color = "black", ls = :dash, legend=(0.3, 0.9),size=(300,200),
		 titlefontsize=14)
	
	td = 10
	r0 = 12.0; rcut = 13.6
	rr, Pr = get_Pr(11, 43, r0, rcut, td)

	plt = plot()
	for n = 1:length(Pr)
	    plot!(plt, rr, Pr[n], lw=2, label = "", 
		 	#label = "P_$n", 
			xlabel = L"r[\AA]", ylabel = "RBF", legend=false, xlims=(0,endg),
			color= "#C64756")
	end
	vline!(r0*[1.0,], label = "r0", lw=1,
		color = "black", ls = :dash, legend=(0.3, 0.9))
	vline!(rcut*[1.0,], label = "rcut", lw=2,
		color = "#C64756",
		ls = :dash,
		legend=(0.3, 0.9),
		yticks = [0,1,2,-1,-2],
		xlims=(1,endg), 
		title = "RBF with "*L"p=11",
		titlefontsize=14
		)
		#F0BB62

	r0 = 12.0; rcut = 13.6
	rr2, Pr2 = get_Pr(25, 43, r0, rcut, td)	 
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
		 color = "#C64756",
		 ls = :dash, 
		 legend=(0.3, 0.9),
		 yticks = [0,1,2,-1,-2], 
		 xlims=(1,endg), 
		 title = "RBF with "*L"p=25",
		 titlefontsize=14
		 )

	#vline!([r0*s2,], lw=2, ls = :dash, c= "#F0BB62", label = L"r_{1}", legend=true)
	plot!(plt_rdf_1, 
	plt_rdf_2,
	plt, 
	plt1, layout=(4,1), size=(1200,800))
	#savefig("_RDF_RBF_Si.png") 
	savefig("_RDF_p.png") 
end



#Bead

od = [2, 3, 4, 5]

RMSE_bead_od = [ 12.629,  10.732, 10.543, 10.295].*4 ./20
MAE_bead_od = [ 8.967,  7.435 , 7.237, 7.074].*4 ./20

model_s_od = [338, 1138, 2106, 2858 ]

RMSE_bead_td = [ 40.660,  12.629, 10.755, 10.793, 11.195, 14.471].*4 ./20
MAE_bead_td = [ 34.639,  8.967, 7.308, 7.407, 7.438, 8.304 ].*4 ./20

RMSE_bead_td = [ 40.660,  12.629, 10.755, 10.793, 11.195, 14.471].*4 ./20
MAE_bead_td = [ 34.639,  8.967, 7.308, 7.407, 7.438, 8.304 ].*4 ./20

#NH
td = [5, 10, 15, 20, 30, 40]
RMSE_NH_od = [ 1.775,  1.626, 1.572, 1.552]
MAE_NH_od = [ 1.236,  1.093 , 1.033 , 1.053]

RMSE_NH_td = [ 4.360,  1.775, 1.645, 1.700, 2.282, 1.775]
MAE_NH_td = [ 2.887,  1.24, 1.236, 1.241, 1.49, 1.236]


begin
	gr(dpi=1200)
	plot(layout=(1,1))
	plot!(size=(850,700))
	plot!(od,
			[RMSE_bead_od MAE_bead_od RMSE_NH_od MAE_NH_od], 
			color = ["#F0BB62" "#F0BB62" "#6E9A50" "#6E9A50" ],
			linestyle = [:solid :dashdot :solid  :dashdot],
			labels = ["RMSE"*L"_{\mathbf{Bead}}" "MAE"*L"_{\mathbf{Bead}}" "RMSE"*L"_{\mathbf{HA}}" "MAE"*L"_{\mathbf{HA}}"],
			label="Line with Markers",
			marker=[:circle :utriangle :circle :utriangle],
			linewidth=[3.5 3.5 3.5 3.5],
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
				color = ["#F0BB62" "#F0BB62" "#6E9A50" "#6E9A50" ],
				linestyle = [:solid :dashdot :solid  :dashdot],
				labels = ["RMSE"*L"_{\mathbf{Bead}}" "MAE"*L"_{\mathbf{Bead}}" "RMSE"*L"_{\mathbf{HA}}" "MAE"*L"_{\mathbf{HA}}"],
				label="Line with Markers",
				marker=[:circle :utriangle :circle :utriangle],
				linewidth=[3.5 3.5 3.5 3.5],
				xlabel = "Total Degrees of Polynomials" ,
				xlims=(0., 45), ylims=(0., 9), ylabel = "RMSE & MAE meV/atom",
				left_margin = 4Plots.mm,
				title="Fine-tuning with Total Degrees",
				titlefontsize=20,
				legendfontsize=12,
				xguidefontsize=14,
				yguidefontsize=14,
		)
	savefig("td_error.png") 
end

rcut_bead = [ 13.0, 13.30, 13.45, 13.60, 13.90, 14.20]

RMSE_bead_rc = ([ 40.997, 13.253, 13.760, 12.629, 12.563, 40.037].-12.563).*4 ./20
MAE_bead_rc = ([ 35.100, 8.819, 9.694, 8.967, 9.095, 33.681].-8.819).*4 ./20

rcut_NH = [ 9.5, 11.0, 12.1, 12.4, 12.7, 13.0, 13.30, 13.45, 13.60, 13.90, 14.20, 14.50, 14.80, 16., 17.5]
RMSE_NH_rc = ([3.669 ,2.085, 1.926, 1.871,1.779, 1.937, 1.776, 1.881, 1.775, 1.777, 1.754, 1.694, 1.724, 1.932, 1.938].-1.754)
MAE_NH_rc = ([2.690 ,1.487, 1.332, 1.309,1.227, 1.347, 1.237, 1.354, 1.236, 1.272, 1.270, 1.192, 1.236, 1.368, 1.356].-1.270)
rcut_NH = [ 13.0, 13.30, 13.45, 13.60, 13.90, 14.20, ]
RMSE_NH_rc = ([ 1.937, 1.776, 1.881, 1.775, 1.777, 1.754,].-1.754)
MAE_NH_rc = ([ 1.347, 1.237, 1.354, 1.236, 1.272, 1.270, ].-1.270)

begin
	gr(dpi=1200)
	plot(layout=(1,1))
	plot!(size=(850,700));
	plot!([rcut_bead ],
			[RMSE_bead_rc MAE_bead_rc], 
			color = ["#F0BB62" "#F0BB62"],
			linestyle = [:solid :dashdot],
			labels = [L"\Delta"*"RMSE"*L"_{\mathbf{Bead}}" L"\Delta"*"MAE"*L"_{\mathbf{Bead}}"],
			label="Line with Markers",
			marker=[:circle :utriangle],
			linewidth=[3.5 3.5 3.5 3.5],
			xlabel = "Cutoff Radias "*L"\AA" ,
			xlims=(12.85, 14.35), xticks = [13.0, 13.30, 13.60, 13.90, 14.20],
			ylims=(-0.2, 6.2), ylabel = L"\Delta"*"RMSE & "*L"\Delta"*"MAE meV/atom",
			left_margin = 4Plots.mm,
			title="Fine-tuning with Cutoff Radias",
			titlefontsize=20,
			legendfontsize=12,
			xguidefontsize=14,
			yguidefontsize=14,
	)
	plot!(rcut_NH ,
	[RMSE_NH_rc MAE_NH_rc], 
	color = ["#6E9A50" "#6E9A50" ],
	linestyle = [ :solid  :dashdot],
	labels = [L"\Delta"*"RMSE"*L"_{\mathbf{NH}}" L"\Delta"*"MAE"*L"_{\mathbf{NH}}"],
	label="Line with Markers",
	marker=[:circle :utriangle],
	linewidth=[3.5 3.5 3.5 3.5],
	xlabel = "Cutoff Radias "*L"\AA" ,
	xlims=(12.85, 14.35), xticks = [13.0, 13.30, 13.60, 13.90, 14.20],
	ylims=(-0.2, 6.2), ylabel = L"\Delta"*"RMSE & "*L"\Delta"*"MAE meV/atom",
	left_margin = 4Plots.mm,
	title="Fine-tuning with Cutoff Radias",
	titlefontsize=20,
	legendfontsize=12,
	xguidefontsize=14,
	yguidefontsize=14, legend=(0.5, 0.9)
)

savefig("rc_error.png") 
end



p = [9, 11, 13, 15, 17, 19, 21, 23, 25]

RMSE_bead_p = ([12.598, 12.445, 12.668, 12.650, 12.632, 12.629, 12.660, 12.707, 12.752].-12.445).*4 ./20
MAE_bead_p = ([8.890, 8.772, 8.994, 9.046, 9.020, 8.967, 8.957, 8.981, 9.011].-8.772).*4 ./20

RMSE_NH_p = ([ 1.768, 1.755, 1.741, 1.756, 1.773, 1.775, 1.773, 1.774, 1.778].-1.741)
MAE_NH_p = ([ 1.224, 1.216, 1.207, 1.215, 1.231, 1.236, 1.237, 1.236, 1.234].-1.207)

begin
	gr(dpi=1200)
	plot(layout=(1,1))
	plot!(size=(850,700));
	plot!(p,
			[RMSE_bead_p MAE_bead_p RMSE_NH_p MAE_NH_p], 
			color = ["#F0BB62" "#F0BB62" "#6E9A50" "#6E9A50" ],
			linestyle = [:solid :dashdot :solid  :dashdot],
			labels = [L"\Delta"*"RMSE"*L"_{\mathbf{Bead}}" L"\Delta"*"MAE"*L"_{\mathbf{Bead}}" L"\Delta"*"RMSE"*L"_{\mathbf{NH}}" L"\Delta"*"MAE"*L"_{\mathbf{NH}}"],
			label="Line with Markers",
			marker=[:circle :utriangle :circle :utriangle],
			linewidth=[3.5 3.5 3.5 3.5],
			xlabel = L"p"*" parameter" ,
			xlims=(7, 27), xticks = [9, 11, 13, 15, 17, 19, 21, 23, 25],
			ylims=(-0.02, 0.08), ylabel = L"\Delta"*"RMSE & "*L"\Delta"*"MAE meV/atom",
			left_margin = 4Plots.mm,
			title="Fine-tuning with "*L"p"*" parameter" ,
			titlefontsize=20,
			legendfontsize=12,
			xguidefontsize=14,
			yguidefontsize=14,
	)
savefig("p_error.png") 
end


q = [37, 39, 41, 43, 45, 47]
RMSE_bead_q = ([40.329 ,40.162, 12.514, 12.629, 41.226, 41.358 ].-12.514).*4 ./20
MAE_bead_q = ([ 34.205,33.940, 8.893, 8.967, 35.387, 35.556].-8.893).*4 ./20

RMSE_NH_q = ([1.778, 1.777, 1.775, 1.775, 1.776, 1.777].-1.775)
MAE_NH_q = ([ 1.250, 1.247, 1.241, 1.236, 1.233, 1.231].-1.231)

begin
	gr(dpi=1200)
	plot(layout=(1,1))
	plot!(size=(850,700));
	plot!(q,
			[RMSE_bead_q MAE_bead_q RMSE_NH_q MAE_NH_q], 
			color = ["#F0BB62" "#F0BB62" "#6E9A50" "#6E9A50" ],
			linestyle = [:solid :dashdot :solid  :dashdot],
			labels = [L"\Delta"*"RMSE"*L"_{\mathbf{Bead}}" L"\Delta"*"MAE"*L"_{\mathbf{Bead}}" L"\Delta"*"RMSE"*L"_{\mathbf{NH}}" L"\Delta"*"MAE"*L"_{\mathbf{NH}}"],
			label="Line with Markers",
			marker=[:circle :utriangle :circle :utriangle],
			linewidth=[3.5 3.5 3.5 3.5],
			xlabel = L"q"*" parameter" ,
			xlims=(35, 48), xticks = [35, 37, 39, 41, 43, 45, 47, 49],
			ylims=(-0.1, 6.2), ylabel = L"\Delta"*"RMSE & "*L"\Delta"*"MAE meV/atom",
			left_margin = 4Plots.mm,
			title="Fine-tuning with "*L"q"*" parameter" ,
			titlefontsize=20,
			legendfontsize=12,
			xguidefontsize=14,
			yguidefontsize=14, legend=(0.5, 0.9)
	)
savefig("q_error.png") 
end