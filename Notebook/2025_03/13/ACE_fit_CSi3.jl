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

begin
	p = 11
	q = 30
	model1 = acemodel(elements = [:Si, :C],
			order = 5,
			totaldegree = 10,
			r0 = 12.0,
	        rcut = 17.5,
			#transform = (:agnesi, p, q),
			pair_transform = (:agnesi, p, q),
			pair_envelope = (:x, 2, 2),
			#envelope = (:x, 2, 2),
	        Eref = [:Si => 0.0, :C => 0.0]
	)
	@show length(model1.basis);
end

begin
	p = 1
	q = 4
	model2 = acemodel(elements = [:Si,:C,],
	        order = 3,
	        totaldegree = [25, 20, 15], #{2,3,4}-body 
			r0 = 12.0,
	        rcut = 16.0,
			#transform = (:agnesi, p, q),
			pair_transform = (:agnesi, p, q),
			pair_envelope = (:x, 2, 2),
			#envelope = (:x, 2, 2),
	        Eref = [:Si => 0.0, :C => 0.0]
	)
	@show length(model2.basis);
end

function get_Pr(p, q, r0, rcut)
    basis = ACE1x.ace_basis(; elements = [:Si,:C,],
	order = 4,
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
	plt_rdf_1 = stephist(rdf_tiny[(:Si, :C)], 
		bins=150, label = "rdf",
		title="Bead_tiny_dataset", titlefontsize=10,
		xlabel = L"r[\AA]", ylabel = "RDF",
		yticks = [0,25,50,75,100], xlims=(1.5,20), ylims=(0.0, 100),
		size=(400,200), left_margin = 2Plots.mm, fill=true, fillalpha=0.5, color= "#F0BB62")
	#vline!(rnn(:Si)*[1.0, 4.65, 5.55], label = "r1, r2, ...", lw=2, color = "black")
	
	rdf = ACEpotentials.get_rdf(bead_dataset, r_cut; rescale = true);
	plt_rdf_2 = stephist(rdf[(:Si, :C)],
		bins=150, label = "rdf",
		title="Bead_dataset", titlefontsize=10,
		xlabel = L"r[\AA]", ylabel = "RDF",
		yticks = [0,100,200,300,400], xlims=(1.5,20), ylims=(0.0, 400),
		size=(400,200), left_margin = 2Plots.mm, fill=true, fillalpha=0.5, color= "#F0BB62")
	#vline!(rnn(:Si)*[1.0, 4.65, 5.55], label = "r1, r2, ...", lw=2, color = "black")
	
	r0 = 12.0; rcut = 14.5
	rr, Pr = get_Pr(11, 30, r0, rcut)
	plt = plot()
	for n = 1:length(Pr) 
	    plot!(plt, rr, Pr[n], lw=2, label = "P_$n", 
			xlabel = L"r[\AA]", ylabel = "RBF", legend=false, xlims=(1.5,20))
	end
	vline!([r0,], lw=2, ls = :dash, c= :black, label = "r0")

	r0 = 12.0; rcut = 18.0
	rr, Pr = get_Pr(1, 4, r0, rcut)
	plt1 = plot()
	for n = 1:length(Pr) 
	    plot!(plt1, rr, Pr[n], lw=2, label = "P_$n", 
			xlabel = L"r[\AA]", ylabel = "RBF", legend=false, xlims=(1.5,20))
	end
	vline!([r0,], lw=2, ls = :dash, c= :black, label = "r0")
	
	
	plot!(plt_rdf_1, plt_rdf_2, plt, plt1, layout=(4,1), size=(1200,600))
	savefig("_RDF_RBF_carbon.png") 
	#savefig("RDF_bead.png") 
end


data_keys1 = (energy_key = "energy", force_key = "", virial_key = "")
solver = ACEfit.BLR(committee_size=500, factorization=:svd)
acefit!(model1, bead_tiny_dataset; solver=solver, data_keys1...);


begin
	@info("Training Errors")
	ACEpotentials.linear_errors(bead_tiny_dataset, model1; data_keys1...);
	
	@info("Test Error")
	ACEpotentials.linear_errors(bead_dataset, model1; data_keys1...);
end

begin
	function extract_energies(dataset)
	    energies = []
	    for atoms in dataset
	        for key in keys(atoms.data)
	            if lowercase(key) == "energy"
	                push!(energies, atoms.data[key].data/length(atoms))
	            end
	        end
	    end 
	    return energies
	end;

	bead_tiny_energies = extract_energies(bead_tiny_dataset)
	bead_energies = extract_energies(bead_dataset)
	
	GC.gc()
end

function assess_model(model, train_dataset)

    plot([-0.5,0.5], [-0.5,0.5]; lc=:black, label="")

    model_energies = []
    model_std = []
    for atoms in bead_dataset
        ene, co_ene = ACE1.co_energy(model.potential, atoms)
        push!(model_energies, ene/length(atoms))
        push!(model_std, std(co_ene/length(atoms)))
    end
	
    rmse = sqrt(sum((model_energies - bead_energies).^2)/length(bead_dataset))
    mae = sum(abs.(model_energies - bead_energies))/length(bead_dataset)

#println(bead_energies)
println("_____________________________")	
println(model_energies)
    scatter!(bead_energies, model_energies;
             label="full dataset",
             title = @sprintf("Structures Used In Training:  %i out of %i\n", length(train_dataset), length(bead_dataset)) *
                     @sprintf("RMSE (MAE) For Entire Dataset:  %.0f (%.0f) meV/atom", 1000*rmse, 1000*mae),
             titlefontsize = 8,
             yerror = model_std,
             xlabel="Energy [eV/atom]", xlims=(-0.03,0.03),
             ylabel="Model Energy [eV/atom]", ylims=(-0.03,0.03),
             aspect_ratio = :equal, color="#F0BB62")

    model_energies = [energy(model.potential,atoms)/length(atoms) for atoms in train_dataset]
    scatter!(extract_energies(train_dataset), model_energies;
             label="training set", color="#064635")
	savefig("fitting_carbon_order5_p11_q30.png") 
end;

assess_model(model1, bead_tiny_dataset)
