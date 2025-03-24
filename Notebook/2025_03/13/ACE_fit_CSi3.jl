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
	savefig("fitting_carbon_order2_p11_p30_rcut14.png") 
end;

assess_model(model1, bead_tiny_dataset)
