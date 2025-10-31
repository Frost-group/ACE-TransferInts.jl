#cd("/Users/kk419/Documents/GitHub/ACE-TransferInts.jl/Notebook/2025_06/03/Fitting/")
cd("C:/Users/Keerati/Documents/GitHub/ACE-TransferInts.jl/Notebook/2025_04/10/Fitting/")


begin
	# --- Load all frames ---
	frames = read_extxyz("all_CSi_modified_rand1000_Num8000.extxyz")
	#frames = read_extxyz("all_CH_modified_rand1000_Num50000_all.extxyz")
	# Total number of frames
	n_total = length(frames)

	# Number of frames to select for training
	n_train =1500

	# Uniformly spaced indices
	train_indices = round.(Int, range(1, stop=n_total, length=n_train))

	# Create training and testing sets
	bead_tiny_dataset = frames[train_indices]
	bead_dataset = [frames[i] for i in 1:n_total if i ∉ train_indices]
end


begin
	p = 19 #19
	q = 45 #45
	od = 2
	td = 10
	rcut = 9.0
	model1 = acemodel(elements = [:Si,:C,],
	        order = od,
	        totaldegree =td, 
			r0 = 7.5,
	        rcut = rcut,
			#wL = 1.0,
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
maximum(bead_energies)

function assess_model(model, train_dataset)
	gr(dpi=1200)
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
    scatter!(bead_energies*4, model_energies*4;
             label="full dataset",
             title = @sprintf("Structures Used In Training:  %i out of %i\n", length(train_dataset), length(bead_dataset)+length(train_dataset)) *
                     @sprintf("RMSE (MAE) For Entire Dataset:  %.0f (%.0f) μeV/molecule", 1e+6*rmse*4, 1e+6*mae*4),
             titlefontsize = 8,
             yerror = model_std,
           	 xlabel="Energy [μeV/molecule]", xlims=(-7e-5, 11e-5),
             ylabel="Model Energy [μeV/molecule]", ylims=(-7e-5, 11e-5),
             aspect_ratio = :equal, color="#F0BB62")

    model_energies = [energy(model.potential,atoms)/length(atoms) for atoms in train_dataset]
    scatter!(extract_energies(train_dataset)*4, model_energies*4;
            label="training set", color="#064635",
			size=(500,500), fontsize=12,
			xticks = ([-6e-5, -4e-5, -2e-5, 0, 2e-5, 4e-5, 6e-5, 8e-5, 10e-5],
            		["-60", "-40","-20", "0", "20", "40", "60", "80" , "100"]),
			yticks = ([-6e-5, -4e-5, -2e-5, 0, 2e-5, 4e-5, 6e-5, 8e-5, 10e-5],
            		["-60", "-40","-20", "0", "20", "40", "60", "80", "100"]),
			)
#savefig("fitting_NH_CH_od$od"*"_td$td"*"_p$p"*"_q$q"*"_rcut75_N3000_2.png") 
#savefig("fitting_NH_CSi_od$od"*"_td$td"*"_p$p"*"_q$q"*"_rcut$rcut"*"_N$n_train"*"_rand_mol.png") 
	println(1e+6*rmse*4)
	println(1e+6*mae*4)
end;

assess_model(model1, bead_tiny_dataset)
 