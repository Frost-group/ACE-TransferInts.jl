cd("C:/Users/Keerati/Documents/GitHub/ACE-TransferInts.jl/Notebook/2025_04/13/Fitting/")

begin
	# --- Load all frames ---
	frames = read_extxyz("all_CSi_modified_N.extxyz")

	# Total number of frames
	n_total = length(frames)

	# Number of frames to select for training
	n_train = 1350
	# Uniformly spaced indices
	train_indices = round.(Int, range(1, stop=n_total, length=n_train))

	# Create training and testing sets
	bead_tiny_dataset = frames[train_indices]
	bead_dataset = [frames[i] for i in 1:n_total if i ∉ train_indices]

end
length(bead_tiny_dataset)+length(bead_dataset)

begin
	p = 19
	q = 43
	od = 2
	td = 10
	rc = 13.6
	model1 = acemodel(elements = [:Si,:C,],
	        order = od,
	        totaldegree =td, 
			r0 = 12.0,
	        rcut = rc,
			#transform = (:agnesi, p, q),
			pair_transform = (:agnesi, p, q),
			pair_envelope = (:x, 2, 2),
			#envelope = (:x, 2, 2),
	        Eref = [:Si => 0.0, :C => 0.0]
	)
	@show length(model1.basis);
end
rC = Int(rc*10)


data_keys1 = (energy_key = "energy", force_key = "", virial_key = "")
solver = ACEfit.BLR(committee_size=500, factorization=:svd)
@time acefit!(model1, bead_tiny_dataset; solver=solver, data_keys1...);

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


# --- Sign prediction accuracy ---
begin
	function sign_accuracy(model, dataset, ref_energies; scale=10)
		# Get predicted energies for all structures in dataset
		pred_energies = [energy(model.potential, atoms)/length(atoms) 
		                 for atoms in dataset]
		
		# Convert to array and scale to per-molecule
		ref  = collect(ref_energies)  .* scale
		pred = collect(pred_energies) .* scale

		# Compute sign agreement
		n_correct = sum(sign.(pred) .== sign.(ref))
		n_total   = length(ref)
		acc = n_correct / n_total * 100  # in percent

		println("_____________________________")
		println("Sign Prediction Accuracy:")
		println("  Correct signs : $n_correct / $n_total")
		@printf("  Accuracy      : %.2f %%\n", acc)
		println("_____________________________")

		return acc
	end

	# Compute sign accuracy on test set
	@info("Sign Prediction Accuracy - Test Set")
	acc_test = sign_accuracy(model1, bead_dataset, bead_energies; scale=10)

	# Compute sign accuracy on training set
	@info("Sign Prediction Accuracy - Training Set")
	acc_train = sign_accuracy(model1, bead_tiny_dataset, bead_tiny_energies; scale=10)
end


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
	
    rmse = sqrt(sum((model_energies - bead_energies).^2)/length(bead_dataset));
    mae = sum(abs.(model_energies - bead_energies))/length(bead_dataset);

	# Compute sign accuracy for annotation
	n_correct = sum(sign.(collect(model_energies) .* 10) .== sign.(collect(bead_energies) .* 10))
	sign_acc  = n_correct / length(bead_dataset) * 100

println("_____________________________")	
println(model_energies)
    scatter!(bead_energies*10, model_energies*10;
             label="full dataset",
             title = @sprintf("Structures Used In Training:  %i out of %i\n", length(train_dataset), length(bead_dataset)+length(train_dataset)) *
                     @sprintf("RMSE (MAE) For Entire Dataset:  %.0f (%.0f) meV/atom\n", 1000*rmse*10, 1000*mae*10) *
                     @sprintf("Sign Accuracy: %.1f %%", sign_acc),
             titlefontsize = 8,
             yerror = model_std,
             xlabel="Energy [eV/atom]", xlims=(-0.3,0.3),
             ylabel="Model Energy [eV/atom]", ylims=(-0.3,0.3),
             aspect_ratio = :equal, color="#F0BB62")

    model_energies = [energy(model.potential,atoms)/length(atoms) for atoms in train_dataset];
    scatter!(extract_energies(train_dataset)*10, model_energies*10;
            label="training set", color="#064635",
			size=(500,500), fontsize=12,
			)
savefig("fitting_NH_CSi_od$od"*"_td$td"*"_p$p"*"_q$q"*"_rcut$rC"*"_N$n_train"*"_mol.png") 
end;

assess_model(model1, bead_tiny_dataset)