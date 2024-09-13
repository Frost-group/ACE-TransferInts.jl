### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ e212233e-6c4c-11ef-115e-476408bfabc1
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

# ╔═╡ 018914e2-b683-4ce1-b38d-c0a3f1a418b0
begin
	Pkg.activate(".")
		Pkg.Registry.add("General")  # only needed when installing Julia for the first time
		Pkg.Registry.add(RegistrySpec(url="https://github.com/ACEsuit/ACEregistry"))
		Pkg.add("ACEpotentials")
		using ACEpotentials
end

# ╔═╡ 953fd3da-8994-4c78-baa9-f1588ad3e5bd
Pkg.status()

# ╔═╡ eea2c9b4-94bd-4947-a4c3-54023c90e6d9
bead_tiny_dataset = read_extxyz("train_bead.extxyz")

# ╔═╡ d54bfe21-dc89-41e6-8fa2-f1f1138ef006
bead_dataset = read_extxyz("test_bead.extxyz")

# ╔═╡ fb701de3-ddd4-4d43-8aac-cdb58e7c3a32
begin
	r_cut = 20.0
	rdf_tiny = ACEpotentials.get_rdf(bead_tiny_dataset, r_cut; rescale = true)
	plt_rdf_1 = histogram(rdf_tiny[(:Si, :Si)], 
		bins=150, label = "rdf",
		title="Bead_tiny_dataset", titlefontsize=10,
		xlabel = L"r[\AA]", ylabel = "RDF",
		yticks = [0,5,10], xlims=(1.5,20), ylims=(0.0, 10),
		size=(400,200), left_margin = 2Plots.mm)
	vline!(rnn(:Si)*[1.0, 4.65, 5.55], label = "r1, r2, ...", lw=2, color = "black")
	
	rdf = ACEpotentials.get_rdf(bead_dataset, r_cut; rescale = true);
	plt_rdf_2 = histogram(rdf[(:Si, :Si)],
		bins=150, label = "rdf",
		title="Bead_dataset", titlefontsize=10,
		xlabel = L"r[\AA]", ylabel = "RDF",
		yticks = [0,50,100], xlims=(1.5,20), ylims=(0.0, 100),
		size=(400,200), left_margin = 2Plots.mm)
	vline!(rnn(:Si)*[1.0, 4.65, 5.55], label = "r1, r2, ...", lw=2, color = "black")
	
	plot!(plt_rdf_1, plt_rdf_2, layout=(2,1), size=(600,400))
	#savefig("RDF_bead.png") 
end

# ╔═╡ 899acd5f-ebde-4e34-804d-3d805274e10e
rdf_tiny

# ╔═╡ d69886e4-daf3-493a-b1b0-f99b9130d2a1
begin
	r_cut_adf = 4.85 * rnn(:Si)
	eq_angle = 1.91 # radians
	adf_tiny = ACEpotentials.get_adf(bead_tiny_dataset, r_cut_adf);
	plt_adf_1 = histogram(adf_tiny, bins=50, label = "adf", yticks = [0,25,50], c = 3,
	                    title = "Bead_tiny_dataset", titlefontsize = 10,
	                    xlabel = L"\theta", ylabel = "ADF",
	                    xlims = (0, π), ylims=(0.0, 50), size=(400,200), left_margin = 2Plots.mm)
	vline!([ eq_angle,], label = "109.5˚", lw=5, color = "black")
	
	adf = ACEpotentials.get_adf(bead_dataset, r_cut_adf);
	plt_adf_2 = histogram(adf, bins=50, label = "adf", yticks = [0,250,500], c = 3,
	                    title = "Bead_dataset", titlefontsize = 10,
	                    xlabel = L"\theta", ylabel = "ADF",
	                    xlims = (0, π), ylims=(0.0, 500), size=(400,200), left_margin = 2Plots.mm)
	vline!([ eq_angle,], label = "109.5˚", lw=5, color = "black")
	
	plot(plt_adf_1, plt_adf_2, layout=(2,1), size=(600,400))
end

# ╔═╡ c697733c-a0ee-45b2-824d-d092ec81a00e
begin
	p = 1
	q = 2
	model1 = acemodel(elements = [:Si,],
	        order = 3,
	        totaldegree = 9, 
			r0 = 0.1,
	        rcut = 30.0,
			pair_transform = (:agnesi, p, q),
			pair_envelope = (:r, 1) 
	        #Eref = [:Si => 0.0]
	)
	@show length(model1.basis);
end

# ╔═╡ 94ee83f4-262f-4476-b7a7-340ee806e874
#solver = ACEfit.QR(lambda = 1e-6)
solver = ACEfit.LSQR(damp = 1e-6, atol = 1e-6)
#solver = ACEfit.RRQR(rtol = 1e-6)
#solver = ACEfit.BLR(committee_size=300, factorization=:svd)

# ╔═╡ a6e36e28-404e-4e97-b205-9aeb89117bef
data_keys1 = (energy_key = "energy", force_key = "force")

# ╔═╡ 76290b18-7c23-43d0-a567-e1e7b47fde60
acefit!(model1, bead_tiny_dataset; solver=solver, data_keys1...);

# ╔═╡ 077e0f84-aee9-4429-a61c-4b3dc887cce3
begin
	@info("Training Errors")
	ACEpotentials.linear_errors(bead_tiny_dataset, model1; data_keys1...);
	
	@info("Test Error")
	ACEpotentials.linear_errors(bead_dataset, model1; data_keys1...);
end

# ╔═╡ 9a394aaf-45c0-41d8-bfcd-a545eba618e4
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

# ╔═╡ e6f17e14-f929-4224-9547-0a48a7ebefc5
function assess_model(model, train_dataset)

    plot([0.001,0.01], [0.001,0.01]; lc=:black, label="")

    model_energies = []
    model_std = []
    for atoms in bead_dataset
        ene, co_ene = ACE1.co_energy(model.potential, atoms)
        push!(model_energies, ene/length(atoms))
        push!(model_std, std(co_ene/length(atoms)))
    end
	
    rmse = sqrt(sum((model_energies - bead_energies).^2)/length(bead_dataset))
    mae = sum(abs.(model_energies - bead_energies))/length(bead_dataset)

println(bead_energies)
println("_____________________________")	
println(model_energies)
    scatter!(bead_energies, model_energies;
             label="full dataset",
             title = @sprintf("Structures Used In Training:  %i out of %i\n", length(train_dataset), length(bead_dataset)) *
                     @sprintf("RMSE (MAE) For Entire Dataset:  %.0f (%.0f) meV/atom", 1000*rmse, 1000*mae),
             titlefontsize = 8,
             yerror = model_std,
             xlabel="Energy [eV/atom]", xlims=(0.001,0.01),
             ylabel="Model Energy [eV/atom]", ylims=(0.001,0.01),
             aspect_ratio = :equal, color=1)

    model_energies = [energy(model.potential,atoms)/length(atoms) for atoms in train_dataset]
    scatter!(extract_energies(train_dataset), model_energies;
             label="training set", color=2)
#savefig("fitting_Bead.png") 
end;

# ╔═╡ 4bc28c08-a51f-475a-8d7c-a06ee96992f0
assess_model(model1, bead_tiny_dataset)

# ╔═╡ Cell order:
# ╠═e212233e-6c4c-11ef-115e-476408bfabc1
# ╠═018914e2-b683-4ce1-b38d-c0a3f1a418b0
# ╠═953fd3da-8994-4c78-baa9-f1588ad3e5bd
# ╠═eea2c9b4-94bd-4947-a4c3-54023c90e6d9
# ╠═d54bfe21-dc89-41e6-8fa2-f1f1138ef006
# ╠═899acd5f-ebde-4e34-804d-3d805274e10e
# ╠═fb701de3-ddd4-4d43-8aac-cdb58e7c3a32
# ╠═d69886e4-daf3-493a-b1b0-f99b9130d2a1
# ╠═c697733c-a0ee-45b2-824d-d092ec81a00e
# ╠═94ee83f4-262f-4476-b7a7-340ee806e874
# ╠═a6e36e28-404e-4e97-b205-9aeb89117bef
# ╠═76290b18-7c23-43d0-a567-e1e7b47fde60
# ╠═077e0f84-aee9-4429-a61c-4b3dc887cce3
# ╠═9a394aaf-45c0-41d8-bfcd-a545eba618e4
# ╠═e6f17e14-f929-4224-9547-0a48a7ebefc5
# ╠═4bc28c08-a51f-475a-8d7c-a06ee96992f0
