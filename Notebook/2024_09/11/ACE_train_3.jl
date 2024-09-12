### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ a4ac8dbe-6f7a-11ef-3b9e-3d92a7a83913
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

# ╔═╡ 623495a8-df2a-423c-ae0e-46cf73ce8522
begin
	Pkg.activate(".")
	Pkg.Registry.add("General")  # only needed when installing Julia for the first time
	Pkg.Registry.add(RegistrySpec(url="https://github.com/ACEsuit/ACEregistry"))
	Pkg.add("ACEpotentials")
end

# ╔═╡ b8f884c6-88f6-4e7a-8519-871f0529c93a
begin
	pkg"activate pathtoproject"
	Pkg.instantiate()
end

# ╔═╡ 79cf027a-82d0-46a4-a400-862731b0716f
using ACEpotentials

# ╔═╡ 84ec872f-59de-4ca7-8255-faeae0fa5800


# ╔═╡ d36203a1-a897-4806-a964-629cd275bba6
bead_tiny_dataset = read_extxyz("train_d_Car.extxyz")

# ╔═╡ d81f1cc1-1086-42fb-b4df-bf1beb4459c0
bead_dataset = read_extxyz("test_d_Car.extxyz")

# ╔═╡ 5e4ed684-8c9c-4c26-ac3c-8fa8be88fa23
begin
	r_cut = 20.0
	rdf_tiny = ACEpotentials.get_rdf(bead_tiny_dataset, r_cut; rescale = true)
	plt_rdf_1 = histogram(rdf_tiny[(:C, :C)], 
		bins=150, label = "rdf",
		title="Bead_tiny_dataset", titlefontsize=10,
		xlabel = L"r[\AA]", ylabel = "RDF",
		yticks = [0,5,10], xlims=(0,20), ylims=(0.0, 100),
		size=(400,200), left_margin = 2Plots.mm)
	vline!(rnn(:C)*[1.0, 5.55, 7.1, 8.7], label = "r1, r2, ...", lw=2, color = "black")
	
	rdf = ACEpotentials.get_rdf(bead_dataset, r_cut; rescale = true);
	plt_rdf_2 = histogram(rdf[(:C, :C)],
		bins=150, label = "rdf",
		title="Bead_dataset", titlefontsize=10,
		xlabel = L"r[\AA]", ylabel = "RDF",
		yticks = [0,50,100], xlims=(0,20), ylims=(0.0, 250),
		size=(400,200), left_margin = 2Plots.mm)
	vline!(rnn(:C)*[1.0, 5.55, 7.1, 8.7], label = "r1, r2, ...", lw=2, color = "black")
	
	plot!(plt_rdf_1, plt_rdf_2, layout=(2,1), size=(600,400))
end

# ╔═╡ cfe6ac36-76ab-4950-8c72-c9bb5f50021a
begin
	r_cut_adf = 3 * rnn(:C)
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

# ╔═╡ a5bf0da5-49c4-42b9-9c6e-b4b3e2d1f21c
begin
	model1 = acemodel(elements = [:C,],
	                 order = 3,
	                 totaldegree = 5,
	                 rcut = 16.0,
	                 Eref = [:C => 0.0])
	@show length(model1.basis);
end

# ╔═╡ 8b3e78c1-26b6-4e80-8c22-1f47c0ed9b20
data_keys1 = (energy_key = "energy", force_key = "force")

# ╔═╡ d47dcd41-6ca8-4f23-8a6a-a6b450e97389
solver = ACEfit.BLR(committee_size=200, factorization=:svd)

# ╔═╡ 55e28fde-00eb-4c73-b75a-e73d45e2f3fd
acefit!(model1, bead_tiny_dataset; solver=solver, data_keys1...);

# ╔═╡ 317713d0-c834-4f22-a776-545e7f341a19
begin
	@info("Training Errors")
	ACEpotentials.linear_errors(bead_tiny_dataset, model1; data_keys1...);
	
	@info("Test Error")
	ACEpotentials.linear_errors(bead_dataset, model1; data_keys1...);
end

# ╔═╡ 8f847bcd-4d15-4b3d-b7ce-1865604d86ff
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

# ╔═╡ 0c25abd6-ffb5-4f1f-8223-2ef1116f6b11
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

end;

# ╔═╡ 84d449ab-d4eb-4999-9122-a642e344eb57
assess_model(model1, bead_tiny_dataset)

# ╔═╡ Cell order:
# ╠═a4ac8dbe-6f7a-11ef-3b9e-3d92a7a83913
# ╠═623495a8-df2a-423c-ae0e-46cf73ce8522
# ╠═b8f884c6-88f6-4e7a-8519-871f0529c93a
# ╠═79cf027a-82d0-46a4-a400-862731b0716f
# ╠═84ec872f-59de-4ca7-8255-faeae0fa5800
# ╠═d36203a1-a897-4806-a964-629cd275bba6
# ╠═d81f1cc1-1086-42fb-b4df-bf1beb4459c0
# ╠═5e4ed684-8c9c-4c26-ac3c-8fa8be88fa23
# ╠═cfe6ac36-76ab-4950-8c72-c9bb5f50021a
# ╠═a5bf0da5-49c4-42b9-9c6e-b4b3e2d1f21c
# ╠═8b3e78c1-26b6-4e80-8c22-1f47c0ed9b20
# ╠═d47dcd41-6ca8-4f23-8a6a-a6b450e97389
# ╠═55e28fde-00eb-4c73-b75a-e73d45e2f3fd
# ╠═317713d0-c834-4f22-a776-545e7f341a19
# ╠═8f847bcd-4d15-4b3d-b7ce-1865604d86ff
# ╠═0c25abd6-ffb5-4f1f-8223-2ef1116f6b11
# ╠═84d449ab-d4eb-4999-9122-a642e344eb57
