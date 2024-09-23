### A Pluto.jl notebook ###
# v0.19.43

using Markdown
using InteractiveUtils

# ╔═╡ fa702ca2-75b4-11ef-04dd-59012347ffc7
begin
	using Pkg
	# Uncomment the next line if installing Julia for the first time
	# Pkg.Registry.add("General")
	Pkg.activate(".")
	Pkg.add(["LaTeXStrings", "MultivariateStats", "Plots", "PrettyTables",
	         "Suppressor", "ExtXYZ", "Unitful", "Distributed", "AtomsCalculators",])
	# Add the ACE registry, which stores the ACEpotentials package information
	Pkg.Registry.add(RegistrySpec(url="https://github.com/ACEsuit/ACEregistry"))
	Pkg.add("ACEpotentials")
end

# ╔═╡ 9d65703b-efef-4212-8d06-c4aae2204a45
begin
	Pkg.activate(".")
	Pkg.status()
end

# ╔═╡ 87a731f7-2a3b-4fc2-b504-82678ec4d3b9
begin
	using LaTeXStrings, MultivariateStats, Plots, PrettyTables, Printf,
	      Statistics, Suppressor, ExtXYZ, Unitful
	
	using ACEpotentials
end

# ╔═╡ 1524818b-9422-4cd1-a60a-54b3a55f0868
 Pkg.instantiate()

# ╔═╡ 90679dc9-7643-4d03-8de8-41aa5265c5c6
begin
	Si_tiny_dataset, _, _ = ACEpotentials.example_dataset("Si_tiny");
	if !isfile("Si_dataset.xyz")
	    download("https://www.dropbox.com/scl/fi/z6lvcpx3djp775zenz032/Si-PRX-2018.xyz?rlkey=ja5e9z99c3ta1ugra5ayq5lcv&st=cs6g7vbu&dl=1",
	         "Si_dataset.xyz");
	end
	
	Si_dataset = ExtXYZ.load("Si_dataset.xyz");
	
	# The last command generates a warning referring to missing pbc in the
	# first structure in the dataset, the isolated atom. We can safely remove this.
	deleteat!(Si_dataset, 1);
end

# ╔═╡ 90b67ec0-7899-42aa-9a23-b95c688b6a9c
begin
	using AtomsCalculators: potential_energy
	
	function assess_model(model, train_dataset)
	
	    plot([-164,-157], [-164,-157]; lc=:black, label="")
	
	    model_energies = []
	    model_std = []
	    for atoms in Si_dataset
	        ene, co_ene = @committee potential_energy(atoms, model)
	        push!(model_energies, ustrip(ene/length(atoms)))
	        push!(model_std, ustrip(std(co_ene/length(atoms))))
	    end
	    rmse = sqrt(sum((model_energies-Si_dataset_energies).^2)/length(Si_dataset))
	    mae = sum(abs.(model_energies-Si_dataset_energies))/length(Si_dataset)
	    scatter!(Si_dataset_energies, model_energies;
	             label="full dataset",
	             title = @sprintf("Structures Used In Training:  %i out of %i\n", length(train_dataset), length(Si_dataset)) *
	                     @sprintf("RMSE (MAE) For Entire Dataset:  %.0f (%.0f) meV/atom", 1000*rmse, 1000*mae),
	             titlefontsize = 8,
	             yerror = model_std,
	             xlabel="Energy [eV/atom]", xlims=(-164,-157),
	             ylabel="Model Energy [eV/atom]", ylims=(-164,-157),
	             aspect_ratio = :equal, color=1)
	
	    model_energies = [ustrip(potential_energy(atoms, model)/length(atoms)) for atoms in train_dataset]
	    scatter!(extract_energies(train_dataset), model_energies;
	             label="training set", color=2)
	end;
end

# ╔═╡ cd045e66-efc2-4f45-9d85-95506460af16
data_keys1 = (energy_key = "dft_energy", force_key = "dft_force", virial_key = "dft_virial")

# ╔═╡ d7b849ea-24a2-4265-a8a8-2d4243b28c3b
begin
	model = acemodel(elements = [:Si,],
	                   Eref = [:Si =>-158.54496821],
	                   order = 3,
	                   totaldegree = 12);
	
	acefit!(model, Si_tiny_dataset;
	        solver = ACEfit.BLR(committee_size = 30, factorization = :svd),
	        data_keys1...);
end

# ╔═╡ Cell order:
# ╠═fa702ca2-75b4-11ef-04dd-59012347ffc7
# ╠═9d65703b-efef-4212-8d06-c4aae2204a45
# ╠═87a731f7-2a3b-4fc2-b504-82678ec4d3b9
# ╠═1524818b-9422-4cd1-a60a-54b3a55f0868
# ╠═90679dc9-7643-4d03-8de8-41aa5265c5c6
# ╠═cd045e66-efc2-4f45-9d85-95506460af16
# ╠═d7b849ea-24a2-4265-a8a8-2d4243b28c3b
# ╠═90b67ec0-7899-42aa-9a23-b95c688b6a9c
