using Pkg
Pkg.activate("/Notebook/2024_11/19")
Pkg.instantiate()
Pkg.add(["LaTeXStrings", "MultivariateStats", "Plots", "PrettyTables", 
         "Suppressor", "ExtXYZ", "Unitful", "Distributed", "AtomsCalculators",
        "AtomsBuilder", "AtomsBase"
         ])

using LaTeXStrings, MultivariateStats, Plots, PrettyTables, Printf, 
		 Statistics, Suppressor, ExtXYZ, Unitful, AtomsBuilder, AtomsBase


Pkg.Registry.add(RegistrySpec(url="https://github.com/ACEsuit/ACEregistry"))
Pkg.add("ACEpotentials")

using ACEpotentials

Pkg.status()
Pkg.update()

begin
	bead_tiny_dataset = ExtXYZ.load("Notebook/2024_11/19/train_carbon_CSi3.extxyz")
	bead_dataset = ExtXYZ.load("Notebook/2024_11/19/test_carbon_CSi3.extxyz")
end

begin
r_cut = 20.0 * u"Å"
rdf = ACEpotentials.get_rdf(bead_tiny_dataset, r_cut; rescale = true)
end


using ACEpotentials, Plots, LaTeXStrings, Unitful
using AtomsBuilder: bulk, rattle!
using AtomsBuilder.Chemistry: rnn    



begin
	r_cut = 20.0 * u"Å"
	rdf_tiny = ACEpotentials.get_rdf(bead_tiny_dataset, r_cut; rescale = true)
	plt_rdf_1 = stephist(rdf_tiny[(:Si, :C)], 
		bins=150, label = "rdf",
		title="Bead_tiny_dataset", titlefontsize=10,
		xlabel = L"r[\AA]", ylabel = "RDF",
		yticks = [0,300,600,900,1200], xlims=(1.5,20), ylims=(0.0, 1200),
		size=(400,200), left_margin = 2Plots.mm, fill=true, fillalpha=0.5, color= "#F0BB62")
	#vline!(rnn(:Si)*[1.0, 4.65, 5.55], label = "r1, r2, ...", lw=2, color = "black")
	
	rdf = ACEpotentials.get_rdf(bead_dataset, r_cut; rescale = true);
	plt_rdf_2 = stephist(rdf[(:Si, :C)],
		bins=150, label = "rdf",
		title="Bead_dataset", titlefontsize=10,
		xlabel = L"r[\AA]", ylabel = "RDF",
		yticks = [0,800,1600,2400,3200], xlims=(1.5,20), ylims=(0.0, 3200),
		size=(400,200), left_margin = 2Plots.mm, fill=true, fillalpha=0.5, color= "#F0BB62")
	#vline!(rnn(:Si)*[1.0, 4.65, 5.55], label = "r1, r2, ...", lw=2, color = "black")
	
	plot!(plt_rdf_1, plt_rdf_2, layout=(2,1), size=(600,400))
	#savefig("RDF_bead.png") 
end

begin
	p = 1
	q = 4
	model1 = ace1_model(elements = [:Si,:C,],
	        order = 3,
	        totaldegree = [25, 20, 15], #{2,3,4} 
			r0 = 12.0,
	        rcut = 16.0,
			#transform = (:agnesi, p, q),
			pair_transform = (:agnesi, p, q),
			pair_envelope = (:x, 2, 2),
			#envelope = (:x, 2, 2),
	        Eref = [:Si => 0.0, :C => 0.0]
	)
end
@show length_basis(model1);


data_keys1 = (energy_key = "energy", force_key = "", virial_key = "")
solver = ACEfit.BLR(committee_size=500, factorization=:svd)
#Sparse solver
acefit!(bead_tiny_dataset, model1; solver=solver, data_keys1...);

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
        for key in keys(atoms)
            if lowercase(String(key)) == "energy"
                push!(energies, atoms[key] / length(atoms))
            end
        end
    end
    return energies
end;

bead_tiny_energies = extract_energies(bead_tiny_dataset)
bead_energies = extract_energies(bead_dataset)
	
	GC.gc()
end

using AtomsCalculators: potential_energy

function assess_model(model, train_dataset)

    plot([-0.04,0.04], [-0.04,0.04]; lc=:black, label="")

    model_energies = []
    model_std = []
    for atoms in bead_dataset
        ene, co_ene = @committee potential_energy(atoms, model)
        #ene, co_ene = ACE1.co_energy(model.potential, atoms)
        push!(model_energies, ustrip(ene/length(atoms)))
        push!(model_std, ustrip(std(co_ene/length(atoms))))
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
             xlabel="Energy [eV/atom]", xlims=(-.02,0.02),
             ylabel="Model Energy [eV/atom]", ylims=(-0.02,0.02),
             aspect_ratio = :equal, color="#F0BB62")

    model_energies = [ustrip(potential_energy(atoms, model)/length(atoms)) for atoms in train_dataset]
    scatter!(extract_energies(train_dataset), model_energies;
             label="training set", color="#064635")
	savefig("C:/Users/Keerati/Documents/GitHub/ACE-TransferInts.jl/Notebook/2024_11/19/fitting_carbon_0_8_5.png") 
#savefig("/Users/kk419/Documents/GitHub/ACE-TransferInts.jl/Notebook/2024_11/19/fitting_Bead_SiO2.png") 
end;

assess_model(model1, bead_tiny_dataset)