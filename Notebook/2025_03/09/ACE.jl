using Pkg
Pkg.activate(".")
# uncomment the next line if the General registry is not yet installed,
# e.g. if using Julia for the first time.
# Pkg.Registry.add("General")
Pkg.activate(".")
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


ACEpotentials.list_example_datasets()

Si_tiny_dataset, _, _ = ACEpotentials.example_dataset("Si_tiny");

if !isfile("Si_dataset.xyz")
    download("https://www.dropbox.com/scl/fi/z6lvcpx3djp775zenz032/Si-PRX-2018.xyz?rlkey=ja5e9z99c3ta1ugra5ayq5lcv&st=cs6g7vbu&dl=1",
         "Si_dataset.xyz");
end

Si_dataset = ExtXYZ.load("Si_dataset.xyz");



config_types_tiny = [at[:config_type] for at in Si_tiny_dataset]
config_types = [ at[:config_type] for at in Si_dataset]

function count_configs(config_types)
    config_counts = [sum(config_types .== ct) for ct in unique(config_types)]
    config_dict = Dict([ct=>cc for (ct,cc) in zip(unique(config_types), config_counts)])
end;

println("There are ", length(unique(config_types_tiny)),
        " unique config_types in the tiny dataset:")
display(count_configs(config_types_tiny))

println("There are ", length(unique(config_types)),
        " unique config_types in the full dataset:")
display(count_configs(config_types))



using ACEpotentials, Plots, LaTeXStrings, Unitful
using AtomsBuilder: bulk, rattle!
using AtomsBuilder.Chemistry: rnn
using AtomsBase: bounding_box   


function extract_energies(dataset)
    energies = []
    for atoms in dataset
        for key in keys(atoms)
            if lowercase(String(key)) == "dft_energy"
                push!(energies, atoms[key] / length(atoms))
            end
        end
    end
    return energies
end;

Si_dataset_energies = extract_energies(Si_dataset)
;  # the ; is just to suppress the ouput


r_cut = 6.0u"Å"

rdf_tiny = ACEpotentials.get_rdf(Si_tiny_dataset, r_cut; rescale = true)
plt_rdf_1 = histogram(rdf_tiny[(:Si, :Si)], bins=150, label = "rdf",
                      title="Si_tiny_dataset", titlefontsize=10,
                      xlabel = L"r[\AA]", ylabel = "RDF", yticks = [],
                      xlims=(1.5,6), size=(400,200), left_margin = 2Plots.mm)
vline!(rnn * [1.0, 1.633, 1.915, 2.3, 2.5], label = "r1, r2, ...", lw=3)

rdf = ACEpotentials.get_rdf(Si_dataset, r_cut; rescale = true);
plt_rdf_2 = histogram(rdf[(:Si, :Si)], bins=150, label = "rdf",
                      title="Si_dataset", titlefontsize=10,
                      xlabel = L"r[\AA]", ylabel = "RDF", yticks = [],
                      xlims=(1.5,6), size=(400,200), left_margin = 2Plots.mm)
vline!(rnn * [1.0, 1.633, 1.915, 2.3, 2.5], label = "r1, r2, ...", lw=3)

plot(plt_rdf_1, plt_rdf_2, layout=(2,1), size=(400,400))

model = ace1_model(elements = [:Si],
                   rcut = 5.5,
                   order = 3,        # body-order - 1
                   totaldegree = 8 );

                   
descriptors = []
for system in Si_tiny_dataset
    struct_descriptor = sum(site_descriptors(system, model)) / length(system)
    push!(descriptors, struct_descriptor)
end


descriptors = reduce(hcat, descriptors)  # convert to matrix
M = fit(PCA, descriptors; maxoutdim=3, pratio=1)
descriptors_trans = transform(M, descriptors)
p = scatter(
     descriptors_trans[1,:], descriptors_trans[2,:], descriptors_trans[3,:],
     marker=:circle, linewidth=0, group=config_types_tiny, legend=:right)
plot!(p, xlabel="PC1", ylabel="PC2", zlabel="PC3", camera=(20,10))



model = ace1_model(
              elements = [:Si,],
              order = 3,
              totaldegree = 8,
              rcut = 5.0,
              Eref = Dict(:Si => -158.54496821))

# `ace1_model` specifies a linear model `model`; because it is linear
# it is implicitly defined by a basis. In `ACEpotentials`, the size of
# this basis (= number of parameters) can be checked as follows
@show length_basis(model);

solver = ACEfit.QR(lambda=1e-1)
data_keys = (energy_key = "dft_energy", force_key = "dft_force", virial_key = "dft_virial")
acefit!(Si_tiny_dataset, model;
        solver=solver, data_keys...);

@info("Training Errors")
compute_errors(Si_tiny_dataset, model; data_keys...);

@info("Test Error")
compute_errors(Si_dataset, model; data_keys...);


function extract_energies(dataset)
    energies = []
    for atoms in dataset
        for key in keys(atoms)
            if lowercase(String(key)) == "dft_energy"
                push!(energies, atoms[key] / length(atoms))
            end
        end
    end
    return energies
end;

Si_dataset_energies = extract_energies(Si_dataset)