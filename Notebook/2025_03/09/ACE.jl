using Pkg
Pkg.activate(".")
# uncomment the next line if the General registry is not yet installed,
# e.g. if using Julia for the first time.
# Pkg.Registry.add("General")
Pkg.activate(".")
Pkg.add(["LaTeXStrings", "MultivariateStats", "Plots", "PrettyTables", 
         "Suppressor", "ExtXYZ", "Unitful", "Distributed", "AtomsCalculators",
        "AtomsBuilder", "AtomsBase", "Lux", "Zygote", "Optimisers", "Folds"
         ])

using LaTeXStrings, MultivariateStats, Plots, PrettyTables, Printf, 
		 Statistics, Suppressor, ExtXYZ, Unitful, AtomsBuilder, AtomsBase

Pkg.Registry.add(RegistrySpec(url="https://github.com/ACEsuit/ACEregistry"))
Pkg.add("ACEpotentials")

using ACEpotentials

Pkg.status()
Pkg.update()
#   ## Part 1: Basic dataset analysis
# 
#   `ACEpotentials` provides quick access to several example datasets, which can
#   be useful for testing. The following command lists these datasets. (We
#   expect to expand this list signifcantly; please feel free to suggest
#   additions.)

Si_tiny_dataset, _, _ = ACEpotentials.example_dataset("Si_tiny");

if !isfile("Si_dataset.xyz")
   download("https://www.dropbox.com/scl/fi/z6lvcpx3djp775zenz032/Si-PRX-2018.xyz?rlkey=ja5e9z99c3ta1ugra5ayq5lcv&st=cs6g7vbu&dl=1",
        "Si_dataset.xyz");
end

Si_dataset = ExtXYZ.load("Si_dataset.xyz");

# The last command generates a warning referring to missing pbc in the
# first structure in the dataset, the isolated atom. We can safely remove this.
deleteat!(Si_dataset, 1);


println("The tiny dataset has ", length(Si_tiny_dataset), " structures.")
println("The large dataset has ", length(Si_dataset), " structures.")

	

r_cut = 6.0u"Å"
rnn = 2.35

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
     