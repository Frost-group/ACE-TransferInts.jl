using Pkg

#cd("/Users/kk419/Documents/GitHub/ACE-TransferInts.jl/Notebook/2025_05/20/Fitting/")
#cd("C:/Users/Keerati/Documents/GitHub/ACE-TransferInts.jl/Notebook/2025_05/20/Fitting/")
cd("C:/Users/Keerati/OneDrive - Imperial College London/Documents/GitHub/ACE-TransferInts.jl/Notebook/2025_05/20/Fitting/")

Pkg.activate(".")
Base.active_project()
Pkg.instantiate()
begin# add and load general packages used in this notebook.
	Pkg.add(["LaTeXStrings", "MultivariateStats", "Plots", "Suppressor", "StatsBase"
         ])
	using LaTeXStrings, MultivariateStats, Plots, Printf, Statistics, Suppressor, StatsBase
end;
begin
	#Pkg.Registry.add("General")  # only needed when installing Julia for the first time
	Pkg.Registry.add(RegistrySpec(url="https://github.com/ACEsuit/ACEregistry"))
	Pkg.add(PackageSpec(name="ACEpotentials", version="0.6.7"))
end

using LaTeXStrings, MultivariateStats, Plots, Printf,
	Statistics, Suppressor, ACEpotentials

println(pwd())
