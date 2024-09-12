using Pkg
Pkg.activate(".")
Pkg.add("LaTeXStrings")
Pkg.add("MultivariateStats")
Pkg.add("Plots")
Pkg.add("Suppressor")
using LaTeXStrings, MultivariateStats, Plots, Printf, Statistics, Suppressor

Pkg.activate(".")
Pkg.Registry.add("General")  # only needed when installing Julia for the first time
Pkg.Registry.add(RegistrySpec(url="https://github.com/ACEsuit/ACEregistry"))
Pkg.add("ACEpotentials")
using ACEpotentials

bead_tiny_dataset = read_extxyz("train_bead.extxyz")
bead_dataset = read_extxyz("test_bead.extxyz")

model1 = acemodel(elements = [:Si,],
	                 order = 3,
	                 totaldegree = 5,
	                 rcut = 15.0,
	                 Eref = [:Si => 0.0])
@show length(model1.basis);

solver1 = ACEfit.BLR(committee_size=100, factorization=:svd)
data_keys1 = (energy_key = "energy")

acefit!(model1, bead_tiny_dataset; solver=solver1, data_keys1...);

@info("Training Errors")
ACEpotentials.linear_errors(bead_tiny_dataset, model1; data_keys1...);

@info("Test Error")
ACEpotentials.linear_errors(bead_dataset, model1; data_keys1...);
