using Pkg
cd(@__DIR__)
Pkg.activate(".")
Pkg.instantiate()

using LaTeXStrings, MultivariateStats, Plots, Printf, Statistics, Suppressor, ACEpotentials, ACE1

# ── Load all frames ────────────────────────────────────────────────────────────
frames  = read_extxyz("all_CSi_modified_N.extxyz")
n_total = length(frames)
n_train = 450
train_indices    = round.(Int, range(1, stop=n_total, length=n_train))
bead_tiny_dataset = frames[train_indices]
bead_dataset      = [frames[i] for i in 1:n_total if i ∉ train_indices]
println("Train: $(length(bead_tiny_dataset))  Test: $(length(bead_dataset))  Total: $(length(bead_tiny_dataset)+length(bead_dataset))")

# ── ACE model ─────────────────────────────────────────────────────────────────
p = 19; q = 43; od = 2; td = 10; rc = 13.6
model1 = acemodel(
    elements     = [:Si, :C],
    order        = od,
    totaldegree  = td,
    r0           = 12.0,
    rcut         = rc,
    pair_transform = (:agnesi, p, q),
    pair_envelope  = (:x, 2, 2),
    Eref         = [:Si => 0.0, :C => 0.0],
)
@show length(model1.basis)
rC = Int(rc * 10)

# ── Fit ───────────────────────────────────────────────────────────────────────
data_keys1 = (energy_key = "energy", force_key = "", virial_key = "")
solver     = ACEfit.BLR(committee_size=500, factorization=:svd)
@time acefit!(model1, bead_tiny_dataset; solver=solver, data_keys1...)

@info "Training Errors"
ACEpotentials.linear_errors(bead_tiny_dataset, model1; data_keys1...)
@info "Test Errors"
ACEpotentials.linear_errors(bead_dataset,      model1; data_keys1...)

# ── Extract energies ──────────────────────────────────────────────────────────
function extract_energies(dataset)
    energies = Float64[]
    for atoms in dataset
        for key in keys(atoms.data)
            if lowercase(key) == "energy"
                push!(energies, atoms.data[key].data / length(atoms))
            end
        end
    end
    return energies
end

bead_tiny_energies = extract_energies(bead_tiny_dataset)
bead_energies      = extract_energies(bead_dataset)
GC.gc()

# ── Assess & plot ─────────────────────────────────────────────────────────────
function assess_model(model, train_dataset)
    gr(dpi=1200)
    plot([-0.5, 0.5], [-0.5, 0.5]; lc=:black, label="")

    model_energies = Float64[]
    model_std      = Float64[]
    for atoms in bead_dataset
        ene, co_ene = ACE1.co_energy(model.potential, atoms)
        push!(model_energies, ene / length(atoms))
        push!(model_std,      std(co_ene / length(atoms)))
    end

    err_rmse = sqrt(sum((model_energies .- bead_energies).^2) / length(bead_dataset))
    err_mae  = sum(abs.(model_energies .- bead_energies))     / length(bead_dataset)

    scatter!(bead_energies .* 10, model_energies .* 10;
             label      = "full dataset",
             title      = @sprintf("Structures Used In Training:  %i out of %i\n",
                                   length(train_dataset),
                                   length(bead_dataset) + length(train_dataset)) *
                          @sprintf("RMSE (MAE) For Entire Dataset:  %.0f (%.0f) meV/atom",
                                   1000 * err_rmse * 10, 1000 * err_mae * 10),
             titlefontsize = 8,
             yerror     = model_std,
             xlabel     = "Energy [eV/atom]", xlims = (-0.3, 0.3),
             ylabel     = "Model Energy [eV/atom]", ylims = (-0.3, 0.3),
             aspect_ratio = :equal, color = "#F0BB62")

    train_pred = [energy(model.potential, atoms) / length(atoms) for atoms in train_dataset]
    scatter!(extract_energies(train_dataset) .* 10, train_pred .* 10;
             label = "training set", color = "#064635",
             size = (500, 500), fontsize = 12)

    outfile = "fitting_NH_CSi_od$(od)_td$(td)_p$(p)_q$(q)_rcut$(rC)_N$(n_train)_mol.png"
    savefig(outfile)
    println("Saved: $outfile")
end

assess_model(model1, bead_tiny_dataset)
