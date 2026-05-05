cd("C:/Users/Keerati/Documents/GitHub/ACE-TransferInts.jl/Notebook/2025_06/10/Fitting/")

# --- Load all frames once ---
frames = read_extxyz("all_CSi_modified_rand1000_Num8000.extxyz")
n_total = length(frames)

# --- Helpers defined once outside loop ---
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
end

function sign_accuracy(model, dataset, ref_energies; scale=4)
    pred_energies = [energy(model.potential, atoms)/length(atoms)
                     for atoms in dataset]
    ref  = collect(ref_energies)  .* scale
    pred = collect(pred_energies) .* scale
    n_correct = sum(sign.(pred) .== sign.(ref))
    n_total   = length(ref)
    acc = n_correct / n_total * 100
    println("_____________________________")
    println("Sign Prediction Accuracy:")
    println("  Correct signs : $n_correct / $n_total")
    @printf("  Accuracy      : %.2f %%\n", acc)
    println("_____________________________")
    return acc
end

# --- Model parameters ---
p    = 19
q    = 45
od   = 2
td   = 10
rcut = 9.0
rC   = Int(rcut*10)

data_keys1 = (energy_key = "energy", force_key = "", virial_key = "")

# --- Storage arrays for learning curve ---
n_train_values  = [100, 250, 500, 1000, 1500, 2000, 3000, 5000]
rmse_values     = Float64[]
mae_values      = Float64[]
sign_acc_values = Float64[]

# --- Loop over training sizes ---
for n_train in n_train_values

    println("\n")
    println("==============================================")
    println("  Training with n_train = $n_train")
    println("==============================================")

    # --- Split dataset ---
    train_indices     = round.(Int, range(1, stop=n_total, length=n_train))
    bead_tiny_dataset = frames[train_indices]
    bead_dataset      = [frames[i] for i in 1:n_total if i ∉ train_indices]

    # --- Build model ---
    model1 = acemodel(elements = [:Si, :C],
        order       = od,
        totaldegree = td,
        r0          = 7.5,
        rcut        = rcut,
        pair_transform = (:agnesi, p, q),
        pair_envelope  = (:x, 2, 2),
        Eref = [:Si => 0.0, :C => 0.0]
    )
    @show length(model1.basis)

    # --- Fit ---
    solver = ACEfit.BLR(committee_size=500, factorization=:svd)
    @time acefit!(model1, bead_tiny_dataset; solver=solver, data_keys1...)

    # --- Errors ---
    @info("Training Errors")
    ACEpotentials.linear_errors(bead_tiny_dataset, model1; data_keys1...)

    @info("Test Errors")
    ACEpotentials.linear_errors(bead_dataset, model1; data_keys1...)

    # --- Extract energies ---
    bead_tiny_energies = extract_energies(bead_tiny_dataset)
    bead_energies      = extract_energies(bead_dataset)
    GC.gc()

    println("Maximum reference energy: $(maximum(bead_energies))")

    # --- Sign accuracy ---
    @info("Sign Prediction Accuracy - Test Set")
    acc_test  = sign_accuracy(model1, bead_dataset,      bead_energies;      scale=4)

    @info("Sign Prediction Accuracy - Training Set")
    acc_train = sign_accuracy(model1, bead_tiny_dataset, bead_tiny_energies; scale=4)

    # --- Compute RMSE and MAE on test set ---
    model_energies = []
    model_std      = []
    for atoms in bead_dataset
        ene, co_ene = ACE1.co_energy(model1.potential, atoms)
        push!(model_energies, ene/length(atoms))
        push!(model_std, std(co_ene/length(atoms)))
    end

    rmse = sqrt(sum((model_energies - bead_energies).^2)/length(bead_dataset))
    mae  = sum(abs.(model_energies - bead_energies))/length(bead_dataset)

    # Scale to μeV/molecule (scale=4) and store
    push!(rmse_values,      1e6*rmse*4)
    push!(mae_values,       1e6*mae*4)
    push!(sign_acc_values,  acc_test)

    println("RMSE: $(1e6*rmse*4) μeV/molecule")
    println("MAE:  $(1e6*mae*4) μeV/molecule")

    # --- Sign accuracy for plot title ---
    n_correct = sum(sign.(collect(model_energies) .* 4) .== sign.(collect(bead_energies) .* 4))
    sign_acc  = n_correct / length(bead_dataset) * 100

    # --- Plot individual fit ---
    gr(dpi=1200)
    plot([-0.5, 0.5], [-0.5, 0.5]; lc=:black, label="")

    scatter!(bead_energies*4, model_energies*4;
             label="full dataset",
             title = @sprintf("Structures Used In Training:  %i out of %i\n", n_train, n_total) *
                     @sprintf("RMSE (MAE):  %.0f (%.0f) μeV/molecule\n", 1e6*rmse*4, 1e6*mae*4) *
                     @sprintf("Sign Accuracy: %.1f %%", sign_acc),
             titlefontsize = 8,
             yerror = model_std,
             xlabel="Energy [μeV/molecule]", xlims=(-7e-5, 11e-5),
             ylabel="Model Energy [μeV/molecule]", ylims=(-7e-5, 11e-5),
             aspect_ratio = :equal, color="#F0BB62")

    model_energies_train = [energy(model1.potential, atoms)/length(atoms) for atoms in bead_tiny_dataset]
    scatter!(extract_energies(bead_tiny_dataset)*4, model_energies_train*4;
             label="training set", color="#064635",
             size=(500,500), fontsize=12,
             xticks=([-6e-5, -4e-5, -2e-5, 0, 2e-5, 4e-5, 6e-5, 8e-5, 10e-5],
                     ["-60", "-40", "-20", "0", "20", "40", "60", "80", "100"]),
             yticks=([-6e-5, -4e-5, -2e-5, 0, 2e-5, 4e-5, 6e-5, 8e-5, 10e-5],
                     ["-60", "-40", "-20", "0", "20", "40", "60", "80", "100"]),
             )

    savefig("LC_fitting_Thi_CSi_od$od"*"_td$td"*"_p$p"*"_q$q"*"_rcut$rC"*"_N$n_train"*"_mol.png")

end # end loop

# --- Save learning curve data to CSV ---
println("\nSaving learning curve data to CSV...")

open("LC_Thi_CSi_od$od"*"_td$td"*"_p$p"*"_q$q"*"_rcut$rC"*".csv", "w") do f
    write(f, "n_train,"  * join(string.(n_train_values), ",")                        * "\n")
    write(f, "rmse,"     * join([@sprintf("%.4f", v) for v in rmse_values],    ",")  * "\n")
    write(f, "mae,"      * join([@sprintf("%.4f", v) for v in mae_values],     ",")  * "\n")
    write(f, "sign_acc," * join([@sprintf("%.2f", v) for v in sign_acc_values], ",") * "\n")
end

println("CSV saved.")
println("\nLearning curve summary:")
println("n_train  : ", n_train_values)
println("RMSE     : ", rmse_values)
println("MAE      : ", mae_values)
println("Sign Acc : ", sign_acc_values)
println("\nAll training runs complete.")