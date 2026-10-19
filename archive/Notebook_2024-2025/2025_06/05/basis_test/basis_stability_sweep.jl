# ============================================================================
#  Basis-stability demonstration sweep  (Referee 2, Point 3 / Referee 1, Point 4)
# ----------------------------------------------------------------------------
#  Take one trained ethylene ACE model (BLR committee, M = 500) and push a
#  cofacial dimer INWARD past the sampled separation range (training band is
#  r in [6.25, 7.00] Angstrom).  At each separation we evaluate the committee
#  mean transfer integral J and the committee standard deviation sigma.
#
#  Expected: J stays finite and smooth as r -> 0 (the radial envelope
#  f_env = (y - 1)^2 (y - y_cut)^2 forces the radial basis smoothly to zero,
#  so the prediction cannot diverge), while sigma grows once the configuration
#  leaves the training domain -- demonstrating both basis stability at high
#  molecular overlap AND honest extrapolation flagging by the committee UQ.
#
#  Run:  julia "<this file>"
#  (Activates the sibling Fitting/ project so the ACEpotentials Manifest is reused.)
# ============================================================================

using Pkg

# --- Reuse the existing, instantiated Fitting environment -------------------
const ROOT     = "C:/Users/Keerati/Documents/GitHub/ACE-TransferInts.jl/archive/Notebook_2024-2025/2025_06/05"
const FITDIR   = joinpath(ROOT, "Fitting")
const DGDIR    = joinpath(ROOT, "Datagen")
const OUTDIR   = joinpath(ROOT, "basis_test")
# Paper figure destination
const PAPERFIG = "C:/Users/Keerati/Documents/GitHub/2025-04-ACEtransfer/Figures/BasisStability"

Pkg.activate(FITDIR)
Base.active_project()
Pkg.instantiate()
Pkg.add(["JuLIP", "StaticArrays", "LinearAlgebra", "Plots","LaTeXStrings", "Rotations", "Random", "Printf"])

using JuLIP
using StaticArrays
using LinearAlgebra
using Plots
using LaTeXStrings
using Rotations  # for Euler rotation matrices
using Random
using Printf

using ACEpotentials, Statistics, Printf, Plots, LinearAlgebra
using JuLIP: Atoms, AtomicNumber, set_cell!, set_pbc!, JVecF, read_extxyz, chemical_symbol

isdir(OUTDIR)   || mkpath(OUTDIR)
isdir(PAPERFIG) || mkpath(PAPERFIG)

# ============================================================================
#  1.  Fit one model  (identical hyperparameters to the paper / ACE_fit_NH_mol.jl)
# ============================================================================
const P, Q, OD, TD, RCUT = 19, 45, 2, 10, 7.5
const N_TRAIN  = 2000
const SCALE    = 2          # per-molecule scale used throughout the paper
const TRAIN_LO = 6.25       # sampled separation range (Angstrom)
const TRAIN_HI = 7.00

data_keys = (energy_key = "energy", force_key = "", virial_key = "")

println("Loading training frames ...")
frames  = read_extxyz(joinpath(FITDIR, "all_CSi_modified_rand1000_Num8000.extxyz"))
n_total = length(frames)
train_indices = round.(Int, range(1, stop = n_total, length = N_TRAIN))
train_data    = frames[train_indices]
println("  total frames = $n_total ; training on N = $N_TRAIN")

model = acemodel(elements = [:Si, :C],
    order          = OD,
    totaldegree    = TD,
    r0             = 6.625,
    rcut           = RCUT,
    pair_transform = (:agnesi, P, Q),
    pair_envelope  = (:x, 2, 2),
    Eref           = [:Si => 0.0, :C => 0.0])
@show length(model.basis)

println("Fitting BLR committee (M = 500) ...")
solver = ACEfit.BLR(committee_size = 500, factorization = :svd)
@time acefit!(model, train_data; solver = solver, data_keys...)

@info "Training errors"
ACEpotentials.linear_errors(train_data, model; data_keys...)

# ============================================================================
#  2.  Build the cofacial ethylene dimer at separation r
# ============================================================================
#  Heavy-atom (carbon-only) representation, pseudo-species labelling:
#    monomer A carbons -> C (Z = 6),  monomer B carbons -> Si (Z = 14).
#  Monomer B is the same C=C unit translated by r along the stacking axis (z),
#  giving a clean cofacial / sandwich geometry -- the high-overlap regime.
mono   = read_extxyz(joinpath(DGDIR, "ethylene.extxyz"))[1]
c_idx  = findall(z -> z == AtomicNumber(:C), mono.Z)
cpos   = [mono.X[i] for i in c_idx]                # two carbon positions of monomer A
@assert length(cpos) == 2 "expected 2 carbons in the ethylene monomer, got $(length(cpos))"

const ZC  = AtomicNumber(:C)
const ZSi = AtomicNumber(:Si)

function cofacial_dimer(r)
    shift = JVecF(0.0, 0.0, r)
    posB  = [x + shift for x in cpos]              # monomer B, stacked along z
    X     = vcat(cpos, posB)
    Z     = [ZC, ZC, ZSi, ZSi]
    at    = Atoms(Z, X)
    set_cell!(at, 20.0 * I(3))
    set_pbc!(at, false)
    return at
end

# ============================================================================
#  3.  Sweep separation inward past the sampled range
# ============================================================================
r_vals = collect(range(2.5, 7.5; step = 0.05))     # Angstrom
J_ueV  = Float64[]                                  # committee-mean J   [ueV/molecule]
s_ueV  = Float64[]                                  # committee-std sigma [ueV/molecule]

println("Sweeping separation r = $(first(r_vals)) ... $(last(r_vals)) Angstrom ...")
for r in r_vals
    at      = cofacial_dimer(r)
    n       = length(at)
    ene, co = ACE1.co_energy(model.potential, at)   # mean energy, committee vector
    push!(J_ueV, ene / n * SCALE * 1e6)
    push!(s_ueV, std(co ./ n) * SCALE * 1e6)
end

# ============================================================================
#  4.  Save data + figure
# ============================================================================
csv_path = joinpath(OUTDIR, "basis_stability_cofacial_ethylene_N$(N_TRAIN).csv")
open(csv_path, "w") do f
    write(f, "r_angstrom,J_ueV_per_mol,sigma_ueV_per_mol\n")
    for (r, J, s) in zip(r_vals, J_ueV, s_ueV)
        write(f, @sprintf("%.4f,%.6f,%.6f\n", r, J, s))
    end
end
println("Wrote $csv_path")

gr(dpi = 1200)
plt = plot(r_vals, J_ueV;
    ribbon     = s_ueV,
    fillalpha  = 0.25,
    lw         = 2,
    lc         = "#064635",
    fillcolor  = "#F0BB62",
    label      = "committee mean J  ±  σ",
    xlabel     = "intermolecular separation  r  (Å)",
    ylabel     = "transfer integral J  (μeV / molecule)",
    title      = "Basis stability of the ACE model at high molecular overlap\n" *
                 @sprintf("cofacial ethylene dimer, N = %d, M = 500 committee", N_TRAIN),
    titlefontsize = 9,
    legend     = :topright,
    size       = (640, 460))

# shade the sampled (training) separation band
vspan!(plt, [TRAIN_LO, TRAIN_HI]; color = :gray, alpha = 0.15, label = "sampled range [6.25, 7.0] Å")

png_local = joinpath(OUTDIR,   "basis_stability_cofacial_ethylene_N$(N_TRAIN).png")
png_paper = joinpath(PAPERFIG, "basis_stability_cofacial_ethylene_N$(N_TRAIN).png")
savefig(plt, png_local)
savefig(plt, png_paper)
println("Wrote $png_local")
println("Wrote $png_paper")

println("\nDone.  Quick summary:")
@printf("  σ in sampled band (r∈[6.25,7.0]) : mean = %.3g μeV\n",
        mean(s_ueV[(r_vals .>= TRAIN_LO) .& (r_vals .<= TRAIN_HI)]))
@printf("  σ deep inside overlap (r ≤ 4 Å)  : mean = %.3g μeV\n",
        mean(s_ueV[r_vals .<= 4.0]))
@printf("  |J| max over sweep               : %.3g μeV  (finite ⇒ no divergence)\n",
        maximum(abs.(J_ueV)))
