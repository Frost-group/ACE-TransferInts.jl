# ============================================================================
#  Basis-stability SANITY CHECK for thiophene (heteroatom case)
# ----------------------------------------------------------------------------
#  Generality check for the ethylene basis-stability figure: confirm that the
#  same behaviour holds for thiophene, the one system with a heteroatom (S).
#  C+S "CSi_SO" representation, pseudo-species labelling:
#     monomer A = 4 C + 1 S,  monomer B = 4 Si + 1 O.
#  Drive a cofacial dimer inward past the sampled band r in [7.25, 7.75] A and
#  check that the committee-mean J stays finite/smooth while sigma grows.
#  No paper figure is produced; this just validates the one-sentence
#  generality claim.  Hyperparameters match ACE_fit_NH_2_mol.jl.
# ============================================================================

using Pkg

const ROOT   = "C:/Users/Keerati/Documents/GitHub/ACE-TransferInts.jl/archive/Notebook_2024-2025/2025_06/10"
const FITDIR = joinpath(ROOT, "Fitting")
const DGDIR  = joinpath(ROOT, "Datagen")
const OUTDIR = joinpath(ROOT, "basis_test")
# The 06/10 Fitting env has a broken MbedTLS_jll artifact, so reuse the
# known-good, already-precompiled 06/05 env (same ACEpotentials 0.6.7).
const ENVDIR = "C:/Users/Keerati/Documents/GitHub/ACE-TransferInts.jl/archive/Notebook_2024-2025/2025_06/05/Fitting"

Pkg.activate(ENVDIR)
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

isdir(OUTDIR) || mkpath(OUTDIR)

# ---- model (identical to ACE_fit_NH_2_mol.jl) ------------------------------
const P, Q, OD, TD, RCUT = 19, 45, 2, 10, 9.0
const N_TRAIN  = 2000
const SCALE    = 5          # thiophene per-molecule scale (5 heavy atoms/monomer)
const TRAIN_LO = 7.25
const TRAIN_HI = 7.75

data_keys = (energy_key = "energy", force_key = "", virial_key = "")

println("Loading thiophene C+S training frames ...")
frames  = read_extxyz(joinpath(FITDIR, "all_CSi_SO_modified_rand1000_Num8000.extxyz"))
n_total = length(frames)
train_indices = round.(Int, range(1, stop = n_total, length = N_TRAIN))
train_data    = frames[train_indices]
println("  total = $n_total ; training N = $N_TRAIN")

model = acemodel(elements = [:Si, :C, :S, :O],
    order          = OD,
    totaldegree    = TD,
    r0             = 7.5,
    rcut           = RCUT,
    pair_transform = (:agnesi, P, Q),
    pair_envelope  = (:x, 2, 2),
    Eref           = [:Si => 0.0, :C => 0.0, :S => 0.0, :O => 0.0])
@show length(model.basis)

println("Fitting BLR committee (M = 500) ...")
solver = ACEfit.BLR(committee_size = 500, factorization = :svd)
@time acefit!(model, train_data; solver = solver, data_keys...)

@info "Training errors"
ACEpotentials.linear_errors(train_data, model; data_keys...)

# ---- cofacial thiophene dimer ----------------------------------------------
mono   = read_extxyz(joinpath(DGDIR, "thiophene.extxyz"))[1]
ZCc, ZSc = AtomicNumber(:C), AtomicNumber(:S)
heavy  = findall(z -> z == ZCc || z == ZSc, mono.Z)        # 4 C + 1 S, in file order
hpos0  = [mono.X[i] for i in heavy]
hZ     = [mono.Z[i] for i in heavy]
cen    = sum(hpos0) / length(hpos0)
hpos   = [p - cen for p in hpos0]                           # centre the heavy-atom monomer
@assert length(hpos) == 5 "expected 5 heavy atoms (4 C + 1 S), got $(length(hpos))"

relabel(z) = z == ZCc ? AtomicNumber(:Si) : (z == ZSc ? AtomicNumber(:O) : z)
ZB = relabel.(hZ)                                          # monomer B: 4 Si + 1 O

function cofacial_dimer(r)
    shift = JVecF(0.0, 0.0, r)
    posB  = [x + shift for x in hpos]
    at    = Atoms(vcat(hZ, ZB), vcat(hpos, posB))
    set_cell!(at, 20.0 * I(3)); set_pbc!(at, false)
    return at
end

# ---- sweep inward past the sampled band ------------------------------------
r_vals = collect(range(3.0, 9.0; step = 0.05))
J_ueV  = Float64[]; s_ueV = Float64[]
println("Sweeping r = $(first(r_vals)) ... $(last(r_vals)) A ...")
for r in r_vals
    at      = cofacial_dimer(r)
    n       = length(at)
    ene, co = ACE1.co_energy(model.potential, at)
    push!(J_ueV, ene / n * SCALE * 1e6)
    push!(s_ueV, std(co ./ n) * SCALE * 1e6)
end

# ---- save + verdict --------------------------------------------------------
csv = joinpath(OUTDIR, "basis_stability_cofacial_thiophene_N$(N_TRAIN).csv")
open(csv, "w") do f
    write(f, "r_angstrom,J_ueV_per_mol,sigma_ueV_per_mol\n")
    for (r, J, s) in zip(r_vals, J_ueV, s_ueV)
        write(f, @sprintf("%.4f,%.6f,%.6f\n", r, J, s))
    end
end

gr(dpi = 600)
plt = plot(r_vals, J_ueV; ribbon = s_ueV, fillalpha = 0.25, lw = 2,
    lc = "#064635", fillcolor = "#F0BB62", label = "committee mean J  ±  σ",
    xlabel = "intermolecular separation  r  (Å)",
    ylabel = "transfer integral J  (μeV / molecule)",
    title = @sprintf("Thiophene (C+S) basis-stability sanity check, N=%d, M=500", N_TRAIN),
    titlefontsize = 9, legend = :topright, size = (640, 460))
vspan!(plt, [TRAIN_LO, TRAIN_HI]; color = :gray, alpha = 0.15, label = "sampled range [7.25, 7.75] Å")
savefig(plt, joinpath(OUTDIR, "basis_stability_cofacial_thiophene_N$(N_TRAIN).png"))

s_in  = mean(s_ueV[(r_vals .>= TRAIN_LO) .& (r_vals .<= TRAIN_HI)])
s_deep = mean(s_ueV[r_vals .<= 4.5])
Jmax  = maximum(abs.(J_ueV))
println("\n==== THIOPHENE SANITY-CHECK VERDICT ====")
@printf("  σ in sampled band [7.25,7.75] : %.3g μeV\n", s_in)
@printf("  σ deep inside overlap (r≤4.5) : %.3g μeV   (ratio %.0fx)\n", s_deep, s_deep/s_in)
@printf("  |J| max over sweep            : %.3g μeV   (finite ⇒ no divergence)\n", Jmax)
@printf("  J finite everywhere?          : %s\n", all(isfinite, J_ueV) ? "YES" : "NO")
println("  Wrote $csv")
println("========================================")
