using Plots, Printf

# ============================================================
# CSV reader — reads our custom row-based CSV format:
#   n_train,  v1, v2, ...
#   rmse,     v1, v2, ...
#   mae,      v1, v2, ...
#   sign_acc, v1, v2, ...
# ============================================================
function read_lc_csv(filepath)
    data = Dict{String, Vector{Float64}}()
    open(filepath, "r") do f
        for line in eachline(f)
            parts = split(strip(line), ",")
            key   = strip(parts[1])
            vals  = parse.(Float64, strip.(parts[2:end]))
            data[key] = vals
        end
    end
    n        = Int.(data["n_train"])
    rmse     = data["rmse"]
    mae      = data["mae"]
    sign_acc = haskey(data, "sign_acc") ? data["sign_acc"] : fill(NaN, length(n))
    return n, rmse, mae, sign_acc
end

# ============================================================
# Load CSV files — update paths to match your directory
# ============================================================
eth_n,  eth_rmse,  eth_mae,  eth_sign  = read_lc_csv("LC_Ethylene_CSi_od2_td10_p19_q45_rcut75.csv")
thi_n,  thi_rmse,  thi_mae,  thi_sign  = read_lc_csv("LC_Thiophene_CSi_od2_td10_p19_q45_rcut90.csv")
thi2_n, thi2_rmse, thi2_mae, thi2_sign = read_lc_csv("LC_Thiophene_CSi_SO_od2_td10_p19_q45_rcut90.csv")
nph_n,  nph_rmse,  nph_mae,  nph_sign  = read_lc_csv("LC_HeavyAtom_naphthalene_CSi_od2_td10_p19_q43_rcut136.csv")
bead_n, bead_rmse, bead_mae, bead_sign = read_lc_csv("LC_Bead_naphthalene_CSi_od2_td10_p19_q43_rcut136.csv")

# Print summary to verify loading
println("Loaded data summary:")
println("  Ethylene        n_train = ", eth_n,  "  final RMSE = ", eth_rmse[end])
println("  Thiophene C,Si  n_train = ", thi_n,  "  final RMSE = ", thi_rmse[end])
println("  Thiophene CSSO  n_train = ", thi2_n, "  final RMSE = ", thi2_rmse[end])
println("  Naph Heavy-atom n_train = ", nph_n,  "  final RMSE = ", nph_rmse[end])
println("  Naph Bead       n_train = ", bead_n, "  final RMSE = ", bead_rmse[end])

# ============================================================
# Colours
# ============================================================
c_eth   = "#F0BB62"   # ethylene
c_thi1  = "#6E9A50"   # thiophene C,Si
c_thi2  = "#184D47"   # thiophene C,Si,S,O
c_nph   = "#C64756"   # naphthalene heavy-atom
c_bead  = "#7A2D36"   # naphthalene bead

# ============================================================
# PLOT 1 — Ethylene & Thiophene (μeV/molecule, log-log scale)
# ============================================================
gr(dpi=1200)

p1 = plot(
    xlabel = "Number of training structures",
    ylabel = "RMSE (μeV/molecule)",
    title  = "Learning curves: Ethylene and Thiophene",
    xscale = :log10,
    yscale = :log10,
    legend = :topright,
    legendfontsize = 8,
    titlefontsize  = 10,
    size = (600, 450),
    grid = true,
    minorgrid = true,
    framestyle = :box,
)

# RMSE lines (solid)
plot!(p1, eth_n,  eth_rmse;
      label="Ethylene (C,Si)",
      color=c_eth,  lw=2, marker=:circle, markersize=5)

plot!(p1, thi_n,  thi_rmse;
      label="Thiophene (C,Si)",
      color=c_thi1, lw=2, marker=:circle, markersize=5)

plot!(p1, thi2_n, thi2_rmse;
      label="Thiophene (C,Si,S,O)",
      color=c_thi2, lw=2, marker=:circle, markersize=5)

# MAE lines (dashed)
plot!(p1, eth_n,  eth_mae;
      label="", color=c_eth,  lw=2, linestyle=:dash, marker=:diamond, markersize=4)

plot!(p1, thi_n,  thi_mae;
      label="", color=c_thi1, lw=2, linestyle=:dash, marker=:diamond, markersize=4)

plot!(p1, thi2_n, thi2_mae;
      label="", color=c_thi2, lw=2, linestyle=:dash, marker=:diamond, markersize=4)

# Legend proxies for line styles
plot!(p1, [NaN], [NaN]; label="RMSE (solid)", color=:black, lw=2, linestyle=:solid)
plot!(p1, [NaN], [NaN]; label="MAE (dashed)", color=:black, lw=2, linestyle=:dash)

# Vertical reference line at n=2000 (training size reported in paper)
vline!(p1, [2000];
       label="n = 2000",
       color=:gray40, lw=1.5, linestyle=:dot)

savefig(p1, "LC_Ethylene_Thiophene.png")
println("Plot 1 saved: LC_Ethylene_Thiophene.png")


# ============================================================
# PLOT 2 — Naphthalene (meV/molecule, log-log scale)
# Note: Bead n=100 is a massive outlier (RMSE ~26572 meV)
#       plotted on log scale so it appears but does not
#       destroy the y-axis range for the other points.
# ============================================================

p2 = plot(
    xlabel = "Number of training structures",
    ylabel = "RMSE (meV/molecule)",
    title  = "Learning curves: Naphthalene",
    xscale = :log10,
    yscale = :log10,
    legend = :topright,
    legendfontsize = 8,
    titlefontsize  = 10,
    size = (600, 450),
    grid = true,
    minorgrid = true,
    framestyle = :box,
)

# Heavy-atom RMSE and MAE
plot!(p2, nph_n, nph_rmse;
      label="Heavy-atom (C,Si)",
      color=c_nph, lw=2, marker=:circle, markersize=5)

plot!(p2, nph_n, nph_mae;
      label="", color=c_nph, lw=2, linestyle=:dash, marker=:diamond, markersize=4)

# Bead RMSE and MAE — full range including outlier at n=100
plot!(p2, bead_n, bead_rmse;
      label="Bead (C,Si)",
      color=c_bead, lw=2, marker=:circle, markersize=5)

plot!(p2, bead_n, bead_mae;
      label="", color=c_bead, lw=2, linestyle=:dash, marker=:diamond, markersize=4)

# Mark the n=100 bead outlier explicitly
#scatter!(p2, [bead_n[1]], [bead_rmse[1]];
#         label="Bead n=100 outlier",
#         color=c_bead, marker=:star5, markersize=8)

# Legend proxies for line styles
plot!(p2, [NaN], [NaN]; label="RMSE (solid)", color=:black, lw=2, linestyle=:solid)
plot!(p2, [NaN], [NaN]; label="MAE (dashed)", color=:black, lw=2, linestyle=:dash)

# Vertical reference line at n=1350 (training size reported in paper)
vline!(p2, [1350];
       label="n = 1350",
       color=:gray40, lw=1.5, linestyle=:dot)

savefig(p2, "LC_Naphthalene.png")
println("Plot 2 saved: LC_Naphthalene.png")

println("\nDone. Summary of final RMSE values:")
println("  Ethylene        (n=$(eth_n[end])):  $(eth_rmse[end]) μeV/molecule")
println("  Thiophene C,Si  (n=$(thi_n[end])):  $(thi_rmse[end]) μeV/molecule")
println("  Thiophene CSSO  (n=$(thi2_n[end])): $(thi2_rmse[end]) μeV/molecule")
println("  Naph Heavy-atom (n=$(nph_n[end])):  $(nph_rmse[end]) meV/molecule")
println("  Naph Bead       (n=$(bead_n[end])): $(bead_rmse[end]) meV/molecule")