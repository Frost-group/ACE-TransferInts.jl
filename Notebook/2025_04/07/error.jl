using Pkg
cd("C:/Users/Keerati/Documents/GitHub/ACE-TransferInts.jl/Notebook/2025_04/07")
Pkg.add(["XLSX", "Plots", "DataFrames"])

using XLSX
using Plots
using DataFrames
# Load the file
xlsx_path = "ACEFIT_20250401.xlsx"

# Function to load sheet and plot
function plot_rmse_vs_parameter(sheet_number::Int, parameter::String)
    data = XLSX.readtable(xlsx_path, sheet_number) |> DataFrame
    x = data[:, parameter]
    y = data[:, "test_RMSE"]  # assuming column name is exactly "test_RMSE"
    
    plt = scatter(x, y,
        xlabel=parameter,
        ylabel="Test RMSE",
        title="Test RMSE vs $parameter",
        legend=false,
        markersize=6)
    return plt
end

# Plot settings
default(markerstrokecolor = :black, markerstrokewidth = 0.5)
gr(dpi=1200) 
# Plot for each parameter and sheet
plt = plot_rmse_vs_parameter(1, "rcut")
savefig(plt,"rcut.png") 
plt = plot_rmse_vs_parameter(3, "p")
savefig(plt,"p.png")
plt = plot_rmse_vs_parameter(4, "q")
savefig(plt,"q.png")
plt = plot_rmse_vs_parameter(5, "totaldegree")
savefig(plt,"totaldegree.png") 
