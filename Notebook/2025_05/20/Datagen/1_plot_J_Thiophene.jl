
dimer_set = read_extxyz("output750.extxyz")

θ_list = range(0.01, stop=π, length=49);
φ_list = range(0.01, stop=2π, length=98);

J_eff = zeros(size(θ_list,1),size(φ_list,1))
get_data(dimer_set[1], "theta")

for i in range(1,size(θ_list,1))
    for j in range(1,size(φ_list,1))

    J_eff[i,j] = dimer_set[ (i-1) * 98 + j]["energy_eff"]

    end
end 
gr(dpi=1200) 
heatmap(
    φ_list, θ_list, ((J_eff*10e+3)),
    
    size=(1200, 1200),
    xlabel= L"$\phi$",
	ylabel=L"$\theta$",

	yticks=([0,pi/2,pi],
            ["0", L"$\pi/2$", L"$\pi$"]),
	xticks=([0,pi/2,pi,3pi/2,2pi],
            ["0", L"$\pi/2$", L"$\pi$", L"$3\pi/2$", L"$2\pi$"]),

	zlabel=L"J"*" (meV)",

	c=cgrad(["#064635","#519259", "#96BB7C", "#F0BB62", "#FAD586","#F4EEA9"]),
    #colorbar = false, st = :surface,

    title = "Thiophene Transfer Integrals "*L"J"*" (meV)",
    titlefont = 24,
    guidefont = 20,  # axes labels
    tickfont  = 16,
)
savefig("new_Thiophene_J_surface.png")
savefig("new_Thiophene_J_heatmap.png")

heatmap(
    φ_list, θ_list, ((J_eff*10e+3)).^2,
    #st = :surface,
    size=(1200, 1200),
    xlabel= L"$\phi$",
	ylabel=L"$\theta$",

	yticks=([0,pi/2,pi],
            ["0", L"$\pi/2$", L"$\pi$"]),
	xticks=([0,pi/2,pi,3pi/2,2pi],
            ["0", L"$\pi/2$", L"$\pi$", L"$3\pi/2$", L"$2\pi$"]),

	zlabel=L"J^2"*" (meV)"*L"^2",

	c=cgrad(["#064635","#519259", "#96BB7C", "#F0BB62", "#FAD586","#F4EEA9"]),
    #colorbar = false, st = :surface,

    title = "Squared Thiophene Transfer Integrals "*L"J^2"*" (meV)"*L"^2",
    titlefont = 24,
    guidefont = 20,  # axes labels
    tickfont  = 16,
)
savefig("new_Thiophene_J^2_surface.png")
savefig("new_Thiophene_J^2_heatmap.png")

heatmap(
    φ_list, θ_list, log.((J_eff).^2),
    #st = :surface,
    size=(1200, 1200),
    xlabel= L"$\phi$",
	ylabel=L"$\theta$",

	yticks=([0, pi/2, pi],
            ["0", L"$\pi/2$", L"$\pi$"]),
	xticks=([0,pi/2,pi,3pi/2,2pi],
            ["0", L"$\pi/2$", L"$\pi$", L"$3\pi/2$", L"$2\pi$"]),

	zlabel=L"\log(\mathbf{J}^2)",
    #zticks=([-40, -20, 0, 10, 20], ["-40", "-20", "0", "10", "20"]),

	c=cgrad(["#064635","#519259", "#96BB7C", "#F0BB62", "#FAD586","#F4EEA9"]),
    #colorbar = false, st = :surface,

    title = "Log Squared Thiophene Transfer Integrals "*L"\log(\mathbf{J}^2)",
    titlefont = 24,
    guidefont = 20,  # axes labels
    tickfont  = 16,
)
savefig("new_Thiophene_logJ^2_surface.png")
savefig("new_Thiophene_logJ^2_heatmap.png")



heatmap(
    φ_list, θ_list, J_eff*10e+3,
    #st = :surface,
    size=(900, 960),
    xlabel= L"$\phi$",
	ylabel=L"$\theta$",
	yticks=([0,pi/2,pi],
            ["0", L"$\pi/2$", L"$\pi$"]),
	xticks=([0,pi/2,pi,3pi/2,2pi],
            ["0", L"$\pi/2$", L"$\pi$", L"$3\pi/2$", L"$2\pi$"]),
	zlabel=L"J",
	c=cgrad(["#064635","#519259", "#96BB7C", "#F0BB62", "#FAD586","#F4EEA9"]),
    title = "Ethylene Transfer Integrals "*L"\mathbf{J}"*" (meV)",
    titlefont = 20,
    guidefont = 18,  # axes labels
    tickfont  = 14,
    #legendtitle = L"J"
    #colorbar = false,
    #colorbar_width = 0.4     # relative height scaling (default 1.0)
)
savefig("ethylene_J_HM.png")

# Load the structures
atoms_list = read_extxyz("output_r8_t8_p8_rx5_ry5_rz5.extxyz")
# Prepare arrays to store (r, θ, energy_eff) data
rs = Float64[]
thetas = Float64[]
energies = Float64[]