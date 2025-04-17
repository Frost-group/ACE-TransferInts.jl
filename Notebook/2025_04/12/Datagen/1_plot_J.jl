using Pkg
cd("C:/Users/Keerati/Documents/GitHub/ACE-TransferInts.jl/Notebook/2025_04/12/Datagen")
Pkg.activate(".")
Base.active_project()
Pkg.instantiate()
Pkg.add(["JuLIP", "StaticArrays", "LinearAlgebra", "Plots","LaTeXStrings","FFTW"])

using JuLIP
using StaticArrays
using LinearAlgebra
using Plots
using LaTeXStrings
using FFTW

dimer_set = read_extxyz("output.extxyz")

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
plot(
    φ_list, θ_list, (J_eff*192.8),
    st = :surface,
    size=(700, 700),
    xlabel= L"$\phi$",
	ylabel=L"$\theta$",
	yticks=([0,pi/2,pi],
            ["0", L"$\pi/2$", L"$\pi$"]),
	xticks=([0,pi/2,pi,3pi/2,2pi],
            ["0", L"$\pi/2$", L"$\pi$", L"$3\pi/2$", L"$2\pi$"]),
	zlabel=L"J",
	c=cgrad(["#064635","#519259", "#96BB7C", "#F0BB62", "#FAD586","#F4EEA9"]),
    title = "Surface Plot of " * L"\mathbf{J}~(\theta, \phi)"
)
savefig("DFethylene_J_surface.png")

plot(
    φ_list, θ_list, ((J_eff*192.8)).^2,
    st = :surface,
    size=(700, 700),
    xlabel= L"$\phi$",
	ylabel=L"$\theta$",
	yticks=([0,pi/2,pi],
            ["0", L"$\pi/2$", L"$\pi$"]),
	xticks=([0,pi/2,pi,3pi/2,2pi],
            ["0", L"$\pi/2$", L"$\pi$", L"$3\pi/2$", L"$2\pi$"]),
	zlabel=L"J",
	c=cgrad(["#064635","#519259", "#96BB7C", "#F0BB62", "#FAD586","#F4EEA9"]),
    title = "Surface Plot of " * L"\mathbf{J}~(\theta, \phi)"
)
savefig("DFethylene_J^2_surface.png")

plot(
    φ_list, θ_list, log.((J_eff*192.8).^2),
    st = :surface,
    size=(700, 700),
    xlabel= L"$\phi$",
	ylabel=L"$\theta$",
	yticks=([0,pi/2,pi],
            ["0", L"$\pi/2$", L"$\pi$"]),
	xticks=([0,pi/2,pi,3pi/2,2pi],
            ["0", L"$\pi/2$", L"$\pi$", L"$3\pi/2$", L"$2\pi$"]),
	zlabel=L"J",
	c=cgrad(["#064635","#519259", "#96BB7C", "#F0BB62", "#FAD586","#F4EEA9"]),
    title = "Surface Plot of " * L"\mathbf{J}~(\theta, \phi)"
)
savefig("DFethylene_logJ^2_surface.png")





heatmap(
    φ_list, θ_list, J_eff,
    size=(700, 700),
    xlabel= L"$\phi$",
	ylabel=L"$\theta$",
	yticks=([0,pi/2,pi],
            ["0", L"$\pi/2$", L"$\pi$"]),
	xticks=([0,pi/2,pi,3pi/2,2pi],
            ["0", L"$\pi/2$", L"$\pi$", L"$3\pi/2$", L"$2\pi$"]),
	zlabel=L"J",
	c=cgrad(["#064635","#519259", "#96BB7C", "#F0BB62", "#FAD586","#F4EEA9"]),
    title = "Surface Plot of " * L"\mathbf{J}_{eff}(\theta, \phi)"
)


0.1 .+ J_eff

fft(J_eff)

plot(
    φ_list, θ_list, abs.(fft(J_eff)),
    st = :surface,
    size=(700, 700),
    xlabel= L"$\phi$",
	ylabel=L"$\theta$",
	yticks=([0,pi/2,pi],
            ["0", L"$\pi/2$", L"$\pi$"]),
	xticks=([0,pi/2,pi,3pi/2,2pi],
            ["0", L"$\pi/2$", L"$\pi$", L"$3\pi/2$", L"$2\pi$"]),
	zlabel=L"J",
	c=cgrad(["#064635","#519259", "#96BB7C", "#F0BB62", "#FAD586","#F4EEA9"]),
    title = "Surface Plot of " * L"\mathbf{J}_{eff}(\theta, \phi)"
)
