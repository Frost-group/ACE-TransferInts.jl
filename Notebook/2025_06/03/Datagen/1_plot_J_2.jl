
dimer_set = read_extxyz("output625_rot1.extxyz")

θ_list = range(0.01, stop=π, length=17);
φ_list = range(0.01, stop=2π, length=33);

J_eff = zeros(size(θ_list,1),size(φ_list,1))
get_data(dimer_set[1], "theta")

for i in range(1,size(θ_list,1))
    for j in range(1,size(φ_list,1))

    J_eff[i,j] = dimer_set[ (i-1) * 33 + j]["energy_eff"]

    end
end    
gr(dpi=1200) 
plot(
    φ_list, θ_list, ( J_eff)
    ,
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
savefig("ethylene_J_surface625_rot1.png")


plot(
    φ_list, θ_list, ((J_eff)).^2,
    st = :surface,
    size=(700, 700),
    xlabel= L"$\phi$",
	ylabel=L"$\theta$",
	yticks=([0,pi/2,pi],
            ["0", L"$\pi/2$", L"$\pi$"]),
	xticks=([0,pi/2,pi,3pi/2,2pi],
            ["0", L"$\pi/2$", L"$\pi$", L"$3\pi/2$", L"$2\pi$"]),
	zlabel=L"J^2",
	c=cgrad(["#064635","#519259", "#96BB7C", "#F0BB62", "#FAD586","#F4EEA9"]),
    title = "Surface Plot of " * L"\mathbf{J}^2~(\theta, \phi)"
)
savefig("ethylene_J^2_surface70.png")

plot(
    φ_list, θ_list, log.((J_eff).^2),
    st = :surface,
    size=(700, 700),
    xlabel= L"$\phi$",
	ylabel=L"$\theta$",
	yticks=([0,pi/2,pi],
            ["0", L"$\pi/2$", L"$\pi$"]),
	xticks=([0,pi/2,pi,3pi/2,2pi],
            ["0", L"$\pi/2$", L"$\pi$", L"$3\pi/2$", L"$2\pi$"]),
	zlabel=L"\log(\mathbf{J}^2)",
	c=cgrad(["#064635","#519259", "#96BB7C", "#F0BB62", "#FAD586","#F4EEA9"]),
    title = "Surface Plot of " * L"\log(\mathbf{J}^2)~(\theta, \phi)"
)
savefig("ethylene_logJ^2_surface70.png")


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

# Load the structures
atoms_list = read_extxyz("output_r8_t8_p8_rx5_ry5_rz5_2.extxyz")
# Prepare arrays to store (r, θ, energy_eff) data

# Define your grid values explicitly to fix ordering
Nr, Nθ, Nϕ, Nx, Ny, Nz = 8, 8, 8, 5, 5, 5
r_vals  = collect(range(6.25, stop=7.00, length=Nr))
θ_vals  = collect(range(0.1, stop=π, length=Nθ))
ϕ_vals  = collect(range(0.1, stop=2π, length=Nϕ))
θx_vals = collect(range(0.1, stop=π, length=Nx))
θy_vals = collect(range(0.1, stop=π, length=Ny))
θz_vals = collect(range(π/6, stop=5*π/6, length=Nz))
#θx_vals = collect(range(0.1, stop=2π, length=Nx))
#θy_vals = collect(range(0.1, stop=2π, length=Ny))
#θz_vals = collect(range(0.1, stop=2π, length=Nz))

# Create mapping from value to index (for precise placement)
function make_index_map(values)
    return Dict(v => i for (i, v) in enumerate(values))
end

r_idx   = make_index_map(r_vals)
θ_idx   = make_index_map(θ_vals)
ϕ_idx   = make_index_map(ϕ_vals)
θx_idx  = make_index_map(θx_vals)
θy_idx  = make_index_map(θy_vals)
θz_idx  = make_index_map(θz_vals)

# Initialize tensor
energy_tensor = fill(NaN, Nr, Nθ, Nϕ,Nx, Ny, Nz);
function closest_index(x, grid; atol=1e-6)
    diffs = abs.(grid .- x)
    minval, idx = findmin(diffs)
    if minval > atol
        error("No close match found for $x (min diff = $minval)")
    end
    return idx
end
# Fill tensor from atoms list
for atoms in atoms_list
    r = get_data(atoms, "r")
    θ = get_data(atoms, "theta")
    ϕ = get_data(atoms, "phi")
    θx = get_data(atoms, "θx")
    θy = get_data(atoms, "θy")
    θz = get_data(atoms, "θz")
    eff   = get_data(atoms, "energy_eff")

    # Use mapping dictionaries to get indices
    # Use approximate matching instead of exact dictionary lookup
    i = closest_index(r, r_vals)
    j = closest_index(θ, θ_vals)
    k = closest_index(ϕ, ϕ_vals)
    l = closest_index(θx, θx_vals)
    m = closest_index(θy, θy_vals)
    n = closest_index(θz, θz_vals)

    energy_tensor[i, j, k, l, m, n] = eff
end

@show size(energy_tensor)
energy_tensor[:, 1, 1, 1, 1, :]

#######################################################
plot(
    θz_vals, r_vals, energy_tensor[:, 1, 1, 1, 1, :],
    st = :surface,
    size=(700, 700),
    xlabel= L"$\theta_{z} (deg)$",
	ylabel=L"$r$ (Å)",
	yticks=([6.25, 6.5, 6.75, 7.00],
            ["6.25", "6.5", "6.75", "7.00"]),
	xticks=([pi/6, pi/3, pi/2, 2pi/3, 5pi/6, pi],
            [L"$\pi/6$", L"$\pi/3$", L"$\pi/2$", L"$2\pi/3$", L"$5\pi/6$"]),
	zlabel=L"J",
	c=cgrad(["#064635","#519259", "#96BB7C", "#F0BB62", "#FAD586","#F4EEA9"]),
    title = "Surface Plot of " * L"\mathbf{J}~(\theta_{z}, r)"
)
savefig("ethylene_J_r_thetaZ.png")

anim_rz_θ = @animate for t in 1:8
    plot(
    θz_vals, r_vals, energy_tensor[:, t, 1, 1, 1, :],
    st = :surface,
    size=(700, 700),
    xlabel= L"$\theta_{z} (rad)$",
	ylabel=L"$r$ (Å)",
	yticks=([6.25, 6.5, 6.75, 7.00],
            ["6.25", "6.5", "6.75", "7.00"]),
	xticks=([pi/6, pi/3, pi/2, 2pi/3, 5pi/6, pi],
            [L"$\pi/6$", L"$\pi/3$", L"$\pi/2$", L"$2\pi/3$", L"$5\pi/6$"]),
	zlabel=L"J",
	c=cgrad(["#064635","#519259", "#96BB7C", "#F0BB62", "#FAD586","#F4EEA9"]),
    title ="Surface Plot of " * L"\mathbf{J}~(\theta_{z}, r)~;"* L"\theta="*string(round(θ_vals[t]/π, digits=2))*L"\pi"
    )
end
gif(anim_rz_θ, "thylene_J_r_thetaZ_theta_2.gif", fps=2)

anim_rz_θ_hm = @animate for t in 1:8
    heatmap(
    θz_vals, r_vals, energy_tensor[:, t, 1, 1, 1, :],
    size=(700, 700),
    xlabel= L"$\theta_{z} (rad)$",
	ylabel=L"$r$ (Å)",
	yticks=([6.25, 6.5, 6.75, 7.00],
            ["6.25", "6.5", "6.75", "7.00"]),
	xticks=([pi/6, pi/3, pi/2, 2pi/3, 5pi/6, pi],
            [L"$\pi/6$", L"$\pi/3$", L"$\pi/2$", L"$2\pi/3$", L"$5\pi/6$"]),
	zlabel=L"J",
	c=cgrad(["#064635","#519259", "#96BB7C", "#F0BB62", "#FAD586","#F4EEA9"]),
    title = "Surface Plot of " * L"\mathbf{J}~(\theta_{z}, r)~;"* L"\theta="*string(round(θ_vals[t]/π, digits=2))*L"\pi"
    )
end
gif(anim_rz_θ_hm, "thylene_J_r_thetaZ_theta_hm_2.gif", fps=2)

#######################################################

anim_rz_ϕ = @animate for t in 1:8
    plot(
    θz_vals, r_vals, energy_tensor[:, 1, t, 1, 1, :],
    st = :surface,
    size=(700, 700),
    xlabel= L"$\theta_{z} (rad)$",
	ylabel=L"$r$ (Å)",
	yticks=([6.25, 6.5, 6.75, 7.00],
            ["6.25", "6.5", "6.75", "7.00"]),
	xticks=([pi/6, pi/3, pi/2, 2pi/3, 5pi/6, pi],
            [L"$\pi/6$", L"$\pi/3$", L"$\pi/2$", L"$2\pi/3$", L"$5\pi/6$"]),
	zlabel=L"J",
	c=cgrad(["#064635","#519259", "#96BB7C", "#F0BB62", "#FAD586","#F4EEA9"]),
    title = "Surface Plot of " * L"\mathbf{J}~(\theta_{z}, r)~;"* L"\phi="*string(round(ϕ_vals[t]/π, digits=2))*L"\pi"
    )
end
gif(anim_rz_ϕ, "thylene_J_r_thetaZ_phi_2.gif", fps=2)

anim_rz_ϕ_hm = @animate for t in 1:8
    heatmap(
    θz_vals, r_vals, energy_tensor[:, 1, t, 1, 1, :],
    size=(700, 700),
    xlabel= L"$\theta_{z} (rad)$",
	ylabel=L"$r$ (Å)",
	yticks=([6.25, 6.5, 6.75, 7.00],
            ["6.25", "6.5", "6.75", "7.00"]),
	xticks=([pi/6, pi/3, pi/2, 2pi/3, 5pi/6, pi],
            [L"$\pi/6$", L"$\pi/3$", L"$\pi/2$", L"$2\pi/3$", L"$5\pi/6$"]),
	zlabel=L"J",
	c=cgrad(["#064635","#519259", "#96BB7C", "#F0BB62", "#FAD586","#F4EEA9"]),
    title = "Surface Plot of " * L"\mathbf{J}~(\theta_{z}, r)~;"* L"\phi="*string(round(ϕ_vals[t]/π, digits=2))*L"\pi"
    )
end
gif(anim_rz_ϕ_hm, "thylene_J_r_thetaZ_phi_hm_2.gif", fps=2)

anim_rz_ϕ = @animate for t in 1:8
    plot(
    θz_vals, r_vals, energy_tensor[:, 1, t, 1, 1, :],
    st = :surface,
    size=(700, 700),
    xlabel= L"$\theta_{z} (rad)$",
	ylabel=L"$r$ (Å)",
	yticks=([6.25, 6.5, 6.75, 7.00],
            ["6.25", "6.5", "6.75", "7.00"]),
	xticks=([pi/6, pi/3, pi/2, 2pi/3, 5pi/6, pi],
            [L"$\pi/6$", L"$\pi/3$", L"$\pi/2$", L"$2\pi/3$", L"$5\pi/6$"]),
	zlabel=L"J",
	c=cgrad(["#064635","#519259", "#96BB7C", "#F0BB62", "#FAD586","#F4EEA9"]),
    title = "Surface Plot of " * L"\mathbf{J}~(\theta_{z}, r)~;"* L"\phi="*string(round(ϕ_vals[t]/π, digits=2))*L"\pi"
    )
    end
    gif(anim_rz_ϕ, "thylene_J_r_thetaZ_phi_2.gif", fps=2)

    #######################################################

    anim_rz_x = @animate for t in 1:5
        plot(
        θz_vals, r_vals, energy_tensor[:, 1, 1, t, 1, :],
        st = :surface,
        size=(700, 700),
        xlabel= L"$\theta_{z} (rad)$",
        ylabel=L"$r$ (Å)",
        yticks=([6.25, 6.5, 6.75, 7.00],
                ["6.25", "6.5", "6.75", "7.00"]),
        xticks=([pi/6, pi/3, pi/2, 2pi/3, 5pi/6, pi],
                [L"$\pi/6$", L"$\pi/3$", L"$\pi/2$", L"$2\pi/3$", L"$5\pi/6$"]),
        zlabel=L"J",
        c=cgrad(["#064635","#519259", "#96BB7C", "#F0BB62", "#FAD586","#F4EEA9"]),
        title = "Surface Plot of " * L"\mathbf{J}~(\theta_{z}, r)~;"* L"\theta_{x}="*string(round(θx_vals[t]/π, digits=2))*L"\pi"
        )
    end
    gif(anim_rz_x, "thylene_J_r_thetaZ_x_2.gif", fps=2)

    anim_rz_x_hm = @animate for t in 1:5
        heatmap(
        θz_vals, r_vals, energy_tensor[:, 1, t, 1, 1, :],
        size=(700, 700),
        xlabel= L"$\theta_{z} (rad)$",
        ylabel=L"$r$ (Å)",
        yticks=([6.25, 6.5, 6.75, 7.00],
                ["6.25", "6.5", "6.75", "7.00"]),
        xticks=([pi/6, pi/3, pi/2, 2pi/3, 5pi/6, pi],
                [L"$\pi/6$", L"$\pi/3$", L"$\pi/2$", L"$2\pi/3$", L"$5\pi/6$"]),
        zlabel=L"J",
        c=cgrad(["#064635","#519259", "#96BB7C", "#F0BB62", "#FAD586","#F4EEA9"]),
        title = "Surface Plot of " * L"\mathbf{J}~(\theta_{z}, r)~;"* L"\theta_{x}="*string(round(θx_vals[t]/π, digits=2))*L"\pi"
        )
    end
    gif(anim_rz_x_hm, "thylene_J_r_thetaZ_x_hm_2.gif", fps=2)

    #######################################################

    anim_rz_y = @animate for t in 1:5
        plot(
        θz_vals, r_vals, energy_tensor[:, 1, 1, t, 1, :],
        st = :surface,
        size=(700, 700),
        xlabel= L"$\theta_{z} (rad)$",
        ylabel=L"$r$ (Å)",
        yticks=([6.25, 6.5, 6.75, 7.00],
                ["6.25", "6.5", "6.75", "7.00"]),
        xticks=([pi/6, pi/3, pi/2, 2pi/3, 5pi/6, pi],
                [L"$\pi/6$", L"$\pi/3$", L"$\pi/2$", L"$2\pi/3$", L"$5\pi/6$"]),
        zlabel=L"J",
        c=cgrad(["#064635","#519259", "#96BB7C", "#F0BB62", "#FAD586","#F4EEA9"]),
        title = "Surface Plot of " * L"\mathbf{J}~(\theta_{z}, r)~;"* L"\theta_{y}="*string(round(θy_vals[t]/π, digits=2))*L"\pi"
        )
    end
    gif(anim_rz_y, "thylene_J_r_thetaZ_y_2.gif", fps=2)

    anim_rz_y_hm = @animate for t in 1:5
        heatmap(
        θz_vals, r_vals, energy_tensor[:, 1, t, 1, 1, :],
        size=(700, 700),
        xlabel= L"$\theta_{z} (rad)$",
        ylabel=L"$r$ (Å)",
        yticks=([6.25, 6.5, 6.75, 7.00],
                ["6.25", "6.5", "6.75", "7.00"]),
        xticks=([pi/6, pi/3, pi/2, 2pi/3, 5pi/6, pi],
                [L"$\pi/6$", L"$\pi/3$", L"$\pi/2$", L"$2\pi/3$", L"$5\pi/6$"]),
        zlabel=L"J",
        c=cgrad(["#064635","#519259", "#96BB7C", "#F0BB62", "#FAD586","#F4EEA9"]),
        title = "Surface Plot of " * L"\mathbf{J}~(\theta_{z}, r)~;"* L"\theta_{y}="*string(round(θy_vals[t]/π, digits=2))*L"\pi"
        )
    end
    gif(anim_rz_y_hm, "thylene_J_r_thetaZ_y_hm_2.gif", fps=2)