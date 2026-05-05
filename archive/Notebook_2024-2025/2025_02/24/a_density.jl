using Plots

# Define the Gaussian density function
function density_function(x, y, neighbors, sigma=1.0)
    density = zeros(size(x))
    for (xj, yj) in neighbors
        density .+= exp.(-((x .- xj).^2 .+ (y .- yj).^2) / (2 * sigma^2))
    end
    return density
end

# Define the grid
x = range(-3, stop=3, length=100)
y = range(-3, stop=3, length=100)
X, Y = [x for x in x, y in y], [y for x in x, y in y]  # Meshgrid

# Define neighboring points
neighbors = [(-1, -1), (1, -1), (-1, 1), (1, 1)]  # Example positions of neighbors

# Compute density
Z = density_function(X, Y, neighbors, sigma=1.0)

# Plot the density function as a 3D surface
plot(
    X, Y, Z,
    st=:surface,
    xlabel="X",
    ylabel="Y",
    zlabel=L"\rho_i(\mathbf{r})",
    title="Density Function " * L"\mathbf{d}_i(\mathbf{r})",
    c=:viridis
)
