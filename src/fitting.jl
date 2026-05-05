"""
Model assessment utilities: energy extraction, RMSE, MAE, sign accuracy.
"""

"""
    extract_energies(dataset) -> Vector

Pull the "energy" field (per atom) from every frame in `dataset`.
"""
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

"""
    rmse(predicted, reference) -> Float64
"""
function rmse(predicted, reference)
    return sqrt(sum((predicted .- reference).^2) / length(reference))
end

"""
    mae(predicted, reference) -> Float64
"""
function mae(predicted, reference)
    return sum(abs.(predicted .- reference)) / length(reference)
end

"""
    sign_accuracy(predicted, reference) -> Float64

Fraction of structures whose predicted coupling sign matches the reference.
"""
function sign_accuracy(predicted, reference)
    return mean(sign.(predicted) .== sign.(reference))
end
