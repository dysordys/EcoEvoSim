using EcoEvoSim
using Plots
using DifferentialEquations
using Distributions
using Random


function patchGrowth(community, d)
    # spTraits is a vector of scalar trait values (one per species)
    # Returns a matrix with species in rows and patches in columns
    spTraits = traits(community)
    nSpecies = numSpecies(community)
    y = [d/2, -d/2]  # Optimal traits for each patch
    growth = zeros(nSpecies, length(y))
    for i in 1:nSpecies
        growth[i, :] = pdf.(Normal(0, 1), spTraits[i] .- y)
    end
    return growth
end


function densDep(densityMatrix, alpha)
    return alpha .* dropdims(sum(densityMatrix; dims = 1), dims = 1)
end


function multipatch(
        community::Community{T, AuxClasses};
        d::Float64 = 1.2, mu::Float64 = 0.1, alpha::Float64 = 1.0
    ) where {T<:Real, AuxClasses}
    nSpecies = numSpecies(community)
    growth = patchGrowth(community, d)  # species x patches
    nPatches = size(growth, 2)
    function odeFn(u::Vector{T}, p, t) where {T}
        # Reconstruct state into density matrix (species x patches)
        densityMatrix = zeros(T, nSpecies, nPatches)
        for i in 1:nSpecies
            for j in 1:nPatches
                idx = nPatches * (i - 1) + j
                densityMatrix[i, j] = u[idx]
            end
        end
        # Compute density dependence per patch based on current densities
        dens = densDep(densityMatrix, alpha)
        # Build dynamical matrix with current densities
        nEqs = nSpecies * nPatches
        M = zeros(T, nEqs, nEqs)
        for i in 1:nSpecies
            for j in 1:nPatches
                idx_j = nPatches * (i - 1) + j
                # Diagonal terms: growth minus density dependence
                M[idx_j, idx_j] = growth[i, j] - dens[j] - mu
                # Off-diagonal terms: migration to other patches
                for k in 1:nPatches
                    if k != j
                        idx_k = nPatches * (i - 1) + k
                        M[idx_j, idx_k] = mu
                    end
                end
            end
        end
        return M * u
    end
    return odeFn
end


config = EcoEvoConfig(
    ecoDyn = comm -> multipatch(comm; d = 1.0, mu = 0.1, alpha = 1.0),
    mutationGenerator = c -> generateMutant(c; invaderPopsize=0.001, variance=0.003^2),
    integrationParams = IntegrationParams(
        maxTime = Inf,
        algorithm = DynamicSS(),
        abstol = 1e-10,
        reltol = 1e-8
    ),
    extThreshold = 0.003
)


Random.seed!(54321)

lineage = Community([1.0 1.0;], [-0.2])
lineage = ecoDyn(lineage, config)
@time lineage = evolve!(lineage, config, 1000);

p = plotEvo(lineage)
