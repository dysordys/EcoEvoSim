using EcoEvoSim
using Plots
using DifferentialEquations
using Distributions


function patchGrowth(community, d)
    # spTraits is a vector of scalar trait values (one per species)
    # Returns a matrix with species in rows and patches in columns
    spTraits = traits(community)
    nSpecies = numSpecies(community)
    y = [d/2, -d/2]  # Optimal traits for patches 1 and 2
    growth = zeros(nSpecies, 2)
    for i in 1:nSpecies
        growth[i, :] = pdf.(Normal(0, 1), spTraits[i] .- y)
    end
    return growth
end


function densDep(densityMatrix, alpha)
    return alpha .* dropdims(sum(densityMatrix; dims = 1), dims = 1)
end


function twopatch(
        community::Community{T, AuxClasses};
        d::Float64 = 1.2, mu::Float64 = 0.1, alpha::Float64 = 1.0
    ) where {T<:Real, AuxClasses}
    nSpecies = numSpecies(community)
    growth = patchGrowth(community, d)  # species x patches
    function odeFn(u::Vector{T}, p, t) where {T}
        # Reconstruct state into density matrix (species x patches)
        densityMatrix = zeros(T, nSpecies, 2)
        for i in 1:nSpecies
            idx1 = 2 * (i - 1) + 1
            idx2 = 2 * (i - 1) + 2
            densityMatrix[i, 1] = u[idx1]
            densityMatrix[i, 2] = u[idx2]
        end
        # Compute density dependence per patch based on current densities
        dens = densDep(densityMatrix, alpha)
        # Build dynamical matrix with current densities
        nEqs = nSpecies * 2
        M = zeros(T, nEqs, nEqs)
        for i in 1:nSpecies
            idx1 = 2 * (i - 1) + 1
            idx2 = 2 * (i - 1) + 2
            # Diagonal terms: growth minus density dependence
            M[idx1, idx1] = growth[i, 1] - dens[1] - mu
            M[idx2, idx2] = growth[i, 2] - dens[2] - mu
            # Off-diagonal terms: migration
            M[idx1, idx2] = mu
            M[idx2, idx1] = mu
        end
        return M * u
    end
    return odeFn
end


config = EcoEvoConfig(
    ecoDyn = comm -> twopatch(comm; d = 1.0, mu = 0.1, alpha = 1.0),
    mutationGenerator = c -> generateMutant(c; invaderPopsize=0.001, variance=0.003^2),
    integrationParams = IntegrationParams(
        maxTime = Inf,
        algorithm = DynamicSS(),
        abstol = 1e-10,
        reltol = 1e-8
    ),
    extThreshold = 0.003
)


lineage = Community([1.0 1.0;], [-0.2])
lineage = ecoDyn(lineage, config)
@time lineage = evolve!(lineage, config, 1000);

p = plotEvo(lineage)
