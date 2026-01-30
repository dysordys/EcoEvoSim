using EcoEvoSim
using Plots
using DifferentialEquations
using DataFrames
using CSV
using Random
using Distributions
using StatsBase
using LinearAlgebra


function patchGrowth(spTraits, d::Float64 = 1.5)
    # spTraits is a vector of scalar trait values (one per species)
    # Returns a matrix with species in rows and patches in columns
    nSpecies = length(spTraits)
    y = [d/2, -d/2]  # Optimal traits for patches 1 and 2
    growth = zeros(nSpecies, 2)
    for i in 1:nSpecies
        growth[i, :] = pdf.(Normal(0, 1), spTraits[i] .- y)
    end
    return growth
end


function densDep(community, alpha::Float64 = 1.0)
    densityMatrix = popsizesToMatrix(community)
    return alpha .* dropdims(sum(densityMatrix; dims = 1), dims = 1)
end


function twopatch(community::Community{T, AuxClasses}, mu::Float64 = 0.1) where {T<:Real, AuxClasses}
    nSpecies = numSpecies(community)
    # Extract traits (one per species, scalar)
    spTraits = [traits(community, i)[1] for i in 1:nSpecies]
    # Compute growth rates for each patch (based on optimal trait distances)
    growth = patchGrowth(spTraits)  # nSpecies x 2 matrix (species x patches)
    # Compute density dependence per patch
    dens = densDep(community)  # 2-element vector (one per patch)
    # Build the migration matrix (2 patches x 2 patches per species)
    # This is the intraspecific migration for a single species
    migrationMat = [-mu mu; mu -mu]
    # Build full dynamical matrix
    # We have nSpecies * 2 equations (each species in 2 patches)
    nEqs = nSpecies * 2
    A = zeros(T, nEqs, nEqs)
    # For each species, construct the block
    for i in 1:nSpecies
        # Row indices for this species (patches 1 and 2)
        idx1 = 2 * (i - 1) + 1  # Patch 1
        idx2 = 2 * (i - 1) + 2  # Patch 2

        # Diagonal terms: growth minus density dependence
        A[idx1, idx1] = growth[i, 1] - dens[1]
        A[idx2, idx2] = growth[i, 2] - dens[2]

        # Off-diagonal terms: migration
        A[idx1, idx2] = migrationMat[1, 2]  # migration from patch 2 to patch 1
        A[idx2, idx1] = migrationMat[2, 1]  # migration from patch 1 to patch 2
    end
    return (u, p, t) -> A * u
end


function generateMutantTwoPatch(community::Community{T, AuxClasses}, config, variance) where {T<:Real, AuxClasses}
    # Generate a mutant with 2 stage classes (one for each patch)
    # Select parent species with probability proportional to total population size
    nParents = numSpecies(community)
    nParents > 0 || throw(ArgumentError("Cannot generate mutant from empty community"))

    # Get total population sizes for each species (sum across patches)
    popMatrix = popsizesToMatrix(community)
    totalPops = dropdims(sum(popMatrix; dims = 2), dims = 2)

    # Select parent with probability proportional to population size
    probabilities = totalPops ./ sum(totalPops)
    parentIdx = sample(1:nParents, Weights(probabilities))
    parentTrait = traits(community, parentIdx)

    # Generate new trait by adding Gaussian noise
    newTrait = parentTrait .+ randn(length(parentTrait)) .* sqrt(variance)

    # Create mutant with 2 stage classes (split invader popsize equally between patches)
    mutantPopsize = config.invaderPopsize / 2  # Equal distribution between patches

    # Create Species with 2 stage classes using the matrix constructor
    mutantSpecies = Species{T}(
        [PopulationSize{T}([mutantPopsize, mutantPopsize])],
        [Phenotype{T}(newTrait)]
    )

    # Create new community with mutant appended (like generateMutant does)
    newSpeciesList = vcat(speciesList(community), mutantSpecies)
    Community(newSpeciesList, community.aux, community.time)
end


config = EcoEvoConfig(
    ecoDyn = twopatch,
    mutationGenerator = (c, cfg) -> generateMutantTwoPatch(c, cfg, 0.002^2),
    integrationParams = IntegrationParams(
        maxTime = Inf,
        algorithm = DynamicSS(),
        abstol = 1e-14,
        reltol = 1e-8
    ),
    invaderPopsize = 0.01,
    extThreshold = 0.0001
)

lineage = Community([1.0 1.0], [0], Float64[])

@time lineage = evolve!(lineage, config, 8000);

p = plotEvo(lineage)


μ = 0.0
σ = 1.0
dist = Normal(μ, σ)

x = range(-4, 4; length = 400)   # points for a smooth curve
y = pdf.(dist, x)                # broadcast pdf over the vector

plot(x, y,
     label = "N(0,1) PDF",
     linewidth = 2,
     color = :steelblue,
     title = "Standard Normal Distribution",
     xlabel = "x",
     ylabel = "Density"
)

