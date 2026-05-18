# Two-patch model with explicit resources using structuredModel helper.
# Compare to two-patch-resource.jl - much simpler.

using EcoEvoSim
using Plots
using OrdinaryDiffEq
using Distributions
using Random


d = 1.0; mu = 0.1
eta = 1.0; chi = 1.0; gamma = 1.0
beta = 1.0; theta = 1.0; sigma = 1.0
y = [d/2, -d/2]
m_base = (beta * eta) / chi

ecology = structuredModel(
    auxDynamics = (aux, n, z, nSpecies, nPatches) ->
        [aux[k] * (eta - chi * aux[k]) - gamma * sum(n[i, k] for i in 1:nSpecies) * aux[k]
         for k in 1:nPatches]
) do i, j, n, z, aux, nSpecies, nPatches
    mort = m_base - theta * pdf(Normal(y[j], sigma), z[i][1])
    (beta * aux[j] - mort - mu) * n[i, j] + mu * sum(n[i, k] for k in 1:nPatches if k != j)
end


config = EcoEvoConfig(
    ecoDyn = ecology,
    mutationGenerator =
        generateMutantSpatial(invaderPopsize = 0.001, variance = 0.003^2),
    integrationParams = IntegrationParams(
        maxTime = Inf,
        algorithm = DynamicSS(),
        abstol = 1e-10,
        reltol = 1e-10
    ),
    extThreshold = 0.003
)


Random.seed!(54321)

lineage = Community([1.0 1.0;], [-0.2], [1.0, 1.0])
lineage = ecoDyn(lineage, config)
@time lineage = evolve(lineage, config, 1000);

p = plotEvo(lineage)
