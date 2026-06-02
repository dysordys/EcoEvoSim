# Two-patch model using structuredModel helper.
# Compare to two-patch.jl - no manual index arithmetic needed.

using EcoEvoSim
using Plots
using Random


d = 1.0
y = [d/2, -d/2]
mu = 0.1

ecology = structuredModel() do i, k, n, z, aux, nSpecies, nPatches
    localGrowth = 1 - (z[i][1] - y[k])^2 - sum(n[i,k] for i in 1:nSpecies)
    l = k == 1 ? 2 : 1 # Non-focal patch index
    (localGrowth - mu) * n[i,k] + mu * n[i,l]
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

lineage = Community([1.0 1.0;], [-0.2])
lineage = ecoDyn(lineage, config)
@time lineage = evolve(lineage, config, 1500);

p = plotEvo(lineage)
