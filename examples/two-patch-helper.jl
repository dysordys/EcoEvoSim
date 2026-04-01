# Two-patch model using structuredModel helper.
# Compare to two-patch.jl - no manual index arithmetic needed.

using EcoEvoSim
using Plots
using DifferentialEquations
using Distributions
using Random


d = 1.0; mu = 0.1; alpha = 1.0
y = [d/2, -d/2]

ecology = structuredModel() do i, j, N, z, R, nSpecies, nPatches
    growth = pdf(Normal(0, 1), z[i][1] - y[j])
    dd = alpha * sum(N[k, j] for k in 1:nSpecies) # density dependence
    (growth - dd - mu) * N[i, j] + mu * sum(N[i, k] for k in 1:nPatches if k != j)
end


config = EcoEvoConfig(
    ecoDyn = ecology,
    mutationGenerator =
        c -> generateMutantSpatial(c; invaderPopsize=0.001, variance=0.003^2),
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
