using EcoEvoSim
using Plots


growthFn = z -> 1.0 - sum(z.^2) / 0.5^2
kernelFn = (zi, zj) -> -exp(-sum((zi .- zj).^2) / 0.25^2)

config = EcoEvoConfig(
    ecoDyn = lotkaVolterra(growthFn, kernelFn),
    mutationGenerator = (c) -> generateMutant(c, 0.001, 0.002^2),
    integrationParams = IntegrationParams(maxTime = 1.0e8),
    extThreshold = 0.003
)

lineage = Community([Species(1.0, [0.3, -0.3])], PopulationSize{Float64}[], 0.0)
lineage = ecoDyn(lineage, config)
lineage = evolve!(lineage, config, 30000)

plotEvo(lineage)
plotEvo(lineage, traitDim=1)
plotEvo(lineage, traitDim=2)
plotEvoTwoTrait(lineage, camera=(60, 25))
