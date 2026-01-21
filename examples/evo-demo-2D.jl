using EcoEvoSim
using Plots


growthFn = z -> 1.0 - z^2
kernelFn = (zi, zj) -> -exp(-((zi - zj) / 0.25)^2)

community = Community([Species(1.0, [0.3, -0.1])], PopulationSize{Float64}[], 0.0)

config = EcoEvoConfig(
    ecoDyn = lotkaVolterra(growthFn, kernelFn),
    mutationGenerator = (c, cfg) -> generateMutant(c, cfg, 0.002^2),
    integrationParams = IntegrationParams(maxTime = 1.0e8),
    invaderPopsize = 0.001,
    extThreshold = 0.003
)

history = evolve!(community, config, 15000)

plotEvo(history)
plotEvo(history, traitDim=1)
plotEvo(history, traitDim=2)
plotEvoTwoTrait(history, camera=(60, 25))
