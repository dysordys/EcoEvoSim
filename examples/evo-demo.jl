using EcoEvoSim
using Plots


# Create a simple evolutionary simulation
growthFn = z -> z^2 / (3 + z^2) #1.0 - z^2
kernelFn = (zi, zj) -> -(tanh((zi - zj) / 0.3) + 1) / 2 #-exp(-((zi - zj) / 0.25)^2)

community = Community([1.0], [0.3], Float64[])

config = EcoEvoConfig(
    ecoDyn = lotkaVolterra(growthFn, kernelFn),
    mutationGenerator = (c, cfg) -> generateMutant(c, cfg, 0.002^2),
    integrationParams = IntegrationParams(maxTime = 1.0e8),
    invaderPopsize = 0.001,
    extThreshold = 0.003
)

history = evolve!(community, config, 15000)

p = plotEvo(history)
savefig(p, "test_plotEvo.png")
