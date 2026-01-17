using EcoEvoSim
using Plots


# Create a simple evolutionary simulation
growthFn = (z) -> @. (tanh((z - 0.5) / 0.2) + 1) / 2 - 0.006692851 # traits -> 1.0 .- traits.^2
interactionFn = (z_i, z_j) -> (tanh((z_i - z_j) / 0.1) + 1) / 2 #exp(-((z_i - z_j) / 0.15)^2)

community = Community([1.0], [0.3], Float64[])

config = EcoEvoConfig(
    ecoDyn = lotkaVolterra(growthFn, interactionFn),
    mutationGenerator = (c, cfg) -> generateMutant(c, cfg, 4.0e-6),
    integrationParams = IntegrationParams(maxTime = 100000.0),
    invaderPopsize = 0.001,
    extThreshold = 0.003
)

history = evolve!(community, config, 15000)

p = plotEvo(history)
savefig(p, "test_plotEvo.png")
