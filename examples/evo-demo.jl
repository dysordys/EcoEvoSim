using EcoEvoSim
using Plots
using DataFrames
using CSV


# Create a simple evolutionary simulation
growthFn = z -> sum(z.^2) / (3 + sum(z.^2)) #1.0 - sum(z.^2)
kernelFn = (zi, zj) -> -(tanh(sum(zi .- zj) / 0.3) + 1) / 2 #-exp(-sum((zi .- zj).^2) / 0.25^2)

community = Community([1.0], [0.3], Float64[])

config = EcoEvoConfig(
    ecoDyn = lotkaVolterra(growthFn, kernelFn),
    mutationGenerator = (c, cfg) -> generateMutant(c, cfg, 0.002^2),
    integrationParams = IntegrationParams(maxTime = 1.0e10),
    invaderPopsize = 0.001,
    extThreshold = 0.003
)

history = evolve!(community, config, 1000)

# Resume the simulation - both styles work:
# 1. Without reassignment (history is modified in place):
#    evolve!(history, config, 500)
# 2. With reassignment (returns the modified history):
history = evolve!(history, config, 500)

p = plotEvo(history)
# savefig(p, "test_plotEvo.png")

# Test historyToTable
# table = historyToTable(history)
# CSV.write("evo.csv", DataFrame(table))
