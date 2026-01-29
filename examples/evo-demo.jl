using EcoEvoSim
using Plots
using DataFrames
using CSV


# Create a simple evolutionary simulation
growthFn = (z) -> (tanh(sum(z .- 0.5) / 0.2) + 1) / 2 - 0.006692851
kernelFn = (zi, zj) -> -(tanh(sum(zi .- zj) / 0.15) + 1) / 2

config = EcoEvoConfig(
    ecoDyn = lotkaVolterra(growthFn, kernelFn),
    mutationGenerator = (c, cfg) -> generateMutant(c, cfg, 0.002^2),
    integrationParams = IntegrationParams(
        maxTime = 1e10,
        abstol = 1e-10,
        reltol = 1e-8
    ),
    invaderPopsize = 0.001,
    extThreshold = 0.003
)

lineage = Community([1.0], [0.3], Float64[])

lineage = evolve!(lineage, config, 1500)

p = plotEvo(lineage)
# savefig(p, "test_plotEvo.png")

# Test historyToTable
table = historyToTable(lineage)
# CSV.write("evo.csv", DataFrame(table))
