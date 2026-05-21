using EcoEvoSim
using Plots
using Random


growthFn(z) = 1 - sum(z.^2) / 0.5^2
kernelFn(zi, zj) = -exp(-sum((zi .- zj).^2) / 0.15^2)


config = EcoEvoConfig(
    ecoDyn = lotkaVolterra(growthFn, kernelFn),
    mutationGenerator = generateMutant(invaderPopsize = 0.001, variance = 0.002^2),
    integrationParams = IntegrationParams(
        maxTime = 1e12,
        abstol = 1e-8,
        reltol = 1e-8
    ),
    extThreshold = 0.003
)


Random.seed!(54321)

lineage = Community([1.0], [0.3], Float64[])
lineage = ecoDyn(lineage, config)
@time lineage = evolve(lineage, config, 1500);

p = plotEvo(lineage)
# savefig(p, "examples/plot.png")

table = historyToTable(lineage)
# CSV.write("examples/evo.csv", table)
