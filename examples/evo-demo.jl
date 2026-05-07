using EcoEvoSim
using Plots
using Random
using DataFrames
using CSV


Q(z) = (tanh(z) + 1) / 2
growthFn(z) = Q((z[1] - 0.5) / 0.2) - Q(-0.5 / 0.2)
kernelFn(zi, zj) = -Q((zi[1] - zj[1]) / 0.15)


config = EcoEvoConfig(
    ecoDyn = lotkaVolterra(growthFn, kernelFn),
    mutationGenerator = generateMutant(invaderPopsize = 0.001, variance = 0.002^2),
    integrationParams = IntegrationParams(
        maxTime = 1e12,
        abstol = 1e-14,
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
# CSV.write("examples/evo.csv", DataFrame(table))
