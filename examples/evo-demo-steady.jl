using EcoEvoSim
using Plots
using DifferentialEquations
using DataFrames
using CSV
using Random


growthFn = (z) -> (tanh(sum(z .- 0.5) / 0.2) + 1) / 2 - 0.006692851
kernelFn = (zi, zj) -> -(tanh(sum(zi .- zj) / 0.15) + 1) / 2

config = EcoEvoConfig(
    ecoDyn = lotkaVolterra(growthFn, kernelFn),
    mutationGenerator = (c, cfg) -> generateMutant(c, cfg, 0.002^2),
    integrationParams = IntegrationParams(
        maxTime = Inf,
        algorithm = DynamicSS(),
        abstol = 1e-14,
        reltol = 1e-8
    ),
    invaderPopsize = 0.001,
    extThreshold = 0.003
)

Random.seed!(54321)

lineage = Community([1.0], [0.3], Float64[])

@time lineage = evolve!(lineage, config, 8000);

p = plotEvo(lineage)
# savefig(p, "examples/plot-steady.png")

table = historyToTable(lineage)
CSV.write("examples/evo.csv", DataFrame(table))
