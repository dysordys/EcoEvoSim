using EcoEvoSim
using Plots
using DataFrames
using CSV
using Random
using DifferentialEquations


growthFn(z) = (tanh(sum(z .- 0.5) / 0.2) + 1) / 2 - 0.006692851
kernelFn(zi, zj) = -(tanh(sum(zi .- zj) / 0.15) + 1) / 2

config = EcoEvoConfig(
    ecoDyn = lotkaVolterra(growthFn, kernelFn),
    mutationGenerator = generateMutant(invaderPopsize = 0.001, variance = 0.002^2),
    integrationParams = IntegrationParams(
        maxTime = Inf,
        algorithm = DynamicSS(),
        abstol = 1e-16,
        reltol = 1e-8
    ),
    extThreshold = 0.003
)

Random.seed!(54321)

lineage = Community([1.0], [0.3], Float64[])
lineage = ecoDyn(lineage, config)
@time lineage = evolve(lineage, config, 1500);

p = plotEvo(lineage)
# savefig(p, "examples/plot-steady.png")

table = historyToTable(lineage)
# CSV.write("examples/evo.csv", DataFrame(table))
