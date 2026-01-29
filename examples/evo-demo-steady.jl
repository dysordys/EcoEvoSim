using EcoEvoSim
using Plots
using DifferentialEquations



# Create a simple evolutionary simulation
growthFn = (z) -> (tanh(sum(z .- 0.5) / 0.2) + 1) / 2 - 0.006692851
kernelFn = (zi, zj) -> -(tanh(sum(zi .- zj) / 0.15) + 1) / 2

config = EcoEvoConfig(
    ecoDyn = lotkaVolterra(growthFn, kernelFn),
    mutationGenerator = (c, cfg) -> generateMutant(c, cfg, 0.002^2),
    integrationParams = IntegrationParams(
        maxTime = Inf,   # No timeout for steady-state solver
        algorithm = DynamicSS(),
        abstol = 1e-14,  # Absolute tolerance for steady state
        reltol = 1e-8    # Relative tolerance for steady state
    ),
    invaderPopsize = 0.001,
    extThreshold = 0.003
)

lineage = Community([1.0], [0.3], Float64[])

lineage = evolve!(lineage, config, 1500)

p = plotEvo(lineage)
# savefig(p, "test_plotEvo_steady.png")
