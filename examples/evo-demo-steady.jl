using EcoEvoSim
using Plots
using DifferentialEquations

# Create a simple evolutionary simulation
growthFn = (z) -> sum(z.^2) / (3 + sum(z.^2))
kernelFn = (zi, zj) -> -(tanh(sum(zi .- zj) / 0.3) + 1) / 2

community = Community([1.0], [0.3], Float64[])

config = EcoEvoConfig(
    ecoDyn = lotkaVolterra(growthFn, kernelFn),
    mutationGenerator = (c, cfg) -> generateMutant(c, cfg, 0.002^2),
    integrationParams = IntegrationParams(
        maxTime = Inf,   # No timeout for steady-state solver
        algorithm = DynamicSS(),
        abstol = 1e-14,  # Absolute tolerance for steady state
        reltol = 1e-8    # Relative tolerance for steady state
    ),
    invaderPopsize = 1e-3,
    extThreshold = 3e-3
)

history = evolve!(community, config, 8000)

p = plotEvo(history)
# savefig(p, "test_plotEvo_steady.png")
