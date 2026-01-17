using EcoEvoSim
using Plots
using DifferentialEquations

# Create a simple evolutionary simulation
growthFn = (z) -> @. (tanh((z - 0.5) / 0.2) + 1) / 2 - 0.006692851
interactionFn = (z_i, z_j) -> (tanh((z_i - z_j) / 0.1) + 1) / 2

community = Community([1.0], [0.3], Float64[])

config = EcoEvoConfig(
    ecoDyn = lotkaVolterra(growthFn, interactionFn),
    mutationGenerator = (c, cfg) -> generateMutant(c, cfg, 4.0e-6),
    integrationParams = IntegrationParams(
        maxTime = Inf,   # No timeout for steady-state solver
        algorithm = DynamicSS(),
        abstol = 1e-14,  # Absolute tolerance for steady state
        reltol = 1e-8    # Relative tolerance for steady state
    ),
    invaderPopsize = 1e-3,
    extThreshold = 1e-3
)

history = evolve!(community, config, 15000)

p = plotEvo(history)
savefig(p, "test_plotEvo_steady.png")
