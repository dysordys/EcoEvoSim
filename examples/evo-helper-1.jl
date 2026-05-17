# Equivalent of evo-demo.jl, but using the function-based helper to define
# the model instead of writing out the full ODE function. Compare to evo-demo.jl.

using EcoEvoSim
using Plots
using Random


growthFn(z) = 1 - sum(z.^2) / 0.5^2
kernelFn(zi, zj) = -exp(-sum((zi .- zj).^2) / 0.15^2)


# Create model using unstructuredModel helper function:
ecology = unstructuredModel() do i, n, z, aux, nSpecies
    n[i] * (growthFn(z[i]) + sum(kernelFn(z[i], z[j]) * n[j] for j in 1:nSpecies))
end

config = EcoEvoConfig(
    ecoDyn = ecology,
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
