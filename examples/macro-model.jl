using EcoEvoSim
using Plots
using Random


growthFn(z) = (tanh(sum(z .- 0.5) / 0.2) + 1) / 2 - 0.006692851
kernelFn(zi, zj) = -(tanh(sum(zi .- zj) / 0.15) + 1) / 2

# Create the model using the new @unstructuredModel macro
# Users simply specify what dn[i]/dt is, and the macro handles the rest
ecology = @unstructuredModel begin
    dn[i] = n[i] * (growthFn(z[i]) + sum(kernelFn(z[i], z[j]) * n[j] for j in 1:nSpecies))
end

config = EcoEvoConfig(
    ecoDyn = ecology,
    mutationGenerator = c -> generateMutant(c; invaderPopsize=0.001, variance=0.002^2),
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
@time lineage = evolve!(lineage, config, 1500);

p = plotEvo(lineage)
# savefig(p, "examples/plot.png")
