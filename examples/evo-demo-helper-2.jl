# Equivalent of evo-demo.jl, but using the function-based helper to define
# the model instead of writing out the full ODE function. It also demonstrates
# the use of precomputed values for efficiency. Compare to evo-demo.jl and
# evo-demo-helper-1.jl.

using EcoEvoSim
using Plots
using Random


growthFn(z) = (tanh(sum(z .- 0.5) / 0.2) + 1) / 2 - 0.006692851
kernelFn(zi, zj) = -(tanh(sum(zi .- zj) / 0.15) + 1) / 2

# Create model using unstructuredModel helper function, using precomputed
# intrinsic growth rates and interaction matrices for efficiency. The precompute
# function is called once per community, and the results are passed to the
# equation function as `pre`.
ecology = unstructuredModel(
    precompute = (z, nSpecies) -> (
        b = [growthFn(z[i]) for i in 1:nSpecies],
        A = [kernelFn(z[i], z[j]) for i in 1:nSpecies, j in 1:nSpecies]
    )
) do i, n, z, nSpecies, pre
    n[i] * (pre.b[i] + sum(pre.A[i, j] * n[j] for j in 1:nSpecies))
end

config = EcoEvoConfig(
    ecoDyn = ecology,
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
