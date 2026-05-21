using EcoEvoSim
using Plots
using Random


growthFn(z) = 1 - sum(z.^2) / 0.5^2
kernelFn(zi, zj) = exp(-sum((zi .- zj).^2) / 0.15^2)

multispeciesBevertonHolt = unstructuredModel(
    precompute = (z, nSpecies) -> (
        b = [growthFn(z[i]) for i in 1:nSpecies],
        A = [kernelFn(z[i], z[j]) for i in 1:nSpecies, j in 1:nSpecies]
    )
) do i, n, z, aux, nSpecies, pre
    n[i] * ((pre.b[i] + 1) / (1 + sum(pre.A[i, j] * n[j] for j in 1:nSpecies)))
end


config = EcoEvoConfig(
    ecoDyn = multispeciesBevertonHolt,
    mutationGenerator = generateMutant(invaderPopsize = 0.001, variance = 0.002^2),
    integrationParams = IntegrationParams(
        maxTime = Inf,
        algorithm = DiscreteSS(),
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
