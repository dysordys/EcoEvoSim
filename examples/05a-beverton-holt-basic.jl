# Multispecies Beverton-Holt model: cancer cell competition in discrete time.
#
# A discrete-time analog of the cancer cell competition model
# (03-cancer-cell-competition.jl). Instead of a differential equation, the
# population is advanced by a difference equation that gives n_i(t+1) directly,
# the multispecies Beverton-Holt map:
#
#      n_i(t + 1) = n_i(t) * (b_i + 1) / (1 + sum_j A_ij n_j(t))
#
# As before, `growthFn` is the quadratic intrinsic growth rate (peaking at z = 0)
# and `kernelFn` is a Gaussian competition kernel in which similar phenotypes
# compete more strongly. Both are evaluated once per ecological phase via
# `precompute`: the per-clone growth rates `b` and the full interaction matrix `A`.
# Here this makes a substantial reduction in the execution time. A discrete-time
# solver is selected by the `FunctionMap()` algorithm; `maxTime = 20000` is then
# a number of iteration steps rather than a continuous duration.

using EcoEvoSim, Plots


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
    mutationGenerator = generateMutant(
        invaderPopsize = 0.001,
        variance = 0.002^2
    ),
    integrationParams = IntegrationParams(
        maxTime = 20000,
        algorithm = FunctionMap()
    ),
    extThreshold = 0.003
)

initCommunity = Community([1.0], [0.3])
initCommunity = ecoDyn(initCommunity, config)
lineage = evolve(initCommunity, config, 1500);
p = plotEvo(lineage)
