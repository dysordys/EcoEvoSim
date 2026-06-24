# Two-patch model: evolutionary branching driven by spatial heterogeneity.
#
# This is a spatially-structured extension of the basic cancer growth model. The
# two patches differ in their optimal phenotype, so clone i's local growth in
# patch k peaks at the patch optimum y[k]. Local competition within a patch is
# density-dependent (the `-sum(n[i,k] ...)` term), and clones disperse symmetrically
# between the patches at rate `mu`. The `structuredModel` helper provides the extra
# patch index `k`, and `n[i,k]` is the density of clone i in patch k. Starting from
# a single clone, the population branches into two lineages, each adapting to one
# patch. A large enough `mu` makes a single generalist best, and will suppress
# the branching.

using EcoEvoSim, Plots


mu = 0.1 # Migration rate
y = [-1/2, 1/2] # The two patch optima

ecology = structuredModel() do i, k, n, z, aux, S, K
    localGrowth = 1 - (sum(z[i]) - y[k])^2 - sum(n[i,k] for i in 1:S)
    (localGrowth - mu) * n[i,k] + mu * sum(n[i,l] for l in 1:K if l != k)
end

config = EcoEvoConfig(
    ecoDyn = ecology,
    mutationGenerator = generateMutant(
        invaderPopsize = 0.001,
        variance = 0.002^2
    ),
    integrationParams = IntegrationParams(
        maxTime = 1e10,
        algorithm = Rodas5(),
        abstol = 1e-8,
        reltol = 1e-8
    ),
    extThreshold = 0.003
)


initCommunity = Community([1.0 0], [0.15])
initCommunity = ecoDyn(initCommunity, config)
lineage = evolve(initCommunity, config, 2500);
p = plotEvo(lineage)
