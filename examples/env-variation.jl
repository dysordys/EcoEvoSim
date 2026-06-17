using EcoEvoSim, OrdinaryDiffEq, Plots


ecology = unstructuredModel() do i, n, z, aux, S, t
    R = 1 - sum(n[j] for j in 1:S)
    m = [1 + cos(sum(z[i])) / 2 for i in 1:S]
    g = [2 + cos(2*pi*t / 10 + sum(z[i])) for i in 1:S]
    n[i] * (g[i] * R - m[i])
end

config = EcoEvoConfig(
    ecoDyn = ecology,
    mutationGenerator = generateMutant(
        invaderPopsize = 0.001,
        variance = 0.002^2
    ),
    integrationParams = IntegrationParams(
        maxTime = 1e4,
        algorithm = Vern7(),
        abstol = 1e-8,
        reltol = 1e-8
    ),
    extThreshold = 0.003
)

initCommunity = Community([1.0], [0.0])
initCommunity = ecoDyn(initCommunity, config)
lineage = evolve(initCommunity, config, 1500);
p = plotEvo(lineage)
