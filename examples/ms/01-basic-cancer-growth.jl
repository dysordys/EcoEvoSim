using EcoEvoSim, OrdinaryDiffEq, Plots


ecology = unstructuredModel() do i, n, z, aux, S
    n[i] * (1 - sum(z[i])^2 - sum(n[j] for j in 1:S))
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

initCommunity = Community([1.0], [0.95])
initCommunity = ecoDyn(initCommunity, config)
lineage = evolve(initCommunity, config, 1500);
p = plotEvo(lineage)
