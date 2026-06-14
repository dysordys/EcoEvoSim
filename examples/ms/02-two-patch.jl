using EcoEvoSim, OrdinaryDiffEq, Plots


mu = 0.1 # Example value for migration rate
y = [-1/2, 1/2] # Example patch optima

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
