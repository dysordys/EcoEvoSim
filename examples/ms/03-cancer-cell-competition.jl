using EcoEvoSim, OrdinaryDiffEq, Plots, GLMakie


growthFn(z) = 1 - sum(z.^2) / 0.5^2
kernelFn(zi, zj) = -exp(-sum((zi .- zj).^2) / 0.2^2)


config = EcoEvoConfig(
    ecoDyn = lotkaVolterra(growthFn, kernelFn),
    mutationGenerator = generateMutant(
        invaderPopsize = 0.001,
        variance = 0.002^2
    ),
    integrationParams = IntegrationParams(
        maxTime = 1e10,
        algorithm = RadauIIA5(),
        abstol = 1e-8,
        reltol = 1e-8
    ),
    extThreshold = 0.003
)


# 1D trait space:
initCommunity1D = Community([1.0], [0.3])
initCommunity1D = ecoDyn(initCommunity1D, config)
lineage1D = evolve(initCommunity1D, config, 1800);
p1D = plotEvo(lineage1D)


# 2D trait space:
initCommunity2D = Community([1.0], [0.3 -0.3])
initCommunity2D = ecoDyn(initCommunity2D, config)
lineage2D = evolve(initCommunity2D, config, 10000);
p2D = plotEvoTwoTrait(lineage2D, camera = (70, 20))
plotEvoTwoTraitInteractive(lineage2D)
