using EcoEvoSim, OrdinaryDiffEq, Plots


Q(x) = (tanh(x) + 1) / 2
growthFn(z) = Q(sum(z .- 0.5) / 0.2) - Q(-0.5 / 0.2)
kernelFn(zi, zj) = -Q(sum(zi .- zj) / 0.1)


config = EcoEvoConfig(
    ecoDyn = lotkaVolterra(growthFn, kernelFn),
    mutationGenerator = generateMutantWeighted(
        invaderPopsize = 0.001,
        variance = 0.002^2
    ),
    integrationParams = IntegrationParams(
        maxTime = 1e10,
        algorithm = DynamicSS(RadauIIA5()),
        abstol = 1e-8,
        reltol = 1e-8
    ),
    extThreshold = 0.003
)


lineage = Community([1.0], [0.2])
lineage = ecoDyn(lineage, config)
lineage = evolve(lineage, config, 2500);
p = plotEvo(lineage)
