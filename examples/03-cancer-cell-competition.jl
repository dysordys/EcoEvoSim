# Cancer cell competition: niche differentiation from similarity-based competition.
#
# A generalized Lotka-Volterra system, built with the `lotkaVolterra` helper.
# `growthFn` is the quadratic intrinsic growth rate b(z) = 1 - z^2/z0^2
# (peaking at z = 0, with z0 = 0.5), and `kernelFn` is a Gaussian competition
# kernel: more similar phenotypes compete more strongly, with the kernel width
# set by the denominator in the exponent. Biologically, the trait expresses cells'
# nutrient requirements, so similar clones consume similar nutrients and live at
# each other's expense. A single ancestor undergoes repeated branching into a
# stable community of coexisting clones spread across trait space.
#
# The same specification works in any number of trait dimensions; only the initial
# trait vector changes. Below it is implemented in both a 1D and 2D trait space
# (each clone then has two independent trait axes, e.g. two metabolic pathways).
# The final two lines showcase the package's interactive plotting option.

using EcoEvoSim, OrdinaryDiffEq, Plots


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

# For an interactive plot of the 2D simulation:
using GLMakie
plotEvoTwoTraitInteractive(lineage2D)
