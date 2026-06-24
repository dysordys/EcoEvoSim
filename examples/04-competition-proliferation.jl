# Competition-proliferation tradeoff: hierarchical coexistence.
#
# A generalized Lotka-Volterra model encoding the competition-proliferation
# (a.k.a. competition-colonization) tradeoff model. The evolving trait `z` is
# proliferation capacity: high z means fast growth, low z means a superior
# competitor. The sigmoid `Q` rises smoothly from 0 to 1 around its inflection
# point. `growthFn` makes the intrinsic growth rate increase monotonically with
# z while pinning b(0) = 0 (so clones with z < 0 cannot grow even without
# competition, making z = 0 the most competitive phenotype that can arise).
# `kernelFn` encodes a competitive hierarchy: a slower-growing (lower z) clone
# outcompetes a faster-growing one.
#
# In the configuration block, `generateMutantWeighted` picks the parent clone with
# probability proportional to its density, mimicking the per capita nature of
# mutation. `DynamicSS(RadauIIA5())` integrates each ecological phase to steady
# state rather than for a set number of time units. The `abstol` and `reltol`
# parameters then control how tightly the steady state must be approached. In the
# simulation, a single ancestor diversifies into a stable, hierarchically structured
# set of clones.

using EcoEvoSim, Plots


Q(x) = (tanh(x) + 1) / 2
growthFn(z) = Q(sum(z .- 0.5) / 0.2) - Q(-0.5 / 0.2)
kernelFn(zi, zj) = -Q(sum(zi .- zj) / 0.15)

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

initCommunity = Community([1.0], [0.2])
initCommunity = ecoDyn(initCommunity, config)
lineage = evolve(initCommunity, config, 2500);
p = plotEvo(lineage)
