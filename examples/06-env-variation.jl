# Evolutionary branching driven by a fluctuating environment (storage effect).
#
# A density-dependent version of the two-season model of Miller & Klausmeier
# (2016 Theoretical Ecology). The environment is a temperature that alternates
# between two values T1 and T2, spending a fraction `phi` of each period `tau` in
# season 1 (a square wave in time). Each phenotype `z` grows best when it matches
# the current temperature: its per-resource birth rate is a Gaussian of width `w`
# centered on the current temperature. All phenotypes draw on a single shared
# resource R = Rtot - sum(n) (so growth is reduced by the total density) and share
# the same mortality `m`. This is an unstructured model with explicit time
# dependence -- note the extra `t` argument, whose presence `unstructuredModel`
# detects automatically.
#
# Coexistence and branching here arise from the storage effect, which requires
# fluctuations that are SLOW relative to the lifespan. The relevant dimensionless
# quantity is the period measured in lifespans, tau*m; Miller & Klausmeier find
# the branching region vanishes around tau*m ~ 10 and is robust by tau*m ~ 100
# (their Figs. 7-8), so we use tau*m = 200. The disruptive selection that drives
# branching is governed by the width of the tolerance curve relative to the
# season separation: in the slow-fluctuation limit the invasion-fitness curvature
# at z = 0 is m*(1 - w^2)/w^4, which is exactly zero at w = 1 (a degenerate,
# non-branching case) and grows as w shrinks. We therefore use a sharper tolerance
# curve, w = 0.5, to obtain a genuine fitness valley at z = 0 -- a branching point
# at which the lineage splits into two specialists, one adapted to each season.
#
# Numerical notes. Because the environment is periodic the dynamics settle onto a
# limit cycle rather than a fixed point, so we cannot use a steady-state solver:
# `maxTime` is instead set to many full periods, long enough for an invading
# mutant to reach its coexistence density (the periodic analogue of integrating
# to equilibrium). Selection near a branching point is second-order weak, so
# divergence is slow and many mutation events are needed. A stiff solver (Rodas5)
# with `tstops` on every season switch is used because the discontinuous square
# wave and the sharp tolerance curve otherwise destabilize the integration. The
# extinction threshold sits ABOVE the invader population size, so a mutant must
# actively grow across it to survive -- the invasion filter that keeps the
# community from accumulating a cloud of unselected near-neutral phenotypes.

using EcoEvoSim, OrdinaryDiffEq, Plots


w = 0.5 # Width of temperature tolerance curve; < 1 for a fitness valley at z = 0
tau = 200.0 # Period of the temperature fluctuations
m = 1.0 # Mortality rate, shared by all phenotypes
phi = 0.5 # Fraction of each period spent in season 1
Rtot = 10.0 # Total resource in the closed system (keeps z = 0 viable at this w and m)
T1, T2 = -1.0, 1.0 # Temperatures in season 1 and 2
g(z, T, w) = exp(-(z - T)^2 / (2 * w^2)) # Birth rate on one unit of resource
temp(t) = (t % tau) < (phi * tau) ? T1 : T2 # Temperature as a function of time

ecology = unstructuredModel() do i, n, z, aux, S, t
    R = Rtot - sum(n[j] for j in 1:S)            # Available shared resource
    n[i] * (g(sum(z[i]), temp(t), w) * R - m)
end

nMut = 1000 # Number of mutation events
maxTime = 20.3 * tau  # Non-integer multiple de-synchronizes introductions
config = EcoEvoConfig(
    ecoDyn = ecology,
    mutationGenerator = generateMutant(
        invaderPopsize = 0.001,
        variance = 0.006^2
    ),
    integrationParams = IntegrationParams(
        maxTime = maxTime,
        algorithm = Rodas5(),
        abstol = 1e-8,
        reltol = 1e-8,
        tstops = 0.0:(tau/2):((nMut + 2) * maxTime) # To help the integrator
    ),
    extThreshold = 0.003
)

initCommunity = Community([1.0], [0.0])
initCommunity = ecoDyn(initCommunity, config)
lineage = evolve(initCommunity, config, nMut);
p = plotEvo(lineage)
