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
# dependence (note the extra `t` argument, whose presence `unstructuredModel`
# detects automatically). A single ancestor at z = 0 branches into two specialists,
# each adapted to one season.
#
# Numerical note 1: The environment is periodic, so the ecological dynamics settle
# onto a limit cycle rather than a fixed point; we cannot use a steady-state solver.
# Instead, we integrate each phase over many periods (`maxTime`), long enough for
# an invading mutant to reach its stable trajectory. `tstops` on every season switch
# handles the discontinuous square wave. Finally, the extinction threshold is set
# BELOW the invader population size: near a branching point a fresh mutant is
# essentially neutral (second-order invasion fitness) and could never grow across a
# higher threshold within a finite phase, so it must be allowed to persist.
#
# Numerical note 2: Densities cannot go negative in the exact dynamics, since dn_i/dt
# is proportional to n_i. While the troughs of the oscillations stay well above the
# solver tolerance, the integration is stable. But if one pushes into a deep-trough
# regime (much larger tau*m; troughs many orders of magnitude below `abstol`), the
# solver can numerically overshoot and produce n_i < 0. If this is a problem, pass
# `isoutofdomain = (u, p, t) -> any(x -> x < 0, u)` to `IntegrationParams`.

using EcoEvoSim, OrdinaryDiffEq, Plots


w = 0.7 # Width of temperature tolerance curve; branching requires w < 1 in this model
tau = 10.0 # Period of temperature fluctuations
m = 1.0 # Mortality rate, shared by all phenotypes
phi = 0.5 # Fraction of each period spent in season 1
Rtot = 10.0 # Total resource in the closed system
T1, T2 = -1.0, 1.0 # Temperatures in season 1 and 2
g(z, T, w) = exp(-(z - T)^2 / (2 * w^2)) # Birth rate on one unit of resource
temp(t) = (t % tau) < (phi * tau) ? T1 : T2 # Temperature as a function of time

ecology = unstructuredModel() do i, n, z, aux, S, t
    R = Rtot - sum(n[j] for j in 1:S) # Available shared resource
    n[i] * (g(sum(z[i]), temp(t), w) * R - m)
end

nMut = 500 # Number of mutation events
maxTime = 10 * tau # Many periods
config = EcoEvoConfig(
    ecoDyn = ecology,
    mutationGenerator = generateMutant(
        invaderPopsize = 0.001,
        variance = 0.015^2
    ),
    integrationParams = IntegrationParams(
        maxTime = maxTime,
        algorithm = AutoVern7(Rodas5()), # Fast, but falls back to stiff Rodas5 if needed
        abstol = 1e-8,
        reltol = 1e-8,
        tstops = 0.0:(tau/2):((nMut + 2) * maxTime) # Help solver at every season switch
    ),
    extThreshold = 0.0005 # Below invaderPopsize
)

initCommunity = Community([1.0], [-0.3])
initCommunity = ecoDyn(initCommunity, config)
lineage = evolve(initCommunity, config, nMut);
p = plotEvo(lineage)
