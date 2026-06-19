# Basic cancer growth, with precompute (same model as 01a-basic-cancer-growth.jl).
#
# Identical dynamics to 01a, but the intrinsic growth rate `1 - z_i^2`, which
# stays constant throughout an ecological phase, is computed once per phase via
# the `precompute` helper instead of being re-evaluated at every integration
# step. The `precompute` function must return a named tuple (hence the trailing
# comma after the closing bracket). This tuple holds the per-clone growth rates `b`
# here; the equation function then reads them as `pre.b[i]`. The result is the
# same trajectory as 01a, but cheaper to compute. Note: the gain in speed will be
# marginal, since 1 - sum(z[i]) is cheap to evaluate. But it does illustrate how
# `precompute` could be used in more serious examples as well.

using EcoEvoSim, OrdinaryDiffEq, Plots


ecology = unstructuredModel(
    precompute = (z, S) -> (b = [1 - sum(z[i])^2 for i in 1:S],)
) do i, n, z, aux, S, pre
    n[i] * (pre.b[i] - sum(n[j] for j in 1:S))
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
