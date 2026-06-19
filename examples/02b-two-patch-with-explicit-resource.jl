# Two-patch model with explicit resource dynamics.
#
# This is the same model as 02a-two-patch.jl, but with the resource that limits
# growth represented explicitly (one resource variable per patch, tracked
# via `auxDynamics`) instead of being folded into a logistic growth term. With
# the parameters below the resources equilibrate to make this model identical to
# the implicit one: at the quasi-steady state of the resource,
# `beta*aux[k] - m0` reduces to `-sum(n[i,k])`, recovering the `localGrowth`
# expression of 02a-two-patch.jl.

using EcoEvoSim, OrdinaryDiffEq, Plots

mu = 0.1 # Example value for migration rate
y = [-1/2, 1/2] # Example patch optima
eta = 1.0; chi = 1.0; gamma = 1.0 # Resource supply, self-limitation, consumption
beta = 1.0; theta = 1.0 # Resource conversion rate and patch-matching strength
m0 = beta * eta / chi # Baseline mortality (equals beta * resource carrying capacity)

ecology = structuredModel(
    auxDynamics = (aux, n, z, S, K) ->
        [aux[k] * (eta - chi * aux[k]) - gamma * aux[k] * sum(n[i,k] for i in 1:S)
         for k in 1:K]
) do i, k, n, z, aux, S, K
    localGrowth = beta * aux[k] + theta * (1 - (sum(z[i]) - y[k])^2) - m0
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


initCommunity = Community([1.0 0], [0.15], [1.0, 1.0])
initCommunity = ecoDyn(initCommunity, config)
lineage = evolve(initCommunity, config, 2500);
p = plotEvo(lineage)
