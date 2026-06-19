# Basic cancer growth: directional evolution toward an environmental optimum.
#
# A single evolving trait `z` measures environmental adaptation, with `z = 0` the
# best-adapted phenotype. Every clone is limited by the same nutrient in the same
# way, so its growth is reduced in proportion to the total community density. The
# per capita growth rate is therefore a downward parabola in `z` (peaking at the
# optimum) minus a term summing all densities. After setting various constants of
# proportionality to 1, the model reads:
#
#     dn_i/dt = n_i * (1 - z_i^2 - sum_j n_j).
#
# Starting with a single clone far from the optimum (z = 0.95), the trait evolves
# monotonically toward z = 0 as fitter phenotypes successively replace less fit ones,
# with no evolutionary branching.

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
