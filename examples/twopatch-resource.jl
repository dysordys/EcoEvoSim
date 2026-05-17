using EcoEvoSim
using Plots
using Distributions
using Random


function patchMortality(community, d, theta, sigma, beta, eta, chi)
    # Compute niche-dependent mortality: m_k = m_base - niche_benefit
    # where m_base = beta*eta/chi and niche_benefit = theta * pdf(Normal(d_k, sigma), z_i)
    # Returns a matrix with species in rows and patches in columns
    spTraits = traits(community)
    nSpecies = numSpecies(community)
    y = [d/2, -d/2]  # Optimal traits for each patch
    m_base = (beta * eta) / chi  # Baseline mortality rate
    mort = zeros(nSpecies, length(y))
    for i in 1:nSpecies
        niche_benefit = theta .* pdf.(Normal.(y, sigma), spTraits[i][1])
        mort[i, :] = m_base .- niche_benefit
    end
    return mort
end


function multipatch(
        community::Community{T, AuxClasses}; d, mu, eta, chi, gamma, beta, theta, sigma
    ) where {T<:Real, AuxClasses}
    nSpecies = numSpecies(community)
    mort = patchMortality(community, d, theta, sigma, beta, eta, chi)  # species × patches
    nPatches = size(mort, 2)
    function odeFn(u::Vector{T}, p, t) where {T}
        # Reconstruct state into density matrix (species x patches)
        densityMatrix = zeros(T, nSpecies, nPatches)
        for i in 1:nSpecies
            for j in 1:nPatches
                idx = nPatches * (i - 1) + j
                densityMatrix[i, j] = u[idx]
            end
        end
        # Build dynamical matrix with current resource state
        nEqs = nSpecies * nPatches
        M = zeros(T, nEqs, nEqs)
        for i in 1:nSpecies
            for j in 1:nPatches
                idx_j = nPatches * (i - 1) + j
                # Current resource in patch j
                R_j = u[nSpecies*nPatches + j]
                # Per-capita growth rate: r_k = beta*R_k - m_k
                r_k = beta * R_j - mort[i, j]
                # Diagonal terms: within-patch growth
                M[idx_j, idx_j] = r_k - mu
                # Off-diagonal terms: migration to other patches
                for k in 1:nPatches
                    if k != j
                        idx_k = nPatches * (i - 1) + k
                        M[idx_j, idx_k] = mu
                    end
                end
            end
        end

        # Species dynamics
        du_species = M * u[1:nSpecies*nPatches]

        # Resource dynamics: dR_k/dt = R_k(eta - chi*R_k) - gamma*sum_j N_j(k)*R_k
        # Resources occupy indices nSpecies*nPatches+1 to nSpecies*nPatches+nPatches
        du_res = zeros(T, nPatches)
        for k in 1:nPatches
            R_k = u[nSpecies*nPatches + k]
            # Total consumer abundance in patch k
            total_consumers = sum(densityMatrix[i, k] for i in 1:nSpecies)
            # Resource dynamics
            du_res[k] = R_k * (eta - chi * R_k) - gamma * total_consumers * R_k
        end

        # Combine species and resource dynamics
        return [du_species; du_res]
    end
    return odeFn
end


config = EcoEvoConfig(
    ecoDyn = (comm::Community) -> multipatch(comm; d = 1.0, mu = 0.1,
                                eta = 1.0, chi = 1.0, gamma = 1.0,
                                beta = 1.0, theta = 1.0, sigma = 1.0),
    mutationGenerator =
        generateMutantSpatial(invaderPopsize = 0.001, variance = 0.003^2),
    integrationParams = IntegrationParams(
        maxTime = Inf,
        algorithm = DynamicSS(),
        abstol = 1e-10,
        reltol = 1e-10
    ),
    extThreshold = 0.003
)


Random.seed!(54321)

lineage = Community([1.0 1.0;], [-0.2], [1.0, 1.0])
lineage = ecoDyn(lineage, config)
@time lineage = evolve(lineage, config, 1000);

p = plotEvo(lineage)
