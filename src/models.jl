# Pre-defined ecological models


"""
    lotkaVolterra(growthFn, kernelFn)

Factory function for creating a multispecies Lotka-Volterra competition model.

Returns a function that takes a `Community` and produces the ecological dynamics
function (`ecoDyn`) for that community.

# Arguments
- `growthFn`: Function that takes a trait vector (for one species) and returns
  a scalar intrinsic growth rate. Signature: `Vector{T} -> T`
- `kernelFn`: Function that takes two trait vectors and returns the
  interaction coefficient. Signature: `(Vector{T}, Vector{T}) -> T`

# Returns
A function with signature `Community -> Function` that takes a community and
returns the dynamics function for the ODE solver.

# Model
The Lotka-Volterra model has the form:
```
dN_i/dt = N_i * (b_i - sum_j a_ij * N_j)
```
where:
- `N_i` is the population size of species i
- `b_i` is the intrinsic growth rate (from `growthFn`)
- `a_ij` is the interaction coefficient between species i and j (from `kernelFn`)

# Example
```julia
using EcoEvoSim

# Define growth rate as a function of trait
growthFn = traits -> 1.0 - sum(traits.^2)

# Define competition based on trait distance with Gaussian kernel
kernelFn = (z_i, z_j) -> -exp(-sum((z_i .- z_j).^2) / 0.15^2)

# Create a community
community = Community([1.0, 1.0, 1.0], [-0.2, 0.0, 0.3], Float64[])

# Create configuration with desired settings
config = EcoEvoConfig(
    ecoDyn = lotkaVolterra(growthFn, kernelFn),
    mutationGenerator = (x) -> x .+ 0.01,
    integrationParams = IntegrationParams(maxTime = 50.0),
    invaderPopsize = 0.001,
    extThreshold = 0.003
)

# Run dynamics
finalCommunity = ecoDyn(community, config)
```
"""


function lotkaVolterra(growthFn, kernelFn)
    # Return a function that creates the dynamics function for a given community
    function createDynamics(community::Community{T, AuxClasses}) where {T<:Real, AuxClasses}
        # Extract traits from the community
        nSpecies = numSpecies(community)
        spTraits = [traits(community, i) for i in 1:nSpecies]

        # Compute intrinsic growth rates using the provided function
        b = [growthFn(spTraits[i]) for i in 1:nSpecies]

        # Precompute interaction coefficient matrix using the provided function
        a = zeros(T, nSpecies, nSpecies)
        for i in 1:nSpecies
            for j in 1:nSpecies
                a[i, j] = kernelFn(spTraits[i], spTraits[j])
            end
        end

        # Return the dynamics function with proper signature for ODEProblem
        # Note: We can safely close over a and b because this function is recreated
        # for each new community in ecoDyn (which calls config.ecoDyn(community))
        return (u, p, t) -> u .* (b .+ (a * u))
    end

    return createDynamics
end


# ─── Function-based model helpers ────────────────────────────────────────────

"""
    unstructuredModel(eqn_fn; precompute=nothing)

Create an unstructured ecological model from a per-species growth rate function.

The user-supplied function `eqn_fn(i, n, z, nSpecies)` should return `dn[i]/dt`
for species `i`, where:
- `i::Int`: current species index
- `n`: population size vector (`n[j]` = population of species `j`)
- `z`: vector of trait vectors (`z[j]` = trait vector of species `j`)
- `nSpecies::Int`: total number of species

# Keyword Arguments
- `precompute`: optional function `(z, nSpecies) -> pre` that is called once per
  community to precompute trait-dependent quantities (e.g., interaction matrices).
  When provided, the equation function receives an extra argument:
  `eqn_fn(i, n, z, nSpecies, pre)`.

Returns a factory function `Community -> (u, p, t) -> du` suitable for
`EcoEvoConfig.ecoDyn`.

# Examples

Basic usage (no precomputation):
```julia
r(z) = 1 - sum(z.^2)
α(zi, zj) = exp(-sum((zi .- zj).^2) / 0.04)

ecology = unstructuredModel() do i, n, z, nSpecies
    n[i] * (r(z[i]) - sum(α(z[i], z[j]) * n[j] for j in 1:nSpecies))
end
```

With precomputation for efficiency (avoids recomputing the interaction matrix
at every ODE timestep):
```julia
ecology = unstructuredModel(
    precompute = (z, nSpecies) -> (
        b = [r(z[i]) for i in 1:nSpecies],
        A = [α(z[i], z[j]) for i in 1:nSpecies, j in 1:nSpecies]
    )
) do i, n, z, nSpecies, pre
    n[i] * (pre.b[i] - sum(pre.A[i, j] * n[j] for j in 1:nSpecies))
end
```
"""
function unstructuredModel(eqn_fn; precompute=nothing)
    function factory(community::Community{T, AuxClasses}) where {T<:Real, AuxClasses}
        nSp = numSpecies(community)
        z = [traits(community, i) for i in 1:nSp]

        if precompute !== nothing
            pre = precompute(z, nSp)
            return function (u, p, t)
                n = @view u[1:nSp]
                du = similar(u)
                for i in 1:nSp
                    du[i] = eqn_fn(i, n, z, nSp, pre)
                end
                return du
            end
        else
            return function (u, p, t)
                n = @view u[1:nSp]
                du = similar(u)
                for i in 1:nSp
                    du[i] = eqn_fn(i, n, z, nSp)
                end
                return du
            end
        end
    end
    return factory
end


"""
    structuredModel(eqn_fn; auxDynamics=nothing, precompute=nothing)

Create a structured ecological model (spatial patches or stage classes) from
per-species-per-patch dynamics.

The user-supplied function `eqn_fn(i, j, N, z, R, nSpecies, nPatches)` should
return `dN[i,j]/dt` for species `i` in patch/stage `j`, where:
- `i::Int`: species index
- `j::Int`: patch/stage index
- `N`: density matrix (`N[i,j]` = density of species `i` in patch `j`)
- `z`: vector of trait vectors (`z[i]` = trait vector of species `i`)
- `R`: auxiliary variable vector (e.g., resource levels per patch)
- `nSpecies::Int`: total number of species
- `nPatches::Int`: number of patches/stages (inferred from community)

# Keyword Arguments
- `auxDynamics`: optional function `(R, N, z, nSpecies, nPatches) -> Vector`
  returning time derivatives for auxiliary variables (e.g., resource dynamics).
- `precompute`: optional function `(z, nSpecies, nPatches) -> pre` that is called
  once per community to precompute trait-dependent quantities. When provided,
  the equation function receives an extra argument:
  `eqn_fn(i, j, N, z, R, nSpecies, nPatches, pre)`.

Returns a factory function `Community -> (u, p, t) -> du` suitable for
`EcoEvoConfig.ecoDyn`.

# Examples

Basic usage:
```julia
using EcoEvoSim, Distributions

d = 1.0; mu = 0.1; alpha = 1.0
y = [d/2, -d/2]

ecology = structuredModel() do i, j, N, z, R, nSpecies, nPatches
    growth = pdf(Normal(0, 1), z[i][1] - y[j])
    dd = alpha * sum(N[k, j] for k in 1:nSpecies)
    (growth - dd - mu) * N[i, j] + mu * sum(N[i, k] for k in 1:nPatches if k != j)
end
```

With precomputation:
```julia
ecology = structuredModel(
    precompute = (z, nSpecies, nPatches) ->
        [pdf(Normal(0, 1), z[i][1] - y[j]) for i in 1:nSpecies, j in 1:nPatches]
) do i, j, N, z, R, nSpecies, nPatches, pre
    dd = alpha * sum(N[k, j] for k in 1:nSpecies)
    (pre[i, j] - dd - mu) * N[i, j] + mu * sum(N[i, k] for k in 1:nPatches if k != j)
end
```
"""
function structuredModel(eqn_fn; auxDynamics=nothing, precompute=nothing)
    function factory(community::Community{T, AuxClasses}) where {T<:Real, AuxClasses}
        nSp = numSpecies(community)
        nPatch = numStages(community)
        z = [traits(community, i) for i in 1:nSp]
        nAux = sum(length(a.popsize) for a in community.aux; init=0)

        pre = precompute !== nothing ? precompute(z, nSp, nPatch) : nothing

        function ode_fn(u, p, t)
            # Reshape species state: u[1:nSp*nPatch] → species(rows) × patches(cols)
            N = reshape(@view(u[1:nSp*nPatch]), nPatch, nSp)'
            R = @view u[nSp*nPatch+1:nSp*nPatch+nAux]

            du = similar(u)

            # Species dynamics
            for i in 1:nSp
                for j in 1:nPatch
                    du[nPatch*(i-1)+j] = pre !== nothing ?
                        eqn_fn(i, j, N, z, R, nSp, nPatch, pre) :
                        eqn_fn(i, j, N, z, R, nSp, nPatch)
                end
            end

            # Auxiliary dynamics
            if auxDynamics !== nothing
                du_aux = auxDynamics(R, N, z, nSp, nPatch)
                du[nSp*nPatch+1:end] .= du_aux
            end

            return du
        end

        return ode_fn
    end
    return factory
end
