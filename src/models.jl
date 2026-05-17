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

        # Analytical Jacobian (out-of-place): J[i,j] = δ(i,j)*(b[i] + (Au)[i]) + u[i]*a[i,j]
        # i.e. J = Diagonal(b .+ a*u) + Diagonal(u) * a
        function jacFn(u, p, t)
            au = a * u
            J = Diagonal(u) * a  # J[i,j] = u[i] * a[i,j]
            @inbounds for i in 1:nSpecies
                J[i, i] += b[i] + au[i]
            end
            return J
        end

        # Return an ODEFunction with the analytical Jacobian attached.
        # Both f and jacFn are out-of-place, satisfying SciMLBase's convention.
        # Note: We can safely close over a and b because this function is
        # recreated for each new community in ecoDyn (which calls
        # config.ecoDyn(community))
        return ODEFunction((u, p, t) -> u .* (b .+ (a * u)); jac = jacFn)
    end

    return createDynamics
end


# ─── Function-based model helpers ────────────────────────────────────────────

# Return the number of positional arguments the first method of fn expects.
# Used to detect whether the user supplied a time-dependent function by checking
# whether it accepts one extra argument beyond the standard signature.
_fn_nargs(fn) = first(methods(fn)).nargs - 1

"""
    unstructuredModel(eqnFn; auxDynamics=nothing, precompute=nothing)

Create an unstructured ecological model from a per-species growth rate function.

The user-supplied function `eqnFn(i, n, z, nSpecies)` should return `dn[i]/dt`
for species `i`, where:
- `i::Int`: current species index
- `n`: population size vector (`n[j]` = population of species `j`)
- `z`: vector of trait vectors (`z[j]` = trait vector of species `j`)
- `nSpecies::Int`: total number of species

**Time dependence:** if `eqnFn` accepts one extra argument beyond the standard
signature, the current integration time `t` is passed as the final argument.
The same applies to `auxDynamics`. This allows explicit time dependence (e.g.
seasonal forcing), without any additional keyword:
- without `precompute`: `eqnFn(i, n, z, nSpecies, t)`
- with `precompute`: `eqnFn(i, n, z, nSpecies, pre, t)`

# Keyword Arguments
- `auxDynamics`: optional function `(R, n, z, nSpecies) -> Vector`
  returning time derivatives for auxiliary variables (e.g., resource dynamics).
  When `precompute` is also provided, `auxDynamics` receives an extra argument:
  `auxDynamics(R, n, z, nSpecies, pre)`.
  Time dependence is detected by the same arity rule described above.
- `precompute`: optional function `(z, nSpecies) -> pre` that is called once per
  community to precompute trait-dependent quantities (e.g., interaction matrices).
  When provided, both the equation function and `auxDynamics` receive an extra
  argument: `eqnFn(i, n, z, nSpecies, pre)`.

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

With explicit time dependence (seasonal forcing):
```julia
ecology = unstructuredModel() do i, n, z, nSpecies, t
    r_t = 1.0 + 0.5 * sin(2π * t)   # time-varying growth rate
    n[i] * (r_t - sum(n[j] for j in 1:nSpecies))
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

With precomputation and time dependence:
```julia
ecology = unstructuredModel(
    precompute = (z, nSpecies) -> (
        A = [α(z[i], z[j]) for i in 1:nSpecies, j in 1:nSpecies],
    )
) do i, n, z, nSpecies, pre, t
    r_t = 1.0 + 0.5 * sin(2π * t)
    n[i] * (r_t - sum(pre.A[i, j] * n[j] for j in 1:nSpecies))
end
```
"""
function unstructuredModel(eqnFn; auxDynamics=nothing, precompute=nothing)
    function factory(community::Community{T, AuxClasses}) where {T<:Real, AuxClasses}
        nSp = numSpecies(community)
        z = [traits(community, i) for i in 1:nSp]
        nAux = sum(length(a.popsize) for a in community.aux; init=0)

        pre = precompute !== nothing ? precompute(z, nSp) : nothing

        # Detect time dependence: one extra arg beyond the standard signature
        # means the user wants t passed as the final argument.
        base_eqn_nargs = pre !== nothing ? 5 : 4
        eqn_uses_time  = _fn_nargs(eqnFn) > base_eqn_nargs
        aux_uses_time  = auxDynamics !== nothing && _fn_nargs(auxDynamics) > base_eqn_nargs

        function ode_fn(u, p, t)
            n = @view u[1:nSp]
            R = @view u[nSp+1:nSp+nAux]

            du = similar(u)

            # Species dynamics
            for i in 1:nSp
                du[i] = if pre !== nothing
                    eqn_uses_time ? eqnFn(i, n, z, nSp, pre, t) : eqnFn(i, n, z, nSp, pre)
                else
                    eqn_uses_time ? eqnFn(i, n, z, nSp, t) : eqnFn(i, n, z, nSp)
                end
            end

            # Auxiliary dynamics
            if auxDynamics !== nothing
                du_aux = if pre !== nothing
                    aux_uses_time ? auxDynamics(R, n, z, nSp, pre, t) : auxDynamics(R, n, z, nSp, pre)
                else
                    aux_uses_time ? auxDynamics(R, n, z, nSp, t) : auxDynamics(R, n, z, nSp)
                end
                du[nSp+1:end] .= du_aux
            end

            return du
        end

        return ode_fn
    end
    return factory
end


"""
    structuredModel(eqnFn; auxDynamics=nothing, precompute=nothing)

Create a structured ecological model (spatial patches or stage classes) from
per-species-per-patch dynamics.

The user-supplied function `eqnFn(i, j, N, z, R, nSpecies, nPatches)` should
return `dN[i,j]/dt` for species `i` in patch/stage `j`, where:
- `i::Int`: species index
- `j::Int`: patch/stage index
- `N`: density matrix (`N[i,j]` = density of species `i` in patch `j`)
- `z`: vector of trait vectors (`z[i]` = trait vector of species `i`)
- `R`: auxiliary variable vector (e.g., resource levels per patch)
- `nSpecies::Int`: total number of species
- `nPatches::Int`: number of patches/stages (inferred from community)

**Time dependence:** if `eqnFn` accepts one extra argument beyond the standard
signature, the current integration time `t` is passed as the final argument.
The same applies to `auxDynamics`. This allows explicit time dependence (e.g.
seasonal forcing), without any additional keyword:
- without `precompute`: `eqnFn(i, j, N, z, R, nSpecies, nPatches, t)`
- with `precompute`: `eqnFn(i, j, N, z, R, nSpecies, nPatches, pre, t)`

# Keyword Arguments
- `auxDynamics`: optional function `(R, N, z, nSpecies, nPatches) -> Vector`
  returning time derivatives for auxiliary variables (e.g., resource dynamics).
  When `precompute` is also provided, `auxDynamics` receives an extra argument:
  `auxDynamics(R, N, z, nSpecies, nPatches, pre)`.
  Time dependence is detected by the same arity rule described above.
- `precompute`: optional function `(z, nSpecies, nPatches) -> pre` that is called
  once per community to precompute trait-dependent quantities. When provided,
  both the equation function and `auxDynamics` receive an extra argument:
  `eqnFn(i, j, N, z, R, nSpecies, nPatches, pre)`.

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

With explicit time dependence (seasonally varying carrying capacity):
```julia
ecology = structuredModel() do i, j, N, z, R, nSpecies, nPatches, t
    K_t = 1.0 + 0.3 * sin(2π * t)   # seasonal forcing
    (K_t - sum(N[k, j] for k in 1:nSpecies) - 0.1) * N[i, j]
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

With precomputation and time dependence:
```julia
ecology = structuredModel(
    precompute = (z, nSpecies, nPatches) ->
        [pdf(Normal(0, 1), z[i][1] - y[j]) for i in 1:nSpecies, j in 1:nPatches]
) do i, j, N, z, R, nSpecies, nPatches, pre, t
    K_t = 1.0 + 0.3 * sin(2π * t)
    dd = alpha * sum(N[k, j] for k in 1:nSpecies)
    (pre[i, j] * K_t - dd - mu) * N[i, j] + mu * sum(N[i, k] for k in 1:nPatches if k != j)
end
```
"""
function structuredModel(eqnFn; auxDynamics=nothing, precompute=nothing)
    function factory(community::Community{T, AuxClasses}) where {T<:Real, AuxClasses}
        nSp = numSpecies(community)
        nPatch = numStages(community)
        z = [traits(community, i) for i in 1:nSp]
        nAux = sum(length(a.popsize) for a in community.aux; init=0)

        pre = precompute !== nothing ? precompute(z, nSp, nPatch) : nothing

        # Detect time dependence: one extra arg beyond the standard signature
        # means the user wants t passed as the final argument.
        base_eqn_nargs = pre !== nothing ? 8 : 7
        eqn_uses_time  = _fn_nargs(eqnFn) > base_eqn_nargs
        base_aux_nargs = pre !== nothing ? 6 : 5
        aux_uses_time  = auxDynamics !== nothing && _fn_nargs(auxDynamics) > base_aux_nargs

        function ode_fn(u, p, t)
            # Reshape species state: u[1:nSp*nPatch] → species(rows) × patches(cols)
            N = reshape(@view(u[1:nSp*nPatch]), nPatch, nSp)'
            R = @view u[nSp*nPatch+1:nSp*nPatch+nAux]

            du = similar(u)

            # Species dynamics
            for i in 1:nSp
                for j in 1:nPatch
                    du[nPatch*(i-1)+j] = if pre !== nothing
                        eqn_uses_time ? eqnFn(i, j, N, z, R, nSp, nPatch, pre, t) :
                                        eqnFn(i, j, N, z, R, nSp, nPatch, pre)
                    else
                        eqn_uses_time ? eqnFn(i, j, N, z, R, nSp, nPatch, t) :
                                        eqnFn(i, j, N, z, R, nSp, nPatch)
                    end
                end
            end

            # Auxiliary dynamics
            if auxDynamics !== nothing
                du_aux = if pre !== nothing
                    aux_uses_time ? auxDynamics(R, N, z, nSp, nPatch, pre, t) :
                                    auxDynamics(R, N, z, nSp, nPatch, pre)
                else
                    aux_uses_time ? auxDynamics(R, N, z, nSp, nPatch, t) :
                                    auxDynamics(R, N, z, nSp, nPatch)
                end
                du[nSp*nPatch+1:end] .= du_aux
            end

            return du
        end

        return ode_fn
    end
    return factory
end
