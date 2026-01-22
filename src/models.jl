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
