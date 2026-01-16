# Pre-defined ecological models


"""
    lotkaVolterra(growthRateFn, interactionFn)

Factory function for creating a multispecies Lotka-Volterra competition model.

Returns a function that takes a `Community` and produces the ecological dynamics
function (`ecoDyn`) for that community.

# Arguments
- `growthRateFn`: Function that takes a vector of trait values and returns
  a vector of intrinsic growth rates. Signature: `Vector{T} -> Vector{T}`
- `interactionFn`: Function that takes two trait values and returns the
  competition coefficient. Signature: `(T, T) -> T`

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
- `b_i` is the intrinsic growth rate (from `growthRateFn`)
- `a_ij` is the competition coefficient between species i and j (from `interactionFn`)

# Example
```julia
using EcoEvoSim

# Define growth rate as a function of trait
growthFn = traits -> 1.0 .- traits.^2

# Define competition based on trait distance with Gaussian kernel
interactionFn = (z_i, z_j) -> exp(-((z_i - z_j) / 0.15)^2)

# Create a community
community = Community([1.0, 1.0, 1.0], [-0.2, 0.0, 0.3], Float64[])

# Create configuration with desired settings
config = EcoEvoConfig(
    ecoDyn = lotkaVolterra(growthFn, interactionFn),
    mutationGenerator = (x) -> x .+ 0.01,
    integrationParams = IntegrationParams(maxTime = 50.0),
    invaderPopsize = 0.001,
    extThreshold = 0.003
)

# Run dynamics
finalCommunity = ecoDyn(community, config)
```
"""


function lotkaVolterra(growthRateFn, interactionFn)
    # Return a function that creates the dynamics function for a given community
    function createDynamics(community::Community{T, AuxClasses}) where {T<:Real, AuxClasses}
        # Extract traits from the community
        nSpecies = numSpecies(community)
        spTraits = [traits(community, i)[1] for i in 1:nSpecies]

        # Compute intrinsic growth rates using the provided function
        b = growthRateFn(spTraits)

        # Validate that growth rates match species count
        length(b) == nSpecies || throw(ArgumentError(
            "growthRateFn must return a vector of length $(nSpecies), got $(length(b))"
        ))

        # Precompute interaction coefficient matrix using the provided function
        A = zeros(T, nSpecies, nSpecies)
        for i in 1:nSpecies
            for j in 1:nSpecies
                A[i, j] = interactionFn(spTraits[i], spTraits[j])
            end
        end

        # Return the dynamics function with proper signature for ODEProblem
        # Note: We can safely close over A and b because this function is recreated
        # for each new community in ecoDyn (which calls config.ecoDyn(community))
        return (u, p, t) -> u .* (b .- A * u)
    end

    return createDynamics
end
