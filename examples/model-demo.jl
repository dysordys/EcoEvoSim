# Demonstration: Lotka-Volterra competition with trait-dependent parameters
#
# This example shows how to set up and run ecological dynamics for a community
# competing under Lotka-Volterra dynamics with trait-dependent growth rates
# and competition coefficients.
#
# Model:
#   dN_i/dt = N_i * (b_i - sum_j a_ij * N_j)
#
# where:
#   - N_i: population size of species i
#   - b_i = 1 - z_i: trait-dependent intrinsic growth rate
#   - a_ij = exp(-(d_ij/w)^2): competition coefficient
#   - d_ij = |z_i - z_j|: trait distance between species i and j
#   - w: competition width parameter

using EcoEvoSim
using DifferentialEquations
using Plots


# Create initial community with 3 species at different trait values
println("Setting up initial community...")
initialPopsizes = [1.0, 1.0, 1.0]  # Equal initial densities
initialTraits = [-0.2, 0.0, 0.3]  # Three species with different traits
community = Community(initialPopsizes, initialTraits, Float64[])
println(community)


# Wrapper function that extracts traits from the community state
# The dynamics function needs signature (u, p, t) for ODEProblem
# where p will contain the community information (traits and parameters)
function makeLVDynamics(community::Community, w::Float64)
    # Extract traits from the community at setup time
    spTraits = [EcoEvoSim.traits(community, i)[1] for i in 1:numSpecies(community)]
    # Precompute competition coefficient matrix: a_ij = exp(-(d_ij/w)^2)
    # Broadcasting for outer difference: spTraits .- spTraits' creates distance matrix
    d = @. abs(spTraits - spTraits')
    A = @. exp(-(d / w)^2)
    # Precompute intrinsic growth rates: b_i = 1 - z_i^2
    b = @. 1.0 - spTraits^2
    # Return dynamics function with proper signature for ODEProblem
    return (u, p, t) -> u .* (b .- A * u)
end


# Create configuration
config = EcoEvoConfig(
    ecoDyn = makeLVDynamics(community, 0.15),
    mutationGenerator = (x) -> x .+ 0.01, # Dummy mutation generator (not used)
    integrationParams = IntegrationParams(
        maxTime = 50.0,     # Integrate for 50 time units
        abstol = 1e-8,      # Absolute tolerance
        reltol = 1e-6       # Relative tolerance
        # Uses default Rodas5() stiff solver
        # Can add any other solver options: maxiters=10000, save_everystep=false, etc.
    ),
    invaderPopsize = 0.001,  # Not used in this ecology-only example
    extThreshold = 0.003
)

# Run ecological dynamics
finalCommunity = ecoDyn(community, config)
println(finalCommunity)

# Calculate final growth rates
finalTraits = collect(Iterators.flatten(traits(finalCommunity)))
finalPops = collect(Iterators.flatten(popsizes(finalCommunity)))


# Plot trait-fitness relationship
println("\nGenerating plot of trait vs. growth rate...")
traitRange = range(-0.6, 0.6, length=201)
growthRates = 1.0 .- traitRange.^2

p = plot(traitRange, growthRates,
         label="Growth rate b(z) = 1 - z^2",
         xlabel="Trait value z",
         ylabel="Intrinsic growth rate b",
         linewidth=2,
         legend=:bottom)

# Mark initial and final species positions
scatter!(p, initialTraits, 1.0 .- initialTraits.^2,
         label="Initial species",
         markersize=8,
         color=:blue)

scatter!(p, finalTraits, 1.0 .- finalTraits.^2,
         label="Final species (size ∝ N)",
         markersize=max.(5, finalPops .* 10),  # Size proportional to population
         color=:red)

display(p)
savefig(p, "examples/lv-demo-fig.png")
println("Plot saved to examples/lv-demo-fig.png")

println("\n✓ Demonstration complete!")
