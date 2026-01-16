# Demonstration: Lotka-Volterra competition with eco-evolutionary dynamics
#
# This example shows how to use the lotkaVolterra factory function to simulate
# ecological dynamics with trait-dependent growth rates and competition.

using EcoEvoSim

# Define growth rate as a function of trait
growthFn = traits -> 1.0 .- traits.^2

# Define competition based on trait distance with Gaussian kernel
interactionFn = (z_i, z_j) -> exp(-((z_i - z_j) / 0.15)^2)

# Create a community
comm = Community([1.0, 1.0, 1.0], [-0.2, 0.0, 0.3], Float64[])

# Create configuration with desired settings
config = EcoEvoConfig(
    ecoDyn = lotkaVolterra(growthFn, interactionFn),
    mutationGenerator = (x) -> x .+ 0.01,
    integrationParams = IntegrationParams(maxTime = 50.0),
    invaderPopsize = 0.001,
    extThreshold = 0.003
)

# Run dynamics
finalCommunity = ecoDyn(comm, config)
println("\nFinal community:")
println(finalCommunity)

println("\n✓ Example complete!")
