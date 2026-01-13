using EcoEvoSim
using StaticArrays


# PopulationSize for a species
popsize = PopulationSize([10.0])

trait = Phenotype([0.5])
species = Species(10.0, 0.5)

# Auxiliary variables (e.g., resources)
aux1 = PopulationSize([100.0])
aux2 = PopulationSize([200.0])

# Create Community
community = Community{Float64, 1, 1, 2}([species], [aux1, aux2])

# Print it (this will use the custom show method)
println(community)


Community([1.2, 2.6], [1.1, -1.1], [400.0])
