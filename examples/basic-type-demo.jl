using EcoEvoSim
using StaticArrays


# PopulationSize for a species
popsize = PopulationSize([10.0, 11.3])

trait = Phenotype([0.5, 1.4])
species = Species(popsize, trait)

# Auxiliary variables (e.g., resources)
aux1 = PopulationSize(100.0)
aux2 = PopulationSize(200.0)

# Create Community
community = Community([species], [aux1, aux2])

# Print it (this will use the custom show method)
println(community)


Community([1.2, 2.6], [1.1, -1.1], [400.0])
