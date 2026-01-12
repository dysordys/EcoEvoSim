# Sample script to demonstrate the custom Base.show for Community

using EcoEvoSim
using StaticArrays

# Create sample data
T = Float64
stage_classes = 1
trait_dim = 1
aux_classes = 2

# PopulationSize for a species
popsize = PopulationSize{T, stage_classes}(SVector{stage_classes, T}(10.0))
trait = Phenotype{T, trait_dim}(SVector{trait_dim, T}(0.5))
species = Species{T, stage_classes, trait_dim}([popsize], [trait])

# Auxiliary variables (e.g., resources)
aux1 = PopulationSize{T, 1}(SVector{1, T}(100.0))
aux2 = PopulationSize{T, 1}(SVector{1, T}(200.0))
aux = [aux1, aux2]

# Create Community
community = Community{T, stage_classes, trait_dim, aux_classes}([species], aux)

# Print it (this will use the custom show method)
println(community)
