module EcoEvoSim

using StaticArrays

include("types/basic-types.jl")
include("types/basic-constructors.jl")
include("types/basic-selectors.jl")
include("types/show-methods.jl")

export PopulationSize, Phenotype, Species, Community, EvoHistory, species, allSpecies, allAux


end # module EcoEvoSim
