module EcoEvoSim

using StaticArrays

include("types/basic-types.jl")
include("types/basic-constructors.jl")
include("types/basic-selectors.jl")
include("types/show-methods.jl")
include("types/basic-utils.jl")

export PopulationSize, Phenotype, Species, Community, EvoHistory,
       speciesList, popsizes, traits, auxs, numSpecies, speciesIndices,
       addSpecies, emptyCommunity, emptyEvoHistory

end # module EcoEvoSim
