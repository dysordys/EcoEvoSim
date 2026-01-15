module EcoEvoSim

using Random

include("types/basic-types.jl")
include("types/basic-constructors.jl")
include("types/basic-selectors.jl")
include("types/show-methods.jl")
include("types/basic-utils.jl")

export PopulationSize, Phenotype, Species, Community, EvoHistory,
       speciesList, popsizes, traits, auxs, numSpecies, speciesIndices,
       randomSpecies, weightedRandomSpecies, speciesBelowThreshold,
       traitSpaceDim, addSpecies, removeSpecies, changePopsizes, changeTraits,
       emptyCommunity, emptyEvoHistory

end # module EcoEvoSim
