module EcoEvoSim

using Random
using Distributions
using LinearAlgebra
using DifferentialEquations


include("basic-types.jl")
include("basic-constructors.jl")
include("basic-selectors.jl")
include("show-methods.jl")
include("basic-utils.jl")
include("ecoevo.jl")


export PopulationSize, Phenotype, Species, Community, EvoHistory,
       speciesList, popsizes, traits, auxs, numSpecies, speciesIndices,
       randomSpecies, weightedRandomSpecies, speciesBelowThreshold, removeExtinct,
       traitSpaceDim, addSpecies, removeSpecies, changePopsizes, changeTraits,
       orderByTrait, emptyCommunity, emptyEvoHistory, IntegrationParams, EcoEvoConfig,
       generateMutant, generateMutantWeighted, ecoDyn, evolve!,
       unpackCommunity, packCommunity


end # module EcoEvoSim
