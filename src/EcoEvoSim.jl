module EcoEvoSim

using Random
using Distributions
using LinearAlgebra
using DifferentialEquations
using Plots
using OrderedCollections


include("basic-types.jl")
include("basic-constructors.jl")
include("basic-selectors.jl")
include("show-methods.jl")
include("basic-utils.jl")
include("ecoevo.jl")
include("models.jl")
include("visualize.jl")


export PopulationSize, Phenotype, Species, Community, EvoHistory,
       speciesList, popsizes, traits, auxs, numSpecies, speciesIndices,
       randomSpecies, weightedRandomSpecies, speciesBelowThreshold, removeExtinct,
       traitSpaceDim, addSpecies, removeSpecies, changePopsizes, changeTraits,
       orderByTrait, emptyCommunity, emptyEvoHistory, IntegrationParams, EcoEvoConfig,
       generateMutant, generateMutantWeighted, ecoDyn, singleEvoStep, evolve!,
       unpackCommunity, packCommunity, lotkaVolterra, plotEvo, niceTickInterval,
       historyToTable


end # module EcoEvoSim
