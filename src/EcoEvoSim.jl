module EcoEvoSim

using Random
using Distributions
using LinearAlgebra
using OrdinaryDiffEq
using SteadyStateDiffEq
using NonlinearSolve
using OrderedCollections


include("basic-types.jl")
include("basic-constructors.jl")
include("basic-selectors.jl")
include("basic-utils.jl")
include("ecoevo.jl")
include("models.jl")
include("show-methods.jl")


export PopulationSize, Phenotype, Species, Community, EvoHistory,
       popsize, trait, commTime,
       speciesList, popsizes, traits, auxs, numSpecies, speciesIndices,
       popsizesToMatrix, randomSpecies, weightedRandomSpecies,
       speciesBelowThreshold, removeExtinct, traitSpaceDim, numStages, addSpecies,
       removeSpecies, addAux, removeAux, changePopsizes,
       changeTraits, orderByTrait, selectTraitDim,
       emptyCommunity, emptyEvoHistory, IntegrationParams, EcoEvoConfig,
       generateMutant, generateMutantWeighted, generateMutantSpatial,
       generateMutantSpatialWeighted, noMutation,
       ecoDyn, ecoDynTimeSeries, singleEvoStep,
       evolve, unpackCommunity, packCommunity, lotkaVolterra,
       plotEvo, plotEvoTwoTrait, plotEvoTwoTraitInteractive,
       niceTickInterval, historyToTable, timeSeriesToTable,
       historyList, filterHistory, lastCommunity,
       unstructuredModel, structuredModel,
       DynamicSS, SSRootfind, FunctionMap, DiscreteSS


# Stubs for the Plots extension — implemented when Plots is loaded
function plotEvo end
function plotEvoTwoTrait end

# Stub for the GLMakie extension — implemented when GLMakie is loaded
function plotEvoTwoTraitInteractive end


end # module EcoEvoSim

