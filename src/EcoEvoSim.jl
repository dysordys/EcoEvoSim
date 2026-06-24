module EcoEvoSim

using Random
using Distributions
using LinearAlgebra
using OrdinaryDiffEq
using SteadyStateDiffEq
using NonlinearSolve
using OrderedCollections


include("types.jl")
include("constructors.jl")
include("selectors.jl")
include("utils.jl")
include("mutgen.jl")
include("evoevo.jl")
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
       # Integration algorithms: EcoEvoSim's own plus a curated set re-exported
       # from the SciML solver stack, so common models need only `using EcoEvoSim`.
       DynamicSS, SSRootfind, FunctionMap, DiscreteSS,
       Rodas5, RadauIIA5, Tsit5, Vern7, AutoVern7


# Stubs for the Plots extension — implemented when Plots is loaded
function plotEvo end
function plotEvoTwoTrait end

# Stub for the GLMakie extension — implemented when GLMakie is loaded
function plotEvoTwoTraitInteractive end


end # module EcoEvoSim

