module EcoEvoSim

using StaticArrays

include("EcoEvoTypes.jl")
include("EcoEvoMethods.jl")
include("show.jl")

export PopulationSize, Phenotype, Species, Community


end # module EcoEvoSim
