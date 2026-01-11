# EcoEvoSim/test/runtests.jl
using Test
using EcoEvoSim                     # loads the top‑level package
using EcoEvoSim.EcoEvoCore         # brings the sub‑module into scope

# Load helper utilities (PropCheck generators, etc.)
include("utils.jl")
using .TestUtils: generateCommunity, communityGen   # relative import

# Include the actual test files
#include("addSpecies.jl")
include("pack_unpack_prop.jl")
#include("pack_unpack.jl")
