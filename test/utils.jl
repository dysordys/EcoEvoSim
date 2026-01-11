module TestUtils

using EcoEvoSim.EcoEvoCore
using StaticArrays
using PropCheck
using Random

export generateSystemState, systemStateGen, sampleSystemState


# Helper to convert a PropCheck generator into an SVector
function svector_gen(::Type{T}, len::Int, g) where {T}
    return map(vector(g, len)) do vec
        SVector{len,T}(vec)
    end
end


function systemStateGen(;
        traitRange = (-5.0, 5.0),
        densRange  = (0.0, 20.0),
        nSpecies   = 1:5,
        nAux       = 0:3,
        traitDim   = 1:4,
        stageCls   = 1:3,
        timeRange  = (0.0, 10.0)
    )
    td   = rand(traitDim)
    sc   = rand(stageCls)
    na   = rand(nAux)
    ns   = rand(nSpecies)

    traitGen = PropCheck.ifloat(traitRange...)
    densGen   = PropCheck.ifloat(densRange...)

    traitSV   = svector_gen(Float64, td, traitGen)
    densSV    = svector_gen(Float64, sc, densGen)

    speciesG  = map(traitSV, densSV) do tr, dn
        EcoEvoCore.Species(
            EcoEvoCore.Phenotype{Float64,td}(tr),
            EcoEvoCore.Population{Float64,sc}(dn)
        )
    end

    auxSV = svector_gen(Float64, na, densGen)
    auxG  = map(auxSV) do d
        EcoEvoCore.Population{Float64,na}(d)
    end

    timeG = PropCheck.ifloat(timeRange...)

    return map(vector(speciesG, ns), vector(auxG, na), timeG) do spVec, auxVec, t
        EcoEvoCore.SystemState{Float64,td,sc,na}(spVec, auxVec, t)
    end
end

function generateSystemState(; kwargs...)
    gen = systemStateGen(; kwargs...)
    return generate(gen)
end

# Simple random sampler that does not rely on PropCheck generators. Useful for lightweight
# randomized tests where shrinking is not required.
function sampleSystemState(; 
        traitRange = (-5.0, 5.0),
        densRange  = (0.0, 20.0),
        nSpecies   = 1:5,
        nAux       = 0:3,
        traitDim   = 1:4,
        stageCls   = 1:3,
        timeRange  = (0.0, 10.0)
    )
    td = rand(traitDim)
    sc = rand(stageCls)
    na = rand(nAux)
    ns = rand(nSpecies)

    species = Vector{EcoEvoCore.Species{Float64,td,sc}}(undef, ns)
    for i in 1:ns
        tr = SVector{td,Float64}(rand(Float64, td) .* (traitRange[2]-traitRange[1]) .+ traitRange[1])
        dn = SVector{sc,Float64}(rand(Float64, sc) .* (densRange[2]-densRange[1]) .+ densRange[1])
        species[i] = EcoEvoCore.Species(EcoEvoCore.Phenotype{Float64,td}(tr), EcoEvoCore.Population{Float64,sc}(dn))
    end

    aux = Vector{EcoEvoCore.Population{Float64,na}}(undef, na)
    for j in 1:na
        d = SVector{na,Float64}(rand(Float64, na) .* (densRange[2]-densRange[1]) .+ densRange[1])
        aux[j] = EcoEvoCore.Population{Float64,na}(d)
    end

    t = rand()*(timeRange[2]-timeRange[1]) + timeRange[1]

    EcoEvoCore.SystemState{Float64,td,sc,na}(species, aux, t)
end


end # module TestUtils
