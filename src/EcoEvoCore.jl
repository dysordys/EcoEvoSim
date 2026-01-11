module EcoEvoCore

using StaticArrays

export Phenotype, Population, Species, SystemState, EvoHistory
export trait, traits, density, densities
export nSpecies, nAux, nStageClasses, nTrait
export totalBiomass, addSpecies!, removeSpecies!
export makePhenotype, makePopulation, makeSpecies
export makeSystemState, makeEvoHistory
export packState, unpackState, updateState!



struct Phenotype{T<:Real,TraitDim}
    trait :: SVector{TraitDim,T}
end

struct Population{T<:Real,StageClasses}
    density :: SVector{StageClasses,T}
end

struct Species{T<:Real,TraitDim,StageClasses}
    phenotype :: Phenotype{T,TraitDim}
    pop       :: Population{T,StageClasses}
end

struct SystemState{T<:Real,TraitDim,StageClasses,NAux}
    species :: Vector{Species{T,TraitDim,StageClasses}}
    aux     :: Vector{Population{T,NAux}}
    time    :: T
end

struct EvoHistory{T<:Real,TraitDim,StageClasses,NAux}
    states :: Vector{SystemState{T,TraitDim,StageClasses,NAux}}
    times  :: Vector{T}
end



trait(p::Phenotype) = p.trait

traits(s::Species) = s.phenotype.trait

density(p::Population) = p.density

densities(s::Species) = s.pop.density

nSpecies(st::SystemState) = length(st.species)

nAux(st::SystemState) = length(st.aux)

nStageClasses(::SystemState{T,DT,SC,NA}) where {T,DT,SC,NA} = SC

nTrait(::SystemState{T,DT,SC,NA}) where {T,DT,SC,NA} = DT

totalBiomass(st::SystemState) = sum(sum(densities(s)) for s in st.species)


function addSpecies!(st::SystemState, sp::Species)
    push!(st.species, sp)
    st
end


function removeSpecies!(st::SystemState, idx::Integer)
    splice!(st.species, idx)
end


makePhenotype(traitVec::AbstractVector{T}) where {T<:Real} =
    Phenotype{T,length(traitVec)}(SVector{length(traitVec),T}(traitVec))


makePopulation(densVec::AbstractVector{T}) where {T<:Real} =
    Population{T,length(densVec)}(SVector{length(densVec),T}(densVec))


makeSpecies(traitVec::AbstractVector{T}, densVec::AbstractVector{T}) where {T<:Real} =
    Species(makePhenotype(traitVec), makePopulation(densVec))


function makeSystemState(species::Vector{<:Species},
                         aux::Vector{<:Population},
                         time::Real = zero(eltype(time)))
    T = promote_type(map(s->eltype(s.phenotype.trait), species)...,
                     map(a->eltype(a.density), aux)...,
                     typeof(time))
    TraitDim = length(first(species).phenotype.trait)
    StageCls = length(first(species).pop.density)
    NAux = isempty(aux) ? 0 : length(first(aux).density)
    SystemState{T,TraitDim,StageCls,NAux}(species, aux, T(time))
end


function makeEvoHistory(::Type{T},
                        ::Val{TraitDim},
                        ::Val{StageCls},
                        ::Val{NAux}) where {T,TraitDim,StageCls,NAux}
    SystemT = SystemState{T,TraitDim,StageCls,NAux}
    EvoHistory{T,TraitDim,StageCls,NAux}(Vector{SystemT}(), Vector{T}())
end


function packState(st::SystemState{T,DT,SC,NA}) where {T,DT,SC,NA}
    speciesBlock = isempty(st.species) ? T[] : reduce(vcat, (densities(s) for s in st.species))
    auxBlock = isempty(st.aux) ? T[] : reduce(vcat, (density(a) for a in st.aux))
    u = vcat(speciesBlock, auxBlock)
    traitMat = isempty(st.species) ? Array{T,2}(undef, DT, 0) : hcat([traits(s) for s in st.species]...)
    p = (traits = traitMat,)
    u, p
end


function unpackState(u::AbstractVector{T},
                     old::SystemState{T,DT,SC,NA}) where {T,DT,SC,NA}
    nSp = length(old.species)
    nAux = length(old.aux)
    speciesLen = nSp * SC
    speciesVec = @view u[1:speciesLen]
    auxVec = @view u[speciesLen+1:end]
    newSpecies = Vector{Species{T,DT,SC}}(undef, nSp)
    for i in 1:nSp
        start = (i-1)*SC + 1
        stop = i*SC
        dens_i = SVector{SC,T}(speciesVec[start:stop])
        tr_i = SVector{DT,T}(old.species[i].phenotype.trait)
        newSpecies[i] = Species(Phenotype{T,DT}(tr_i), Population{T,SC}(dens_i))
    end
    newAux = Vector{Population{T,NA}}(undef, nAux)
    for j in 1:nAux
        start = (j-1)*NA + 1
        stop = j*NA
        dens_j = SVector{NA,T}(auxVec[start:stop])
        newAux[j] = Population{T,NA}(dens_j)
    end
    SystemState{T,DT,SC,NA}(newSpecies, newAux, old.time)
end


function updateState!(st::SystemState{T,DT,SC,NA},
                      u::AbstractVector{T}) where {T,DT,SC,NA}
    nSp = length(st.species)
    nAux = length(st.aux)
    speciesLen = nSp * SC
    @assert length(u) == speciesLen + nAux * NA
    for i in 1:nSp
        start = (i-1)*SC + 1
        stop = i*SC
        newdens = SVector{SC,T}(u[start:stop])
        st.species[i] = Species(st.species[i].phenotype,
                                Population{T,SC}(newdens))
    end
    offset = speciesLen
    for j in 1:nAux
        start = offset + (j-1)*NA + 1
        stop = offset + j*NA
        newdens = SVector{NA,T}(u[start:stop])
        st.aux[j] = Population{T,NA}(newdens)
    end
    st
end


import Base: show
function show(io::IO, st::SystemState{T,DT,SC,NA}) where {T,DT,SC,NA}
    println(io, "SystemState (t = $(st.time))")
    for (i, sp) in enumerate(st.species)
        println(io, "  Species $i:")
        println(io, "    traits  = $(sp.phenotype.trait)")
        println(io, "    density = $(sp.pop.density)")
    end
    for (j, aux) in enumerate(st.aux)
        println(io, "  Aux $j density = $(aux.density)")
    end
end


end  # module EcoEvoCore
