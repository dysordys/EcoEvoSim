# Alternative outer constructors for basic types

PopulationSize(popsize::T) where {T<:Real} =
    PopulationSize{T, 1}(SVector{1, T}(popsize))


PopulationSize(popsizeVec::AbstractVector{T}) where {T<:Real} =
    PopulationSize{T, length(popsizeVec)}(SVector{length(popsizeVec), T}(popsizeVec))


Phenotype(trait::T) where {T<:Real} =
    Phenotype{T, 1}(SVector{1, T}(trait))


Phenotype(traitVec::AbstractVector{T}) where {T<:Real} =
    Phenotype{T, length(traitVec)}(SVector{length(traitVec), T}(traitVec))


Species(popsizeVal::T, traitVal::T) where {T<:Real} =
    Species{T, 1, 1}([PopulationSize(popsizeVal)], [Phenotype(traitVal)])


Species(popsizeVec::AbstractVector{T}, traitVal::T) where {T<:Real} = begin
    stageClasses = length(popsizeVec[1].popsize)
    Species{T, stageClasses, 1}(popsizeVec, traitVal)
end


Species(popsizeVal::T, traitVec::AbstractVector{T}) where {T<:Real} = begin
    traitDim = length(traitVec[1].trait)
    Species{T, 1, traitDim}(popsizeVal, traitVec)
end


Species(popsizeVec::AbstractVector{T}, traitVec::AbstractVector{T}) where {T<:Real} = begin
    stageClasses = length(popsizeVec[1].popsize)
    traitDim = length(traitVec[1].trait)
    Species{T, stageClasses, traitDim}(popsizeVec, traitVec)
end


function Community{T, StageClasses, TraitDim, AuxClasses}(
        species::Vector{Species{T, StageClasses, TraitDim}},
        aux::Vector{PopulationSize{T, 1}}
    ) where {T<:Real, StageClasses, TraitDim, AuxClasses}
    Community{T, StageClasses, TraitDim, AuxClasses}(
        species, aux, Community{T, StageClasses, TraitDim, AuxClasses}[]
    )
end


function Community(
        popVals::AbstractVector{T},
        traitVals::AbstractVector{T},
        auxVals::AbstractVector{T}
    ) where {T<:Real}
    length(popVals) == length(traitVals) ||
        throw(ArgumentError("`popVals` and `traitVals` must have the same length"))
    stageClasses = length(popVals)
    traitDim = length(traitVals)
    sp = Species{T, stageClasses, traitDim}(
        [PopulationSize{T, stageClasses}(SVector{stageClasses, T}(popVals))],
        [Phenotype{T, traitDim}(SVector{traitDim, T}(traitVals))]
    )
    aux = [PopulationSize{T, 1}(SVector{1, T}(a)) for a in auxVals]
    AuxClasses = length(aux)
    Community{T, stageClasses, traitDim, AuxClasses}(
        [sp], aux, Community{T, stageClasses, traitDim, AuxClasses}[]
    )
end


# Selectors
species(comm::Community) = comm.species
aux(comm::Community) = comm.aux
history(comm::Community) = comm.history
