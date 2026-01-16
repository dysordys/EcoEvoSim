PopulationSize(mat::AbstractMatrix{T}) where {T<:Real} =
    [PopulationSize(vec) for vec in eachrow(mat)]


Phenotype(mat::AbstractMatrix{T}) where {T<:Real} =
    [Phenotype(vec) for vec in eachrow(mat)]


PopulationSize(popsize::T) where {T<:Real} =
    PopulationSize{T}([popsize])


PopulationSize(popsizeVec::AbstractVector{T}) where {T<:Real} =
    PopulationSize{T}(Vector{T}(popsizeVec))


Phenotype(trait::T) where {T<:Real} =
    Phenotype{T}([trait])


Phenotype(traitVec::AbstractVector{T}) where {T<:Real} =
    Phenotype{T}(Vector{T}(traitVec))


function Species(
        popsize::PopulationSize{T},
        phenotype::Phenotype{T}
    ) where {T<:Real}
    Species{T}([popsize], [phenotype])
end


Species(popsizeVal::T, traitVal::T) where {T<:Real} =
    Species{T}([PopulationSize(popsizeVal)], [Phenotype(traitVal)])


Species(popsizeVec::AbstractVector{T}, traitVal::T) where {T<:Real} =
    Species{T}([PopulationSize(popsizeVec)], [Phenotype(traitVal)])


Species(popsizeVal::T, traitVec::AbstractVector{T}) where {T<:Real} =
    Species{T}([PopulationSize(popsizeVal)], [Phenotype(traitVec)])


Species(popsizeVec::AbstractVector{T}, traitVec::AbstractVector{T}) where {T<:Real} =
    Species{T}([PopulationSize(popsizeVec)], [Phenotype(traitVec)])


function Species(popMat::AbstractMatrix{T}, traitMat::AbstractMatrix{T}) where {T<:Real}
    n_species = size(popMat, 1)
    n_species == size(traitMat, 1) ||
        throw(ArgumentError("popMat and traitMat must have the same number of rows (species)"))
    [Species{T}(
        [PopulationSize{T}(Vector{T}(popMat[i, :]))],
        [Phenotype{T}(Vector{T}(traitMat[i, :]))]
    ) for i in 1:n_species]
end


function Species(
    popVecs::Vector{<:AbstractVector{T}}, traitVecs::Vector{<:AbstractVector{T}}
) where {T<:Real}
    n_species = length(popVecs)
    n_species == length(traitVecs) ||
        throw(ArgumentError("popVecs and traitVecs must have the same length"))
    [
        Species{T}([PopulationSize(popVecs[i])], [Phenotype(traitVecs[i])])
        for i in 1:n_species
    ]
end


function Community(
        species::Vector{Species{T}},
        aux::Vector{PopulationSize{T}},
        time::T = zero(T)
    ) where {T<:Real}
    AuxClasses = length(aux)
    Community{T, AuxClasses}(species, aux, time)
end


function Community(
        popVals::AbstractVector{T},
        traitVals::AbstractVector{T},
        auxVals::AbstractVector{T}
    ) where {T<:Real}
    length(popVals) == length(traitVals) ||
        throw(ArgumentError("`popVals` and `traitVals` must have the same length"))
    # Create vector of species where each species has one stage class and one trait value
    species = [Species{T}(
        [PopulationSize{T}([popVals[i]])],
        [Phenotype{T}([traitVals[i]])]
    ) for i in 1:length(popVals)]
    aux = [PopulationSize{T}([a]) for a in auxVals]
    AuxClasses = length(aux)
    Community{T, AuxClasses}(
        species, aux, zero(T)
    )
end


function EvoHistory(
        comm::Community{T, AuxClasses}
    ) where {T<:Real, AuxClasses}
    EvoHistory{T, AuxClasses}([comm])
end


function EvoHistory(
        comms::Vector{Community{T, AuxClasses}}
    ) where {T<:Real, AuxClasses}
    EvoHistory{T, AuxClasses}(comms)
end


function emptyCommunity(T::Type{<:Real} = Float64)
    Community{T, 0}(Species{T}[], PopulationSize{T}[], zero(T))
end


function emptyEvoHistory(T::Type{<:Real} = Float64)
    EvoHistory{T, 0}([emptyCommunity(T)])
end
