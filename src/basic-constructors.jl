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
        species::Vector{Species{T}},
        aux::AbstractVector,
        time::T = zero(T)
    ) where {T<:Real}
    # Convert an untyped/heterogeneous aux vector to the expected
    # Vector{PopulationSize{T}}. This handles empty vectors produced by
    # comprehensions (which may have element type `Any`).
    n = length(aux)
    aux_typed = Vector{PopulationSize{T}}(undef, n)
    for (i, a) in enumerate(aux)
        if a isa PopulationSize{T}
            aux_typed[i] = a
        elseif a isa PopulationSize
            aux_typed[i] = PopulationSize{T}(Vector{T}(a.popsize))
        elseif a isa Real
            aux_typed[i] = PopulationSize{T}(a)
        elseif a isa AbstractVector
            aux_typed[i] = PopulationSize{T}(Vector{T}(a))
        else
            throw(ArgumentError("Cannot convert auxiliary value to PopulationSize{T}"))
        end
    end
    AuxClasses = n
    Community{T, AuxClasses}(species, aux_typed, time)
end


function Community(
        popVals::AbstractVector{T},
        traitVals::AbstractVector{T}
    ) where {T<:Real}
    length(popVals) == length(traitVals) ||
        throw(ArgumentError("`popVals` and `traitVals` must have the same length"))
    # Create vector of species where each species has one stage class and one trait value
    species = [Species{T}(
        [PopulationSize{T}([popVals[i]])],
        [Phenotype{T}([traitVals[i]])]
    ) for i in eachindex(popVals)]
    Community{T, 0}(species, PopulationSize{T}[], zero(T))
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
    ) for i in eachindex(popVals)]
    aux = [PopulationSize{T}([a]) for a in auxVals]
    AuxClasses = length(aux)
    Community{T, AuxClasses}(
        species, aux, zero(T)
    )
end


function Community(
        popMat::AbstractMatrix{T},
        traitVec::AbstractVector{T}
    ) where {T<:Real}
    n_species = size(popMat, 1)
    n_species == length(traitVec) ||
        throw(ArgumentError("popMat must have the same number of rows as the length of" *
                            "traitVec (one trait value per species)"))
    # Create vector of species where each species has stage classes from a row of popMat
    # and a single trait value from traitVec
    species = [Species{T}(
        [PopulationSize{T}(Vector{T}(popMat[i, :]))],
        [Phenotype{T}([traitVec[i]])]
    ) for i in 1:n_species]
    Community{T, 0}(species, PopulationSize{T}[], zero(T))
end


function Community(
        popMat::AbstractMatrix{T},
        traitVec::AbstractVector{T},
        auxVals::AbstractVector{T}
    ) where {T<:Real}
    n_species = size(popMat, 1)
    n_species == length(traitVec) ||
        throw(ArgumentError("popMat must have the same number of rows as the length of " *
                            "traitVec (one trait value per species)"))
    # Create vector of species where each species has stage classes from a row of popMat
    # and a single trait value from traitVec
    species = [Species{T}(
        [PopulationSize{T}(Vector{T}(popMat[i, :]))],
        [Phenotype{T}([traitVec[i]])]
    ) for i in 1:n_species]
    aux = [PopulationSize{T}([a]) for a in auxVals]
    AuxClasses = length(aux)
    Community{T, AuxClasses}(species, aux, zero(T))
end


function Community(
        popVec::AbstractVector{T},
        traitMat::AbstractMatrix{T}
    ) where {T<:Real}
    n_species = size(traitMat, 1)
    n_species == length(popVec) ||
        throw(ArgumentError("popVec must have the same length as the number of rows " *
                            "of traitMat (one population value per species)"))
    # Create vector of species where each species has a single population value
    # and trait dimensions from a row of traitMat
    species = [Species{T}(
        [PopulationSize{T}([popVec[i]])],
        [Phenotype{T}(Vector{T}(traitMat[i, :]))]
    ) for i in 1:n_species]
    Community{T, 0}(species, PopulationSize{T}[], zero(T))
end


function Community(
        popVec::AbstractVector{T},
        traitMat::AbstractMatrix{T},
        auxVals::AbstractVector{T}
    ) where {T<:Real}
    n_species = size(traitMat, 1)
    n_species == length(popVec) ||
        throw(ArgumentError("popVec must have the same length as the number of rows " *
                            "of traitMat (one population value per species)"))
    # Create vector of species where each species has a single population value
    # and trait dimensions from a row of traitMat
    species = [Species{T}(
        [PopulationSize{T}([popVec[i]])],
        [Phenotype{T}(Vector{T}(traitMat[i, :]))]
    ) for i in 1:n_species]
    aux = [PopulationSize{T}([a]) for a in auxVals]
    AuxClasses = length(aux)
    Community{T, AuxClasses}(species, aux, zero(T))
end


function Community(
        popMat::AbstractMatrix{T},
        traitMat::AbstractMatrix{T}
    ) where {T<:Real}
    n_species = size(popMat, 1)
    n_species == size(traitMat, 1) ||
        throw(ArgumentError("popMat and traitMat must have the same number of rows (species)"))
    # Create vector of species where each species has stage classes from a row of popMat
    # and trait dimensions from a row of traitMat
    species = [Species{T}(
        [PopulationSize{T}(Vector{T}(popMat[i, :]))],
        [Phenotype{T}(Vector{T}(traitMat[i, :]))]
    ) for i in 1:n_species]
    Community{T, 0}(species, PopulationSize{T}[], zero(T))
end


function Community(
        popMat::AbstractMatrix{T},
        traitMat::AbstractMatrix{T},
        auxVals::AbstractVector{T}
    ) where {T<:Real}
    n_species = size(popMat, 1)
    n_species == size(traitMat, 1) ||
        throw(ArgumentError("popMat and traitMat must have the same number of rows (species)"))
    # Create vector of species where each species has stage classes from a row of popMat
    # and trait dimensions from a row of traitMat
    species = [Species{T}(
        [PopulationSize{T}(Vector{T}(popMat[i, :]))],
        [Phenotype{T}(Vector{T}(traitMat[i, :]))]
    ) for i in 1:n_species]
    aux = [PopulationSize{T}([a]) for a in auxVals]
    AuxClasses = length(aux)
    Community{T, AuxClasses}(species, aux, zero(T))
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


"""
    emptyCommunity(T::Type{<:Real} = Float64)

Create an empty community with no species and no auxiliary variables.

# Arguments
- `T::Type{<:Real}`: Numeric type for the community (default: `Float64`)

# Returns
- `Community{T, 0}`: Empty community at time zero

# Example
```julia
comm = emptyCommunity()
comm = emptyCommunity(Float32)
```
"""
function emptyCommunity(T::Type{<:Real} = Float64)
    Community{T, 0}(Species{T}[], PopulationSize{T}[], zero(T))
end


"""
    emptyEvoHistory(T::Type{<:Real} = Float64)

Create an empty evolutionary history containing only an empty community.

# Arguments
- `T::Type{<:Real}`: Numeric type for the history (default: `Float64`)

# Returns
- `EvoHistory{T, 0}`: History with a single empty community

# Example
```julia
history = emptyEvoHistory()
history = emptyEvoHistory(Float32)
```
"""
function emptyEvoHistory(T::Type{<:Real} = Float64)
    EvoHistory{T, 0}([emptyCommunity(T)])
end
