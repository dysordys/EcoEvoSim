# Add a species to an existing community

function addSpecies(comm::Community, sp::Species)
    newSpeciesList = [speciesList(comm); sp]
    Community(newSpeciesList, auxs(comm), comm.time)
end


function addSpecies(sp::Species, comm::Community)
    addSpecies(comm, sp)
end


function addSpecies(comm::Community, species_vec::Vector)
    newSpeciesList = [speciesList(comm); species_vec]
    Community(newSpeciesList, auxs(comm), comm.time)
end


function addSpecies(species_vec::Vector, comm::Community)
    addSpecies(comm, species_vec)
end



# Change population sizes in a community

function changePopsizes(
        comm::Community{T, AuxClasses},
        newPopsizes::AbstractVector{T}
    ) where {T<:Real, AuxClasses}
    # Determine if single stage class by checking first species
    numSp = numSpecies(comm)
    numSp > 0 || throw(ArgumentError("Cannot change popsizes in empty community"))
    stageClasses = length(popsizes(comm, 1))
    stageClasses == 1 || throw(ArgumentError(
        "For multiple stage classes, use a matrix for newPopsizes"
    ))
    length(newPopsizes) == numSp || throw(ArgumentError(
        "newPopsizes must have length $numSp (number of species in community)"
    ))

    # Create new species with updated population sizes
    oldSpecies = speciesList(comm)
    newSpecies = Species{T}[]
    for i in 1:numSp
        newPopsize = PopulationSize(newPopsizes[i])
        newSp = Species{T}([newPopsize], oldSpecies[i].trait)
        push!(newSpecies, newSp)
    end

    Community(newSpecies, auxs(comm), comm.time)
end


function changePopsizes(
        comm::Community{T, AuxClasses},
        newPopsizes::AbstractMatrix{T}
    ) where {T<:Real, AuxClasses}
    numSp = numSpecies(comm)
    size(newPopsizes, 1) == numSp || throw(ArgumentError(
        "newPopsizes must have $numSp rows (number of species in community)"
    ))
    # Determine stage classes from first species
    numSp > 0 || throw(ArgumentError("Cannot change popsizes in empty community"))
    stageClasses = length(popsizes(comm, 1))
    size(newPopsizes, 2) == stageClasses || throw(ArgumentError(
        "newPopsizes must have $stageClasses columns (number of stage classes)"
    ))

    # Create new species with updated population sizes
    oldSpecies = speciesList(comm)
    newSpecies = Species{T}[]
    for i in 1:numSp
        newPopsize = PopulationSize(Vector{T}(newPopsizes[i, :]))
        newSp = Species{T}([newPopsize], oldSpecies[i].trait)
        push!(newSpecies, newSp)
    end

    Community(newSpecies, auxs(comm), comm.time)
end



# Change traits in a community

function changeTraits(
        comm::Community{T, AuxClasses},
        newTraits::AbstractVector{T}
    ) where {T<:Real, AuxClasses}
    # Determine if single trait dimension by checking first species
    numSp = numSpecies(comm)
    numSp > 0 || throw(ArgumentError("Cannot change traits in empty community"))
    traitDim = length(traits(comm, 1))
    traitDim == 1 || throw(ArgumentError(
        "For multiple trait dimensions, use a matrix for newTraits"
    ))
    length(newTraits) == numSp || throw(ArgumentError(
        "newTraits must have length $numSp (number of species in community)"
    ))

    # Create new species with updated traits
    oldSpecies = speciesList(comm)
    newSpecies = Species{T}[]
    for i in 1:numSp
        newTrait = Phenotype(newTraits[i])
        newSp = Species{T}(oldSpecies[i].popsize, [newTrait])
        push!(newSpecies, newSp)
    end

    Community(newSpecies, auxs(comm), comm.time)
end


function changeTraits(
        comm::Community{T, AuxClasses},
        newTraits::AbstractMatrix{T}
    ) where {T<:Real, AuxClasses}
    numSp = numSpecies(comm)
    size(newTraits, 1) == numSp || throw(ArgumentError(
        "newTraits must have $numSp rows (number of species in community)"
    ))
    # Determine trait dimension from first species
    numSp > 0 || throw(ArgumentError("Cannot change traits in empty community"))
    traitDim = length(traits(comm, 1))
    size(newTraits, 2) == traitDim || throw(ArgumentError(
        "newTraits must have $traitDim columns (trait space dimension)"
    ))

    # Create new species with updated traits
    oldSpecies = speciesList(comm)
    newSpecies = Species{T}[]
    for i in 1:numSp
        newTrait = Phenotype(Vector{T}(newTraits[i, :]))
        newSp = Species{T}(oldSpecies[i].popsize, [newTrait])
        push!(newSpecies, newSp)
    end

    Community(newSpecies, auxs(comm), comm.time)
end



# Remove species from a community

function removeSpecies(comm::Community, index::Integer)
    numSp = numSpecies(comm)
    1 <= index <= numSp || throw(ArgumentError(
        "Species index $index out of bounds (community has $numSp species)"
    ))

    oldSpecies = speciesList(comm)
    newSpecies = [oldSpecies[i] for i in 1:numSp if i != index]
    Community(newSpecies, auxs(comm), comm.time)
end


function removeSpecies(comm::Community, indices)
    numSp = numSpecies(comm)
    oldSpecies = speciesList(comm)

    # Convert to Set for efficient lookup and validate indices
    indices_set = Set{Int}()
    for idx in indices
        1 <= idx <= numSp || throw(ArgumentError(
            "Species index $idx out of bounds (community has $numSp species)"
        ))
        push!(indices_set, idx)
    end

    # Keep species not in the removal set
    newSpecies = [oldSpecies[i] for i in 1:numSp if i ∉ indices_set]
    Community(newSpecies, auxs(comm), comm.time)
end



# Find species below extinction threshold

function speciesBelowThreshold(comm::Community, extThreshold::Real)
    extThreshold >= 0.0 || throw(ArgumentError(
        "Extinction threshold must be non-negative"
    ))

    numSp = numSpecies(comm)
    if numSp == 0
        return Int[]  # Empty community, no species below threshold
    end

    # Find species below threshold
    indices_below = Int[]
    for i in 1:numSp
        popsizeVec = popsizes(comm, i)
        total = sum(popsizeVec)
        if total < extThreshold
            push!(indices_below, i)
        end
    end

    return indices_below
end



# Remove extinct species (those below the extinction threshold)

function removeExtinct(comm::Community, extThreshold::Real)
    indices_to_remove = speciesBelowThreshold(comm, extThreshold)

    if isempty(indices_to_remove)
        return comm  # No species to remove
    else
        return removeSpecies(comm, indices_to_remove)
    end
end
