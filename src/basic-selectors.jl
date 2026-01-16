# Functions for extracting species from a community

numSpecies(comm::Community) = length(speciesList(comm))


speciesIndices(comm::Community) = 1:numSpecies(comm)


speciesList(comm::Community) = comm.species


function speciesList(comm::Community, i::Integer)
    numSp = length(speciesList(comm))
    1 <= i <= numSp || throw(ArgumentError(
        "Species index $i out of bounds (community has $numSp species)"
    ))
    speciesList(comm)[i]
end


function speciesList(comm::Community, indices)
    speciesVec = speciesList(comm)
    numSp = length(speciesVec)
    result = eltype(speciesVec)[]
    for i in indices
        1 <= i <= numSp || throw(ArgumentError(
            "Species index $i out of bounds (community has $numSp species)"
        ))
        push!(result, speciesVec[i])
    end
    return result
end



# Functions for extracting species' population sizes from a community

function popsizes(comm::Community, i::Integer)
    numSp = length(speciesList(comm))
    1 <= i <= numSp || throw(ArgumentError(
        "Species index $i out of bounds (community has $numSp species)"
    ))
    speciesList(comm)[i].popsize[1].popsize
end


function popsizes(comm::Community)
    [sp.popsize[1].popsize for sp in speciesList(comm)]
end


function popsizes(comm::Community, indices)
    speciesVec = speciesList(comm)
    numSp = length(speciesVec)
    result = Vector{typeof(speciesVec[1].popsize[1].popsize)}()
    for i in indices
        1 <= i <= numSp || throw(ArgumentError(
            "Species index $i out of bounds (community has $numSp species)"
        ))
        push!(result, speciesVec[i].popsize[1].popsize)
    end
    return result
end



# Functions for extracting species' traits from a community

function traitSpaceDim(comm::Community)
    numSp = numSpecies(comm)
    numSp > 0 ||
        throw(ArgumentError("Cannot determine trait space dimension from empty community"))
    length(traits(comm, 1))
end


function traits(comm::Community, i::Integer)
    numSp = length(speciesList(comm))
    1 <= i <= numSp || throw(ArgumentError(
        "Species index $i out of bounds (community has $numSp species)"
    ))
    speciesList(comm)[i].trait[1].trait
end


function traits(comm::Community)
    [sp.trait[1].trait for sp in speciesList(comm)]
end


function traits(comm::Community, indices)
    speciesVec = speciesList(comm)
    numSp = length(speciesVec)
    result = Vector{typeof(speciesVec[1].trait[1].trait)}()
    for i in indices
        1 <= i <= numSp || throw(ArgumentError(
            "Species index $i out of bounds (community has $numSp species)"
        ))
        push!(result, speciesVec[i].trait[1].trait)
    end
    return result
end



# Functions for extracting auxiliary variables from a community

auxs(comm::Community) = comm.aux


function auxs(comm::Community, i::Integer)
    numAux = length(auxs(comm))
    1 <= i <= numAux || throw(ArgumentError(
        "Auxiliary variable index $i out of bounds (community has $numAux auxiliary variables)"
    ))
    auxs(comm)[i]
end


function auxs(comm::Community, indices)
    auxVec = auxs(comm)
    numAux = length(auxVec)
    result = eltype(auxVec)[]
    for i in indices
        1 <= i <= numAux || throw(ArgumentError(
            "Auxiliary variable index $i out of bounds (community has $numAux auxiliary variables)"
        ))
        push!(result, auxVec[i])
    end
    return result
end



# Random selection of species indices

function randomSpecies(comm::Community, n::Integer)
    numSp = numSpecies(comm)
    n >= 0 || throw(ArgumentError("n must be non-negative"))
    n <= numSp || throw(ArgumentError(
        "Cannot sample $n species from community with only $numSp species"
    ))
    sort(randperm(numSp)[1:n])
end


function randomSpecies(comm::Community)
    numSp = numSpecies(comm)
    numSp > 0 || throw(ArgumentError("Cannot sample from empty community"))
    rand(1:numSp)
end


function weightedRandomSpecies(comm::Community)
    numSp = numSpecies(comm)
    numSp > 0 || throw(ArgumentError("Cannot sample from empty community"))

    # Calculate total population size for each species (sum across stage classes)
    totalPops = Float64[]
    for i in 1:numSp
        popsizeVec = popsizes(comm, i)
        total = sum(popsizeVec)
        total >= 0.0 || throw(ArgumentError(
            "Species $i has negative total population size"
        ))
        push!(totalPops, total)
    end

    # Check that at least one species has positive population
    sumPops = sum(totalPops)
    sumPops > 0.0 || throw(ArgumentError(
        "All species have zero population size; cannot perform weighted sampling"
    ))

    # Sample using cumulative distribution
    r = rand() * sumPops
    cumsum = 0.0
    for i in 1:numSp
        cumsum += totalPops[i]
        if r <= cumsum
            return i
        end
    end

    # Fallback (should never reach here due to numerical precision)
    return numSp
end

