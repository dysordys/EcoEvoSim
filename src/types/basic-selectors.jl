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


function popsizes(comm::Community, i::Integer)
    numSp = length(speciesList(comm))
    1 <= i <= numSp || throw(ArgumentError(
        "Species index $i out of bounds (community has $numSp species)"
    ))
    speciesList(comm)[i].popsize
end


function popsizes(comm::Community)
    [sp.popsize for sp in speciesList(comm)]
end


function popsizes(comm::Community, indices)
    speciesVec = speciesList(comm)
    numSp = length(speciesVec)
    result = Vector{typeof(speciesVec[1].popsize)}()
    for i in indices
        1 <= i <= numSp || throw(ArgumentError(
            "Species index $i out of bounds (community has $numSp species)"
        ))
        push!(result, speciesVec[i].popsize)
    end
    return result
end


function traits(comm::Community, i::Integer)
    numSp = length(speciesList(comm))
    1 <= i <= numSp || throw(ArgumentError(
        "Species index $i out of bounds (community has $numSp species)"
    ))
    speciesList(comm)[i].trait
end


function traits(comm::Community)
    [sp.trait for sp in speciesList(comm)]
end


function traits(comm::Community, indices)
    speciesVec = speciesList(comm)
    numSp = length(speciesVec)
    result = Vector{typeof(speciesVec[1].trait)}()
    for i in indices
        1 <= i <= numSp || throw(ArgumentError(
            "Species index $i out of bounds (community has $numSp species)"
        ))
        push!(result, speciesVec[i].trait)
    end
    return result
end


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
