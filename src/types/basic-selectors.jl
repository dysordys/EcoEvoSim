allSpecies(comm::Community) = comm.species
allAux(comm::Community) = comm.aux


function species(comm::Community, i::Integer)
    1 <= i <= length(allSpecies(comm)) || throw(ArgumentError(
        "Species index $i out of bounds (community has $(length(allSpecies(comm))) species)"
    ))
    allSpecies(comm)[i]
end
