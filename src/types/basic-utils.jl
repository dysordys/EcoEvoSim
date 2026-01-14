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
