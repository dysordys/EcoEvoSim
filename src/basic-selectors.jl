# Functions for extracting species from a community

"""
    numSpecies(comm::Community)

Return the number of species in a community.

# Example
```julia
comm = Community([1.0, 2.0, 3.0], [0.1, 0.2, 0.3], Float64[])
n = numSpecies(comm)  # Returns 3
```
"""
numSpecies(comm::Community) = length(speciesList(comm))


"""
    speciesIndices(comm::Community)

Return a range of all species indices in a community.

# Example
```julia
comm = Community([1.0, 2.0, 3.0], [0.1, 0.2, 0.3], Float64[])
indices = speciesIndices(comm)  # Returns 1:3
```
"""
speciesIndices(comm::Community) = 1:numSpecies(comm)


"""
    speciesList(comm::Community)
    speciesList(comm::Community, i::Integer)
    speciesList(comm::Community, indices)

Extract species from a community.

# Arguments
- `comm::Community`: The community
- `i::Integer`: Single species index (optional)
- `indices`: Collection of species indices (optional)

# Returns
- Without index: Vector of all species
- With single index: Single species
- With multiple indices: Vector of specified species

# Example
```julia
comm = Community([1.0, 2.0, 3.0], [0.1, 0.2, 0.3], Float64[])
all_sp = speciesList(comm)
sp1 = speciesList(comm, 1)
subset = speciesList(comm, [1, 3])
```
"""
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

"""
    popsizes(comm::Community)
    popsizes(comm::Community, i::Integer)
    popsizes(comm::Community, indices)

Extract population sizes from species in a community.

# Arguments
- `comm::Community`: The community
- `i::Integer`: Single species index (optional)
- `indices`: Collection of species indices (optional)

# Returns
- Without index: Vector of population size vectors for all species
- With single index: Population size vector for one species
- With multiple indices: Vector of population size vectors

# Example
```julia
comm = Community([1.0, 2.0, 3.0], [0.1, 0.2, 0.3], Float64[])
all_pops = popsizes(comm)
pop1 = popsizes(comm, 1)
subset = popsizes(comm, [1, 3])
```
"""
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

"""
    traitSpaceDim(comm::Community)

Return the dimension of the trait space (number of trait dimensions).

# Example
```julia
comm = Community([1.0, 2.0], [0.1, 0.2], Float64[])
dim = traitSpaceDim(comm)  # Returns 1
```
"""
function traitSpaceDim(comm::Community)
    numSp = numSpecies(comm)
    numSp > 0 ||
        throw(ArgumentError("Cannot determine trait space dimension from empty community"))
    length(traits(comm, 1))
end


"""
    traits(comm::Community)
    traits(comm::Community, i::Integer)
    traits(comm::Community, indices)

Extract trait values from species in a community.

# Arguments
- `comm::Community`: The community
- `i::Integer`: Single species index (optional)
- `indices`: Collection of species indices (optional)

# Returns
- Without index: Vector of trait vectors for all species
- With single index: Trait vector for one species
- With multiple indices: Vector of trait vectors

# Example
```julia
comm = Community([1.0, 2.0, 3.0], [0.1, 0.2, 0.3], Float64[])
all_traits = traits(comm)
trait1 = traits(comm, 1)
subset = traits(comm, [1, 3])
```
"""
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

"""
    auxs(comm::Community)
    auxs(comm::Community, i::Integer)
    auxs(comm::Community, indices)

Extract auxiliary variables (e.g., resources) from a community.

# Arguments
- `comm::Community`: The community
- `i::Integer`: Single auxiliary variable index (optional)
- `indices`: Collection of indices (optional)

# Returns
- Without index: Vector of all auxiliary variables
- With single index: Single auxiliary variable
- With multiple indices: Vector of specified auxiliary variables

# Example
```julia
# Community with auxiliary variables
species = [Species(1.0, 0.1), Species(2.0, 0.2)]
aux = [PopulationSize(100.0)]
comm = Community(species, aux, 0.0)
all_aux = auxs(comm)
aux1 = auxs(comm, 1)
```
"""
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

"""
    randomSpecies(comm::Community)
    randomSpecies(comm::Community, n::Integer)

Randomly select species indices from a community with uniform probability.

# Arguments
- `comm::Community`: The community
- `n::Integer`: Number of species to select (optional, defaults to 1)

# Returns
- Without n: Single random species index
- With n: Sorted vector of n unique random species indices

# Example
```julia
comm = Community([1.0, 2.0, 3.0], [0.1, 0.2, 0.3], Float64[])
idx = randomSpecies(comm)
indices = randomSpecies(comm, 2)
```
"""
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


"""
    weightedRandomSpecies(comm::Community)

Randomly select a species index with probability proportional to population size.

Species with larger population sizes have higher probability of being selected.
Sums total population across all stage classes for each species.

# Returns
Single species index (weighted by population size)

# Example
```julia
comm = Community([1.0, 10.0, 5.0], [0.1, 0.2, 0.3], Float64[])
idx = weightedRandomSpecies(comm)  # More likely to select species 2
```
"""
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


"""
    speciesBelowThreshold(comm::Community, extThreshold::Real)

Find all species with total population size below the specified threshold.

# Arguments
- `comm::Community`: The community
- `extThreshold::Real`: Extinction threshold (non-negative)

# Returns
Vector of indices for species below the threshold

# Example
```julia
comm = Community([0.001, 2.0, 0.005], [0.1, 0.2, 0.3], Float64[])
extinct_indices = speciesBelowThreshold(comm, 0.01)  # Returns [1, 3]
```
"""
