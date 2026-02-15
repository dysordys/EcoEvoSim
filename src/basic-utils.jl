# Get the length of an EvoHistory

"""
    length(h::EvoHistory)

Return the number of communities stored in the evolutionary history.

# Arguments
- `h::EvoHistory`: The evolutionary history object

# Returns
Number of time points (communities) in the history

# Example
```julia
# Assuming you have an EvoHistory object called 'history'
n = length(history)  # Number of snapshots in history
```
"""
Base.length(h::EvoHistory) = length(historyList(h))



# Filter evolutionary histories

"""
    filterHistory(h::EvoHistory, f::Function)

Filter an evolutionary history by applying a predicate function to each snapshot index.

The predicate function should take a single integer argument (the 1-based snapshot index)
and return `true` if that snapshot should be kept, `false` otherwise.

# Arguments
- `h::EvoHistory`: The evolutionary history to filter
- `f::Function`: A function that takes an integer index and returns a boolean

# Returns
New `EvoHistory` containing only the snapshots where `f(i)` returns `true`

# Example
```julia
# Keep only snapshots at indices 1, 5, 10, and 15
filtered = filterHistory(history, i -> i in [1, 5, 10, 15])

# Keep every other snapshot starting from the first
filtered = filterHistory(history, i -> isodd(i))

# Keep the last 50 snapshots
n = length(history)
filtered = filterHistory(history, i -> i > n - 50)
```
"""
function filterHistory(h::EvoHistory{T, AuxClasses}, f::Function) where {T<:Real, AuxClasses}
    n = length(h)
    selected_indices = [i for i in 1:n if f(i)]

    if length(selected_indices) == 0
        throw(ArgumentError("Filter function selected no snapshots"))
    end

    filtered_communities = historyList(h, selected_indices)
    EvoHistory{T, AuxClasses}(filtered_communities)
end


"""
    filterHistory(h::EvoHistory, indices::AbstractVector{<:Integer})

Filter an evolutionary history by selecting specific snapshot indices.

# Arguments
- `h::EvoHistory`: The evolutionary history to filter
- `indices::AbstractVector{<:Integer}`: Vector of 1-based indices to keep

# Returns
New `EvoHistory` containing only the specified snapshots in the given order

# Example
```julia
# Keep only the first, fifth, and tenth snapshots
filtered = filterHistory(history, [1, 5, 10])

# Keep the first 100 snapshots
filtered = filterHistory(history, 1:100)
```
"""
function filterHistory(h::EvoHistory{T, AuxClasses}, indices::AbstractVector{<:Integer}) where {T<:Real, AuxClasses}
    n = length(h)
    for i in indices
        1 <= i <= n || throw(ArgumentError(
            "Index $i out of bounds (history has $n snapshots)"
        ))
    end

    filtered_communities = historyList(h, indices)
    EvoHistory{T, AuxClasses}(filtered_communities)
end


"""
    filterHistory(h::EvoHistory, range::UnitRange)

Filter an evolutionary history by selecting a contiguous range of snapshots.

# Arguments
- `h::EvoHistory`: The evolutionary history to filter
- `range::UnitRange`: A range of 1-based indices (e.g., `10:50`)

# Returns
New `EvoHistory` containing only the snapshots in the specified range

# Example
```julia
# Keep snapshots from index 10 to 50
filtered = filterHistory(history, 10:50)

# Keep the first 100 snapshots (equivalent to history[1:100])
filtered = filterHistory(history, 1:100)
```
"""
function filterHistory(h::EvoHistory{T, AuxClasses}, range::UnitRange{<:Integer}) where {T<:Real, AuxClasses}
    n = length(h)
    first(range) >= 1 || throw(ArgumentError(
        "Range start $(first(range)) is less than 1"
    ))
    last(range) <= n || throw(ArgumentError(
        "Range end $(last(range)) exceeds history length $n"
    ))

    filtered_communities = historyList(h, range)
    EvoHistory{T, AuxClasses}(filtered_communities)
end


"""
    filterHistory(h::EvoHistory, step::Integer)

Filter an evolutionary history by keeping every `step`-th snapshot.

# Arguments
- `h::EvoHistory`: The evolutionary history to filter
- `step::Integer`: Keep every step-th snapshot (step >= 1)

# Returns
New `EvoHistory` containing every step-th snapshot starting from the first

# Example
```julia
# Keep every other snapshot
filtered = filterHistory(history, 2)

# Keep every 10th snapshot
filtered = filterHistory(history, 10)

# Keep only the first snapshot (step = n where n >= length(history))
filtered = filterHistory(history, length(history))
```
"""
function filterHistory(h::EvoHistory{T, AuxClasses}, step::Integer) where {T<:Real, AuxClasses}
    step >= 1 || throw(ArgumentError("step must be >= 1, got $step"))

    n = length(h)
    selected_indices = 1:step:n
    filtered_communities = historyList(h, selected_indices)
    EvoHistory{T, AuxClasses}(filtered_communities)
end



# Add a species to an existing community

"""
    addSpecies(comm::Community, sp::Species)
    addSpecies(sp::Species, comm::Community)
    addSpecies(comm::Community, species_vec::Vector)
    addSpecies(species_vec::Vector, comm::Community)

Add one or more species to an existing community.

Creates a new community with the added species. Preserves auxiliary variables and time.

# Arguments
- `comm::Community`: The community to add to
- `sp::Species`: Single species to add
- `species_vec::Vector`: Multiple species to add

# Returns
New community with added species

# Example
```julia
comm = Community([1.0, 2.0], [0.1, 0.2], Float64[])
new_sp = Species(1.0, 0.5)
comm2 = addSpecies(comm, new_sp)
```
"""
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

"""
    changePopsizes(comm::Community, newPopsizes::AbstractVector)
    changePopsizes(comm::Community, newPopsizes::AbstractMatrix)

Update population sizes for all species in a community.

# Arguments
- `comm::Community`: The community to update
- `newPopsizes::AbstractVector`: New population sizes (for single stage class)
- `newPopsizes::AbstractMatrix`: New population sizes (rows=species, cols=stage classes)

# Returns
New community with updated population sizes

# Example
```julia
comm = Community([1.0, 2.0, 3.0], [0.1, 0.2, 0.3], Float64[])
# Single stage class
comm2 = changePopsizes(comm, [2.0, 3.0, 4.0])

# Multiple stage classes (2 species, 2 stage classes each)
comm_stages = Community(Species([1.0, 2.0], [0.1, 0.2]), emptyCommunity().aux, 0.0)
comm2 = changePopsizes(comm_stages, [1.5 2.5; 3.5 4.5])
```
"""
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

"""
    changeTraits(comm::Community, newTraits::AbstractVector)
    changeTraits(comm::Community, newTraits::AbstractMatrix)

Update trait values for all species in a community.

# Arguments
- `comm::Community`: The community to update
- `newTraits::AbstractVector`: New trait values (for 1D trait space)
- `newTraits::AbstractMatrix`: New trait values (rows=species, cols=trait dimensions)

# Returns
New community with updated traits

# Example
```julia
comm = Community([1.0, 2.0, 3.0], [0.1, 0.2, 0.3], Float64[])
# 1D trait space
comm2 = changeTraits(comm, [0.15, 0.25, 0.35])

# Multidimensional trait space (2 species, 2 trait dimensions)
comm_multi = Community([1.0, 2.0], [[0.1, 0.5], [0.2, 0.6]])
comm2 = changeTraits(comm_multi, [0.15 0.55; 0.25 0.65])
```
"""
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



# Extract a single trait dimension

"""
    selectTraitDim(comm::Community, dimIndex::Integer)

Extract a single trait dimension from a community with multidimensional traits.

Creates a new community where each species retains only the specified trait dimension.
This is useful for visualizing or analyzing evolution along a single trait axis.

# Arguments
- `comm::Community`: The community to extract from
- `dimIndex::Integer`: The trait dimension to retain (1-indexed)

# Returns
New community with only the specified trait dimension for each species

# Example
```julia
# Community with 2D traits
comm = Community([Species(1.0, [0.3, 0.5]), Species(2.0, [0.1, 0.8])],
                 PopulationSize{Float64}[], 0.0)

# Extract only the first trait dimension
comm1D = selectTraitDim(comm, 1)  # Traits become [0.3] and [0.1]

# Extract only the second trait dimension
comm2D = selectTraitDim(comm, 2)  # Traits become [0.5] and [0.8]
```
"""
function selectTraitDim(
        comm::Community{T, AuxClasses},
        dimIndex::Integer
    ) where {T<:Real, AuxClasses}
    numSp = numSpecies(comm)
    numSp > 0 || throw(ArgumentError("Cannot select trait dimension from empty community"))

    # Validate dimension index
    traitDim = length(traits(comm, 1))
    1 <= dimIndex <= traitDim || throw(ArgumentError(
        "Trait dimension $dimIndex out of bounds (trait space has dimension $traitDim)"
    ))

    # Create new species with only the selected trait dimension
    oldSpecies = speciesList(comm)
    newSpecies = Species{T}[]
    for sp in oldSpecies
        traitVec = sp.trait[1].trait  # Get the trait vector
        selectedTrait = Phenotype(traitVec[dimIndex])  # Extract single dimension
        newSp = Species{T}(sp.popsize, [selectedTrait])
        push!(newSpecies, newSp)
    end

    Community(newSpecies, auxs(comm), comm.time)
end


"""
    selectTraitDim(comm::Community, dimIndices::AbstractVector{<:Integer})

Extract multiple trait dimensions from a community with multidimensional traits.

Creates a new community where each species retains only the specified trait dimensions
in the given order. This is useful for projecting high-dimensional trait spaces onto
lower-dimensional subspaces.

# Arguments
- `comm::Community`: The community to extract from
- `dimIndices::AbstractVector{<:Integer}`: Vector of trait dimensions to retain (1-indexed)

# Returns
New community with only the specified trait dimensions for each species

# Example
```julia
# Community with 4D traits
comm = Community([Species(1.0, [0.1, 0.2, 0.3, 0.4]),
                  Species(2.0, [0.5, 0.6, 0.7, 0.8])],
                 PopulationSize{Float64}[], 0.0)

# Extract dimensions 1 and 3
comm2D = selectTraitDim(comm, [1, 3])  # Traits become [0.1, 0.3] and [0.5, 0.7]

# Reorder dimensions
commReordered = selectTraitDim(comm, [4, 1])  # Traits become [0.4, 0.1] and [0.8, 0.5]
```
"""
function selectTraitDim(
        comm::Community{T, AuxClasses},
        dimIndices::AbstractVector{<:Integer}
    ) where {T<:Real, AuxClasses}
    numSp = numSpecies(comm)
    numSp > 0 || throw(ArgumentError("Cannot select trait dimensions from empty community"))

    length(dimIndices) > 0 || throw(ArgumentError("dimIndices must contain at least one dimension"))

    # Validate all dimension indices
    traitDim = length(traits(comm, 1))
    for idx in dimIndices
        1 <= idx <= traitDim || throw(ArgumentError(
            "Trait dimension $idx out of bounds (trait space has dimension $traitDim)"
        ))
    end

    # Check for duplicate indices
    if length(unique(dimIndices)) != length(dimIndices)
        throw(ArgumentError("dimIndices contains duplicate dimensions"))
    end

    # Create new species with only the selected trait dimensions
    oldSpecies = speciesList(comm)
    newSpecies = Species{T}[]
    for sp in oldSpecies
        traitVec = sp.trait[1].trait  # Get the trait vector
        selectedTraits = [traitVec[i] for i in dimIndices]  # Extract selected dimensions
        newTrait = Phenotype(selectedTraits)
        newSp = Species{T}(sp.popsize, [newTrait])
        push!(newSpecies, newSp)
    end

    Community(newSpecies, auxs(comm), comm.time)
end


# Remove species from a community

"""
    removeSpecies(comm::Community, index::Integer)
    removeSpecies(comm::Community, indices)

Remove one or more species from a community.

# Arguments
- `comm::Community`: The community
- `index::Integer`: Single species index to remove
- `indices`: Collection of species indices to remove

# Returns
New community with specified species removed

# Example
```julia
comm = Community([1.0, 2.0, 3.0], [0.1, 0.2, 0.3], Float64[])
comm2 = removeSpecies(comm, 2)
comm3 = removeSpecies(comm, [1, 3])
```
"""
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

"""
    removeExtinct(comm::Community, extThreshold::Real)

Remove all species with total population size below the extinction threshold.

# Arguments
- `comm::Community`: The community
- `extThreshold::Real`: Extinction threshold (non-negative)

# Returns
New community with extinct species removed

# Example
```julia
comm = Community([0.001, 2.0, 0.005, 3.0], [0.1, 0.2, 0.3, 0.4], Float64[])
comm2 = removeExtinct(comm, 0.01)  # Keeps only species 2 and 4
```
"""
function removeExtinct(comm::Community, extThreshold::Real)
    indices_to_remove = speciesBelowThreshold(comm, extThreshold)

    if isempty(indices_to_remove)
        return comm  # No species to remove
    else
        return removeSpecies(comm, indices_to_remove)
    end
end



# Order species by trait component

"""
    orderByTrait(comm::Community)
    orderByTrait(comm::Community, n::Integer)

Sort species in a community by trait value.

# Arguments
- `comm::Community`: The community to sort
- `n::Integer`: Trait dimension to sort by (default: 1)

# Returns
New community with species sorted in ascending order by the specified trait

# Example
```julia
comm = Community([1.0, 2.0, 3.0], [0.3, 0.1, 0.2], Float64[])
sorted_comm = orderByTrait(comm)  # Sorts by trait: species order becomes [2, 3, 1]

# For multidimensional traits
comm_multi = Community([1.0, 2.0], [[0.5, 0.9], [0.3, 0.7]])
sorted_comm2 = orderByTrait(comm_multi, 2)  # Sort by second trait dimension
```
"""
function orderByTrait(comm::Community, n::Integer)
    numSp = numSpecies(comm)

    if numSp == 0
        return comm  # Empty community, nothing to reorder
    end

    # Validate trait dimension
    traitDim = length(traits(comm, 1))
    1 <= n <= traitDim || throw(ArgumentError(
        "Trait component $n out of bounds (trait space has dimension $traitDim)"
    ))

    # Extract nth trait component for each species
    trait_values = Float64[]
    for i in 1:numSp
        trait_vec = traits(comm, i)
        push!(trait_values, trait_vec[n])
    end

    # Get sorted indices
    sorted_indices = sortperm(trait_values)

    # Reorder species
    oldSpecies = speciesList(comm)
    newSpecies = [oldSpecies[i] for i in sorted_indices]

    Community(newSpecies, auxs(comm), comm.time)
end


function orderByTrait(comm::Community)
    orderByTrait(comm, 1)
end


"""
    historyToTable(history::EvoHistory)

Convert an EvoHistory to a tabular format suitable for saving as CSV/TSV.

Returns an `OrderedDict` where keys are column names and values are vectors of the data.
Each row represents one species at one mutation event, with multiple rows per mutation event.

# Column structure
- `mutNo`: Mutation event number (0 for initial community)
- `time`: Integration time for the given step
- `species`: Species index within the community
- `popsize_i`: Population size in stage class i
- `trait_j`: Trait value in dimension j
- `aux_k`: Auxiliary variable k (same for all species in a mutation event)

The maximum number of stage classes, trait dimensions, and auxiliary variables
across all communities determines the number of columns.

# Example
```julia
history = evolve!(community, config, 100)
table = historyToTable(history)

# Save to CSV (requires DataFrames and CSV packages)
using DataFrames, CSV
df = DataFrame(table)
CSV.write("evolution.csv", df)
```
"""
function historyToTable(history::EvoHistory{T, AuxClasses}) where {T<:Real, AuxClasses}
    nSteps = length(history)
    nSteps > 0 || throw(ArgumentError("Cannot convert empty history to table"))

    # Determine maximum dimensions and count total rows
    dims = _getMaxDimensions(historyList(history))
    nRows = sum(numSpecies(comm) for comm in historyList(history))

    # Initialize table as an ordered dictionary of vectors to preserve column order
    table = OrderedDict{String, Vector{Union{T, Missing, Int}}}()

    # Initialize columns
    table["mutNo"] = Vector{Int}(undef, nRows)
    table["time"] = Vector{T}(undef, nRows)
    table["species"] = Vector{Int}(undef, nRows)

    # Add popsize columns
    for i in 1:dims.maxStageClasses
        table["popsize_$i"] = Vector{Union{T, Missing}}(undef, nRows)
    end

    # Add trait columns
    for j in 1:dims.maxTraitDims
        table["trait_$j"] = Vector{Union{T, Missing}}(undef, nRows)
    end

    # Add auxiliary variable columns
    for k in 1:dims.maxAuxVars
        for c in 1:dims.maxAuxComponents
            table["aux_$(k)_$(c)"] = Vector{Union{T, Missing}}(undef, nRows)
        end
    end

    # Fill in the data row by row
    rowIdx = 1
    for (mutEvent, comm) in enumerate(historyList(history))
        nSpecies = numSpecies(comm)
        for sp in 1:nSpecies
            # Basic info
            table["mutNo"][rowIdx] = mutEvent - 1  # 0-indexed
            table["time"][rowIdx] = comm.time
            table["species"][rowIdx] = sp

            # Popsize data
            popsize_vals = popsizes(comm, sp)
            for i in 1:dims.maxStageClasses
                table["popsize_$i"][rowIdx] = i <= length(popsize_vals) ? popsize_vals[i] : missing
            end

            # Trait data
            trait_vals = traits(comm, sp)
            for j in 1:dims.maxTraitDims
                table["trait_$j"][rowIdx] = j <= length(trait_vals) ? trait_vals[j] : missing
            end

            # Auxiliary variable data (same for all species in this mutation event)
            aux_vars = auxs(comm)
            for k in 1:dims.maxAuxVars
                for c in 1:dims.maxAuxComponents
                    if k <= length(aux_vars)
                        aux_components = aux_vars[k].popsize
                        table["aux_$(k)_$(c)"][rowIdx] = c <= length(aux_components) ?
                                                         aux_components[c] : missing
                    else
                        table["aux_$(k)_$(c)"][rowIdx] = missing
                    end
                end
            end

            rowIdx += 1
        end
    end

    return table
end


"""
Helper function to determine maximum dimensions across all communities in history.
"""
function _getMaxDimensions(communities::Vector{<:Community})
    maxSpecies = maximum(numSpecies(comm) for comm in communities)
    maxStageClasses = 0
    maxTraitDims = 0
    maxAuxVars = 0
    maxAuxComponents = 0

    for comm in communities
        if numSpecies(comm) > 0
            maxStageClasses = max(maxStageClasses, length(popsizes(comm, 1)))
            maxTraitDims = max(maxTraitDims, length(traits(comm, 1)))
        end
        if length(auxs(comm)) > 0
            maxAuxVars = length(auxs(comm))
            maxAuxComponents = max(maxAuxComponents, length(auxs(comm)[1].popsize))
        end
    end

    return (maxSpecies = maxSpecies, maxStageClasses = maxStageClasses,
            maxTraitDims = maxTraitDims, maxAuxVars = maxAuxVars,
            maxAuxComponents = maxAuxComponents)
end


# Utility function for axis tick spacing

"""
    niceTickInterval(maxValue::Real, targetTicks::Int=10)

Calculate a human-readable tick interval for axis labels.

Returns an interval rounded to (1, 2, or 5) × 10^n that produces approximately
`targetTicks` tick marks between 0 and `maxValue`.

# Arguments
- `maxValue::Real`: The maximum value on the axis
- `targetTicks::Int=10`: Approximate number of tick marks desired

# Returns
- `Int`: A nice round interval for tick marks

# Example
```julia
niceTickInterval(1500)  # Returns 200 (gives ticks: 0, 200, 400, ..., 1400, 1600)
niceTickInterval(95000) # Returns 20000 (gives ticks: 0, 20000, 40000, ..., 80000, 100000)
```
"""
function niceTickInterval(maxValue::Real, targetTicks::Int=6)
    targetTicks > 0 || throw(ArgumentError("targetTicks must be positive"))
    maxValue >= 0 || throw(ArgumentError("maxValue must be non-negative"))

    approximate_interval = maxValue / targetTicks
    magnitude = 10.0 ^ floor(log10(max(approximate_interval, 1)))

    # Choose the nicest multiplier (1, 2, or 5) for this magnitude
    ratio = approximate_interval / magnitude
    tick_step = if ratio <= 1.5
        magnitude
    elseif ratio <= 3.5
        2 * magnitude
    elseif ratio <= 7.5
        5 * magnitude
    else
        10 * magnitude
    end

    return Int(round(tick_step))
end
