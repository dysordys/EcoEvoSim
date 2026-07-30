# Mutant generation for eco-evolutionary simulations

using Distributions
using LinearAlgebra


"""
    generateMutant(community::Community{T, AuxClasses}, invaderPopsize::T,
                   covMat::AbstractMatrix{T}, selectionFunc) where {T, AuxClasses}

General mutant generation function with customizable parent selection.
- Selects a parent species using the provided selection function
- Creates a mutant with trait = parent trait + random variate from Normal(0, covMat)
- Sets mutant population size to invaderPopsize
- Returns a new community with the mutant appended to the species list
"""
function generateMutant(
        community::Community{T, AuxClasses},
        invaderPopsize::T,
        covMat::AbstractMatrix{T},
        selectionFunc
    ) where {T<:Real, AuxClasses}
    # Validate covariance matrix
    traitDim = traitSpaceDim(community)
    size(covMat) == (traitDim, traitDim) || throw(ArgumentError(
        "covMat must be $traitDim × $traitDim to match trait space dimension"
    ))

    # Select parent species using provided selection function
    parentIdx = selectionFunc(community)
    parentTrait = traits(community, parentIdx)

    # Get number of stage classes from parent species
    nStages = numStages(community)

    # Generate mutant trait: parent + multivariate normal
    mutantTrait = parentTrait .+ rand(MvNormal(covMat))

    # Create mutant species with invader population size distributed across stage classes
    mutantPopsizeVec = fill(invaderPopsize / nStages, nStages)
    mutantSpecies = Species(mutantPopsizeVec, mutantTrait)

    # Create new community with mutant appended
    newSpeciesList = vcat(speciesList(community), mutantSpecies)
    Community(newSpeciesList, community.aux, commTime(community))
end


# User-facing factory functions (return Community -> Community closures)

# Shared helper — reduces the factory bodies to one place
function _makeMutantFactory(workerFn;
    invaderPopsize::Real,
    covMat::Union{AbstractMatrix{<:Real}, Nothing} = nothing,
    variance::Union{Real, Nothing} = nothing,
    parentSelection = weightedRandomSpecies
)
    if covMat !== nothing && variance !== nothing
        throw(ArgumentError("Specify either covMat or variance, not both"))
    elseif covMat === nothing && variance === nothing
        throw(ArgumentError("Specify either covMat or variance"))
    end
    if variance !== nothing
        variance > 0 || throw(ArgumentError("variance must be positive"))
    end
    isempty(methods(parentSelection)) && throw(ArgumentError(
        "parentSelection must be a function mapping a Community to a species index " *
        "(e.g. weightedRandomSpecies or randomSpecies); got a $(typeof(parentSelection))"
    ))
    return function(community::Community{T, AuxClasses}) where {T<:Real, AuxClasses}
        cov = covMat !== nothing ?
            Matrix{T}(covMat) :
            Matrix{T}(T(variance) * I(traitSpaceDim(community)))
        workerFn(community, T(invaderPopsize), cov, parentSelection)
    end
end


"""
    generateMutant(; invaderPopsize, covMat=nothing, variance=nothing,
                   parentSelection=weightedRandomSpecies)

Factory: return a mutation generator. Returns a function `Community -> Community`
suitable for the `mutationGenerator` field of `EcoEvoConfig`.

Specify either `covMat` (full covariance matrix) or `variance` (diagonal covariance).
Argument validation (excluding covariance matrix dimension) occurs at factory creation time.

# Arguments
- `invaderPopsize`: Population size the mutant is introduced at
- `covMat`: Full covariance matrix of the mutational step distribution
- `variance`: Scalar variance, giving a diagonal covariance matrix (alternative to `covMat`)
- `parentSelection`: Any function mapping a `Community` to the index of the parent
  species. The default, [`weightedRandomSpecies`](@ref), draws the parent with
  probability proportional to its population size, so that every *individual* is
  equally likely to give rise to the mutant. Pass [`randomSpecies`](@ref) to give
  every *species* an equal chance instead, regardless of abundance, or supply a
  selector of your own.

# Example
```julia
comm = Community([1.0, 10.0], [0.0, 0.3])

# Parent drawn in proportion to density (the default):
gen = generateMutant(invaderPopsize=0.001, variance=0.01^2)
mutantComm = gen(comm)

# Parent drawn uniformly over species instead:
genUniform = generateMutant(invaderPopsize=0.001, variance=0.01^2,
                            parentSelection=randomSpecies)
mutantCommUniform = genUniform(comm)
```
"""
generateMutant(; kw...) = _makeMutantFactory(generateMutant; kw...)


# Spatial mutant generation - for populations with explicit spatial structure

"""
    generateMutantSpatial(community::Community{T, AuxClasses}, invaderPopsize::T,
                          covMat::AbstractMatrix{T}, selectionFunc) where {T, AuxClasses}

Generate a mutant for spatially-structured populations by placing it in a single patch.

For populations where each stage class represents a spatial patch, this function:
- Selects a parent species using the provided selection function
- Creates a mutant with trait = parent trait + random variate from Normal(0, covMat)
- Places the entire invader population in a single patch, chosen probabilistically based
  on the relative population size in each patch across all species
- Returns a new community with the mutant appended to the species list

# Arguments
- `community::Community{T}`: The current community
- `invaderPopsize::T`: Total population size of the invader
- `covMat::AbstractMatrix{T}`: Covariance matrix for trait mutations
- `selectionFunc`: Function selecting a parent species (e.g., `randomSpecies`, `weightedRandomSpecies`)
"""
function generateMutantSpatial(
        community::Community{T, AuxClasses},
        invaderPopsize::T,
        covMat::AbstractMatrix{T},
        selectionFunc
    ) where {T<:Real, AuxClasses}
    # Validate covariance matrix
    traitDim = traitSpaceDim(community)
    size(covMat) == (traitDim, traitDim) || throw(ArgumentError(
        "covMat must be $traitDim × $traitDim to match trait space dimension"
    ))

    # Select parent species using provided selection function
    parentIdx = selectionFunc(community)
    parentTrait = traits(community, parentIdx)

    # Get number of patches (stage classes)
    nPatches = numStages(community)

    # Calculate total population in each patch across all species
    # This is used to determine probability of mutant appearing in each patch
    patchPopulations = zeros(T, nPatches)
    for sp in speciesList(community)
        patchPopulations .+= popsize(sp)
    end

    # Check that at least one patch has positive population
    totalPopulation = sum(patchPopulations)
    totalPopulation > 0.0 || throw(ArgumentError(
        "All patches have zero population size; cannot place mutant"
    ))

    # Select the patch where the mutant will appear using cumulative distribution
    r = rand() * totalPopulation
    cumsum = zero(T)
    selectedPatch = 1  # Default fallback
    for p in 1:nPatches
        cumsum += patchPopulations[p]
        if r <= cumsum
            selectedPatch = p
            break
        end
    end

    # Generate mutant trait: parent + multivariate normal
    mutantTrait = parentTrait .+ rand(MvNormal(covMat))

    # Create mutant species with all population in the selected patch
    mutantPopsizeVec = zeros(T, nPatches)
    mutantPopsizeVec[selectedPatch] = invaderPopsize
    mutantSpecies = Species(mutantPopsizeVec, mutantTrait)

    # Create new community with mutant appended
    newSpeciesList = vcat(speciesList(community), mutantSpecies)
    Community(newSpeciesList, community.aux, commTime(community))
end


"""
    generateMutantSpatial(; invaderPopsize, covMat=nothing, variance=nothing,
                          parentSelection=weightedRandomSpecies)

Factory: return a mutation generator for spatially-structured populations. Returns a
function `Community -> Community` suitable for the `mutationGenerator` field of
`EcoEvoConfig`.

The mutant appears in a single patch, chosen probabilistically based on patch population sizes.
Specify either `covMat` (full covariance matrix) or `variance` (diagonal covariance).
Argument validation (excluding covariance matrix dimension) occurs at factory creation time.

# Arguments
- `invaderPopsize`: Total population size the mutant is introduced at, all of it
  placed in the selected patch
- `covMat`: Full covariance matrix of the mutational step distribution
- `variance`: Scalar variance, giving a diagonal covariance matrix (alternative to `covMat`)
- `parentSelection`: Any function mapping a `Community` to the index of the parent
  species. The default, [`weightedRandomSpecies`](@ref), draws the parent with
  probability proportional to its total population size across patches, so that every
  *individual* is equally likely to give rise to the mutant — matching the way the
  patch itself is chosen. Pass [`randomSpecies`](@ref) to give every *species* an
  equal chance instead, or supply a selector of your own.

# Example
```julia
# Create a spatially-structured community (2 species with densities in 2 patches)
comm = Community([1.0 1.0; 10.0 10.0], [0.0, 0.3])

# Create generator and apply it
gen = generateMutantSpatial(invaderPopsize=0.001, variance=0.01^2)
mutant_comm = gen(comm)

# Give every species an equal chance of being the parent instead
genUniform = generateMutantSpatial(invaderPopsize=0.001, variance=0.01^2,
                                   parentSelection=randomSpecies)
mutant_comm_uniform = genUniform(comm)
```
"""
generateMutantSpatial(; kw...) = _makeMutantFactory(generateMutantSpatial; kw...)


"""
    noMutation(community::Community) -> Community

Dummy mutation generator that returns the community unchanged, without introducing
any new mutant species. Useful for running pure ecological dynamics (no evolution)
within the eco-evolutionary framework.
"""
noMutation(community::Community) = community
