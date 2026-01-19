"""
    PopulationSize{T<:Real}

Wrapper type for population sizes, represented as a vector to support
structured populations (e.g., stage classes or spatial locations). For
unstructured populations, use a single-element vector.

# Fields
- `popsize::Vector{T}`: Vector of population sizes across structure (e.g., stages/locations)

# Example
```julia
# Single stage class
pop = PopulationSize(10.5)

# Multiple stage classes
pop = PopulationSize([5.0, 15.0, 8.0])
```
"""
struct PopulationSize{T<:Real}
    popsize :: Vector{T}
    function PopulationSize{T}(popsize::Vector{T}) where {T<:Real}
        length(popsize) > 0 || throw(ArgumentError("popsize must have positive length"))
        new{T}(popsize)
    end
end


"""
    Phenotype{T<:Real}

Wrapper type for phenotypic trait values, represented as a vector to support
multidimensional trait spaces.

# Fields
- `trait::Vector{T}`: Vector of trait values across trait dimensions

# Example
```julia
# Single trait
trait = Phenotype(0.5)

# Multiple traits
trait = Phenotype([0.5, -0.2, 1.0])
```
"""
struct Phenotype{T<:Real}
    trait :: Vector{T}
    function Phenotype{T}(trait::Vector{T}) where {T<:Real}
        length(trait) > 0 || throw(ArgumentError("trait must have positive length"))
        new{T}(trait)
    end
end


"""
    Species{T<:Real}

Represents a single species with population structure and phenotypic traits.

# Fields
- `popsize::Vector{PopulationSize{T}}`: Population sizes (typically one element)
- `trait::Vector{Phenotype{T}}`: Phenotypic traits (typically one element)

# Example
```julia
# Simple species with single value
sp = Species(10.0, 0.5)

# Species with stage classes and multidimensional traits
sp = Species([5.0, 15.0], [0.5, -0.2])
```
"""
struct Species{T<:Real}
    popsize :: Vector{PopulationSize{T}}
    trait :: Vector{Phenotype{T}}
    function Species{T}(popsize::Vector{PopulationSize{T}},
                        trait::Vector{Phenotype{T}}) where {T<:Real}
        length(popsize) > 0 || throw(ArgumentError("popsize must have positive length"))
        length(trait) > 0 || throw(ArgumentError("trait must have positive length"))
        new{T}(popsize, trait)
    end
end


"""
    Community{T<:Real, AuxClasses}

Represents a community of species with optional auxiliary variables (e.g., resources).

# Type Parameters
- `T<:Real`: Numeric type for population sizes and traits
- `AuxClasses`: Number of auxiliary variables (e.g., resource pools)

# Fields
- `species::Vector{Species{T}}`: Collection of species in the community
- `aux::Vector{PopulationSize{T}}`: Auxiliary variables (e.g., resources)
- `time::T`: Current time point in the simulation

# Example
```julia
# Community with 3 species and no auxiliary variables
comm = Community([1.0, 2.0, 3.0], [0.1, 0.2, 0.3], Float64[])

# Community with auxiliary variables (resources)
comm = Community(species_vec, aux_vec, 0.0)
```
"""
struct Community{T<:Real, AuxClasses}
    species :: Vector{Species{T}}
    aux :: Vector{PopulationSize{T}}  # Auxiliary variables (e.g., resources)
    time :: T
    function Community{T, AuxClasses}(species::Vector{Species{T}},
                                      aux::Vector{PopulationSize{T}},
                                      time::T = zero(T)) where {T <: Real, AuxClasses}
        AuxClasses >= 0 || throw(ArgumentError("AuxClasses must be non-negative"))
        length(aux) == AuxClasses || throw(ArgumentError(
            "aux must have length $AuxClasses"
        ))
        new{T, AuxClasses}(species, aux, time)
    end
end


"""
    EvoHistory{T<:Real, AuxClasses}

Container for storing the evolutionary history of a community as a sequence of
community snapshots after each mutation event.

# Type Parameters
- `T<:Real`: Numeric type for population sizes and traits
- `AuxClasses`: Number of auxiliary state variables

# Fields
- `history::Vector{Community{T, AuxClasses}}`: Sequence of community states

# Example
```julia
# Create history from initial community
history = EvoHistory(initial_community)

# Run evolution
evolve!(history, config, 100)
```
"""
struct EvoHistory{T<:Real, AuxClasses}
    history :: Vector{Community{T, AuxClasses}}
    function EvoHistory{T, AuxClasses}(
                history::Vector{Community{T, AuxClasses}}
            ) where {T<:Real, AuxClasses}
        AuxClasses >= 0 || throw(ArgumentError("AuxClasses must be non-negative"))
        new{T, AuxClasses}(history)
    end
end
