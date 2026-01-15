struct PopulationSize{T<:Real}
    popsize :: Vector{T}
    function PopulationSize{T}(popsize::Vector{T}) where {T<:Real}
        length(popsize) > 0 || throw(ArgumentError("popsize must have positive length"))
        new{T}(popsize)
    end
end


struct Phenotype{T<:Real}
    trait :: Vector{T}
    function Phenotype{T}(trait::Vector{T}) where {T<:Real}
        length(trait) > 0 || throw(ArgumentError("trait must have positive length"))
        new{T}(trait)
    end
end


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


struct EvoHistory{T<:Real, AuxClasses}
    history :: Vector{Community{T, AuxClasses}}
    function EvoHistory{T, AuxClasses}(
                history::Vector{Community{T, AuxClasses}}
            ) where {T<:Real, AuxClasses}
        AuxClasses >= 0 || throw(ArgumentError("AuxClasses must be non-negative"))
        new{T, AuxClasses}(history)
    end
end
