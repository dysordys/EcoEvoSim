struct PopulationSize{T<:Real, StageClasses}
    popsize :: SVector{StageClasses, T}
    function PopulationSize{T, StageClasses}(
            popsize::SVector{StageClasses, T}
        ) where {T<:Real, StageClasses}
        StageClasses > 0 || throw(ArgumentError("StageClasses must be positive"))
        new{T, StageClasses}(popsize)
    end
end


struct Phenotype{T<:Real, TraitDim}
    trait :: SVector{TraitDim, T}
    function Phenotype{T, TraitDim}(
            trait::SVector{TraitDim, T}
        ) where {T<:Real, TraitDim}
        TraitDim > 0 || throw(ArgumentError("TraitDim must be positive"))
        new{T, TraitDim}(trait)
    end
end


struct Species{T<:Real, StageClasses, TraitDim}
    popsize :: Vector{PopulationSize{T, StageClasses}}
    trait :: Vector{Phenotype{T, TraitDim}}
    function Species{T, StageClasses, TraitDim}(
            popsize::Vector{PopulationSize{T, StageClasses}},
            trait::Vector{Phenotype{T, TraitDim}}
        ) where {T<:Real, StageClasses, TraitDim}
        StageClasses > 0 || throw(ArgumentError("StageClasses must be positive"))
        TraitDim > 0 || throw(ArgumentError("TraitDim must be positive"))
        new{T, StageClasses, TraitDim}(popsize, trait)
    end
end


struct Community{T<:Real, StageClasses, TraitDim, AuxClasses}
    species :: Vector{Species{T, StageClasses, TraitDim}}
    aux :: Vector{PopulationSize{T, 1}}  # Auxiliary variables (e.g., resources)
    time :: T
    function Community{T, StageClasses, TraitDim, AuxClasses}(
            species::Vector{Species{T, StageClasses, TraitDim}},
            aux::Vector{PopulationSize{T, 1}},
            time::T = zero(T)
        ) where {T <: Real, StageClasses, TraitDim, AuxClasses}
        StageClasses > 0 || throw(ArgumentError("StageClasses must be positive"))
        TraitDim > 0 || throw(ArgumentError("TraitDim must be positive"))
        AuxClasses >= 0 || throw(ArgumentError("AuxClasses must be non-negative"))
        new{T, StageClasses, TraitDim, AuxClasses}(species, aux, time)
    end
end


struct EvoHistory{T<:Real, StageClasses, TraitDim, AuxClasses}
    history :: Vector{Community{T, StageClasses, TraitDim, AuxClasses}}
    function EvoHistory{T, StageClasses, TraitDim, AuxClasses}(
            history::Vector{Community{T, StageClasses, TraitDim, AuxClasses}}
        ) where {T<:Real, StageClasses, TraitDim, AuxClasses}
        StageClasses > 0 || throw(ArgumentError("StageClasses must be positive"))
        TraitDim > 0 || throw(ArgumentError("TraitDim must be positive"))
        AuxClasses >= 0 || throw(ArgumentError("AuxClasses must be non-negative"))
        new{T, StageClasses, TraitDim, AuxClasses}(history)
    end
end
