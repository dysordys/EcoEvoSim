# Outer constructor for Community with default empty history
function Community{T, StageClasses, TraitDim, AuxClasses}(
        species::Vector{Species{T, StageClasses, TraitDim}},
        aux::Vector{PopulationSize{T, 1}}
    ) where {T<:Real, StageClasses, TraitDim, AuxClasses}
    Community{T, StageClasses, TraitDim, AuxClasses}(
        species, aux, Community{T, StageClasses, TraitDim, AuxClasses}[]
    )
end

# Selectors
species(comm::Community) = comm.species
aux(comm::Community) = comm.aux
history(comm::Community) = comm.history
