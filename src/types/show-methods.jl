function Base.show(io::IO, ps::PopulationSize)
    print(io, ps.popsize)
end


function Base.show(io::IO, ph::Phenotype)
    print(io, ph.trait)
end


function Base.show(io::IO, sp::Species)
    if length(sp.popsize) == 1
        density_vals = Vector(sp.popsize[1].popsize)
    else
        density_vals = [Vector(ps.popsize) for ps in sp.popsize]
    end
    if length(sp.trait) == 1
        trait_vals = Vector(sp.trait[1].trait)
    else
        trait_vals = [Vector(ph.trait) for ph in sp.trait]
    end
    println(io, "- popsize:   $density_vals")
    println(io, "- phenotype: $trait_vals")
end


function Base.show(
        io::IO, comm::Community{T, AuxClasses}
    ) where {T, AuxClasses}
    println(io, "Community (t = $(comm.time))")
    for (i, sp) in enumerate(comm.species)
        println(io, "  Species $i:")
        buf = IOBuffer()
        show(buf, sp)
        str = String(take!(buf))
        indented = join("    " * line for line in split(str, '\n') if !isempty(line), "\n")
        println(io, indented)
    end
    for (j, aux) in enumerate(comm.aux)
        println(io, "  Auxiliary $j: $(aux)")
    end
end
