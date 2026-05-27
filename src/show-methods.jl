function Base.show(io::IO, ps::PopulationSize)
    print(io, ps.popsize)
end


function Base.show(io::IO, ph::Phenotype)
    print(io, ph.trait)
end


function Base.show(io::IO, sp::Species)
    density_vals = Vector(sp.popsize.popsize)
    trait_vals = Vector(sp.trait.trait)
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


function Base.show(io::IO, h::EvoHistory)
    n = length(h)
    if n == 0
        print(io, "EvoHistory (empty)")
        return
    end
    last_comm = h.history[end]
    print(io, "EvoHistory ($n snapshot$(n == 1 ? "" : "s"), " *
              "last: t = $(last_comm.time), $(numSpecies(last_comm)) species)")
end


function Base.show(io::IO, ts::Vector{<:Community})
    n = length(ts)
    if n == 0
        print(io, "EcoDynTimeSeries (empty)")
        return
    end
    t0     = ts[1].time
    tn     = ts[end].time
    ns     = numSpecies(ts[end])
    plural = n == 1 ? "" : "s"
    print(io, "EcoDynTimeSeries ($n snapshot$plural, t = $t0 to $tn, $ns species)")
end

function Base.show(io::IO, ::MIME"text/plain", ts::Vector{<:Community})
    show(io, ts)
end


function Base.show(io::IO, params::IntegrationParams)
    println(io, "IntegrationParams:")
    println(io, "  maxTime:   $(params.maxTime)")
    println(io, "  algorithm: $(nameof(typeof(params.algorithm)))")
    print(io,   "  options:   $(params.solver_options)")
end


function Base.show(io::IO, config::EcoEvoConfig)
    println(io, "EcoEvoConfig:")
    println(io, "  ecoDyn:       $(nameof(typeof(config.ecoDyn)))")
    println(io, "  mutGenerator: $(nameof(typeof(config.mutationGenerator)))")
    println(io, "  extThreshold: $(config.extThreshold)")
    print(io,   "  integration:  ")
    show(io, config.integrationParams)
end
