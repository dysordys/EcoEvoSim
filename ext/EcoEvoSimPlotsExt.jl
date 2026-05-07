module EcoEvoSimPlotsExt

using EcoEvoSim
using Plots


function EcoEvoSim.plotEvo(
        history::EvoHistory{T, AuxClasses};
        traitDim::Integer = 1,
        colormap = cgrad([:lightgray, :navy]),
        markersize::Real = 3.0,
        alpha::Real = 0.7,
        xlabel::String = "Trait value",
        ylabel::String = "Mutation event",
        title::Union{String, Nothing} = nothing,
        size::Tuple = (800, 600),
        legend::Bool = true
    ) where {T<:Real, AuxClasses}

    # Check for empty history
    if length(history.history) == 0
        throw(ArgumentError("Cannot plot empty evolutionary history"))
    end

    # Validate trait dimension
    actualTraitDim = traitSpaceDim(history.history[1])
    traitDim >= 1 || throw(ArgumentError("traitDim must be positive (got $traitDim)"))
    traitDim <= actualTraitDim || throw(ArgumentError(
        "traitDim = $traitDim exceeds trait space dimension $actualTraitDim"
    ))

    # Extract the specified trait dimension from history if multidimensional
    processedHistory = if actualTraitDim == 1
        history
    else
        newComms = [selectTraitDim(comm, traitDim) for comm in history.history]
        EvoHistory(newComms)
    end

    # Extract data from processed history
    nSteps = length(processedHistory.history)
    allTraits = T[]
    allPopsizes = T[]
    allTimes = Int[]

    for (step, comm) in enumerate(processedHistory.history)
        for sp_idx in 1:numSpecies(comm)
            trait_val = traits(comm, sp_idx)[1]
            pop_val = sum(popsizes(comm, sp_idx))
            push!(allTraits, trait_val)
            push!(allPopsizes, pop_val)
            push!(allTimes, step)
        end
    end

    if length(allPopsizes) == 0
        throw(ArgumentError("No species found in evolutionary history"))
    end

    plot_kwargs = (
        marker_z = allPopsizes,
        color = colormap,
        markersize = markersize,
        markerstrokewidth = 0,
        alpha = alpha,
        xlabel = xlabel,
        ylabel = ylabel,
        size = size,
        legend = legend ? :right : false,
        colorbar = legend,
        colorbar_title = legend ? "Population size" : "",
        label = ""
    )

    if title !== nothing
        plot_kwargs = merge(plot_kwargs, (title = title,))
    end

    p = Plots.scatter(allTraits, allTimes; plot_kwargs...)

    tick_step = niceTickInterval(nSteps)
    Plots.plot!(p, yticks = 0:tick_step:nSteps)

    return p
end


function EcoEvoSim.plotEvoTwoTrait(
        history::EvoHistory{T, AuxClasses};
        colormap = cgrad([:lightgray, :navy]),
        markersize::Real = 3.0,
        alpha::Real = 0.7,
        xlabel::String = "Trait 1",
        ylabel::String = "Trait 2",
        zlabel::String = "Mutation event",
        title::Union{String, Nothing} = nothing,
        size::Tuple = (800, 600),
        legend::Bool = true,
        camera::Tuple = (45, 30)
    ) where {T<:Real, AuxClasses}

    if length(history.history) > 0
        traitDim = traitSpaceDim(history.history[1])
        traitDim == 2 || throw(ArgumentError(
            "plotEvoTwoTrait only works for 2-dimensional trait spaces (got dimension $traitDim)"
        ))
    else
        throw(ArgumentError("Cannot plot empty evolutionary history"))
    end

    nSteps = length(history.history)
    allTrait1 = T[]
    allTrait2 = T[]
    allPopsizes = T[]
    allTimes = Int[]

    for (step, comm) in enumerate(history.history)
        for sp_idx in 1:numSpecies(comm)
            traitsVec = traits(comm, sp_idx)
            push!(allTrait1, traitsVec[1])
            push!(allTrait2, traitsVec[2])
            pop_val = sum(popsizes(comm, sp_idx))
            push!(allPopsizes, pop_val)
            push!(allTimes, step)
        end
    end

    if length(allPopsizes) == 0
        throw(ArgumentError("No species found in evolutionary history"))
    end

    plot_kwargs = (
        marker_z = allPopsizes,
        color = colormap,
        markersize = markersize,
        markerstrokewidth = 0,
        alpha = alpha,
        xlabel = xlabel,
        ylabel = ylabel,
        zlabel = zlabel,
        size = size,
        legend = legend ? :right : false,
        colorbar = legend,
        colorbar_title = legend ? "\n\nPopulation size" : "",
        colorbar_title_location = :right,
        colorbar_titleoffsetx = 0,
        label = "",
        camera = camera
    )

    if title !== nothing
        plot_kwargs = merge(plot_kwargs, (title = title,))
    end

    p = Plots.scatter(allTrait1, allTrait2, allTimes; plot_kwargs...)

    tick_step = niceTickInterval(nSteps)
    Plots.plot!(p, zticks = 0:tick_step:nSteps)

    return p
end


end # module EcoEvoSimPlotsExt
