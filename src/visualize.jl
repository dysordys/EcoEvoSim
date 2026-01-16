# Visualization functions for eco-evolutionary dynamics

using Plots


"""
    plotEvo(history::EvoHistory; kwargs...)

Visualize evolutionary dynamics over time for systems with 1-dimensional trait space.

Creates a plot showing trait evolution over mutational steps, with population densities
indicated by color only.

# Arguments
- `history::EvoHistory`: The evolutionary history to visualize

# Keyword Arguments
- `colormap::Symbol = :viridis`: Color scheme for population densities
- `markersize::Real = 5.0`: Size of markers
- `alpha::Real = 0.7`: Transparency of markers
- `xlabel::String = "Trait value"`: Label for x-axis
- `ylabel::String = "Mutation step"`: Label for y-axis
- `title::String = "Evolutionary Dynamics"`: Plot title
- `size::Tuple = (800, 600)`: Figure size
- `legend::Bool = true`: Whether to show legend

# Returns
A Plots.jl plot object

# Example
```julia
using EcoEvoSim
using Plots

# Run evolutionary simulation
history = evolve!(community, config, 100)

# Visualize
p = plotEvo(history)
display(p)
```

# Errors
Throws an `ArgumentError` if the trait space is not 1-dimensional.
"""
function plotEvo(
        history::EvoHistory{T, AuxClasses};
        colormap::Symbol = :viridis,
        markersize::Real = 5.0,
        alpha::Real = 0.7,
        xlabel::String = "Trait value",
        ylabel::String = "Mutation step",
        title::String = "Evolutionary Dynamics",
        size::Tuple = (800, 600),
        legend::Bool = true
    ) where {T<:Real, AuxClasses}

    # Check that trait space is 1-dimensional
    if length(history.history) > 0
        traitDim = traitSpaceDim(history.history[1])
        traitDim == 1 || throw(ArgumentError(
            "plotEvo only works for 1-dimensional trait spaces (got dimension $traitDim)"
        ))
    else
        throw(ArgumentError("Cannot plot empty evolutionary history"))
    end

    # Extract data from history
    nSteps = length(history.history)
    allTraits = T[]
    allPopsizes = T[]
    allTimes = Int[]

    for (step, comm) in enumerate(history.history)
        for sp_idx in 1:numSpecies(comm)
            trait_val = traits(comm, sp_idx)[1]  # Get first (and only) trait dimension
            pop_val = sum(popsizes(comm, sp_idx))  # Sum across stage classes if any
            push!(allTraits, trait_val)
            push!(allPopsizes, pop_val)
            push!(allTimes, step)
        end
    end

    # Normalize population sizes for coloring
    if length(allPopsizes) > 0
        max_pop = maximum(allPopsizes)
        min_pop = minimum(allPopsizes)
        pop_range = max_pop - min_pop

        # Avoid division by zero
        if pop_range > 0
            normalized_pops = (allPopsizes .- min_pop) ./ pop_range
        else
            normalized_pops = ones(length(allPopsizes))
        end
    else
        throw(ArgumentError("No species found in evolutionary history"))
    end

    # Create plot with fixed marker size
    p = scatter(
        allTraits,
        allTimes,
        marker_z = allPopsizes,
        color = colormap,
        markersize = markersize,
        markerstrokewidth = 0,
        alpha = alpha,
        xlabel = xlabel,
        ylabel = ylabel,
        title = title,
        size = size,
        legend = legend ? :right : false,
        colorbar = legend,
        colorbar_title = legend ? "Population\nsize" : "",
        label = ""
    )

    # Set y-axis to show mutation steps
    plot!(p, yticks = 1:max(1, div(nSteps, 10)):nSteps)

    return p
end
