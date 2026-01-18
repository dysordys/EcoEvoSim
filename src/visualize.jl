# Visualization functions for eco-evolutionary dynamics

using Plots


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
niceTickInterval(95000) # Returns 10000 (gives ticks: 0, 10000, 20000, ..., 90000)
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


"""
    plotEvo(history::EvoHistory; kwargs...)

Visualize evolutionary dynamics over time for systems with 1-dimensional trait space.

Creates a plot showing trait evolution over mutational steps, with population densities
indicated by color only.

# Arguments
- `history::EvoHistory`: The evolutionary history to visualize

# Keyword Arguments
- `colormap = cgrad([RGB(0.92, 0.92, 0.92), :navy])`: Color scheme for population sizes
- `markersize::Real = 5.0`: Size of markers
- `alpha::Real = 0.7`: Transparency of markers
- `xlabel::String = "Trait value"`: Label for x-axis
- `ylabel::String = "Mutation event"`: Label for y-axis
- `title::Union{String, Nothing} = nothing`: Plot title
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
        colormap = cgrad([RGB(0.92, 0.92, 0.92), :navy]),
        markersize::Real = 3.0,
        alpha::Real = 0.7,
        xlabel::String = "Trait value",
        ylabel::String = "Mutation event",
        title::Union{String, Nothing} = nothing,
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

    # Add title only if not nothing
    if title !== nothing
        plot_kwargs = merge(plot_kwargs, (title = title,))
    end

    p = scatter(allTraits, allTimes; plot_kwargs...)

    # Set y-axis to show mutation steps with nice round intervals
    tick_step = niceTickInterval(nSteps)
    plot!(p, yticks = 0:tick_step:nSteps)

    return p
end
