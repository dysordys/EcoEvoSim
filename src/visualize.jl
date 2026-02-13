# Visualization functions for eco-evolutionary dynamics

using Plots
using GLMakie


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
niceTickInterval(95000) # Returns 20000 (gives ticks: 0, 20000, 40000, ..., 80000, 100000)
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
    plotEvo(history::EvoHistory; traitDim::Integer=1, kwargs...)

Visualize evolutionary dynamics over time along a single trait dimension.

Creates a plot showing trait evolution over mutational steps, with population densities
indicated by color. For multidimensional trait spaces, specify which dimension to plot.

# Arguments
- `history::EvoHistory`: The evolutionary history to visualize

# Keyword Arguments
- `traitDim::Integer = 1`: Trait dimension to plot (for multidimensional trait spaces)
- `colormap = cgrad([:lightgray, :navy])`: Color scheme for population sizes
- `markersize::Real = 3.0`: Size of markers
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

# Visualize (default: first trait dimension)
p = plotEvo(history)
display(p)

# For multidimensional traits, plot second dimension
p2 = plotEvo(history, traitDim=2)
display(p2)
```
"""
function plotEvo(
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
        # Create new history with only the selected trait dimension
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
            trait_val = traits(comm, sp_idx)[1]  # Get the (now single) trait dimension
            pop_val = sum(popsizes(comm, sp_idx))  # Sum across stage classes if any
            push!(allTraits, trait_val)
            push!(allPopsizes, pop_val)
            push!(allTimes, step)
        end
    end

    # Check if we have data
    if length(allPopsizes) == 0
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

    p = Plots.scatter(allTraits, allTimes; plot_kwargs...)

    # Set y-axis to show mutation steps with nice round intervals
    tick_step = niceTickInterval(nSteps)
    Plots.plot!(p, yticks = 0:tick_step:nSteps)

    return p
end


"""
    plotEvoTwoTrait(history::EvoHistory; kwargs...)

Visualize evolutionary dynamics over time for systems with 2-dimensional trait space.

Creates a 3D plot showing trait evolution over mutational steps, with population densities
indicated by color. The three axes represent the two trait dimensions and time.

# Arguments
- `history::EvoHistory`: The evolutionary history to visualize

# Keyword Arguments
- `colormap = cgrad([:lightgray, :navy])`: Color scheme for population sizes
- `markersize::Real = 3.0`: Size of markers
- `alpha::Real = 0.7`: Transparency of markers
- `xlabel::String = "Trait 1"`: Label for first trait axis
- `ylabel::String = "Trait 2"`: Label for second trait axis
- `zlabel::String = "Mutation event"`: Label for time axis
- `title::Union{String, Nothing} = nothing`: Plot title
- `size::Tuple = (800, 600)`: Figure size
- `legend::Bool = true`: Whether to show colorbar legend
- `camera::Tuple = (45, 30)`: Camera viewing angle (azimuth, elevation)

# Returns
A Plots.jl 3D plot object

# Example
```julia
using EcoEvoSim
using Plots

# Run evolutionary simulation with 2D traits
history = evolve!(community, config, 100)

# Visualize
p = plotEvoTwoTrait(history)
display(p)
```

# Errors
Throws an `ArgumentError` if the trait space is not 2-dimensional.
"""
function plotEvoTwoTrait(
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

    # Check that trait space is 2-dimensional
    if length(history.history) > 0
        traitDim = traitSpaceDim(history.history[1])
        traitDim == 2 || throw(ArgumentError(
            "plotEvoTwoTrait only works for 2-dimensional trait spaces (got dimension $traitDim)"
        ))
    else
        throw(ArgumentError("Cannot plot empty evolutionary history"))
    end

    # Extract data from history
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
            pop_val = sum(popsizes(comm, sp_idx))  # Sum across stage classes if any
            push!(allPopsizes, pop_val)
            push!(allTimes, step)
        end
    end

    # Check if we have data
    if length(allPopsizes) == 0
        throw(ArgumentError("No species found in evolutionary history"))
    end

    # Create 3D plot with fixed marker size
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

    # Add title only if not nothing
    if title !== nothing
        plot_kwargs = merge(plot_kwargs, (title = title,))
    end

    p = Plots.scatter(allTrait1, allTrait2, allTimes; plot_kwargs...)

    # Set z-axis to show mutation steps with nice round intervals
    tick_step = niceTickInterval(nSteps)
    Plots.plot!(p, zticks = 0:tick_step:nSteps)

    return p
end


function plotEvoTwoTraitInteractive(
        history::EvoHistory{T, AuxClasses};
        colormap = cgrad([:lightgray, :navy]),
        markersize::Real = 8.0,
        xlabel::String = "Trait 1",
        ylabel::String = "Trait 2",
        zlabel::String = "Mutation event"
    ) where {T<:Real, AuxClasses}

    # Validation
    if length(history.history) > 0
        traitDim = traitSpaceDim(history.history[1])
        traitDim == 2 || throw(ArgumentError(
            "Function only works for 2-dimensional trait spaces (got dimension $traitDim)"
        ))
    else
        throw(ArgumentError("Cannot plot empty evolutionary history"))
    end

    # Extract data
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

    # Normalize population sizes for coloring
    pop_min, pop_max = extrema(allPopsizes)
    normalized_sizes = (allPopsizes .- pop_min) ./ (pop_max - pop_min)

    # Create figure and axis
    fig = Figure(size = (1000, 800))
    ax = Axis3(fig[1, 1];
        xlabel = xlabel,
        ylabel = ylabel,
        zlabel = zlabel
    )

    # Create scatter plot
    Makie.scatter!(ax, allTrait1, allTrait2, allTimes;
        color = normalized_sizes,
        colormap = colormap,
        markersize = markersize,
        alpha = 0.8
    )

    # Add colorbar
    Colorbar(fig[1, 2]; limits = (pop_min, pop_max), colormap = colormap,
        label = "Population size")

    return fig
end
