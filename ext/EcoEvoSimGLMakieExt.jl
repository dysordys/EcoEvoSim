module EcoEvoSimGLMakieExt

using EcoEvoSim
using GLMakie


"""
    plotEvoTwoTraitInteractive(history::EvoHistory; colormap, markersize, xlabel, ylabel, zlabel)

Create an interactive 3D visualization of evolutionary dynamics for 2D trait spaces.

Uses GLMakie to create an interactive plot that can be rotated and zoomed, showing
trait evolution over time with population densities indicated by color.

# Arguments
- `history::EvoHistory`: The evolutionary history to visualize

# Keyword Arguments
- `colormap = [:lightgray, :navy]`: Color scheme for population sizes
- `markersize::Real = 8.0`: Size of markers
- `xlabel::String = "Trait 1"`: Label for first trait axis
- `ylabel::String = "Trait 2"`: Label for second trait axis
- `zlabel::String = "Mutation event"`: Label for time (mutation step) axis

# Returns
A GLMakie Figure object (interactive 3D plot)

# Example
```julia
using EcoEvoSim
using GLMakie

# Run evolutionary simulation with 2D traits
history = evolve!(community, config, 100)

# Create interactive visualization
fig = plotEvoTwoTraitInteractive(history)
display(fig)
```

# Errors
Throws an `ArgumentError` if the trait space is not 2-dimensional or history is empty.
"""
function EcoEvoSim.plotEvoTwoTraitInteractive(
        history::EvoHistory{T, AuxClasses};
        colormap = [:lightgray, :navy],
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
    normalized_sizes = pop_max > pop_min ?
        (allPopsizes .- pop_min) ./ (pop_max - pop_min) :
        zeros(T, length(allPopsizes))

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


end # module EcoEvoSimGLMakieExt
