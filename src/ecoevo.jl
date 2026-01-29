# Eco-evolutionary dynamics simulation configuration and functions

using Distributions
using LinearAlgebra
using DifferentialEquations
using Pipe


"""
    IntegrationParams{T<:Real, Alg}

Parameters for ODE integration during ecological dynamics.

# Fields
- `maxTime::T`: Maximum integration time
- `algorithm::Alg`: DifferentialEquations.jl solver algorithm (e.g., `Rodas5()`, `DynamicSS()`)
- `solver_options::NamedTuple`: Keyword arguments for `solve()` (e.g., `abstol`, `reltol`)

# Constructors
```julia
# Default constructor with Rodas5 solver
params = IntegrationParams(maxTime = 50.0)

# Custom solver and tolerances
params = IntegrationParams(maxTime = 100.0, algorithm = Tsit5(),
                          abstol = 1e-10, reltol = 1e-8)

# Steady-state solver
params = IntegrationParams(maxTime = 100.0, algorithm = DynamicSS())
```
"""
struct IntegrationParams{T<:Real, Alg}
    maxTime :: T
    algorithm :: Alg  # ODE solver algorithm
    solver_options :: NamedTuple  # Keyword arguments for solve() (abstol, reltol, etc.)
    function IntegrationParams{T, Alg}(
            maxTime::T, algorithm::Alg, solver_options::NamedTuple
        ) where {T<:Real, Alg}
        maxTime > 0 || throw(ArgumentError("maxTime must be positive"))
        new{T, Alg}(maxTime, algorithm, solver_options)
    end
end

IntegrationParams(maxTime::T, algorithm::Alg, solver_options::NamedTuple) where {T<:Real, Alg} =
    IntegrationParams{T, Alg}(maxTime, algorithm, solver_options)

# Convenience constructor with default Rodas5 and common tolerances
IntegrationParams(maxTime::T; algorithm::Alg = Rodas5(),
                  solver_options::NamedTuple = (abstol=1e-8, reltol=1e-6)) where {T<:Real, Alg} =
    IntegrationParams(maxTime, algorithm, solver_options)

# Keyword argument constructor with explicit abstol/reltol (captures any additional kwargs)
function IntegrationParams(; maxTime::T,
                            algorithm::Alg = Rodas5(),
                            abstol = 1e-8,
                            reltol = 1e-6,
                            kwargs...) where {T<:Real, Alg}
    solver_options = merge((abstol=convert(T, abstol), reltol=convert(T, reltol)), kwargs)
    IntegrationParams(maxTime, algorithm, solver_options)
end


"""
    EcoEvoConfig{T<:Real, EcoDynamics, MutGenerator, Alg}

Configuration for eco-evolutionary simulations.

# Fields
- `ecoDyn::EcoDynamics`: Factory function taking `Community` and returning ODE function
- `mutationGenerator::MutGenerator`: Function generating mutants from community and config
- `integrationParams::IntegrationParams`: ODE integration parameters
- `invaderPopsize::T`: Initial population size for new mutants
- `extThreshold::T`: Population size below which species go extinct

# Example
```julia
config = EcoEvoConfig(
    ecoDyn = lotkaVolterra(growthFn, interactionFn),
    mutationGenerator = (comm, cfg) -> generateMutant(comm, cfg, 0.01),
    integrationParams = IntegrationParams(maxTime = 50.0),
    invaderPopsize = 0.001,
    extThreshold = 0.003
)
```
"""
struct EcoEvoConfig{T<:Real, EcoDynamics, MutGenerator, Alg}
    ecoDyn :: EcoDynamics  # Factory function: Community -> (u, p, t) -> du
    mutationGenerator :: MutGenerator   # Function generating mutant traits
    integrationParams :: IntegrationParams{T, Alg}
    invaderPopsize :: T # Population size at which new mutants are introduced
    extThreshold :: T  # Population size below which species are considered extinct
    function EcoEvoConfig{T, ED, MG, Alg}(
            ecoDyn::ED,
            mutationGenerator::MG,
            integrationParams::IntegrationParams{T, Alg},
            invaderPopsize::T,
            extThreshold::T
        ) where {T<:Real, ED, MG, Alg}
        invaderPopsize > 0 || throw(ArgumentError("invaderPopsize must be positive"))
        extThreshold > 0 || throw(ArgumentError("extThreshold must be positive"))
        new{T, ED, MG, Alg}(ecoDyn, mutationGenerator, integrationParams,
                            invaderPopsize, extThreshold)
    end
end

function EcoEvoConfig(
        ecoDyn::ED,
        mutationGenerator::MG,
        integrationParams::IntegrationParams{T, Alg},
        invaderPopsize::T,
        extThreshold::T
    ) where {T<:Real, ED, MG, Alg}
    EcoEvoConfig{T, ED, MG, Alg}(
        ecoDyn, mutationGenerator, integrationParams, invaderPopsize, extThreshold
    )
end

# Keyword argument constructor
function EcoEvoConfig(; ecoDyn::ED,
                        mutationGenerator::MG,
                        integrationParams::IntegrationParams{T, Alg},
                        invaderPopsize::Real,
                        extThreshold::Real) where {T<:Real, ED, MG, Alg}
    EcoEvoConfig{T, ED, MG, Alg}(ecoDyn, mutationGenerator, integrationParams,
                 convert(T, invaderPopsize), convert(T, extThreshold))
end


"""
    unpackCommunity(community::Community{T, AuxClasses}) where {T, AuxClasses}

Extract a flat vector of state variables from a Community.
Concatenates all population sizes (across species and stage classes)
followed by auxiliary variables.
"""
function unpackCommunity(community::Community{T, AuxClasses}) where {T<:Real, AuxClasses}

    u0 = T[]

    # Add all population sizes (flattened across species and stage classes)
    for sp in speciesList(community)
        append!(u0, sp.popsize[1].popsize)
    end

    # Add auxiliary variables
    for aux in community.aux
        append!(u0, aux.popsize)
    end

    return u0
end


"""
    packCommunity(u::Vector{T}, community::Community{T, AuxClasses}, time::T) where {T, AuxClasses}

Reconstruct a Community from a flat state vector.
Uses the structure of the original community to determine how to partition
the state vector back into species and auxiliary variables.
"""
function packCommunity(
        u::Vector{T},
        community::Community{T, AuxClasses},
        time::T
    ) where {T<:Real, AuxClasses}

    idx = 1
    newSpecies = Species{T}[]

    # Reconstruct each species
    for sp in speciesList(community)
        stageClasses = length(sp.popsize[1].popsize)
        newPopsize = PopulationSize(u[idx:idx+stageClasses-1])
        newSpecies_i = Species(newPopsize, sp.trait[1])
        push!(newSpecies, newSpecies_i)
        idx += stageClasses
    end

    # Reconstruct auxiliary variables
    newAux = PopulationSize{T}[]
    for aux in community.aux
        auxClasses = length(aux.popsize)
        newAuxVar = PopulationSize(u[idx:idx+auxClasses-1])
        push!(newAux, newAuxVar)
        idx += auxClasses
    end

    return Community(newSpecies, newAux, time)
end


"""
    ecoDyn(community::Community{T, AuxClasses}, config::EcoEvoConfig{T}) where {T, AuxClasses}

Integrate ecological dynamics for the specified time period.
Extracts population sizes and auxiliary variables, integrates the ODE system,
and returns a new community with updated values.

The config.ecoDyn field should be a factory function that takes a Community
and returns an ODE function (u, p, t) -> du.
"""
function ecoDyn(
        community::Community{T, AuxClasses},
        config::EcoEvoConfig{T}
    ) where {T<:Real, AuxClasses}

    # Unpack community into state vector
    u0 = unpackCommunity(community)

    # Generate the ODE function for this specific community
    ode_fn = config.ecoDyn(community)

    # Check if we're using a steady-state solver
    alg = config.integrationParams.algorithm
    if alg isa DynamicSS
        # For steady-state integration, use ODEProblem with a callback that terminates
        # when the solution reaches steady state (i.e., when derivatives are near zero)
        tspan = (community.time, community.time + config.integrationParams.maxTime)
        prob = ODEProblem(ode_fn, u0, tspan)

        # Get the inner ODE solver and tolerances from the algorithm
        # If no inner algorithm is specified, use a robust default
        inner_alg = alg.alg === nothing ? Rodas5() : alg.alg

        # Get tolerances for steady-state detection
        abstol = get(config.integrationParams.solver_options, :abstol, 1e-8)
        reltol = get(config.integrationParams.solver_options, :reltol, 1e-6)

        # Use the ODE solver with steady-state termination callback
        sol = solve(prob, inner_alg; config.integrationParams.solver_options...,
                    callback = TerminateSteadyState(abstol, reltol))
        u_final = sol.u[end]
        t_final = sol.t[end]
    elseif alg isa SSRootfind
        # SSRootfind uses initial conditions as a starting guess for rootfinding
        prob = SteadyStateProblem(ode_fn, u0)
        sol = solve(prob, config.integrationParams.algorithm;
                    config.integrationParams.solver_options...)
        u_final = sol.u
        # For rootfinding, time is meaningless - keep original time
        t_final = community.time
    else
        # Use ODEProblem for time-integration solvers
        tspan = (community.time, community.time + config.integrationParams.maxTime)
        prob = ODEProblem(ode_fn, u0, tspan)
        sol = solve(prob, config.integrationParams.algorithm;
                    config.integrationParams.solver_options...)
        u_final = sol.u[end]
        t_final = sol.t[end]
    end

    return packCommunity(u_final, community, t_final)
end


"""
    generateMutant(community::Community{T, AuxClasses}, config::EcoEvoConfig{T}, covMat::AbstractMatrix{T}, selectionFunc) where {T, AuxClasses}

General mutant generation function with customizable parent selection.
- Selects a parent species using the provided selection function
- Creates a mutant with trait = parent trait + random variate from Normal(0, covMat)
- Sets mutant population size to config.invaderPopsize
- Returns a new community with the mutant appended to the species list
"""
function generateMutant(
        community::Community{T, AuxClasses},
        config::EcoEvoConfig{T},
        covMat::AbstractMatrix{T},
        selectionFunc
    ) where {T<:Real, AuxClasses}
    # Validate covariance matrix
    traitDim = traitSpaceDim(community)
    size(covMat) == (traitDim, traitDim) || throw(ArgumentError(
        "covMat must be $traitDim × $traitDim to match trait space dimension"
    ))

    # Select parent species using provided selection function
    parentIdx = selectionFunc(community)
    parentTrait = traits(community, parentIdx)

    # Generate mutant trait: parent + multivariate normal
    mutantTrait = parentTrait .+ rand(MvNormal(covMat))

    # Create mutant species with invader population size
    mutantSpecies = Species(config.invaderPopsize, mutantTrait)

    # Create new community with mutant appended
    newSpeciesList = vcat(speciesList(community), mutantSpecies)
    Community(newSpeciesList, community.aux, community.time)
end


# User-facing methods with uniform random selection

"""
    generateMutant(community::Community{T, AuxClasses}, config::EcoEvoConfig{T}, covMat::AbstractMatrix{T}) where {T, AuxClasses}

Generate a mutant with uniform random parent selection (full covariance matrix).
"""
generateMutant(
    community::Community{T, AuxClasses},
    config::EcoEvoConfig{T},
    covMat::AbstractMatrix{T}
) where {T<:Real, AuxClasses} = generateMutant(community, config, covMat, randomSpecies)


"""
    generateMutant(community::Community{T, AuxClasses}, config::EcoEvoConfig{T}, variance::T) where {T, AuxClasses}

Generate a mutant with uniform random parent selection (diagonal covariance).
"""
function generateMutant(
        community::Community{T, AuxClasses},
        config::EcoEvoConfig{T},
        variance::T
    ) where {T<:Real, AuxClasses}
    variance > 0 || throw(ArgumentError("variance must be positive"))
    covMat = Matrix{T}(variance * I(traitSpaceDim(community)))
    generateMutant(community, config, covMat, randomSpecies)
end


# User-facing methods with population-weighted random selection

"""
    generateMutantWeighted(community::Community{T, AuxClasses}, config::EcoEvoConfig{T}, covMat::AbstractMatrix{T}) where {T, AuxClasses}

Generate a mutant with population-weighted parent selection (full covariance matrix).
"""
generateMutantWeighted(
    community::Community{T, AuxClasses},
    config::EcoEvoConfig{T},
    covMat::AbstractMatrix{T}
) where {T<:Real, AuxClasses} =
    generateMutant(community, config, covMat, weightedRandomSpecies)


"""
    generateMutantWeighted(community::Community{T, AuxClasses}, config::EcoEvoConfig{T}, variance::T) where {T, AuxClasses}

Generate a mutant with population-weighted parent selection (diagonal covariance).
"""
function generateMutantWeighted(
        community::Community{T, AuxClasses},
        config::EcoEvoConfig{T},
        variance::T
    ) where {T<:Real, AuxClasses}
    variance > 0 || throw(ArgumentError("variance must be positive"))
    covMat = Matrix{T}(variance * I(traitSpaceDim(community)))
    generateMutant(community, config, covMat, weightedRandomSpecies)
end


"""
    singleEvoStep(community::Community{T, AuxClasses}, config::EcoEvoConfig{T}) where {T, AuxClasses}

Perform one eco-evolutionary step:
1. Add a mutant species using the mutation generator from config
2. Integrate ecological dynamics using ecoDyn
3. Remove species below extinction threshold
4. Return the resulting community
"""
function singleEvoStep(
        community::Community{T, AuxClasses},
        config::EcoEvoConfig{T}
    ) where {T<:Real, AuxClasses}

    @pipe community |>
        config.mutationGenerator(_, config) |>
        ecoDyn(_, config) |>
        removeExtinct(_, config.extThreshold)
end


"""
    evolve!(history::EvoHistory{T, AuxClasses}, config::EcoEvoConfig{T}, nMutEvents::Int) where {T, AuxClasses}

Run the main eco-evolutionary simulation loop:
1. Start with the last community in history
2. Integrate ecological dynamics
3. Introduce a mutant
4. Remove extinct species
5. Repeat for nMutEvents mutation events

Modifies history in place by appending new communities.
"""
function evolve!(
        history::EvoHistory{T, AuxClasses},
        config::EcoEvoConfig{T},
        nMutEvents::Int;
        showProgress::Bool = true
    ) where {T<:Real, AuxClasses}
    nMutEvents >= 0 || throw(ArgumentError("nMutEvents must be non-negative"))

    # Get the current community (last in history)
    length(history.history) > 0 ||
        throw(ArgumentError("History must contain at least one community"))
    currentComm = history.history[end]

    # Perform nMutEvents evolutionary steps
    if showProgress
        print_interval = max(1, div(nMutEvents, 100))  # Update every ~1%
    end

    for i in 1:nMutEvents
        # Apply one evolutionary step (mutant addition → dynamics → extinction removal)
        currentComm = singleEvoStep(currentComm, config)

        # Append to history
        push!(history.history, currentComm)

        # Update progress indicator
        if showProgress && (i % print_interval == 0 || i == nMutEvents)
            pct = round(100 * i / nMutEvents, digits=1)
            print("\rProgress: $i / $nMutEvents ($pct%)")
            flush(stdout)
        end
    end
    showProgress && println()  # New line after completion

    return history
end


"""
    evolve!(community::Community{T, AuxClasses}, config::EcoEvoConfig{T}, nMutEvents::Int; showProgress::Bool=true) where {T, AuxClasses}

Convenience method that creates an EvoHistory from a single initial community
and runs the evolutionary simulation.

Returns the EvoHistory containing the initial community plus all evolved communities.
"""
function evolve!(
        community::Community{T, AuxClasses},
        config::EcoEvoConfig{T},
        nMutEvents::Int;
        showProgress::Bool = true
    ) where {T<:Real, AuxClasses}
    # Create history with initial community
    history = EvoHistory(community)

    # Run evolution
    evolve!(history, config, nMutEvents; showProgress=showProgress)

    # Return the history for user access
    return history
end


"""
    evolve(history::EvoHistory{T, AuxClasses}, config::EcoEvoConfig{T}, nMutEvents::Int; showProgress::Bool=true) where {T, AuxClasses}

Non-mutating version of `evolve!`. Creates a deep copy of the history and runs
the evolutionary simulation on the copy, leaving the original unchanged.

Returns a new EvoHistory with the additional evolved communities.

# Example
```jldoctest
julia> using EcoEvoSim

julia> comm = Community([1.0], [0.3], Float64[]);

julia> growthFn = z -> sum(z.^2) / (3 + sum(z.^2));

julia> kernelFn = (zi, zj) -> -(tanh(sum(zi .- zj) / 0.3) + 1) / 2;

julia> config = EcoEvoConfig(
           ecoDyn = lotkaVolterra(growthFn, kernelFn),
           mutationGenerator = (c, cfg) -> generateMutant(c, cfg, 0.002^2),
           integrationParams = IntegrationParams(maxTime = 1.0e8),
           invaderPopsize = 0.001,
           extThreshold = 0.003
       );

julia> history1 = evolve!(comm, config, 10);

julia> history2 = evolve(history1, config, 5);

julia> length(history1)  # Original unchanged
11

julia> length(history2)  # New history extended
16
```
"""
function evolve(
        history::EvoHistory{T, AuxClasses},
        config::EcoEvoConfig{T},
        nMutEvents::Int;
        showProgress::Bool = true
    ) where {T<:Real, AuxClasses}
    # Create a deep copy to avoid mutating the original
    history_copy = EvoHistory{T, AuxClasses}(deepcopy(history.history))

    # Run evolution on the copy
    evolve!(history_copy, config, nMutEvents; showProgress=showProgress)

    return history_copy
end


"""
    evolve(community::Community{T, AuxClasses}, config::EcoEvoConfig{T}, nMutEvents::Int; showProgress::Bool=true) where {T, AuxClasses}

Non-mutating version of `evolve!` that starts from a Community.

Creates an EvoHistory from the initial community and runs the evolutionary simulation.
This is functionally equivalent to `evolve!` for a Community since a new history is
always created, but provided for API consistency.

Returns the EvoHistory containing the initial community plus all evolved communities.

# Example
```jldoctest
julia> using EcoEvoSim

julia> comm = Community([1.0], [0.3], Float64[]);

julia> growthFn = z -> sum(z.^2) / (3 + sum(z.^2));

julia> kernelFn = (zi, zj) -> -(tanh(sum(zi .- zj) / 0.3) + 1) / 2;

julia> config = EcoEvoConfig(
           ecoDyn = lotkaVolterra(growthFn, kernelFn),
           mutationGenerator = (c, cfg) -> generateMutant(c, cfg, 0.002^2),
           integrationParams = IntegrationParams(maxTime = 1.0e8),
           invaderPopsize = 0.001,
           extThreshold = 0.003
       );

julia> history = evolve(comm, config, 10);

julia> length(history)
11
```
"""
function evolve(
        community::Community{T, AuxClasses},
        config::EcoEvoConfig{T},
        nMutEvents::Int;
        showProgress::Bool = true
    ) where {T<:Real, AuxClasses}
    # For Community, this is equivalent to evolve! since we always create new history
    # But we provide it for consistency
    return evolve!(community, config, nMutEvents; showProgress=showProgress)
end
