# Eco-evolutionary dynamics: configuration, integration, and the evolution loop

using LinearAlgebra
using OrdinaryDiffEq
using SteadyStateDiffEq


"""
    DiscreteSS()

Pseudo-algorithm for discrete-time fixed-point iteration.

When used as the `algorithm` in `IntegrationParams`, `ecoDyn` iterates the map
`u ← f(u, nothing, t)` for up to `maxTime` steps, stopping early when
`‖u_new − u_old‖ < abstol + reltol * ‖u_old‖`.

The `abstol` and `reltol` tolerances are taken from `solver_options` (defaulting
to `1e-8` and `1e-6` respectively).  The community time is advanced by one unit
per step.  Setting `maxTime = Inf` (the natural default) imposes no cap; a finite
value acts as a safety limit on the number of iterations.

Analogous to `DynamicSS()` for continuous-time systems.
"""
struct DiscreteSS end


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

# Discrete-time recursion (maxTime = number of generations)
    params = IntegrationParams(maxTime = 100.0, algorithm = FunctionMap())

# Discrete fixed-point iteration (stop when the map converges; maxTime is a safety cap)
    params = IntegrationParams(maxTime = Inf, algorithm = DiscreteSS())
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
IntegrationParams(maxTime::Real; algorithm::Alg = Rodas5(),
                  solver_options::NamedTuple = (abstol=1e-8, reltol=1e-6)) where {Alg} =
    IntegrationParams(float(maxTime), algorithm, solver_options)

# Keyword argument constructor with explicit abstol/reltol (captures any additional kwargs)
function IntegrationParams(; maxTime::Real,
                            algorithm::Alg = Rodas5(),
                            abstol = 1e-8,
                            reltol = 1e-6,
                            kwargs...) where {Alg}
    T = typeof(float(maxTime))
    solver_options = merge((abstol=convert(T, abstol), reltol=convert(T, reltol)), kwargs)
    IntegrationParams(convert(T, maxTime), algorithm, solver_options)
end


"""
    EcoEvoConfig{T<:Real, EcoDynamics, MutGenerator, Alg}

Configuration for eco-evolutionary simulations.

# Fields
- `ecoDyn::EcoDynamics`: Factory function taking `Community` and returning ODE function
- `mutationGenerator::MutGenerator`: Function generating mutants from community
  (takes only community as input; returned by a factory such as `generateMutant`)
- `integrationParams::IntegrationParams`: ODE integration parameters
- `extThreshold::T`: Population size below which species go extinct

# Example
```julia
config = EcoEvoConfig(
    ecoDyn = lotkaVolterra(growthFn, interactionFn),
    mutationGenerator = generateMutant(invaderPopsize=0.001, variance=0.01),
    integrationParams = IntegrationParams(maxTime = 50.0),
    extThreshold = 0.003
)
```
"""
struct EcoEvoConfig{T<:Real, EcoDynamics, MutGenerator, Alg}
    ecoDyn :: EcoDynamics  # Factory function: Community -> (u, p, t) -> du
    mutationGenerator :: MutGenerator   # Function generating mutant traits
    integrationParams :: IntegrationParams{T, Alg}
    extThreshold :: T  # Population size below which species are considered extinct
    function EcoEvoConfig{T, ED, MG, Alg}(
            ecoDyn::ED,
            mutationGenerator::MG,
            integrationParams::IntegrationParams{T, Alg},
            extThreshold::T
        ) where {T<:Real, ED, MG, Alg}
        extThreshold > 0 || throw(ArgumentError("extThreshold must be positive"))
        new{T, ED, MG, Alg}(ecoDyn, mutationGenerator, integrationParams, extThreshold)
    end
end

function EcoEvoConfig(
        ecoDyn::ED,
        mutationGenerator::MG,
        integrationParams::IntegrationParams{T, Alg},
        extThreshold::T
    ) where {T<:Real, ED, MG, Alg}
    EcoEvoConfig{T, ED, MG, Alg}(
        ecoDyn, mutationGenerator, integrationParams, extThreshold
    )
end

# Keyword argument constructor
function EcoEvoConfig(; ecoDyn::ED,
                        mutationGenerator::MG,
                        integrationParams::IntegrationParams{T, Alg},
                        extThreshold::Real) where {T<:Real, ED, MG, Alg}
    EcoEvoConfig{T, ED, MG, Alg}(ecoDyn, mutationGenerator, integrationParams,
                 convert(T, extThreshold))
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
        append!(u0, popsize(sp))
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
        stageClasses = length(popsize(sp))
        newPopsize = PopulationSize(u[idx:idx+stageClasses-1])
        newSpecies_i = Species(newPopsize, sp.trait)
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
    ecoDyn(community::Community{TC, AuxClasses}, config::EcoEvoConfig{TE}) where {TC<:Real, TE<:Real, AuxClasses}

Integrate ecological dynamics for the specified time period.
Extracts population sizes and auxiliary variables, integrates the system,
and returns a new community with updated values.

The `config.ecoDyn` field should be a factory function `Community -> Function`.
Depending on the algorithm in `config.integrationParams`:
- **Continuous time** (default, e.g. `Rodas5()`, `Tsit5()`): the returned
  function has signature `(u, p, t) -> du` and is used as an ODE right-hand side.
- **Steady state** (`DynamicSS()`, `SSRootfind()`): same signature; solved via
  `SteadyStateProblem`.
- **Discrete time** (`FunctionMap()`): the returned function has signature
  `(u, p, t) -> u_next` and is applied as a recursion map.  `maxTime` is then
  the number of discrete generations to iterate.
- **Discrete fixed-point** (`DiscreteSS()`): same map semantics, but iteration
  stops early once the state has converged within `abstol`/`reltol`.  `maxTime`
  is a safety cap on the number of steps.
"""
function ecoDyn(
        community::Community{TC, AuxClasses},
        config::EcoEvoConfig{TE}
    ) where {TC<:Real, TE<:Real, AuxClasses}

    # Unpack community into state vector
    u0 = unpackCommunity(community)

    # Generate the ODE function for this specific community
    ode_fn = config.ecoDyn(community)

    # Check if we're using a steady-state solver
    alg = config.integrationParams.algorithm
    if alg isa DynamicSS
        # For steady-state integration, use SteadyStateProblem with DynamicSS.
        # Build an effective DynamicSS with an explicit inner ODE solver; if none
        # was provided by the user, fall back to Rodas5().
        inner_alg = alg.alg === nothing ? Rodas5() : alg.alg
        maxTime = config.integrationParams.maxTime
        effective_alg = isfinite(maxTime) ? DynamicSS(inner_alg; tspan = maxTime) :
                                            DynamicSS(inner_alg)
        prob = SteadyStateProblem(ode_fn, u0)
        sol = solve(prob, effective_alg; config.integrationParams.solver_options...)
        u_final = sol.u
        t_final = commTime(community)
    elseif alg isa SSRootfind
        # SSRootfind uses initial conditions as a starting guess for rootfinding
        prob = SteadyStateProblem(ode_fn, u0)
        sol = solve(prob, config.integrationParams.algorithm;
                    config.integrationParams.solver_options...)
        u_final = sol.u
        # For rootfinding, time is meaningless - keep original time
        t_final = commTime(community)
    elseif alg isa FunctionMap
        # Discrete-time recursion: interpret the factory output as u(t+1) = f(u(t), p, t)
        # rather than a derivative.  maxTime is the number of discrete steps.
        tspan = (commTime(community), commTime(community) + config.integrationParams.maxTime)
        prob = DiscreteProblem(ode_fn, u0, tspan)
        sol = solve(prob, FunctionMap(); config.integrationParams.solver_options...)
        u_final = sol.u[end]
        t_final = convert(TC, sol.t[end])
    elseif alg isa DiscreteSS
        # Discrete fixed-point iteration: iterate u ← f(u, nothing, t) until
        # ‖u_new − u‖ < abstol + reltol * ‖u‖, or until maxTime steps.
        opts = config.integrationParams.solver_options
        abstol_val = TC(get(opts, :abstol, 1e-8))
        reltol_val = TC(get(opts, :reltol, 1e-6))
        maxSteps   = isfinite(config.integrationParams.maxTime) ?
                         round(Int, config.integrationParams.maxTime) : typemax(Int)
        u   = copy(u0)
        t   = commTime(community)
        for _ in 1:maxSteps
            u_new = ode_fn(u, nothing, t)
            if norm(u_new .- u) < abstol_val + reltol_val * norm(u)
                u = u_new
                break
            end
            u = u_new
            t += one(TC)
        end
        u_final = u
        t_final = t
    else
        # Use ODEProblem for time-integration solvers
        tspan = (commTime(community), commTime(community) + config.integrationParams.maxTime)
        prob = ODEProblem(ode_fn, u0, tspan)
        sol = solve(prob, config.integrationParams.algorithm;
                    config.integrationParams.solver_options...)
        u_final = sol.u[end]
        t_final = sol.t[end]
    end

    return packCommunity(u_final, community, t_final)
end


"""
    ecoDyn(community::Community, configs::AbstractVector)

Multi-stage variant of [`ecoDyn`](@ref) that applies each config's ecological
dynamics in sequence, with extinction removal after each stage.  The first
config's `ecoDyn` is applied to `community`, the second to the result of the
first, and so on.

Useful for a two-step equilibration strategy where a transient ODE integration
(e.g. `DynamicSS`) is followed by a Newton polish (e.g. `SSRootfind`):

```julia
ecoDyn(community, [configTransient, configPolish])
```
"""
function ecoDyn(
        community::Community{TC, AuxClasses},
        configs::AbstractVector
    ) where {TC<:Real, AuxClasses}
    isempty(configs) && throw(ArgumentError("configs must be non-empty"))
    c = community
    for cfg in configs
        c = ecoDyn(c, cfg)
        c = removeExtinct(c, cfg.extThreshold)
    end
    return c
end


"""
    ecoDynTimeSeries(community::Community{TC, AuxClasses}, config::EcoEvoConfig{TE};
                     stepsize=nothing) where {TC<:Real, TE<:Real, AuxClasses}

Integrate ecological dynamics and return the full trajectory as a
`Vector{Community}`, one entry per saved time point (including the initial state).

Behaves like `ecoDyn` but captures intermediate states:
- **Continuous time** (e.g. `Rodas5()`, `Tsit5()`): `stepsize` is the gap in
  time between successive saved states (starting from `community.time`). Omit to
  use the solver's adaptive output points.
- **Discrete time** (`FunctionMap()`): `stepsize` is passed to `solve()` as the
  generation gap between saves; omit to save every generation.
- **Discrete fixed-point** (`DiscreteSS()`): every iteration is saved; `stepsize`
  is ignored.
- **Steady-state** (`DynamicSS()`, `SSRootfind()`): not supported — these
  compute a fixed point directly with no intermediate trajectory. Use `ecoDyn`
  instead.

# Example
```julia
# Save trajectory every 10 time units
ts = ecoDynTimeSeries(community, config; stepsize = 10)
popsizes(ts[end], 1)  # population sizes of species 1 at final saved time
```
"""
function ecoDynTimeSeries(
        community::Community{TC, AuxClasses},
        config::EcoEvoConfig{TE};
        stepsize = nothing
    ) where {TC<:Real, TE<:Real, AuxClasses}

    alg = config.integrationParams.algorithm

    alg isa DynamicSS && throw(ArgumentError(
        "ecoDynTimeSeries does not support DynamicSS (no intermediate trajectory). " *
        "Use ecoDyn instead."
    ))
    alg isa SSRootfind && throw(ArgumentError(
        "ecoDynTimeSeries does not support SSRootfind (no intermediate trajectory). " *
        "Use ecoDyn instead."
    ))

    u0     = unpackCommunity(community)
    ode_fn = config.ecoDyn(community)
    opts   = config.integrationParams.solver_options

    if alg isa FunctionMap
        tspan = (commTime(community), commTime(community) + config.integrationParams.maxTime)
        prob  = DiscreteProblem(ode_fn, u0, tspan)
        sol   = stepsize !== nothing ?
            solve(prob, FunctionMap(); opts..., saveat = stepsize) :
            solve(prob, FunctionMap(); opts...)
        return [packCommunity(sol.u[i], community, sol.t[i])
                for i in eachindex(sol.t)]

    elseif alg isa DiscreteSS
        abstol_val = TC(get(opts, :abstol, 1e-8))
        reltol_val = TC(get(opts, :reltol, 1e-6))
        maxSteps   = isfinite(config.integrationParams.maxTime) ?
                         round(Int, config.integrationParams.maxTime) : typemax(Int)
        u      = copy(u0)
        t      = commTime(community)
        states = [packCommunity(u, community, t)]
        for _ in 1:maxSteps
            u_new     = ode_fn(u, nothing, t)
            converged = norm(u_new .- u) < abstol_val + reltol_val * norm(u)
            u  = u_new
            t += one(TC)
            push!(states, packCommunity(u, community, t))
            converged && break
        end
        return states

    else
        # Continuous-time ODE
        tspan = (commTime(community), commTime(community) + config.integrationParams.maxTime)
        prob  = ODEProblem(ode_fn, u0, tspan)
        sol   = stepsize !== nothing ?
            solve(prob, alg; opts..., saveat = stepsize) :
            solve(prob, alg; opts...)
        return [packCommunity(sol.u[i], community, sol.t[i])
                for i in eachindex(sol.t)]
    end
end


"""
    singleEvoStep(community::Community{TC, AuxClasses},
                  config::EcoEvoConfig{TE}) where {TC<:Real, TE<:Real, AuxClasses}

Perform one eco-evolutionary step:
1. Add a mutant species using the mutation generator from config
2. Integrate ecological dynamics using ecoDyn
3. Remove species below extinction threshold
4. Return the resulting community
"""
function singleEvoStep(
        community::Community{TC, AuxClasses},
        config::EcoEvoConfig{TE}
    ) where {TC<:Real, TE<:Real, AuxClasses}

    community |>
        config.mutationGenerator |>
        c -> ecoDyn(c, config) |>
        c -> removeExtinct(c, config.extThreshold)
end


"""
    singleEvoStep(community, configs::AbstractVector)

Multi-stage variant of [`singleEvoStep`](@ref) that chains several integration
stages within a single evolutionary step.

The **first** element of `configs` supplies the mutation generator; every
element supplies its own `ecoDyn` + `extThreshold` stage, applied in sequence.
This is useful for a two-step equilibration strategy: integrate away fast
transients with one solver, then polish with a Newton-based rootfinder:

```julia
configEvo    = EcoEvoConfig(lotkaVolterra(...), generateMutant(...),
                            IntegrationParams(2000.0, Rodas5()), 0.003)
configPolish = EcoEvoConfig(lotkaVolterra(...), noMutation,
                            IntegrationParams(1e8, SSRootfind();
                                             abstol=1e-12, reltol=1e-12), 0.003)
history = evolve(community, [configEvo, configPolish], nSteps)
```
"""
function singleEvoStep(
        community::Community{TC, AuxClasses},
        configs::AbstractVector
    ) where {TC<:Real, AuxClasses}
    isempty(configs) && throw(ArgumentError("configs must be non-empty"))
    # Mutation from the first config only
    c = configs[1].mutationGenerator(community)
    # Apply each config's ecoDyn + extinction removal in sequence
    for cfg in configs
        c = ecoDyn(c, cfg)
        c = removeExtinct(c, cfg.extThreshold)
    end
    return c
end


"""
    evolve!(history::EvoHistory{TC, AuxClasses}, config::EcoEvoConfig{TE},
            nMutEvents::Int) where {TC<:Real, TE<:Real, AuxClasses}

Run the main eco-evolutionary simulation loop:
1. Start with the last community in history
2. Integrate ecological dynamics
3. Introduce a mutant
4. Remove extinct species
5. Repeat for nMutEvents mutation events

Modifies history in place by appending new communities.
"""
function evolve!(
        history::EvoHistory{TC, AuxClasses},
        config::EcoEvoConfig{TE},
        nMutEvents::Int;
        showProgress::Bool = true
    ) where {TC<:Real, TE<:Real, AuxClasses}
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
    evolve!(community::Community{TC, AuxClasses}, config::EcoEvoConfig{TE}, nMutEvents::Int;
            showProgress::Bool=true) where {TC<:Real, TE<:Real, AuxClasses}

Convenience method that creates an EvoHistory from a single initial community
and runs the evolutionary simulation.

Returns the EvoHistory containing the initial community plus all evolved communities.
"""
function evolve!(
        community::Community{TC, AuxClasses},
        config::EcoEvoConfig{TE},
        nMutEvents::Int;
        showProgress::Bool = true
    ) where {TC<:Real, TE<:Real, AuxClasses}
    # Create history with initial community
    history = EvoHistory(community)

    # Run evolution
    evolve!(history, config, nMutEvents; showProgress=showProgress)

    # Return the history for user access
    return history
end


"""
    evolve(history::EvoHistory{TC, AuxClasses}, config::EcoEvoConfig{TE}, nMutEvents::Int;
           showProgress::Bool=true) where {TC<:Real, TE<:Real, AuxClasses}

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
           mutationGenerator = generateMutant(invaderPopsize=0.001, variance=0.002^2),
           integrationParams = IntegrationParams(maxTime = 1.0e8),
           extThreshold = 0.003
       );

julia> history1 = evolve(comm, config, 10; showProgress=false);

julia> history2 = evolve(history1, config, 5; showProgress=false);

julia> length(history1)  # Original unchanged
11

julia> length(history2)  # New history extended
16
```
"""
function evolve(
        history::EvoHistory{TC, AuxClasses},
        config::EcoEvoConfig{TE},
        nMutEvents::Int;
        showProgress::Bool = true
    ) where {TC<:Real, TE<:Real, AuxClasses}
    # Create a deep copy to avoid mutating the original
    history_copy = EvoHistory{TC, AuxClasses}(deepcopy(history.history))

    # Run evolution on the copy
    evolve!(history_copy, config, nMutEvents; showProgress=showProgress)

    return history_copy
end


"""
    evolve(community::Community{TC, AuxClasses}, config::EcoEvoConfig{TE}, nMutEvents::Int;
           showProgress::Bool=true) where {TC<:Real, TE<:Real, AuxClasses}

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
           mutationGenerator = generateMutant(invaderPopsize=0.001, variance=0.002^2),
           integrationParams = IntegrationParams(maxTime = 1.0e8),
           extThreshold = 0.003
       );

julia> history = evolve(comm, config, 10; showProgress=false);

julia> length(history)
11
```
"""
function evolve(
        community::Community{TC, AuxClasses},
        config::EcoEvoConfig{TE},
        nMutEvents::Int;
        showProgress::Bool = true
    ) where {TC<:Real, TE<:Real, AuxClasses}
    # For Community, this is equivalent to evolve! since we always create new history
    # But we provide it for consistency
    return evolve!(community, config, nMutEvents; showProgress=showProgress)
end


# ── Multi-stage evolve overloads (AbstractVector of configs) ─────────────────

"""
    evolve!(history::EvoHistory, configs::AbstractVector, nMutEvents; showProgress=true)

Multi-stage variant of [`evolve!`](@ref) that uses several integration stages
per evolutionary step.  See [`singleEvoStep(community, configs)`](@ref) for the
per-step semantics (mutation from `configs[1]`; each config's `ecoDyn` +
`extThreshold` applied in sequence).
"""
function evolve!(
        history::EvoHistory{TC, AuxClasses},
        configs::AbstractVector,
        nMutEvents::Int;
        showProgress::Bool = true
    ) where {TC<:Real, AuxClasses}
    nMutEvents >= 0 || throw(ArgumentError("nMutEvents must be non-negative"))
    isempty(configs) && throw(ArgumentError("configs must be non-empty"))
    length(history.history) > 0 ||
        throw(ArgumentError("History must contain at least one community"))
    currentComm = history.history[end]

    if showProgress
        print_interval = max(1, div(nMutEvents, 100))
    end

    for i in 1:nMutEvents
        currentComm = singleEvoStep(currentComm, configs)
        push!(history.history, currentComm)

        if showProgress && (i % print_interval == 0 || i == nMutEvents)
            pct = round(100 * i / nMutEvents, digits=1)
            print("\rProgress: $i / $nMutEvents ($pct%)")
            flush(stdout)
        end
    end
    showProgress && println()

    return history
end

"""
    evolve!(community::Community, configs::AbstractVector, nMutEvents; showProgress=true)

Multi-stage variant of [`evolve!`](@ref) starting from a `Community`.
"""
function evolve!(
        community::Community{TC, AuxClasses},
        configs::AbstractVector,
        nMutEvents::Int;
        showProgress::Bool = true
    ) where {TC<:Real, AuxClasses}
    history = EvoHistory(community)
    evolve!(history, configs, nMutEvents; showProgress=showProgress)
    return history
end

"""
    evolve(history::EvoHistory, configs::AbstractVector, nMutEvents; showProgress=true)

Non-mutating multi-stage variant of [`evolve`](@ref).  See
[`singleEvoStep(community, configs)`](@ref) for the per-step semantics.
"""
function evolve(
        history::EvoHistory{TC, AuxClasses},
        configs::AbstractVector,
        nMutEvents::Int;
        showProgress::Bool = true
    ) where {TC<:Real, AuxClasses}
    history_copy = EvoHistory{TC, AuxClasses}(deepcopy(history.history))
    evolve!(history_copy, configs, nMutEvents; showProgress=showProgress)
    return history_copy
end

"""
    evolve(community::Community, configs::AbstractVector, nMutEvents; showProgress=true)

Non-mutating multi-stage variant of [`evolve`](@ref) starting from a `Community`.
"""
function evolve(
        community::Community{TC, AuxClasses},
        configs::AbstractVector,
        nMutEvents::Int;
        showProgress::Bool = true
    ) where {TC<:Real, AuxClasses}
    return evolve!(community, configs, nMutEvents; showProgress=showProgress)
end
