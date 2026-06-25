# EcoEvoSim

[![Docs (stable)](https://img.shields.io/badge/docs-stable-blue.svg)](https://dysordys.github.io/EcoEvoSim/stable)
[![Docs (dev)](https://img.shields.io/badge/docs-dev-blue.svg)](https://dysordys.github.io/EcoEvoSim/dev)
[![CI](https://github.com/dysordys/EcoEvoSim/actions/workflows/CI.yml/badge.svg)](https://github.com/dysordys/EcoEvoSim/actions/workflows/CI.yml)

**EcoEvoSim** is a Julia package for simulating the eco-evolutionary dynamics of
clonally reproducing replicating entities. The package provides a high-level interface
for studying how ecological interactions and evolutionary processes shape community
composition and trait distributions over time. Broadly speaking, the way the package
works is that the user provides a system of equations governing the ecological dynamics
of a system, plus the traits of the species that can evolve. The system then automates
the following steps:

0. We start with a community of phenotypes.
1. The ecological equations are integrated for some specified number of time units
(or until equilibrium is reached).
2. A new mutant is thrown into the community at a low density, and with a trait value
that is a slight change from a randomly-picked resident's trait.
3. And then we repeat from Step 1, for a specified number of mutation events.



## Features

- **Flexible Ecological Dynamics**: Define custom systems of equations governing
population dynamics. Populations structured by space or age/stage classes are supported.
- **Support for Auxiliary Variables**: It is possible to add auxiliary variables
(such as resources) that are part of the system but do not undergo mutations.
- **Trait Evolution**: Simulate evolution in single or multi-dimensional trait spaces.
- **Mutation and Selection**: Automated mutation-selection processes with configurable
parameters.
- **Community Dynamics**: One can track population sizes, trait values, and community
composition through time.
- **Visualization Tools**: There are built-in plotting functions for visualizing
evolutionary trajectories.
- **Data Export**: Support for exporting simulation results to CSV files for further
analysis.



## Installation

EcoEvoSim is not yet registered in the Julia General Registry. To install it, clone
the repository and activate the project:

```julia
using Pkg
Pkg.activate("/path/to/EcoEvoSim")
Pkg.instantiate()
```

Alternatively, you can add it directly from the repository:

```julia
using Pkg
Pkg.add(url="https://github.com/dysordys/EcoEvoSim.git")
```



## Quick Start

Here is a simple example demonstrating how to set up and run an evolutionary simulation:

```julia
using EcoEvoSim
using Plots

# Define quadratic intrinsic growth, with maximum at z = 0 and roots at +/-0.5:
growthFn(z) = 1 - sum(z.^2) / 0.5^2
# Define a Gaussian interaction kernel with width 0.15:
kernelFn(zi, zj) = -exp(-sum((zi .- zj).^2) / (2 * 0.15^2))

# Configure the simulation
config = EcoEvoConfig(
    # Use built-in function `lotkaVolterra` to create Lotka-Volterra dynamics
    # with specific growth and interaction functions:
    ecoDyn = lotkaVolterra(growthFn, kernelFn),
    # This is how new mutants should be generated - a resident is chosen at random,
    # we add a normally-distributed variate to its trait with mean zero and standard
    # deviation 0.002, and we initialize the new phenotype with population size 0.001:
    mutationGenerator = generateMutant(invaderPopsize = 0.001, variance = 0.002^2),
    # Integrate each ecological step for 1e12 time units:
    integrationParams = IntegrationParams(maxTime = 1.0e12),
    # Species with pop. size below 0.003 are removed after ecological simulation:
    extThreshold = 0.003
)

# Initialize a community with a single species, with density 1 and trait value -0.3:
initCommunity = Community([1.0], [-0.3])

# Make sure the initial community is at its ecological equilibrium:
initCommunity = ecoDyn(initCommunity, config)

# Run the simulation for 1000 mutation events:
evoHistory = evolve(initCommunity, config, 1000)

# Run for another 1000 events - this is as easy as applying evolve() to evoHistory:
evoHistory = evolve(evoHistory, config, 1000)

# Visualize the results for all 2000 mutation events:
plotEvo(evoHistory)
```



## Defining Models

EcoEvoSim provides two helper functions that handle all the boilerplate of creating
the factory functions required by `EcoEvoConfig.ecoDyn`. You write the equations;
the helpers take care of extracting traits, setting up the ODE state vector, and
reshaping structured populations.

### `unstructuredModel` -- scalar populations

Use this when each species has a single population size (no spatial or stage structure).

```julia
unstructuredModel(; auxDynamics = nothing, precompute = nothing) do i, n, z, aux, nSpecies
    # return dn[i]/dt
end
```

| Argument   | Description |
|------------|-------------|
| `i`        | Current species index |
| `n`        | Population size vector (`n[j]` = pop. of species `j`) |
| `z`        | Vector of trait vectors (`z[j]` = traits of species `j`) |
| `aux`      | Auxiliary variable vector (e.g., resource levels) |
| `nSpecies` | Total number of species |

The body should return the time derivative $\text{d}n_i/\text{d}t$. External
parameters are captured via closures as usual. For example, to implement the same
Lotka-Volterra model as in the example above:

```julia
ecology = unstructuredModel() do i, n, z, aux, nSpecies
    n[i] * (growthFn(z[i]) + sum(kernelFn(z[i], z[j]) * n[j] for j in 1:nSpecies))
end
```


#### Precomputation

If your equation involves trait-dependent quantities (e.g., interaction matrices)
that are expensive to recompute at every ODE timestep, pass a `precompute` function.
It runs once per community and its result is available as an extra argument `pre`.
For example, the above version of the Lotka-Volterra model is inefficient, because
it recomputes the growth rates and interaction coefficients every iteration of the
integrator, even though they stay constant until the next mutation event and so
they only need to be computed once. Below is a more efficient way of defining the
same model, using `precompute`:

```julia
ecology = unstructuredModel(
    precompute = (z, nSpecies) -> (
        b = [growthFn(z[i]) for i in 1:nSpecies],
        A = [kernelFn(z[i], z[j]) for i in 1:nSpecies, j in 1:nSpecies]
    )
) do i, n, z, aux, nSpecies, pre
    n[i] * (pre.b[i] + sum(pre.A[i, j] * n[j] for j in 1:nSpecies))
end
```

### `structuredModel` -- spatial patches or stage classes

Use this when populations are structured (e.g., multiple patches or age/stage classes).

```julia
structuredModel(; auxDynamics = nothing, precompute = nothing) do i, k, n, z, aux, nSpecies, nPatches
    # return dn[i,k]/dt
end
```

| Argument   | Description |
|------------|-------------|
| `i`        | Species index |
| `k`        | Patch/stage index |
| `n`        | Density matrix (`n[i,k]` = density of species `i` in patch `k`) |
| `z`        | Vector of trait vectors (`z[i]` = traits of species `i`) |
| `aux`      | Auxiliary variable vector (e.g., resource levels per patch) |
| `nSpecies` | Total number of species |
| `nPatches` | Number of patches/stages |

The optional `auxDynamics` keyword accepts a function
`(aux, n, z, nSpecies, nPatches) -> Vector` that returns time derivatives for
auxiliary variables (e.g., resource dynamics).

The optional `precompute` keyword accepts a function `(z, nSpecies, nPatches) -> pre`
and adds `pre` as an extra argument to the equation function, just like for
`unstructuredModel`.

### Pre-built model: `lotkaVolterra`

For standard Lotka-Volterra interactions, a convenience function is also available:

```julia
ecology = lotkaVolterra(growthFn, kernelFn)
```

where `growthFn(z)` returns the intrinsic growth rate of phenotype `z`, and
`kernelFn(z_i, z_j)` returns the interaction coefficient between phenotypes `z_i`
and `z_j`.



## Core Concepts

### Species and Communities

The system relies heavily on the following custom types:

- `PopulationSize`: Wrapper for population densities (can be vectors, for structured
populations).
- `Phenotype`: Wrapper for trait vectors. If single-dimensional, it represents a single
trait value. Otherwise, it represents trait vectors in multidimensional trait spaces.
- `Species`: Combines population size and phenotype above. It has the two fields
`popsize` and `trait`.
- `Community`: Holds three things: (i) `species`, which is a vector of `Species`,
(ii) `aux`, a vector of auxiliary dynamical variables such as resources that do not
undergo mutation (can be empty, if no such variables are needed), and (iii) `time`,
measuring the number of time units for which the ecological dynamics were integrated.
- `EvoHistory`: its single field, `history`, holds a vector of `Community` states,
and is thus used to represent an evolving community through time.

### Configuration

The `EcoEvoConfig` struct allows one to specify:
- The functions and equations governing the ecological dynamics.
- Methods for generating mutant phenotypes.
- Integration (hyper-)parameters.
- The extinction threshold, below which a phenotype is considered extinct.

### Solver algorithms

How each ecological phase is integrated is controlled by the `algorithm` passed to
`IntegrationParams` (part of the integration parameters in `EcoEvoConfig`). For
convenience, EcoEvoSim re-exports a small, curated set of algorithms, so the most
common models need only `using EcoEvoSim`:

- **Continuous-time ODE solvers:** `Rodas5` (the default; a stiff solver),
`RadauIIA5` (stiff), `Tsit5` (non-stiff), `Vern7` (high-order non-stiff), and
`AutoVern7` (auto-switching between `Vern7` and a stiff solver). These come from
[OrdinaryDiffEq](https://docs.sciml.ai/OrdinaryDiffEq/stable/).
- **Steady-state solvers:** `DynamicSS` (integrate until a steady state is reached)
and `SSRootfind` (Newton-type rootfinding). These come from
[SteadyStateDiffEq](https://docs.sciml.ai/DiffEqDocs/stable/solvers/steady_state_solve/).
- **Discrete-time solvers:** `FunctionMap` (iterate a map for a set number of steps)
and `DiscreteSS` (iterate a map until a fixed point is reached; EcoEvoSim's own).

Because `Rodas5` is the default, specifying it is optional:

```julia
IntegrationParams(maxTime = 1.0e8)                       # uses Rodas5()
IntegrationParams(maxTime = 1.0e8, algorithm = Tsit5())  # an explicit choice
```

This re-exported set is only a convenience subset. EcoEvoSim depends on the full
OrdinaryDiffEq solver suite, so any other solver can be used by loading it
explicitly:

```julia
using EcoEvoSim, OrdinaryDiffEq

config = EcoEvoConfig(
    ecoDyn = lotkaVolterra(growthFn, kernelFn),
    mutationGenerator = generateMutant(invaderPopsize = 0.001, variance = 0.002^2),
    integrationParams = IntegrationParams(maxTime = 100.0, algorithm = KenCarp4()),
    extThreshold = 0.003,
)
```

See the [OrdinaryDiffEq documentation](https://docs.sciml.ai/OrdinaryDiffEq/stable/)
for the full catalogue of available solvers.



## Examples

See the `examples/` directory for more detailed demonstrations:

- `01a-basic-cancer-growth.jl`: Directional evolution of a single trait toward an
environmental optimum, in a logistic cancer-growth model. The trait evolves
monotonically, with no evolutionary branching.
- `01b-basic-cancer-growth-with-precompute.jl`: The same model as `01a`, but using
the `precompute` helper to evaluate the constant intrinsic growth rate once per
ecological phase instead of at every integration step.
- `02a-two-patch.jl`: A spatially structured (two-patch) version in which
heterogeneity between the patches drives evolutionary branching into two
locally adapted lineages. Built with `structuredModel`; strong enough migration
suppresses the branching.
- `02b-two-patch-with-explicit-resource.jl`: The same as `02a`, but with the
limiting resource tracked explicitly (one resource per patch, via `auxDynamics`).
The parameters are chosen so it reproduces the implicit-resource model.
- `03-cancer-cell-competition.jl`: A generalized Lotka-Volterra model with a
Gaussian competition kernel, built with the `lotkaVolterra` helper. A single
ancestor branches repeatedly into a community of coexisting clones. Shown in both
1D and 2D trait spaces, and demonstrates the interactive plotting option.
- `04-competition-proliferation.jl`: A competition-proliferation (competition-
colonization) tradeoff model that yields a hierarchically structured community.
Uses `generateMutantWeighted` for density-proportional parent selection and
steady-state integration via `DynamicSS`.
- `05a-beverton-holt-basic.jl`: A discrete-time analog of `03` using the
multispecies Beverton-Holt map, selected via the `FunctionMap()` algorithm, with
`precompute` for the growth rates and interaction matrix.
- `05b-beverton-holt-steadystate.jl`: The same as `05a`, but iterating each
ecological phase to a fixed point using the `DiscreteSS()` algorithm.
- `06-environmental-variation.jl`: Evolutionary branching driven by a fluctuating,
two-season environment (a storage-effect model). An unstructured model with
explicit time dependence; a single ancestor branches into two season specialists.



## Documentation

The full rendered documentation, including the API reference, is hosted at
[dysordys.github.io/EcoEvoSim](https://dysordys.github.io/EcoEvoSim/).

For more details on the implementation, the source files in `src/` are organized as
follows:

- `EcoEvoSim.jl`: Where everything is brought together (module definition, includes,
and exports).
- `types.jl`: Core type definitions.
- `constructors.jl`: Constructor functions.
- `selectors.jl`: Selector functions (one should favor these over explicit
field access).
- `utils.jl`: Utility functions that help work with the system's custom types.
- `mutgen.jl`: Mutant generation for eco-evolutionary simulations.
- `evoevo.jl`: Eco-evolutionary configuration, ecological integration, and the
evolution loop.
- `models.jl`: Predefined ecological models, plus a syntax for conveniently defining
any model.
- `show-methods.jl`: Custom pretty-printing for the system's types.

Plotting is provided through package extensions in `ext/`, loaded automatically when
the corresponding backend is available: `EcoEvoSimPlotsExt.jl` (static plots via
Plots) and `EcoEvoSimGLMakieExt.jl` (interactive plots via GLMakie).

Most source files have a corresponding test file under `/test/` — for example, the
tests for `selectors.jl` live in `/test/test-selectors.jl`. The `runtests.jl` script
runs all of them.



## Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for
development setup, the conventions the codebase follows, and how to submit changes.



## License

This package is licensed under the MIT License. See the [LICENSE](LICENSE) file
for more details.



## Author

György Barabás (dysordys@protonmail.com)



## Acknowledgments

Development of this package was assisted by Anthropic's Claude models — much of the
codebase was written with Claude Haiku, with later refactoring, documentation, and
project setup carried out with Claude Opus 4.8.
