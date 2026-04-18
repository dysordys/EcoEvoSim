# EcoEvoSim

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
    mutationGenerator = c -> generateMutant(c; invaderPopsize=0.001, variance=0.002^2),
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

# Run for another 500 events - this is as easy as applying evolve() to evoHistory:
evoHistory = evolve(evoHistory, config, 500)

# Visualize the results for all 1500 mutation events:
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
unstructuredModel(; precompute = nothing) do i, n, z, nSpecies
    # return dn[i]/dt
end
```

| Argument   | Description |
|------------|-------------|
| `i`        | Current species index |
| `n`        | Population size vector (`n[j]` = pop. of species `j`) |
| `z`        | Vector of trait vectors (`z[j]` = traits of species `j`) |
| `nSpecies` | Total number of species |

The body should return the time derivative $\text{d}n_i/\text{d}t$. External
parameters are captured via closures as usual. For example, to implement the same
Lotka-Volterra model as in the example above:

```julia
ecology = unstructuredModel() do i, n, z, nSpecies
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
) do i, n, z, nSpecies, pre
    n[i] * (pre.b[i] + sum(pre.A[i, j] * n[j] for j in 1:nSpecies))
end
```

### `structuredModel` -- spatial patches or stage classes

Use this when populations are structured (e.g., multiple patches or age/stage classes).

```julia
structuredModel(; auxDynamics = nothing, precompute = nothing) do i, j, N, z, R, nSpecies, nPatches
    # return dN[i,j]/dt
end
```

| Argument   | Description |
|------------|-------------|
| `i`        | Species index |
| `j`        | Patch/stage index |
| `N`        | Density matrix (`N[i,j]` = density of species `i` in patch `j`) |
| `z`        | Vector of trait vectors (`z[i]` = traits of species `i`) |
| `R`        | Auxiliary variable vector (e.g., resource levels per patch) |
| `nSpecies` | Total number of species |
| `nPatches` | Number of patches/stages |

The optional `auxDynamics` keyword accepts a function
`(R, N, z, nSpecies, nPatches) -> Vector` that returns time derivatives for
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



## Examples

See the `examples/` directory for more detailed demonstrations:

- `evo-demo.jl`: Basic evolutionary simulation using Lotka-Volterra competition.
- `evo-demo-steady.jl`: As `evo-demo.jl`, but explicitly integrating the ecological
dynamics until steady state is reached, instead of integrating for a prescribed
number of time units.
- `evo-demo-2D.jl`: Lotka-Volterra competition in a two-dimensional trait space.
- `evo-demo-helper-1.jl`: Same as `evo-demo.jl`, but defining the model with
`unstructuredModel` instead of the pre-built `lotkaVolterra`.
- `evo-demo-helper-2.jl`: Same as `evo-demo.jl`, but defining the model with
`unstructuredModel` instead of the pre-built `lotkaVolterra`, and also relying on
the `precompute` method to make the running of the model more efficient.
- `two-patch.jl`: A model with two spatial patches that have different environmental
conditions and limited migration in between. Depending on the parameters, either a
single generalist prevails or two specialist species emerge.
- `two-patch-helper.jl`: Same as `two-patch.jl`, but using the `structuredModel`
helper method for a much simpler definition.
- `two-patch-resource.jl`: Two-patch model with explicit resource dynamics -- gives
identical results to the implicit-resource model.
- `two-patch-resource-helper.jl`: Same as `two-patch-resource.jl`, but using
the `structuredModel` helper method together with `auxDynamics`.



## Documentation

For more details on the implementation and API, see the source files in `src/`:

- `EcoEvoSim.jl`: The heart of the system where everything is brought together.
- `basic-types.jl`: Core type definitions.
- `basic-constructors.jl`: Constructor functions.
- `basic-selectors.jl`: Selector functions (one should favor these over explicit
field access).
- `basic-utils.jl`: Simple utility functions that help work with the system's custom
types.
- `show-methods.jl`: Custom pretty-printing for the system's types.
- `ecoevo.jl`: Main eco-evolutionary simulation logic.
- `models.jl`: Predefined ecological models, plus a syntax for conveniently defining
any model.
- `visualize.jl`: Plotting utilities.

These source files also have corresponding test suites, under `/test`. for example,
tests for `ecoevo.jl` are in `/test/test-ecoevo.jl`, and so on. The `runtests.jl`
script runs all tests across these test files.



## License

This package is licensed under the MIT License. See the [LICENSE](LICENSE) file
for more details.



## Author

György Barabás (dysordys@protonmail.com)
