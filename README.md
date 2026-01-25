# EcoEvoSim

**EcoEvoSim** is a Julia package for simulating the eco-evolutionary dynamics of clonally reproducing replicating entities. The package provides a high-level interface for studying how ecological interactions and evolutionary processes shape community composition and trait distributions over time. Broadly speaking, the way the package works is that the user provides a system of equations governing the ecological dynamics of a system, plus the traits of the species that can evolve. The system then automates the following steps:

0. We start with a community of phenotypes.
1. The ecological equations are integrated for some specified number of time units (or until equilibrium is reached).
2. A new mutant is thrown into the community at a low density, and with a trait value that is a slight change from a randomly-picked resident's trait.
3. And then we repeat from Step 1, for a specified number of mutation events.



## Features

- **Flexible Ecological Dynamics**: Define custom systems of equations governing population dynamics. Populations structured by space or age/stage classes are supported.
- **Support for Auxiliary Variables**: It is possible to add auxiliary variables (such as resources) that are part of the system but do not undergo mutations.
- **Trait Evolution**: Simulate evolution in single or multi-dimensional trait spaces.
- **Mutation and Selection**: Automated mutation-selection processes with configurable parameters.
- **Community Dynamics**: One can track population sizes, trait values, and community composition through time.
- **Visualization Tools**: There are built-in plotting functions for visualizing evolutionary trajectories.
- **Data Export**: Support for exporting simulation results to CSV files for further analysis.



## Installation

EcoEvoSim is not yet registered in the Julia General Registry. To install it, clone the repository and activate the project:

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

# Define the ecological dynamics
# Quadratic intrinsic growth, with maximum at z = 0 and roots at +/-0.5:
growthFn = z -> 1 - sum(z.^2) / 0.5^2
# Gaussian interaction kernel with width 0.15:
kernelFn = (zi, zj) -> -exp(-sum((zi .- zj).^2) / (2 * 0.15^2))

# Initialize a community with a single species of pop. size 1 and trait value -0.3:
community = Community([1.0], [-0.3], Float64[])

# Configure the simulation
config = EcoEvoConfig(
    # Use Lotka-Volterra dynamics, with the specified functions growthFn and kernelFn:
    ecoDyn = lotkaVolterra(growthFn, kernelFn),
    # Generate normally distributed mutations with standard deviation 0.002:
    mutationGenerator = (comm, confg) -> generateMutant(comm, confg, 0.002^2),
    # Integrate each ecological step for 1.0e10 time units:
    integrationParams = IntegrationParams(maxTime = 1.0e10),
    # Each new mutant has a population size of 0.001:
    invaderPopsize = 0.001,
    # Species with pop. size below 0.003 after eco. simulation are removed:
    extThreshold = 0.003
)

# Run the simulation for 1000 mutation events
evoHistory = evolve!(community, config, 1000)

# Run for another 500 events - this is as easy as applying evolve!() to evoHistory:
evoHistory = evolve!(evoHistory, config, 500)

# Visualize the results
plotEvo(evoHistory)
```



## Core Concepts

### Species and Communities

The system relies heavily on the following custom types:

- `PopulationSize`: Wrapper for population densities (can be vectors, for structured populations).
- `Phenotype`: Wrapper for trait vectors. If single-dimensional, it represents a single trait value. Otherwise, it represents trait vectors in multidimensional trait spaces.
- `Species`: Combines population size and phenotype above. It has the two fields `popsize` and `trait`.
- `Community`: Holds three things: (i) `species`, which is a vector of `Species`, (ii) `aux`, a vector of auxiliary dynamical variables such as resources that do not undergo mutation (can be empty, if no such variables are needed), and (iii) `time`, measuring the number of time units for which the ecological dynamics were integrated.
- `EvoHistory`: its single field, `history`, holds a vector of `Community` states, and is thus used to represent an evolving community through time.

### Configuration

The `EcoEvoConfig` struct allows one to specify:
- The functions and equations governing the ecological dynamics.
- Methods for generating mutant phenotypes.
- Integration (hyper-)parameters.
- Invasion and extinction thresholds.



## Examples

See the `examples/` directory for more detailed demonstrations:

- `evo-demo.jl`: Basic evolutionary simulation.
- `evo-demo-2D.jl`: Evolution in two-dimensional trait space.
- `evo-demo-steady.jl`: As `evo-demo.jl`, but explicitly integrating the ecological dynamics until steady state.



## Documentation

For more details on the implementation and API, see the source files in `src/`:

- `basic-types.jl`: Core type definitions
- `basic-constructors.jl`: Constructor functions
- `basic-selectors.jl`: Selector functions (one should favor these over explicit field access)
- `basic-utils.jl`: Simple utility functions that help work with the system's custom types
- `show-methods.jl`: Custom pretty-printing for the system's types
- `ecoevo.jl`: Main eco-evolutionary simulation logic
- `models.jl`: Predefined ecological models
- `visualize.jl`: Plotting utilities
- `EcoEvoSim.jl`: The heart of the system where everything is brought together

These source files also have corresponding test suites, under `/test`. for example, tests for `ecoevo.jl` are in `/test/test-ecoevo.jl`, and so on. The `runtests.jl` script runs all tests across these test files.



## License

This package is licensed under the MIT License. See the [LICENSE](LICENSE) file for more details.



## Author

György Barabás (dysordys@protonmail.com)
