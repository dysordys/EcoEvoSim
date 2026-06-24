```@meta
CurrentModule = EcoEvoSim
```

# EcoEvoSim

[EcoEvoSim](https://github.com/dysordys/EcoEvoSim) is a Julia package for
simulating the eco-evolutionary dynamics of clonally reproducing entities. It
provides a high-level interface for studying how ecological interactions and
evolutionary processes shape community composition and trait distributions over
time.

You supply a system of equations governing the ecological dynamics together with
the traits that can evolve, and the package automates the eco-evolutionary loop:

1. Start with a community of phenotypes.
2. Integrate the ecological equations for a specified time (or until equilibrium).
3. Introduce a new mutant at low density, with a trait slightly perturbed from a
   randomly chosen resident.
4. Repeat from step 2 for a specified number of mutation events.

## Installation

EcoEvoSim is not yet registered in the Julia General Registry. Add it directly
from the repository:

```julia
using Pkg
Pkg.add(url = "https://github.com/dysordys/EcoEvoSim.git")
```

## Quick start

```julia
using EcoEvoSim
using Plots

# Quadratic intrinsic growth (maximum at z = 0, roots at ±0.5):
growthFn(z) = 1 - sum(z .^ 2) / 0.5^2
# Gaussian competition kernel (width 0.15):
kernelFn(zi, zj) = -exp(-sum((zi .- zj) .^ 2) / (2 * 0.15^2))

config = EcoEvoConfig(
    ecoDyn = lotkaVolterra(growthFn, kernelFn),
    mutationGenerator = generateMutant(invaderPopsize = 0.001, variance = 0.002^2),
    integrationParams = IntegrationParams(maxTime = 1.0e12),
    extThreshold = 0.003,
)

initCommunity = Community([1.0], [-0.3])         # one species, density 1, trait -0.3
initCommunity = ecoDyn(initCommunity, config)    # relax to ecological equilibrium

evoHistory = evolve(initCommunity, config, 1000) # 1000 mutation events
evoHistory = evolve(evoHistory, config, 1000)    # 1000 more, continuing on

plotEvo(evoHistory)                              # visualise the trajectory
```

## Where to go next

- The [README](https://github.com/dysordys/EcoEvoSim#readme) walks through defining
  models with [`unstructuredModel`](@ref) and [`structuredModel`](@ref), the core
  concepts, and a set of worked examples.
- The [`examples/`](https://github.com/dysordys/EcoEvoSim/tree/main/examples)
  directory contains runnable, self-contained demonstrations.
- The [API Reference](@ref) documents every exported type and function.
