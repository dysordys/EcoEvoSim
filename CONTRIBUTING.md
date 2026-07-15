# Contributing to EcoEvoSim

Thanks for your interest in EcoEvoSim! This document describes how to set up a
development environment, the conventions the codebase follows, and how to get
changes merged. It is aimed at both external contributors and future
maintainers.

EcoEvoSim is an MIT-licensed Julia package for high-level eco-evolutionary
simulation. The canonical repository is
<https://github.com/dysordys/EcoEvoSim>.


## Development environment

### Use a single, shared environment for interactive work

The most important setup tip (and the one most likely to save you time) is to
do all interactive work (REPL, running examples, ad-hoc experimentation) from
one consistent Julia environment, and let the automated test runner be the
only other context.

The recommended setup is a *shared* environment (e.g. your default
`~/.julia/environments/vX.Y`) into which EcoEvoSim is developed alongside your
plotting and development tooling:

```julia
pkg> dev /path/to/EcoEvoSim     # or: dev https://github.com/dysordys/EcoEvoSim
pkg> add Plots GLMakie Revise   # plotting backends (weak deps) + live reload
```

Then always launch the REPL the same way — plain `julia` (which activates that
default environment) — and configure your editor to match. For VS Code, the
Julia extension's `julia.environmentPath` (see `.vscode/settings.json`) should
point at this same environment so the integrated REPL and a terminal REPL share
one configuration.

**Why this matters.** Julia keeps a separate precompiled cache for each
distinct build configuration -- that is, for each resolved dependency set plus
each set of compiler flags. EcoEvoSim's own project environment (`--project=.`)
has a deliberately lean dependency tree (`OrdinaryDiffEq`, no plotting), whereas
a tooling-rich shared environment pulls in much heavier dependencies
(`DifferentialEquations`, `Makie`/`GLMakie`, ...). Alternating your interactive
REPL between these two environments forces Julia to recompile large dependency
trees back and forth. Picking one interactive environment and sticking with it
keeps the cache warm and startup fast.

A few related notes:

- **`Pkg.test()` is expected to recompile a little.** It runs in its own
  sandbox environment with `--check-bounds=yes`, which is a separate cache
  variant from interactive use. A short recompile of EcoEvoSim itself after a
  test run is normal and cheap; it is not the same as a full dependency rebuild.
- **Pin your Julia version while developing.** A `juliaup update` that bumps the
  patch version invalidates the whole precompile cache and triggers a full
  rebuild. Update deliberately, not by accident.
- **Diagnosing a surprise recompile.** Launch with `JULIA_DEBUG=loading julia
  …`; Julia will log the exact reason it rejected a cache file (mismatched
  build flags, changed dependencies, stale mtime, etc.).

### Using Revise

With `Revise` loaded before `EcoEvoSim`, edits to the source are picked up live
without restarting the REPL. Adding `using Revise` to
`~/.julia/config/startup.jl` makes this automatic.


## Repository layout

```
src/
  EcoEvoSim.jl          # module: includes, exports, extension stubs
  types.jl              # core types (PopulationSize, Phenotype, Species, Community, …)
  constructors.jl       # constructors / smart constructors
  selectors.jl          # accessors and queries over the core types
  utils.jl              # community/history manipulation helpers
  mutgen.jl             # mutant generation (generateMutant* factories, noMutation)
  ecoevo.jl             # config types, ecological integration, the evolution loop
  models.jl             # model builders (unstructuredModel, structuredModel, lotkaVolterra)
  show-methods.jl       # pretty-printing
ext/
  EcoEvoSimPlotsExt.jl  # plotting via Plots (weak dependency)
  EcoEvoSimGLMakieExt.jl# interactive plotting via GLMakie (weak dependency)
test/
  runtests.jl           # includes the individual test files
  test-*.jl             # one test file per source file (see below)
examples/               # runnable, self-contained example scripts
```

### Test files mirror source files

Each source file has a corresponding `test/test-<name>.jl`. When you change
`src/foo.jl`, the matching tests live in `test/test-foo.jl`. Keeping this
one-to-one mapping is a deliberate maintenance aid — please preserve it when
adding or moving code.


## Coding conventions

These are descriptive of the existing code; match the surrounding style.

- **Naming.** Types are `PascalCase` (`Community`, `EcoEvoConfig`); functions and
  variables are `camelCase` (`generateMutant`, `ecoDynTimeSeries`). Internal
  helpers not meant for users are prefixed with an underscore
  (e.g. `_makeMutantFactory`).
- **Docstrings.** Every exported type and function has a docstring with a
  signature line, a prose description, and -- where useful -- `# Arguments`,
  `# Returns`, and `# Example` sections. New public API should be documented this
  way. The heavy docstring investment is intentional; it is the package's
  primary user documentation.
- **Doctested examples.** Examples in docstrings are validated. Runnable examples
  are exercised in `test/test-docstring-examples.jl`; examples written as
  ```` ```jldoctest ```` blocks are additionally checked by Documenter once the
  docs build is in place. If you add or change an example, update the
  corresponding check.
- **The factory pattern.** User-facing mutation generators (`generateMutant`,
  `generateMutantWeighted`, etc.) are *factories*: keyword
  constructors that return a `Community -> Community` closure for the
  `mutationGenerator` field of `EcoEvoConfig`. Follow this pattern for new
  generators, and route shared logic through `_makeMutantFactory`.
- **Validation.** Constructors and factories validate their arguments eagerly and
  throw `ArgumentError` with a clear message (see `IntegrationParams`,
  `EcoEvoConfig`, `_makeMutantFactory`).


## Tests

Run the full suite with:

```julia
pkg> test                       # from the EcoEvoSim project environment
```

or, from a shell:

```sh
julia --project=. -e 'using Pkg; Pkg.test()'
```

While iterating, you can run a single test file against an already-loaded
package, e.g. `include("test/test-mutgen.jl")`.

Conventions for tests:

- Top-level group per file: `@testset "tests_of_<area>" begin … end`, with one
  nested `@testset "testing_<specific_behavior>" begin … end` per behavior.
- Tests that rely on randomness seed the RNG and, where they check statistical
  properties, repeat over `numTests` iterations.
- New functionality should come with tests in the mirroring test file, and any
  bug fix should add a regression test.

**All tests must pass before a change is merged**, including the docstring
example tests.


## Submitting changes

1. Branch off `main` (don't commit directly to `main`).
2. Make your change, with tests and docstrings as described above.
3. Run the full test suite and make sure it passes.
4. For user-visible changes, update the `README.md` and/or `examples/` as
   appropriate, and add an entry describing the change for the next release.
5. Open a pull request against `dysordys/EcoEvoSim` with a clear description of
   the motivation and the change.


## Versioning and releases

EcoEvoSim follows [Semantic Versioning](https://semver.org/). The exported names
in `src/EcoEvoSim.jl` constitute the public API; breaking them requires a
major-version bump (or a `0.x` minor bump while pre-1.0).

When preparing a release:

- Bump `version` in `Project.toml`.
- Review the `[compat]` bounds in `Project.toml` and update them for any new or
  upgraded dependencies. Every dependency should have a compat entry.
- Ensure tests pass on the supported Julia versions (see the `julia` compat
  entry).


## Questions

For questions, bug reports, or feature ideas, please open an issue on the
[GitHub repository](https://github.com/dysordys/EcoEvoSim/issues).
