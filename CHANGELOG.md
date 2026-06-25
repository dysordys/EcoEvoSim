# Changelog

All notable changes to EcoEvoSim will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.1.0]

Initial public release.

### Added

- Eco-evolutionary simulation loop (`evolve`, `singleEvoStep`) driving a
  configurable mutation–selection process via `EcoEvoConfig` and
  `IntegrationParams`.
- Core types for communities of evolving phenotypes — `PopulationSize`,
  `Phenotype`, `Species`, `Community`, and `EvoHistory` — together with accessors
  (`selectors`) and manipulation utilities.
- Model builders `unstructuredModel` and `structuredModel`, each with optional
  auxiliary-variable dynamics (`auxDynamics`) and a `precompute` hook, plus the
  pre-built `lotkaVolterra` model.
- Mutant-generation factories: `generateMutant`, `generateMutantWeighted`,
  `generateMutantSpatial`, `generateMutantSpatialWeighted`, and `noMutation`.
- Continuous-time, steady-state, and discrete-time ecological dynamics, including
  the package's own `DiscreteSS` fixed-point iteration.
- Convenience re-exports of common solver algorithms (`Rodas5`, `RadauIIA5`,
  `Tsit5`, `Vern7`, `AutoVern7`).
- Plotting through package extensions for Plots and GLMakie.
- Documentation site, runnable examples, and a comprehensive test suite.

[Unreleased]: https://github.com/dysordys/EcoEvoSim/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/dysordys/EcoEvoSim/releases/tag/v0.1.0
