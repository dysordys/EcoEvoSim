# Copilot / Agent instructions — EcoEvoSim

Purpose: Give AI coding agents the precise, actionable context they need to be immediately useful when changing or refactoring this Julia package.

What this project is (big picture)
- `EcoEvoSim` is a small Julia package providing an eco‑evolutionary simulation core.
- Top-level module: `src/EcoEvoSim.jl` (it `include`s and re‑exports `EcoEvoCore`).
- Core logic lives in `src/EcoEvoCore.jl` — this is the primary file to edit for behavior changes.

Key files to index
- `src/EcoEvoSim.jl` (package entry) and `src/EcoEvoCore.jl` (core types & functions)
- `Project.toml` / `Manifest.toml` (package dependencies & environment)
- `test/runtests.jl`, `test/utils.jl`, and other `test/*.jl` for usage examples and property tests

How to run and verify changes (commands)
- Open REPL in project environment: `julia --project=.`
- Run full test suite: `julia --project -e 'using Pkg; Pkg.test()'` or `julia --project test/runtests.jl`
- Run a single test file by editing `test/runtests.jl` includes and re-running the same command above

Project-specific conventions & patterns (do not assume general Julia defaults)
- Uses `StaticArrays.SVector` extensively: trait vectors and density vectors are fixed-size `SVector`s; many types carry dimensions in their type parameters (e.g., `Phenotype{T, TraitDim}`, `Population{T, StageClasses}`).
  - Preserve fixed-length `SVector` usage to avoid changing performance characteristics.
- Types are heavily parametric (TraitDim, StageClasses, NAux). Be careful when changing signatures — update `packState`, `unpackState`, and `updateState!` accordingly.
- Tests use helper generators in `test/utils.jl` (PropCheck-based). Prefer adding property tests there for randomized coverage.
- Export list in `EcoEvoCore` shows the public API; follow it when making breaking changes.

Integration & dependencies
- Primary deps in `Project.toml`: `StaticArrays`, `DifferentialEquations` (extras: `Test`, `PropCheck`) — ensure changes keep compatibility with these.

Known gotchas & places to inspect
- `makeSystemState(..., time::Real = zero(eltype(time)))` — this default expression refers to `time` itself and is likely a bug; be careful when modifying or creating overloads involving default args for `time`.
- Top-level package `EcoEvoSim` simply re‑exports `EcoEvoCore`: prefer working in `EcoEvoCore` and add tests to `test/`.

How to propose changes for review
- Provide a short summary of the intent, list files edited, and include a minimal test that demonstrates the expected behavior (add or modify tests under `test/`).
- Run the full test suite locally before proposing a patch.

Useful example prompts for the Agent
- "Refactor `packState` and `unpackState` to be allocation-minimal; provide tests that assert correctness and performance guidelines using small sized SVectors." 
- "Suggest a safe fix for the default `time` argument in `makeSystemState` and add tests that cover default and explicit time values." 
- "Add property tests for `updateState!` using generators in `test/utils.jl` to ensure type stability and correctness across random inputs."

If something is unclear
- Check `src/EcoEvoCore.jl` for type definitions and function implementations, and `test/runtests.jl` for concrete usage examples.
- Ask the maintainer (György Barabás) to clarify intended semantics when behavior isn't obvious from code & tests.

---
If you'd like, I can iterate on this file to add more examples (small code snippets or sample PR message templates), or merge it with an existing `.github` guide if you have one elsewhere.