# EcoEvoSim Copilot Notes

This document provides context when asking Copilot/Copilot Chat to work on this repo.



## Purpose and scope

This repository aims to develop a Julia package for simulating eco-evolutonary dynamics of clonally-reproducing replicating entities. Broadly speaking, the way it should work is that the user provides a system of equations governing the ecological dynamics of a system, plus the traits of the species that can evolve. The system then automates the following steps:

1. We start with some community of phenotypes;
2. We integrate the ecological equations until equilibrium (or until some specified number of time units);
3. We introduce a new mutant at a low density, and with a trait value that is a slight change from a randomly-picked resident's trait;
4. We repeat from Step 2 for a specified number of mutation events



## Organization

Source files are under `/src` (organized in subdirectories if and as necessary and logical). Tests are in `/test`. There is also `/examples`, with simple example code for demonstration purposes (mainly to show the developer how the various functionalities work at any one time).



## Basic types

Any one species is defined by two properties: its population size and its phenotype. The former is given via a struct wrapper `PopulationSize`, and may be a vector (this could represent populations structured by spatial location or stage class). The phenotype also uses a wrapper, `Phenotype`, and may also be a vector in case the trait space in which species evolve is multidimensional. `Species` then consists of two fields: a population size and a phenotype. Finally, a community is a collection of species.



## Style

- Keep APIs minimal and typed; avoid broad `Any` or implicit conversions.
- Validate inputs early and throw `ArgumentError` with clear messages.
- Prefer small, pure functions; document non-obvious behavior with brief comments.
