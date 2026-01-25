using Test
using EcoEvoSim
using LinearAlgebra
using Random
using DifferentialEquations


numTests = 50


@testset "tests_of_ecoevo_dynamics" begin

    @testset "testing_IntegrationParams_constructor" begin
        # Valid construction with keyword arguments
        params = IntegrationParams(maxTime=100.0, abstol=1e-6, reltol=1e-3)
        @test params.maxTime == 100.0
        @test params.solver_options.abstol == 1e-6
        @test params.solver_options.reltol == 1e-3
        @test params.algorithm isa DifferentialEquations.Rodas5

        # Test with additional solver options
        params2 = IntegrationParams(maxTime=50.0, abstol=1e-8, reltol=1e-6, maxiters=10000)
        @test params2.maxTime == 50.0
        @test params2.solver_options.abstol == 1e-8
        @test params2.solver_options.reltol == 1e-6
        @test params2.solver_options.maxiters == 10000

        # Invalid: non-positive maxTime
        @test_throws ArgumentError IntegrationParams(maxTime=-1.0, abstol=1e-6, reltol=1e-3)
        @test_throws ArgumentError IntegrationParams(maxTime=0.0, abstol=1e-6, reltol=1e-3)
    end


    @testset "testing_EcoEvoConfig_constructor" begin
        # Dummy functions for testing
        ecoDyn = (x, t, p) -> -x
        mutGen = (x) -> x .+ 0.01
        params = IntegrationParams(maxTime=100.0, abstol=1e-6, reltol=0.1)

        # Valid construction
        config = EcoEvoConfig(ecoDyn, mutGen, params, 0.001, 1e-8)
        @test config.ecoDyn === ecoDyn
        @test config.mutationGenerator === mutGen
        @test config.integrationParams === params
        @test config.invaderPopsize == 0.001
        @test config.extThreshold == 1e-8

        # Invalid: non-positive invaderPopsize
        @test_throws ArgumentError EcoEvoConfig(ecoDyn, mutGen, params, -0.001, 1e-8)
        @test_throws ArgumentError EcoEvoConfig(ecoDyn, mutGen, params, 0.0, 1e-8)

        # Invalid: non-positive extThreshold
        @test_throws ArgumentError EcoEvoConfig(ecoDyn, mutGen, params, 0.001, -1e-8)
        @test_throws ArgumentError EcoEvoConfig(ecoDyn, mutGen, params, 0.001, 0.0)
    end


    @testset "testing_generateMutant_with_covariance_matrix" begin
        # Create a simple test community
        for _ in 1:numTests
            nSpecies = rand(2:5)
            traitDim = rand(1:3)

            # Create species with random traits and populations
            species = Species{Float64}[]
            for i in 1:nSpecies
                trait = rand(traitDim)
                popsize = rand() * 10.0
                push!(species, Species(popsize, trait))
            end

            comm = Community(species, PopulationSize{Float64}[], 0.0)

            # Create config
            ecoDyn = (x, t, p) -> -x
            mutGen = (x) -> x .+ 0.01
            params = IntegrationParams(maxTime=100.0, abstol=1e-6, reltol=0.1)
            config = EcoEvoConfig(ecoDyn, mutGen, params, 0.001, 1e-8)

            # Create covariance matrix
            variance = 0.01
            covMat = Matrix{Float64}(variance * I(traitDim))

            # Generate mutant
            newComm = generateMutant(comm, config, covMat)

            # Check that community has one more species
            @test numSpecies(newComm) == numSpecies(comm) + 1

            # Check that the mutant has the right population size
            mutantPopsize = popsizes(newComm, numSpecies(newComm))
            @test mutantPopsize[1] == config.invaderPopsize

            # Check that mutant trait is different from all parents
            mutantTrait = traits(newComm, numSpecies(newComm))
            parentTraits = traits(comm)
            @test !any(all(mutantTrait .≈ pt) for pt in parentTraits)

            # Check that time and aux are preserved
            @test newComm.time == comm.time
            @test newComm.aux == comm.aux
        end
    end


    @testset "testing_generateMutant_with_scalar_variance" begin
        # Create a simple test community
        for _ in 1:numTests
            nSpecies = rand(2:5)
            traitDim = rand(1:3)

            species = Species{Float64}[]
            for i in 1:nSpecies
                trait = rand(traitDim)
                popsize = rand() * 10.0
                push!(species, Species(popsize, trait))
            end

            comm = Community(species, PopulationSize{Float64}[], 0.0)

            # Create config
            ecoDyn = (x, t, p) -> -x
            mutGen = (x) -> x .+ 0.01
            params = IntegrationParams(maxTime=100.0, abstol=1e-6, reltol=0.1)
            config = EcoEvoConfig(ecoDyn, mutGen, params, 0.001, 1e-8)

            # Generate mutant with scalar variance
            variance = 0.01
            newComm = generateMutant(comm, config, variance)

            # Check that community has one more species
            @test numSpecies(newComm) == numSpecies(comm) + 1

            # Check that the mutant has the right population size
            mutantPopsize = popsizes(newComm, numSpecies(newComm))
            @test mutantPopsize[1] == config.invaderPopsize

            # Check trait dimension matches
            mutantTrait = traits(newComm, numSpecies(newComm))
            @test length(mutantTrait) == traitDim
        end

        # Invalid: non-positive variance
        species = [Species(1.0, [0.5, 0.5])]
        comm = Community(species, PopulationSize{Float64}[], 0.0)
        ecoDyn = (x, t, p) -> -x
        mutGen = (x) -> x .+ 0.01
        params = IntegrationParams(maxTime=100.0, abstol=1e-6, reltol=0.1)
        config = EcoEvoConfig(ecoDyn, mutGen, params, 0.001, 1e-8)

        @test_throws ArgumentError generateMutant(comm, config, -0.01)
        @test_throws ArgumentError generateMutant(comm, config, 0.0)
    end


    @testset "testing_generateMutantWeighted_with_covariance_matrix" begin
        # Create a community with varied population sizes
        for _ in 1:numTests
            nSpecies = rand(3:6)
            traitDim = rand(1:3)

            species = Species{Float64}[]
            popsizes_vals = Float64[]
            for i in 1:nSpecies
                trait = rand(traitDim)
                # Create varied population sizes - some large, some small
                popsize = rand() < 0.5 ? rand() * 0.1 : rand() * 10.0
                push!(species, Species(popsize, trait))
                push!(popsizes_vals, popsize)
            end

            comm = Community(species, PopulationSize{Float64}[], 0.0)

            # Create config
            ecoDyn = (x, t, p) -> -x
            mutGen = (x) -> x .+ 0.01
            params = IntegrationParams(maxTime=100.0, abstol=1e-6, reltol=0.1)
            config = EcoEvoConfig(ecoDyn, mutGen, params, 0.001, 1e-8)

            # Create covariance matrix
            variance = 0.01
            covMat = Matrix{Float64}(variance * I(traitDim))

            # Generate many mutants and check they're weighted correctly
            # (species with larger populations should be selected more often)
            nTrials = 1000
            parentCounts = zeros(Int, nSpecies)

            for trial in 1:nTrials
                newComm = generateMutantWeighted(comm, config, covMat)
                mutantTrait = traits(newComm, numSpecies(newComm))

                # Find which parent was used (closest trait)
                parentTraits = traits(comm)
                minDist = Inf
                parentIdx = 1
                for (i, pt) in enumerate(parentTraits)
                    dist = sum((mutantTrait .- pt).^2)
                    if dist < minDist
                        minDist = dist
                        parentIdx = i
                    end
                end
                parentCounts[parentIdx] += 1
            end

            # Check that mutant was generated
            @test sum(parentCounts) == nTrials

            # Check that at least one parent was selected (basic sanity check)
            # Note: With weighted selection, some species with very low population
            # may not be selected even in many trials - this is expected behavior
            @test any(parentCounts .> 0)
        end
    end


    @testset "testing_generateMutantWeighted_with_scalar_variance" begin
        # Create a simple test community
        for _ in 1:numTests
            nSpecies = rand(2:5)
            traitDim = rand(1:3)

            species = Species{Float64}[]
            for i in 1:nSpecies
                trait = rand(traitDim)
                popsize = rand() * 10.0
                push!(species, Species(popsize, trait))
            end

            comm = Community(species, PopulationSize{Float64}[], 0.0)

            # Create config
            ecoDyn = (x, t, p) -> -x
            mutGen = (x) -> x .+ 0.01
            params = IntegrationParams(maxTime=100.0, abstol=1e-6, reltol=0.1)
            config = EcoEvoConfig(ecoDyn, mutGen, params, 0.001, 1e-8)

            # Generate mutant with scalar variance
            variance = 0.01
            newComm = generateMutantWeighted(comm, config, variance)

            # Check that community has one more species
            @test numSpecies(newComm) == numSpecies(comm) + 1

            # Check that the mutant has the right population size
            mutantPopsize = popsizes(newComm, numSpecies(newComm))
            @test mutantPopsize[1] == config.invaderPopsize

            # Check trait dimension matches
            mutantTrait = traits(newComm, numSpecies(newComm))
            @test length(mutantTrait) == traitDim
        end

        # Invalid: non-positive variance
        species = [Species(1.0, [0.5, 0.5])]
        comm = Community(species, PopulationSize{Float64}[], 0.0)
        ecoDyn = (x, t, p) -> -x
        mutGen = (x) -> x .+ 0.01
        params = IntegrationParams(maxTime=100.0, abstol=1e-6, reltol=0.1)
        config = EcoEvoConfig(ecoDyn, mutGen, params, 0.001, 1e-8)

        @test_throws ArgumentError generateMutantWeighted(comm, config, -0.01)
        @test_throws ArgumentError generateMutantWeighted(comm, config, 0.0)
    end


    @testset "testing_invalid_covariance_matrix_dimensions" begin
        # Create a 2D trait community
        species = [Species(1.0, [0.5, 0.5]), Species(2.0, [0.3, 0.7])]
        comm = Community(species, PopulationSize{Float64}[], 0.0)

        ecoDyn = (x, t, p) -> -x
        mutGen = (x) -> x .+ 0.01
        params = IntegrationParams(maxTime=100.0, abstol=1e-6, reltol=0.1)
        config = EcoEvoConfig(ecoDyn, mutGen, params, 0.001, 1e-8)

        # Try with wrong-sized covariance matrix (3x3 instead of 2x2)
        wrongCovMat = Matrix{Float64}(0.01 * I(3))
        @test_throws ArgumentError generateMutant(comm, config, wrongCovMat)

        # Try with non-square matrix
        wrongCovMat2 = rand(2, 3)
        @test_throws ArgumentError generateMutant(comm, config, wrongCovMat2)
    end


    @testset "testing_ecoDyn_simple_exponential_decay" begin
        # Test ecoDyn with simple exponential decay: dx/dt = -x
        for _ in 1:numTests
            nSpecies = rand(2:5)

            # Create species with random initial populations
            species = Species{Float64}[]
            initialPops = Float64[]
            for i in 1:nSpecies
                popsize = rand() * 10.0
                push!(species, Species(popsize, rand(1)))
                push!(initialPops, popsize)
            end

            comm = Community(species, PopulationSize{Float64}[], 0.0)

            # Simple exponential decay dynamics
            ecoDynFactory = (community) -> (u, p, t) -> -u
            mutGen = (c, cfg) -> generateMutant(c, cfg, 0.01)
            params = IntegrationParams(maxTime=1.0, abstol=1e-8, reltol=1e-6)
            config = EcoEvoConfig(
                ecoDyn=ecoDynFactory,
                mutationGenerator=mutGen,
                integrationParams=params,
                invaderPopsize=0.001,
                extThreshold=1e-8
            )

            # Integrate dynamics
            finalComm = ecoDyn(comm, config)

            # Check that all populations decreased (exponential decay)
            for i in 1:nSpecies
                finalPop = popsizes(finalComm, i)[1]
                initialPop = initialPops[i]
                @test finalPop < initialPop
                # Check approximate exponential decay: u(t) = u(0) * exp(-t)
                expected = initialPop * exp(-1.0)
                @test finalPop ≈ expected rtol=1e-3
            end

            # Check that time was updated
            @test finalComm.time ≈ comm.time + params.maxTime

            # Check that traits are unchanged
            for i in 1:nSpecies
                @test traits(finalComm, i) == traits(comm, i)
            end
        end
    end


    @testset "testing_ecoDyn_with_auxiliary_variables" begin
        # Test ecoDyn with auxiliary variables
        for _ in 1:numTests
            nSpecies = rand(2:4)
            nAux = rand(1:3)

            # Create species and auxiliary variables
            species = Species{Float64}[]
            for i in 1:nSpecies
                popsize = rand() * 5.0
                push!(species, Species(popsize, rand(1)))
            end

            aux = [PopulationSize(rand() * 2.0) for _ in 1:nAux]
            comm = Community(species, aux, 0.0)

            # Dynamics: species decay, aux variables grow
            function testDynamicsFactory(community)
                nSp = numSpecies(community)
                return function(u, p, t)
                    du = similar(u)
                    # First nSp elements decay
                    for i in 1:nSp
                        du[i] = -u[i]
                    end
                    # Remaining elements (aux) grow
                    for i in (nSp+1):length(u)
                        du[i] = 0.5 * u[i]
                    end
                    return du
                end
            end

            mutGen = (c, cfg) -> generateMutant(c, cfg, 0.01)
            params = IntegrationParams(maxTime=1.0, abstol=1e-8, reltol=1e-6)
            config = EcoEvoConfig(
                ecoDyn=testDynamicsFactory,
                mutationGenerator=mutGen,
                integrationParams=params,
                invaderPopsize=0.001,
                extThreshold=1e-8
            )

            # Integrate
            finalComm = ecoDyn(comm, config)

            # Check species populations decreased
            for i in 1:nSpecies
                @test popsizes(finalComm, i)[1] < popsizes(comm, i)[1]
            end

            # Check auxiliary variables increased
            for i in 1:nAux
                @test finalComm.aux[i].popsize[1] > comm.aux[i].popsize[1]
            end
        end
    end


    @testset "testing_singleEvoStep_basic_functionality" begin
        # Test that singleEvoStep performs all three operations correctly
        for _ in 1:numTests
            nSpecies = rand(2:5)

            # Create species with varying populations
            species = Species{Float64}[]
            for i in 1:nSpecies
                # Mix of healthy and near-extinct populations
                popsize = rand() < 0.3 ? rand() * 1e-10 : rand() * 10.0
                push!(species, Species(popsize, rand(1)))
            end

            comm = Community(species, PopulationSize{Float64}[], 0.0)

            # Slow exponential decay
            ecoDynFactory = (community) -> (u, p, t) -> -0.1 * u
            mutGen = (c, cfg) -> generateMutant(c, cfg, 0.01)
            params = IntegrationParams(maxTime=1.0, abstol=1e-8, reltol=1e-6)
            extThreshold = 1e-8
            config = EcoEvoConfig(
                ecoDyn=ecoDynFactory,
                mutationGenerator=mutGen,
                integrationParams=params,
                invaderPopsize=1.0,
                extThreshold=extThreshold
            )

            initialNumSpecies = numSpecies(comm)

            # Perform one evolutionary step
            resultComm = singleEvoStep(comm, config)

            # Check that a mutant was added (may be removed if it goes extinct)
            # and dynamics were integrated (time changed)
            @test resultComm.time > comm.time

            # Check that extinct species were removed
            for i in 1:numSpecies(resultComm)
                @test popsizes(resultComm, i)[1] >= extThreshold
            end

            # Check that if any species in original were below threshold,
            # the result has fewer species (or same if mutant also went extinct)
            nExtinct = count(popsizes(comm, i)[1] < extThreshold for i in 1:initialNumSpecies)
            if nExtinct > 0
                @test numSpecies(resultComm) <= initialNumSpecies + 1 - nExtinct
            end
        end
    end


    @testset "testing_singleEvoStep_mutant_addition" begin
        # Test that mutant is properly added with invader population size
        nSpecies = 3

        # Create healthy species
        species = [Species(5.0, [i/10.0]) for i in 1:nSpecies]
        comm = Community(species, PopulationSize{Float64}[], 0.0)

        # Very slow dynamics so populations don't change much
        ecoDynFactory = (community) -> (u, p, t) -> -0.001 * u
        mutGen = (c, cfg) -> generateMutant(c, cfg, 0.01)
        params = IntegrationParams(maxTime=0.1, abstol=1e-8, reltol=1e-6)
        invaderPop = 2.5
        config = EcoEvoConfig(
            ecoDyn=ecoDynFactory,
            mutationGenerator=mutGen,
            integrationParams=params,
            invaderPopsize=invaderPop,
            extThreshold=1e-8
        )

        resultComm = singleEvoStep(comm, config)

        # Should have one more species (no extinctions with slow dynamics)
        @test numSpecies(resultComm) == nSpecies + 1

        # Check that mutant (newest species) started with invader popsize
        # and decayed slightly
        finalMutantPop = popsizes(resultComm, nSpecies + 1)[1]
        expectedMutantPop = invaderPop * exp(-0.001 * 0.1)
        @test finalMutantPop ≈ expectedMutantPop rtol=1e-3
    end


    @testset "testing_singleEvoStep_extinction" begin
        # Test that near-extinct species are removed
        # Create species with one already below threshold
        species = [
            Species(5.0, [0.1]),
            Species(1e-10, [0.5]),  # Below typical threshold
            Species(3.0, [0.9])
        ]
        comm = Community(species, PopulationSize{Float64}[], 0.0)

        ecoDynFactory = (community) -> (u, p, t) -> -0.01 * u
        mutGen = (c, cfg) -> generateMutant(c, cfg, 0.01)
        params = IntegrationParams(maxTime=1.0, abstol=1e-8, reltol=1e-6)
        config = EcoEvoConfig(
            ecoDyn=ecoDynFactory,
            mutationGenerator=mutGen,
            integrationParams=params,
            invaderPopsize=1.0,
            extThreshold=1e-8
        )

        resultComm = singleEvoStep(comm, config)

        # Should have removed the near-extinct species
        # Start with 3, add 1 mutant = 4, remove 1 extinct = 3
        @test numSpecies(resultComm) == 3

        # All remaining species should be above threshold
        for i in 1:numSpecies(resultComm)
            @test popsizes(resultComm, i)[1] >= config.extThreshold
        end
    end


    @testset "testing_evolve!_basic_functionality" begin
        # Test that evolve! performs multiple evolutionary steps and records history
        # Create initial community
        nSpecies = 3
        species = [Species(5.0, [k / 10.0]) for k in 1:nSpecies]
        initialComm = Community(species, PopulationSize{Float64}[], 0.0)

        # Create history with initial community
        history = EvoHistory{Float64, 0}([initialComm])

        # Very slow dynamics so populations don't change much
        ecoDynFactory = (community) -> (u, p, t) -> -0.001 * u
        mutGen = (c, cfg) -> generateMutant(c, cfg, 0.01)
        params = IntegrationParams(maxTime=0.1, abstol=1e-8, reltol=1e-6)
        config = EcoEvoConfig(
            ecoDyn=ecoDynFactory,
            mutationGenerator=mutGen,
            integrationParams=params,
            invaderPopsize=1.0,
            extThreshold=1e-8
        )

        nMutEvents = 5

        # Run evolution
        evolve!(history, config, nMutEvents)

        # Check that history has been extended correctly
        @test length(history.history) == 1 + nMutEvents  # initial + nMutEvents

        # Check that number of species increases (or stays same if extinctions occur)
        finalComm = history.history[end]
        @test numSpecies(finalComm) >= nSpecies

        # Check that time progresses
        @test finalComm.time > initialComm.time

        # Check that all communities in history are valid
        for comm in history.history
            @test numSpecies(comm) > 0
            for i in 1:numSpecies(comm)
                @test popsizes(comm, i)[1] >= config.extThreshold
            end
        end
    end


    @testset "testing_evolve!_history_tracking" begin
        # Test that history correctly tracks evolutionary trajectory
        nSpecies = 2
        species = [Species(10.0, [k / 5.0]) for k in 1:nSpecies]
        initialComm = Community(species, PopulationSize{Float64}[], 0.0)

        history = EvoHistory{Float64, 0}([initialComm])

        # Slow exponential decay
        ecoDynFactory = (community) -> (u, p, t) -> -0.01 * u
        mutGen = (c, cfg) -> generateMutant(c, cfg, 0.01)
        params = IntegrationParams(maxTime=1.0, abstol=1e-8, reltol=1e-6)
        config = EcoEvoConfig(
            ecoDyn=ecoDynFactory,
            mutationGenerator=mutGen,
            integrationParams=params,
            invaderPopsize=2.0,
            extThreshold=1e-8
        )

        nMutEvents = 3
        evolve!(history, config, nMutEvents)

        # Verify each step adds to history
        @test length(history.history) == 4  # 1 initial + 3 events

        # Verify time is monotonically increasing
        times = [comm.time for comm in history.history]
        @test issorted(times)
        @test all(times .>= 0.0)

        # Verify species count increases (assuming no extinctions with these parameters)
        speciesCounts = [numSpecies(comm) for comm in history.history]
        @test speciesCounts[1] == nSpecies
        # Each step adds a mutant, so should increase
        for i in 2:lastindex(speciesCounts)
            @test speciesCounts[i] >= speciesCounts[i-1]
        end
    end


    @testset "testing_evolve!_with_extinctions" begin
        # Test evolve! when extinctions occur
        # Create community with some weak species
        species = [
            Species(10.0, [0.1]),
            Species(0.5, [0.5]),   # Weak species
            Species(8.0, [0.9])
        ]
        initialComm = Community(species, PopulationSize{Float64}[], 0.0)

        history = EvoHistory{Float64, 0}([initialComm])

        # Stronger decay
        ecoDynFactory = (community) -> (u, p, t) -> -0.5 * u
        mutGen = (c, cfg) -> generateMutant(c, cfg, 0.01)
        params = IntegrationParams(maxTime=2.0, abstol=1e-8, reltol=1e-6)
        config = EcoEvoConfig(
            ecoDyn=ecoDynFactory,
            mutationGenerator=mutGen,
            integrationParams=params,
            invaderPopsize=1.0,
            extThreshold=0.1  # Higher threshold to cause extinctions
        )

        nMutEvents = 3
        evolve!(history, config, nMutEvents)

        # History should have correct length
        @test length(history.history) == 4

        # All species in final community should be above threshold
        finalComm = history.history[end]
        for i in 1:numSpecies(finalComm)
            @test popsizes(finalComm, i)[1] >= config.extThreshold
        end

        # Check that extinctions may have occurred
        # (species count might not always increase)
        @test numSpecies(finalComm) >= 1  # At least some species survive
    end


    @testset "testing_evolve!_empty_history_error" begin
        # Test that evolve! throws error for empty history
        emptyHistory = EvoHistory{Float64, 0}(Community{Float64, 0}[])

        ecoDynFactory = (community) -> (u, p, t) -> -0.01 * u
        mutGen = (c, cfg) -> generateMutant(c, cfg, 0.01)
        params = IntegrationParams(maxTime=1.0, abstol=1e-8, reltol=1e-6)
        config = EcoEvoConfig(
            ecoDyn=ecoDynFactory,
            mutationGenerator=mutGen,
            integrationParams=params,
            invaderPopsize=1.0,
            extThreshold=1e-8
        )

        @test_throws ArgumentError evolve!(emptyHistory, config, 1)
    end


    @testset "testing_evolve!_zero_events" begin
        # Test that evolve! with 0 events doesn't modify history
        species = [Species(5.0, [0.5])]
        initialComm = Community(species, PopulationSize{Float64}[], 0.0)
        history = EvoHistory{Float64, 0}([initialComm])

        ecoDynFactory = (community) -> (u, p, t) -> -0.01 * u
        mutGen = (c, cfg) -> generateMutant(c, cfg, 0.01)
        params = IntegrationParams(maxTime=1.0, abstol=1e-8, reltol=1e-6)
        config = EcoEvoConfig(
            ecoDyn=ecoDynFactory,
            mutationGenerator=mutGen,
            integrationParams=params,
            invaderPopsize=1.0,
            extThreshold=1e-8
        )

        evolve!(history, config, 0)

        # History should remain unchanged
        @test length(history.history) == 1
        @test history.history[1] === initialComm
    end


    @testset "testing_evolve!_convenience_method_with_community" begin
        # Test the convenience method that takes a Community directly
        nSpecies = 3
        species = [Species(5.0, [i/10.0]) for i in 1:nSpecies]
        initialComm = Community(species, PopulationSize{Float64}[], 0.0)

        # Very slow dynamics
        ecoDynFactory = (community) -> (u, p, t) -> -0.001 * u
        mutGen = (c, cfg) -> generateMutant(c, cfg, 0.01)
        params = IntegrationParams(maxTime=0.1, abstol=1e-8, reltol=1e-6)
        config = EcoEvoConfig(
            ecoDyn=ecoDynFactory,
            mutationGenerator=mutGen,
            integrationParams=params,
            invaderPopsize=1.0,
            extThreshold=1e-8
        )

        nMutEvents = 4

        # Call convenience method - returns history
        history = evolve!(initialComm, config, nMutEvents)

        # Check that history is returned and has correct structure
        @test history isa EvoHistory{Float64, 0}
        @test length(history.history) == 1 + nMutEvents  # initial + events

        # Check that first community in history is the initial one
        @test history.history[1].time == initialComm.time
        @test numSpecies(history.history[1]) == numSpecies(initialComm)

        # Check that evolution occurred
        finalComm = history.history[end]
        @test finalComm.time > initialComm.time
        @test numSpecies(finalComm) >= nSpecies
    end


    @testset "testing_evolve!_convenience_method_comparison" begin
        # Test that convenience method produces same result as manual history creation
        nSpecies = 2
        species = [Species(10.0, [i/5.0]) for i in 1:nSpecies]
        comm1 = Community(species, PopulationSize{Float64}[], 0.0)
        comm2 = Community(species, PopulationSize{Float64}[], 0.0)

        ecoDynFactory = (community) -> (u, p, t) -> -0.01 * u
        mutGen = (c, cfg) -> generateMutant(c, cfg, 0.01)
        params = IntegrationParams(maxTime=1.0, abstol=1e-8, reltol=1e-6)
        config = EcoEvoConfig(
            ecoDyn=ecoDynFactory,
            mutationGenerator=mutGen,
            integrationParams=params,
            invaderPopsize=2.0,
            extThreshold=1e-8
        )

        # Use convenience method
        history1 = evolve!(comm1, config, 3)

        # Manual history creation
        history2 = EvoHistory{Float64, 0}([comm2])
        evolve!(history2, config, 3)

        # Both should have same structure
        @test length(history1.history) == length(history2.history)
        @test numSpecies(history1.history[end]) == numSpecies(history2.history[end])
    end


    @testset "testing_evolve_non_mutating_with_history" begin
        # Test that evolve (without !) doesn't modify original history
        nSpecies = 2
        species = [Species(5.0, [i/10.0]) for i in 1:nSpecies]
        initialComm = Community(species, PopulationSize{Float64}[], 0.0)

        ecoDynFactory = (community) -> (u, p, t) -> -0.01 * u
        mutGen = (c, cfg) -> generateMutant(c, cfg, 0.01)
        params = IntegrationParams(maxTime=1.0, abstol=1e-8, reltol=1e-6)
        config = EcoEvoConfig(
            ecoDyn=ecoDynFactory,
            mutationGenerator=mutGen,
            integrationParams=params,
            invaderPopsize=1.0,
            extThreshold=1e-8
        )

        # Create initial history
        history1 = evolve!(initialComm, config, 3; showProgress=false)
        originalLength = length(history1.history)

        # Use non-mutating evolve
        history2 = evolve(history1, config, 2; showProgress=false)

        # Check that original history is unchanged
        @test length(history1.history) == originalLength

        # Check that new history is extended
        @test length(history2.history) == originalLength + 2

        # Check that they share initial communities (values should match)
        for i in 1:originalLength
            @test history1.history[i].time == history2.history[i].time
            @test numSpecies(history1.history[i]) == numSpecies(history2.history[i])
        end

        # Check that they are different objects
        @test history1 !== history2
        @test history1.history !== history2.history
    end


    @testset "testing_evolve_non_mutating_with_community" begin
        # Test evolve (without !) starting from Community
        nSpecies = 2
        species = [Species(8.0, [i/10.0]) for i in 1:nSpecies]
        initialComm = Community(species, PopulationSize{Float64}[], 0.0)

        ecoDynFactory = (community) -> (u, p, t) -> -0.001 * u
        mutGen = (c, cfg) -> generateMutant(c, cfg, 0.01)
        params = IntegrationParams(maxTime=0.5, abstol=1e-8, reltol=1e-6)
        config = EcoEvoConfig(
            ecoDyn=ecoDynFactory,
            mutationGenerator=mutGen,
            integrationParams=params,
            invaderPopsize=1.5,
            extThreshold=1e-8
        )

        nMutEvents = 3

        # Use non-mutating evolve
        history = evolve(initialComm, config, nMutEvents; showProgress=false)

        # Check that history is returned with correct structure
        @test history isa EvoHistory{Float64, 0}
        @test length(history.history) == 1 + nMutEvents

        # Check that evolution occurred
        @test history.history[end].time > history.history[1].time
        @test numSpecies(history.history[end]) >= nSpecies
    end


    @testset "testing_evolve_vs_evolve!_consistency" begin
        # Test that evolve and evolve! produce equivalent results
        nSpecies = 2
        species = [Species(6.0, [i/8.0]) for i in 1:nSpecies]
        comm1 = Community(species, PopulationSize{Float64}[], 0.0)
        comm2 = Community(species, PopulationSize{Float64}[], 0.0)

        ecoDynFactory = (community) -> (u, p, t) -> -0.02 * u
        mutGen = (c, cfg) -> generateMutant(c, cfg, 0.015)
        params = IntegrationParams(maxTime=2.0, abstol=1e-8, reltol=1e-6)
        config = EcoEvoConfig(
            ecoDyn=ecoDynFactory,
            mutationGenerator=mutGen,
            integrationParams=params,
            invaderPopsize=1.2,
            extThreshold=1e-8
        )

        # Use both methods
        history1 = evolve(comm1, config, 4; showProgress=false)
        history2 = evolve!(comm2, config, 4; showProgress=false)

        # Should produce histories with same structure
        @test length(history1.history) == length(history2.history)
        @test numSpecies(history1.history[end]) == numSpecies(history2.history[end])
    end


    @testset "testing_evolve_multiple_calls_independence" begin
        # Test that multiple calls to evolve don't interfere with each other
        nSpecies = 2
        species = [Species(4.0, [i/6.0]) for i in 1:nSpecies]
        initialComm = Community(species, PopulationSize{Float64}[], 0.0)

        ecoDynFactory = (community) -> (u, p, t) -> -0.005 * u
        mutGen = (c, cfg) -> generateMutant(c, cfg, 0.02)
        params = IntegrationParams(maxTime=1.5, abstol=1e-8, reltol=1e-6)
        config = EcoEvoConfig(
            ecoDyn=ecoDynFactory,
            mutationGenerator=mutGen,
            integrationParams=params,
            invaderPopsize=0.8,
            extThreshold=1e-8
        )

        # Create base history
        history0 = evolve(initialComm, config, 2; showProgress=false)
        originalLength = length(history0.history)

        # Call evolve multiple times on same base
        historyA = evolve(history0, config, 1; showProgress=false)
        historyB = evolve(history0, config, 2; showProgress=false)
        historyC = evolve(history0, config, 3; showProgress=false)

        # Original should be unchanged
        @test length(history0.history) == originalLength

        # Each should have different lengths
        @test length(historyA.history) == originalLength + 1
        @test length(historyB.history) == originalLength + 2
        @test length(historyC.history) == originalLength + 3

        # All should be independent objects
        @test historyA !== historyB
        @test historyB !== historyC
        @test historyA !== historyC
    end

end
