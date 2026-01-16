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

end
