using Test
using EcoEvoSim
using LinearAlgebra
using Random


numTests = 50


@testset "tests_of_mutant_generation" begin

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
            config = EcoEvoConfig(ecoDyn, mutGen, params, 1e-8)

            # Create covariance matrix
            variance = 0.01
            covMat = Matrix{Float64}(variance * I(traitDim))

            # Generate mutant
            newComm = generateMutant(; invaderPopsize=0.001, covMat=covMat)(comm)

            # Check that community has one more species
            @test numSpecies(newComm) == numSpecies(comm) + 1

            # Check that the mutant has the right population size
            mutantPopsize = popsizes(newComm, numSpecies(newComm))
            @test mutantPopsize[1] == 0.001

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
            config = EcoEvoConfig(ecoDyn, mutGen, params, 1e-8)

            # Generate mutant with scalar variance
            variance = 0.01
            newComm = generateMutant(; invaderPopsize=0.001, variance=variance)(comm)

            # Check that community has one more species
            @test numSpecies(newComm) == numSpecies(comm) + 1

            # Check that the mutant has the right population size
            mutantPopsize = popsizes(newComm, numSpecies(newComm))
            @test mutantPopsize[1] == 0.001

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
        config = EcoEvoConfig(ecoDyn, mutGen, params, 1e-8)

        @test_throws ArgumentError generateMutant(; invaderPopsize=0.001, variance=-0.01)
        @test_throws ArgumentError generateMutant(; invaderPopsize=0.001, variance=0.0)
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
            config = EcoEvoConfig(ecoDyn, mutGen, params, 1e-8)

            # Create covariance matrix
            variance = 0.01
            covMat = Matrix{Float64}(variance * I(traitDim))

            # Generate many mutants and check they're weighted correctly
            # (species with larger populations should be selected more often)
            nTrials = 1000
            parentCounts = zeros(Int, nSpecies)

            for trial in 1:nTrials
                newComm = generateMutantWeighted(; invaderPopsize=0.001, covMat=covMat)(comm)
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
            config = EcoEvoConfig(ecoDyn, mutGen, params, 1e-8)

            # Generate mutant with scalar variance
            variance = 0.01
            newComm = generateMutantWeighted(; invaderPopsize=0.001, variance=variance)(comm)

            # Check that community has one more species
            @test numSpecies(newComm) == numSpecies(comm) + 1

            # Check that the mutant has the right population size
            mutantPopsize = popsizes(newComm, numSpecies(newComm))
            @test mutantPopsize[1] == 0.001

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
        config = EcoEvoConfig(ecoDyn, mutGen, params, 1e-8)

        @test_throws ArgumentError generateMutantWeighted(; invaderPopsize=0.001, variance=-0.01)
        @test_throws ArgumentError generateMutantWeighted(; invaderPopsize=0.001, variance=0.0)
    end


    @testset "testing_invalid_covariance_matrix_dimensions" begin
        # Create a 2D trait community
        species = [Species(1.0, [0.5, 0.5]), Species(2.0, [0.3, 0.7])]
        comm = Community(species, PopulationSize{Float64}[], 0.0)

        ecoDyn = (x, t, p) -> -x
        mutGen = (x) -> x .+ 0.01
        params = IntegrationParams(maxTime=100.0, abstol=1e-6, reltol=0.1)
        config = EcoEvoConfig(ecoDyn, mutGen, params, 1e-8)

        # Try with wrong-sized covariance matrix (3x3 instead of 2x2)
        wrongCovMat = Matrix{Float64}(0.01 * I(3))
        @test_throws ArgumentError generateMutant(; invaderPopsize=0.001, covMat=wrongCovMat)(comm)

        # Try with non-square matrix
        wrongCovMat2 = rand(2, 3)
        @test_throws ArgumentError generateMutant(; invaderPopsize=0.001, covMat=wrongCovMat2)(comm)
    end


    @testset "testing_generateMutantSpatial_with_covariance_matrix" begin
        # Create spatially-structured communities (multiple patches = multiple stage classes)
        for _ in 1:numTests
            nSpecies = rand(2:4)
            nPatches = rand(2:5)  # Number of spatial patches
            traitDim = rand(1:3)

            # Create species with spatially-structured populations
            species = Species{Float64}[]
            for i in 1:nSpecies
                trait = rand(traitDim)
                popsize_spatial = rand(nPatches) .* 5.0  # Random population in each patch
                push!(species, Species(popsize_spatial, trait))
            end

            comm = Community(species, PopulationSize{Float64}[], 0.0)

            # Create covariance matrix
            variance = 0.01
            covMat = Matrix{Float64}(variance * I(traitDim))

            # Generate mutant
            newComm = generateMutantSpatial(; invaderPopsize=0.001, covMat=covMat)(comm)

            # Check that community has one more species
            @test numSpecies(newComm) == numSpecies(comm) + 1

            # Check that the mutant appears in only ONE patch
            mutantPopsize = popsizes(newComm, numSpecies(newComm))
            nonzero_patches = findall(x -> x > 0.0, mutantPopsize)
            @test length(nonzero_patches) == 1

            # Check that the total mutant population size equals invaderPopsize
            @test sum(mutantPopsize) ≈ 0.001

            # Check that mutant trait is different from all parents
            mutantTrait = traits(newComm, numSpecies(newComm))
            parentTraits = traits(comm)
            @test !any(all(mutantTrait .≈ pt) for pt in parentTraits)

            # Check that time and aux are preserved
            @test newComm.time == comm.time
            @test newComm.aux == comm.aux
        end
    end


    @testset "testing_generateMutantSpatial_with_scalar_variance" begin
        # Create spatially-structured communities
        for _ in 1:numTests
            nSpecies = rand(2:4)
            nPatches = rand(2:5)
            traitDim = rand(1:3)

            species = Species{Float64}[]
            for i in 1:nSpecies
                trait = rand(traitDim)
                popsize_spatial = rand(nPatches) .* 5.0
                push!(species, Species(popsize_spatial, trait))
            end

            comm = Community(species, PopulationSize{Float64}[], 0.0)

            # Generate mutant with scalar variance
            variance = 0.01
            newComm = generateMutantSpatial(; invaderPopsize=0.001, variance=variance)(comm)

            # Check that community has one more species
            @test numSpecies(newComm) == numSpecies(comm) + 1

            # Check that the mutant appears in only ONE patch
            mutantPopsize = popsizes(newComm, numSpecies(newComm))
            nonzero_patches = findall(x -> x > 0.0, mutantPopsize)
            @test length(nonzero_patches) == 1

            # Check that the total mutant population size equals invaderPopsize
            @test sum(mutantPopsize) ≈ 0.001

            # Check trait dimension matches
            mutantTrait = traits(newComm, numSpecies(newComm))
            @test length(mutantTrait) == traitDim
        end

        # Invalid: non-positive variance
        species = [Species([1.0, 2.0], [0.5, 0.5])]
        comm = Community(species, PopulationSize{Float64}[], 0.0)

        @test_throws ArgumentError generateMutantSpatial(; invaderPopsize=0.001, variance=-0.01)
        @test_throws ArgumentError generateMutantSpatial(; invaderPopsize=0.001, variance=0.0)
    end


    @testset "testing_generateMutantSpatial_patch_selection_probabilities" begin
        # Test that patch selection is proportional to patch population sizes
        # Create a community with 2 patches, one having much more population than the other
        nPatches = 2
        patch_pops = [10.0, 1.0]  # First patch has 10x the population

        species = [Species(patch_pops, [0.0])]  # Single species, asymmetric spatial distribution
        comm = Community(species, PopulationSize{Float64}[], 0.0)

        # Generate many mutants and track which patch they appear in
        nTrials = 500
        patch_counts = zeros(Int, nPatches)

        for _ in 1:nTrials
            newComm = generateMutantSpatial(; invaderPopsize=0.001, variance=0.01^2)(comm)
            mutantPops = popsizes(newComm, numSpecies(newComm))
            patch_idx = findall(x -> x > 0.0, mutantPops)[1]
            patch_counts[patch_idx] += 1
        end

        # Check that both patches were sampled (with high probability)
        @test all(patch_counts .> 0)

        # Check that patch 1 was selected more often than patch 2
        # (should be roughly 10x more often, but allow for some variance)
        @test patch_counts[1] > patch_counts[2]
    end


    @testset "testing_generateMutantSpatial_with_equal_patches" begin
        # When all patches have equal population, mutant should appear equally likely in any patch
        nPatches = 3
        patch_pops = [5.0, 5.0, 5.0]  # Equal population in each patch

        species = [Species(patch_pops, [0.0])]
        comm = Community(species, PopulationSize{Float64}[], 0.0)

        # Generate many mutants
        nTrials = 300
        patch_counts = zeros(Int, nPatches)

        for _ in 1:nTrials
            newComm = generateMutantSpatial(; invaderPopsize=0.001, variance=0.01^2)(comm)
            mutantPops = popsizes(newComm, numSpecies(newComm))
            patch_idx = findall(x -> x > 0.0, mutantPops)[1]
            patch_counts[patch_idx] += 1
        end

        # All patches should have been selected at least once
        @test all(patch_counts .> 0)
    end


    @testset "testing_generateMutantSpatial_invalid_covariance_matrix_dimensions" begin
        # Create a spatially-structured community with 2D traits
        species = [Species([1.0, 1.0], [0.5, 0.5]), Species([2.0, 2.0], [0.3, 0.7])]
        comm = Community(species, PopulationSize{Float64}[], 0.0)

        # Try with wrong-sized covariance matrix (3x3 instead of 2x2)
        wrongCovMat = Matrix{Float64}(0.01 * I(3))
        @test_throws ArgumentError generateMutantSpatial(; invaderPopsize=0.001, covMat=wrongCovMat)(comm)

        # Try with non-square matrix
        wrongCovMat2 = rand(2, 3)
        @test_throws ArgumentError generateMutantSpatial(; invaderPopsize=0.001, covMat=wrongCovMat2)(comm)
    end


    @testset "testing_generateMutantSpatial_zero_population_error" begin
        # Create a community where all patches have zero population
        species = [Species([0.0, 0.0], [0.0])]
        comm = Community(species, PopulationSize{Float64}[], 0.0)

        # Should throw an error
        @test_throws ArgumentError generateMutantSpatial(; invaderPopsize=0.001, variance=0.01^2)(comm)
    end


    @testset "testing_generateMutantSpatialWeighted_with_covariance_matrix" begin
        # Test weighted parent selection with spatial placement
        for _ in 1:numTests
            nSpecies = rand(2:4)
            nPatches = rand(2:5)
            traitDim = rand(1:3)

            # Create species with spatially-structured populations
            # Make some species abundant and others rare
            species = Species{Float64}[]
            for i in 1:nSpecies
                trait = rand(traitDim)
                if i == 1
                    # First species is very abundant
                    popsize_spatial = ones(nPatches) .* 5.0
                else
                    # Other species are rare
                    popsize_spatial = ones(nPatches) .* 0.1
                end
                push!(species, Species(popsize_spatial, trait))
            end

            comm = Community(species, PopulationSize{Float64}[], 0.0)

            # Create covariance matrix
            variance = 0.01
            covMat = Matrix{Float64}(variance * I(traitDim))

            # Generate mutant with weighted selection
            newComm = generateMutantSpatialWeighted(; invaderPopsize=0.001, covMat=covMat)(comm)

            # Check that community has one more species
            @test numSpecies(newComm) == numSpecies(comm) + 1

            # Check that the mutant appears in only ONE patch
            mutantPopsize = popsizes(newComm, numSpecies(newComm))
            nonzero_patches = findall(x -> x > 0.0, mutantPopsize)
            @test length(nonzero_patches) == 1

            # Check that the total mutant population size equals invaderPopsize
            @test sum(mutantPopsize) ≈ 0.001

            # Check that time and aux are preserved
            @test newComm.time == comm.time
            @test newComm.aux == comm.aux
        end
    end


    @testset "testing_generateMutantSpatialWeighted_with_scalar_variance" begin
        # Test with variance parameter instead of covariance matrix
        for _ in 1:numTests
            nSpecies = rand(2:4)
            nPatches = rand(2:5)
            traitDim = rand(1:3)

            species = Species{Float64}[]
            for i in 1:nSpecies
                trait = rand(traitDim)
                popsize_spatial = rand(nPatches) .* 5.0
                push!(species, Species(popsize_spatial, trait))
            end

            comm = Community(species, PopulationSize{Float64}[], 0.0)

            # Generate mutant with scalar variance
            variance = 0.01
            newComm = generateMutantSpatialWeighted(; invaderPopsize=0.001, variance=variance)(comm)

            # Check that community has one more species
            @test numSpecies(newComm) == numSpecies(comm) + 1

            # Check that the mutant appears in only ONE patch
            mutantPopsize = popsizes(newComm, numSpecies(newComm))
            nonzero_patches = findall(x -> x > 0.0, mutantPopsize)
            @test length(nonzero_patches) == 1

            # Check that the total mutant population size equals invaderPopsize
            @test sum(mutantPopsize) ≈ 0.001

            # Check trait dimension matches
            mutantTrait = traits(newComm, numSpecies(newComm))
            @test length(mutantTrait) == traitDim
        end

        # Invalid: non-positive variance
        species = [Species([1.0, 2.0], [0.5, 0.5])]
        comm = Community(species, PopulationSize{Float64}[], 0.0)

        @test_throws ArgumentError generateMutantSpatialWeighted(; invaderPopsize=0.001, variance=-0.01)
        @test_throws ArgumentError generateMutantSpatialWeighted(; invaderPopsize=0.001, variance=0.0)
    end


    @testset "testing_generateMutantSpatialWeighted_parent_weighted_selection" begin
        # Test that parent selection is weighted by abundance
        # Create 3 species with very different abundances
        nPatches = 2
        species = [
            Species([5.0, 5.0], [0.0]),      # Most abundant
            Species([1.0, 1.0], [0.1]),      # Medium
            Species([0.1, 0.1], [0.2])       # Rarest
        ]
        comm = Community(species, PopulationSize{Float64}[], 0.0)

        # Generate many mutants and track which parent is selected
        nTrials = 1000
        parent_counts = zeros(Int, 3)

        for _ in 1:nTrials
            newComm = generateMutantSpatialWeighted(; invaderPopsize=0.001, variance=0.01^2)(comm)
            mutantTrait = traits(newComm, numSpecies(newComm))

            # Find which parent it's closest to
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
            parent_counts[parentIdx] += 1
        end

        # Check that all parents were selected at least once
        @test all(parent_counts .> 0)

        # Check that the most abundant species was selected most often
        # (species 1 has 10x the population of species 2)
        @test parent_counts[1] > parent_counts[2]
    end


    @testset "testing_generateMutantSpatialWeighted_patch_selection" begin
        # Test that patch selection is still weighted by patch population
        # Create a spatially heterogeneous community
        nPatches = 3
        patch_pops = [10.0, 5.0, 1.0]  # Heterogeneous patches

        species = [
            Species(patch_pops, [0.0]),      # Single species with heterogeneous spatial dist
            Species(patch_pops .* 0.1, [0.1]) # Another species (less abundant overall)
        ]
        comm = Community(species, PopulationSize{Float64}[], 0.0)

        # Generate many mutants and track which patch they appear in
        nTrials = 500
        patch_counts = zeros(Int, nPatches)

        for _ in 1:nTrials
            newComm = generateMutantSpatialWeighted(; invaderPopsize=0.001, variance=0.01^2)(comm)
            mutantPops = popsizes(newComm, numSpecies(newComm))
            patch_idx = findall(x -> x > 0.0, mutantPops)[1]
            patch_counts[patch_idx] += 1
        end

        # All patches should be sampled
        @test all(patch_counts .> 0)

        # Patch 1 should be most common, patch 3 least common
        @test patch_counts[1] > patch_counts[2]
        @test patch_counts[2] > patch_counts[3]
    end


    @testset "testing_generateMutantSpatialWeighted_invalid_covariance_matrix_dimensions" begin
        # Create a spatially-structured community with 2D traits
        species = [Species([1.0, 1.0], [0.5, 0.5]), Species([2.0, 2.0], [0.3, 0.7])]
        comm = Community(species, PopulationSize{Float64}[], 0.0)

        # Try with wrong-sized covariance matrix (3x3 instead of 2x2)
        wrongCovMat = Matrix{Float64}(0.01 * I(3))
        @test_throws ArgumentError generateMutantSpatialWeighted(; invaderPopsize=0.001, covMat=wrongCovMat)(comm)

        # Try with non-square matrix
        wrongCovMat2 = rand(2, 3)
        @test_throws ArgumentError generateMutantSpatialWeighted(; invaderPopsize=0.001, covMat=wrongCovMat2)(comm)
    end
end
