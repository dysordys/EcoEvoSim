using Test
using EcoEvoSim


numTests = 50


@testset "tests_of_basic_selectors" begin

    @testset "testing_species_selector_by_index" begin
        for _ in 1:numTests
            T = Float64
            auxClasses = rand(0:3)
            numSpecies = rand(1:5)
            sps = Species{T}[]
            for _ in 1:numSpecies
                ps = PopulationSize(rand(T, rand(1:3)))
                ph = Phenotype(rand(T, rand(1:3)))
                push!(sps, Species(ps, ph))
            end
            aux = Vector{PopulationSize{T}}(undef, auxClasses)
            for i in 1:auxClasses
                aux[i] = PopulationSize(rand(T))
            end
            comm = Community(sps, aux)
            for i in 1:numSpecies
                @test speciesList(comm, i) === sps[i]
            end
        end
    end


    @testset "testing_species_selector_by_an_index_that_is_out_of_bounds" begin
        T = Float64
        stageClasses = 2
        traitDim = 2
        auxClasses = 1
        ps = PopulationSize(rand(T, rand(1:3)))
        ph = Phenotype(rand(T, rand(1:3)))
        sp = Species(ps, ph)
        aux = [PopulationSize(rand(T))]
        comm = Community([sp], aux)
        @test_throws ArgumentError speciesList(comm, 0)
        @test_throws ArgumentError speciesList(comm, 2)
    end


    @testset "testing_species_selector_by_collection_of_indices" begin
        for _ in 1:numTests
            T = Float64
            auxClasses = rand(0:3)
            numSpecies = rand(2:5)
            sps = Species{T}[]
            for _ in 1:numSpecies
                ps = PopulationSize(rand(T, rand(1:3)))
                ph = Phenotype(rand(T, rand(1:3)))
                push!(sps, Species(ps, ph))
            end
            aux = Vector{PopulationSize{T}}(undef, auxClasses)
            for i in 1:auxClasses
                aux[i] = PopulationSize(rand(T))
            end
            comm = Community(sps, aux)

            # Test with empty vector
            result = speciesList(comm, [])
            @test result == Species{T}[]
            @test length(result) == 0
            @test eltype(result) == Species{T}

            # Test with vector of indices (only if we have enough species)
            if numSpecies >= 3
                indices = [1, 3]
                result = speciesList(comm, indices)
                @test result == [sps[1], sps[3]]
            end

            # Test with range
            result = speciesList(comm, 1:2)
            @test result == [sps[1], sps[2]]

            # Test with all indices
            result = speciesList(comm, 1:numSpecies)
            @test result == sps
        end
    end


    @testset "testing_species_selector_with_collection_out_of_bounds" begin
        T = Float64
        stageClasses = 2
        traitDim = 2
        auxClasses = 1
        ps = PopulationSize(rand(T, rand(1:3)))
        ph = Phenotype(rand(T, rand(1:3)))
        sp = Species(ps, ph)
        aux = [PopulationSize(rand(T))]
        comm = Community([sp], aux)

        # Test with out of bounds index in collection
        @test_throws ArgumentError speciesList(comm, [0, 1])
        @test_throws ArgumentError speciesList(comm, [1, 2])
        @test_throws ArgumentError speciesList(comm, [2])
    end


    @testset "testing_popsizes_extractor" begin
        for _ in 1:numTests
            T = Float64
            auxClasses = rand(0:3)
            numSpecies = rand(1:5)
            sps = Species{T}[]
            for _ in 1:numSpecies
                ps = PopulationSize(rand(T, rand(1:3)))
                ph = Phenotype(rand(T, rand(1:3)))
                push!(sps, Species(ps, ph))
            end
            aux = [PopulationSize(rand(T)) for _ in 1:auxClasses]
            comm = Community(sps, aux)

            # Test getting all popsizes
            all_pops = popsizes(comm)
            @test length(all_pops) == numSpecies
            for i in 1:numSpecies
                @test all_pops[i] === sps[i].popsize.popsize
            end

            # Test single index
            for i in 1:numSpecies
                @test popsizes(comm, i) === sps[i].popsize.popsize
            end

            # Test collection of indices
            if numSpecies >= 2
                indices = [1, 2]
                result = popsizes(comm, indices)
                @test result == [sps[1].popsize.popsize, sps[2].popsize.popsize]
            end
        end
    end


    @testset "testing_popsizesToMatrix" begin
        for _ in 1:numTests
            T = Float64
            auxClasses = rand(0:3)
            numSpecies = rand(1:5)
            nstage = rand(1:4)
            sps = Species{T}[]
            for _ in 1:numSpecies
                ps = PopulationSize(rand(T, nstage))
                ph = Phenotype(rand(T, rand(1:3)))
                push!(sps, Species(ps, ph))
            end
            aux = [PopulationSize(rand(T)) for _ in 1:auxClasses]
            comm = Community(sps, aux)

            M = popsizesToMatrix(comm)
            @test size(M) == (numSpecies, nstage)
            for i in 1:numSpecies, j in 1:nstage
                @test M[i, j] == sps[i].popsize.popsize[j]
            end
        end

        # Test that differing stage-class lengths throw an error
        T = Float64
        ps1 = PopulationSize(rand(T, 1))
        ps2 = PopulationSize(rand(T, 2))
        sp1 = Species(ps1, Phenotype(rand(T,1)))
        sp2 = Species(ps2, Phenotype(rand(T,1)))
        comm_mismatch = Community([sp1, sp2], PopulationSize{T}[])
        @test_throws ArgumentError popsizesToMatrix(comm_mismatch)
    end


    @testset "testing_traits_extractor" begin
        for _ in 1:numTests
            T = Float64
            auxClasses = rand(0:3)
            numSpecies = rand(1:5)
            sps = Species{T}[]
            for _ in 1:numSpecies
                ps = PopulationSize(rand(T, rand(1:3)))
                ph = Phenotype(rand(T, rand(1:3)))
                push!(sps, Species(ps, ph))
            end
            aux = [PopulationSize(rand(T)) for _ in 1:auxClasses]
            comm = Community(sps, aux)

            # Test getting all traits
            all_traits = traits(comm)
            @test length(all_traits) == numSpecies
            for i in 1:numSpecies
                @test all_traits[i] === sps[i].trait.trait
            end

            # Test single index
            for i in 1:numSpecies
                @test traits(comm, i) === sps[i].trait.trait
            end

            # Test collection of indices
            if numSpecies >= 2
                indices = [1, 2]
                result = traits(comm, indices)
                @test result == [sps[1].trait.trait, sps[2].trait.trait]
            end
        end
    end


    @testset "testing_aux_selector_by_index" begin
        for _ in 1:numTests
            T = Float64
            auxClasses = rand(1:5)
            numSpecies = rand(1:5)
            sps = Species{T}[]
            for _ in 1:numSpecies
                ps = PopulationSize(rand(T, rand(1:3)))
                ph = Phenotype(rand(T, rand(1:3)))
                push!(sps, Species(ps, ph))
            end
            aux = [PopulationSize(rand(T)) for _ in 1:auxClasses]
            comm = Community(sps, aux)
            for i in 1:auxClasses
                @test auxs(comm, i) === aux[i]
            end
        end
    end


    @testset "testing_aux_selector_by_an_index_that_is_out_of_bounds" begin
        T = Float64
        stageClasses = 2
        traitDim = 2
        numSpecies = 1
        ps = PopulationSize(rand(T, rand(1:3)))
        ph = Phenotype(rand(T, rand(1:3)))
        sp = Species(ps, ph)
        aux = [PopulationSize(rand(T))]
        comm = Community([sp], aux)
        @test_throws ArgumentError auxs(comm, 0)
        @test_throws ArgumentError auxs(comm, 2)
    end


    @testset "testing_aux_selector_by_collection_of_indices" begin
        for _ in 1:numTests
            T = Float64
            auxClasses = rand(2:5)
            numSpecies = rand(1:5)
            sps = Species{T}[]
            for _ in 1:numSpecies
                ps = PopulationSize(rand(T, rand(1:3)))
                ph = Phenotype(rand(T, rand(1:3)))
                push!(sps, Species(ps, ph))
            end
            aux = [PopulationSize(rand(T)) for _ in 1:auxClasses]
            comm = Community(sps, aux)

            # Test with empty vector
            result = auxs(comm, [])
            @test result == PopulationSize{T}[]
            @test length(result) == 0
            @test eltype(result) == PopulationSize{T}

            # Test with vector of indices (only if we have enough aux)
            if auxClasses >= 3
                indices = [1, 3]
                result = auxs(comm, indices)
                @test result == [aux[1], aux[3]]
            end

            # Test with range
            result = auxs(comm, 1:2)
            @test result == [aux[1], aux[2]]

            # Test with all indices
            result = auxs(comm, 1:auxClasses)
            @test result == aux
        end
    end


    @testset "testing_aux_selector_with_collection_out_of_bounds" begin
        T = Float64
        stageClasses = 2
        traitDim = 2
        numSpecies = 1
        ps = PopulationSize(rand(T, rand(1:3)))
        ph = Phenotype(rand(T, rand(1:3)))
        sp = Species(ps, ph)
        aux = [PopulationSize(rand(T))]
        comm = Community([sp], aux)

        # Test with out of bounds index in collection
        @test_throws ArgumentError auxs(comm, [0, 1])
        @test_throws ArgumentError auxs(comm, [1, 2])
        @test_throws ArgumentError auxs(comm, [2])
    end


    @testset "testing_numSpecies" begin
        for _ in 1:numTests
            T = Float64
            auxClasses = rand(0:3)
            numSpecies_expected = rand(1:5)
            sps = Species{T}[]
            for _ in 1:numSpecies_expected
                ps = PopulationSize(rand(T, rand(1:3)))
                ph = Phenotype(rand(T, rand(1:3)))
                push!(sps, Species(ps, ph))
            end
            aux = [PopulationSize(rand(T)) for _ in 1:auxClasses]
            comm = Community(sps, aux)

            @test numSpecies(comm) == numSpecies_expected
            @test numSpecies(comm) == length(sps)
        end
    end


    @testset "testing_speciesIndices" begin
        for _ in 1:numTests
            T = Float64
            auxClasses = rand(0:3)
            numSpecies_expected = rand(1:5)
            sps = Species{T}[]
            for _ in 1:numSpecies_expected
                ps = PopulationSize(rand(T, rand(1:3)))
                ph = Phenotype(rand(T, rand(1:3)))
                push!(sps, Species(ps, ph))
            end
            aux = [PopulationSize(rand(T)) for _ in 1:auxClasses]
            comm = Community(sps, aux)

            indices = speciesIndices(comm)
            @test indices == 1:numSpecies_expected
            @test collect(indices) == collect(1:numSpecies_expected)
        end
    end

    @testset "testing_numStages" begin
        for _ in 1:numTests
            T = Float64
            numSpecies = rand(1:5)
            numStages_expected = rand(1:4)
            auxClasses = rand(0:3)

            # Create species with the same number of stages
            sps = Species{T}[]
            for _ in 1:numSpecies
                ps = PopulationSize(rand(T, numStages_expected))
                ph = Phenotype(rand(T, rand(1:3)))
                push!(sps, Species(ps, ph))
            end
            aux = [PopulationSize(rand(T)) for _ in 1:auxClasses]
            comm = Community(sps, aux)

            # Test that numStages returns the correct value
            @test numStages(comm) == numStages_expected
        end
    end

    @testset "testing_numStages_with_empty_community" begin
        T = Float64
        comm = Community(Species{T}[], PopulationSize{T}[])
        @test_throws ArgumentError numStages(comm)
    end

end
