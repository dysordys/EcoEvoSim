using Test
using EcoEvoSim
using StaticArrays


numTests = 50


@testset "tests_of_basic_utils" begin

    @testset "testing_addSpecies_single_species" begin
        for _ in 1:numTests
            T = Float64
            stageClasses = rand(1:3)
            traitDim = rand(1:3)
            auxClasses = rand(0:3)
            numSpeciesInit = rand(1:5)
            sps = Species{T, stageClasses, traitDim}[]
            for _ in 1:numSpeciesInit
                ps = PopulationSize(rand(T, stageClasses))
                ph = Phenotype(rand(T, traitDim))
                push!(sps, Species(ps, ph))
            end
            aux = [PopulationSize(rand(T)) for _ in 1:auxClasses]
            comm = Community(sps, aux)

            # Create a new species to add
            new_ps = PopulationSize(rand(T, stageClasses))
            new_ph = Phenotype(rand(T, traitDim))
            new_sp = Species(new_ps, new_ph)

            # Test addSpecies with different argument orders
            comm_new1 = addSpecies(comm, new_sp)
            comm_new2 = addSpecies(new_sp, comm)

            @test numSpecies(comm_new1) == numSpeciesInit + 1
            @test numSpecies(comm_new2) == numSpeciesInit + 1
            @test speciesList(comm_new1) == [sps; new_sp]
            @test speciesList(comm_new2) == [sps; new_sp]
        end
    end


    @testset "testing_addSpecies_multiple_species" begin
        for _ in 1:numTests
            T = Float64
            stageClasses = rand(1:3)
            traitDim = rand(1:3)
            auxClasses = rand(0:3)
            numSpeciesInit = rand(1:5)
            sps = Species{T, stageClasses, traitDim}[]
            for _ in 1:numSpeciesInit
                ps = PopulationSize(rand(T, stageClasses))
                ph = Phenotype(rand(T, traitDim))
                push!(sps, Species(ps, ph))
            end
            aux = [PopulationSize(rand(T)) for _ in 1:auxClasses]
            comm = Community(sps, aux)

            # Create multiple new species to add
            numNew = rand(1:3)
            new_sps = Species{T, stageClasses, traitDim}[]
            for _ in 1:numNew
                ps = PopulationSize(rand(T, stageClasses))
                ph = Phenotype(rand(T, traitDim))
                push!(new_sps, Species(ps, ph))
            end

            # Test addSpecies with different argument orders
            comm_new1 = addSpecies(comm, new_sps)
            comm_new2 = addSpecies(new_sps, comm)

            @test numSpecies(comm_new1) == numSpeciesInit + numNew
            @test numSpecies(comm_new2) == numSpeciesInit + numNew
            @test speciesList(comm_new1) == [sps; new_sps]
            @test speciesList(comm_new2) == [sps; new_sps]
        end
    end


    @testset "testing_addSpecies_preserves_auxiliaries_and_time" begin
        for _ in 1:numTests
            T = Float64
            stageClasses = rand(1:3)
            traitDim = rand(1:3)
            auxClasses = rand(1:3)
            numSpecies = rand(1:5)
            sps = Species{T, stageClasses, traitDim}[]
            for _ in 1:numSpecies
                ps = PopulationSize(rand(T, stageClasses))
                ph = Phenotype(rand(T, traitDim))
                push!(sps, Species(ps, ph))
            end
            aux = [PopulationSize(rand(T)) for _ in 1:auxClasses]
            time = rand(T)
            comm = Community(sps, aux, time)

            # Create a new species to add
            new_ps = PopulationSize(rand(T, stageClasses))
            new_ph = Phenotype(rand(T, traitDim))
            new_sp = Species(new_ps, new_ph)

            comm_new = addSpecies(comm, new_sp)

            # Check that auxiliaries and time are preserved
            @test auxs(comm_new) == aux
            @test comm_new.time == time
        end
    end

end
