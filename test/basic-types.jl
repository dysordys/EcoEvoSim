using Test
using EcoEvoSim
using StaticArrays


@testset "tests_of_basic_ecoevo_types" begin

    @testset "testing_PopulationSize_constructors" begin
        for _ in 1:100
            T = Float64
            stageClasses = rand(1:10)
            popsizeVals = rand(Float64, stageClasses)
            popsize = SVector{stageClasses, T}(popsizeVals)
            ps = PopulationSize{T, stageClasses}(popsize)
            @test ps.popsize == popsize
        end
        # Test invalid: non-positive StageClasses
        @test_throws ArgumentError PopulationSize{Float64, 0}(SVector{0, Float64}())

        # Outer constructor: accept a regular Vector and convert to SVector
        for _ in 1:100
            T = Float64
            stageClasses = rand(1:10)
            popsizeVals = rand(Float64, stageClasses)
            ps = PopulationSize(popsizeVals)
            @test ps.popsize == SVector{stageClasses, T}(popsizeVals)
        end
    end


    @testset "testing_Phenotype_constructor" begin
        for _ in 1:100
            T = Float64
            traitDim = rand(1:10)
            traitVals = rand(Float64, traitDim)
            trait = SVector{traitDim, T}(traitVals)
            ph = Phenotype{T, traitDim}(trait)
            @test ph.trait == trait
        end
        # Test invalid: non-positive TraitDim
        @test_throws ArgumentError Phenotype{Float64, 0}(SVector{0, Float64}())
    end


    @testset "testing_Species_constructor" begin
        for _ in 1:100
            T = Float64
            stageClasses = rand(1:5)
            traitDim = rand(1:5)
            num_instances = rand(1:3)
            popsize = [PopulationSize{T, stageClasses}(SVector{stageClasses, T}(
                rand(Float64, stageClasses)
            )) for _ in 1:num_instances]
            trait = [Phenotype{T, traitDim}(SVector{traitDim, T}(
                rand(Float64, traitDim)
            )) for _ in 1:num_instances]
            sp = Species{T, stageClasses, traitDim}(popsize, trait)
            @test length(sp.popsize) == num_instances
            @test length(sp.trait) == num_instances
        end
    end


    @testset "testing_Community_constructor" begin
        for _ in 1:100
            T = Float64
            stageClasses = rand(1:3)
            traitDim = rand(1:3)
            auxClasses = rand(0:3)
            sp = Species{T, stageClasses, traitDim}(
                [PopulationSize{T, stageClasses}(SVector{stageClasses, T}(
                    rand(Float64, stageClasses)))],
                [Phenotype{T, traitDim}(SVector{traitDim, T}(
                    rand(Float64, traitDim)))])
            aux = [PopulationSize{T, 1}(SVector{1, T}(rand(Float64)))
                for _ in 1:auxClasses]
            history = Community{T, stageClasses, traitDim, auxClasses}[]
            comm = Community{T, stageClasses, traitDim, auxClasses}([sp], aux, history)
            @test length(comm.species) == 1
            @test length(comm.aux) == auxClasses
        end
    end

    @testset "testing_Community_constructor_with_explicit_time" begin
        for _ in 1:100
            T = Float64
            stageClasses = rand(1:3)
            traitDim = rand(1:3)
            auxClasses = rand(0:3)
            time = rand(Float64)
            sp = Species{T, stageClasses, traitDim}(
                [PopulationSize{T, stageClasses}(SVector{stageClasses, T}(
                    rand(Float64, stageClasses)))],
                [Phenotype{T, traitDim}(SVector{traitDim, T}(
                    rand(Float64, traitDim)))])
            aux = [PopulationSize{T, 1}(SVector{1, T}(rand(Float64)))
                for _ in 1:auxClasses]
            history = Community{T, stageClasses, traitDim, auxClasses}[]
            comm = Community{T, stageClasses, traitDim, auxClasses}(
                [sp], aux, history, time)
            @test length(comm.species) == 1
            @test length(comm.aux) == auxClasses
            @test comm.time == time
        end
    end


    @testset "testing_Community_constructor_from_values" begin
        for _ in 1:100
            T = Float64
            stageClasses = rand(1:3)
            popvals = rand(T, stageClasses)
            traitvals = rand(T, stageClasses) # must match length of popvals
            auxvals = rand(T, rand(0:3))

            comm = Community(popvals, traitvals, auxvals)
            @test length(comm.species) == 1
            sp = comm.species[1]
            @test sp.popsize[1].popsize == SVector{stageClasses, T}(popvals)
            @test sp.trait[1].trait == SVector{stageClasses, T}(traitvals)
            @test length(comm.aux) == length(auxvals)
            for (a,v) in zip(comm.aux, auxvals)
                @test a.popsize[1] == v
            end
        end

        # Invalid: pop and trait must have same length
        @test_throws ArgumentError Community(rand(Float64, 3), rand(Float64, 4), Float64[])
    end

end
