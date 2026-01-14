using Test
using EcoEvoSim
using StaticArrays


@testset "tests_of_basic_selectors" begin

    @testset "testing_species_selector_by_index" begin
        for _ in 1:100
            T = Float64
            stageClasses = rand(1:3)
            traitDim = rand(1:3)
            auxClasses = rand(0:3)
            numSpecies = rand(1:5)
            sps = Species{T, stageClasses, traitDim}[]
            for _ in 1:numSpecies
                ps = PopulationSize(rand(T, stageClasses))
                ph = Phenotype(rand(T, traitDim))
                push!(sps, Species(ps, ph))
            end
            aux = [PopulationSize(rand(T)) for _ in 1:auxClasses]
            comm = Community(sps, aux)
            for i in 1:numSpecies
                @test species(comm, i) === sps[i]
            end
        end
    end


    @testset "testing_species_selector_by_an_index_that_is_out_of_bounds" begin
        T = Float64
        stageClasses = 2
        traitDim = 2
        auxClasses = 1
        ps = PopulationSize(rand(T, stageClasses))
        ph = Phenotype(rand(T, traitDim))
        sp = Species(ps, ph)
        aux = [PopulationSize(rand(T))]
        comm = Community([sp], aux)
        @test_throws ArgumentError species(comm, 0)
        @test_throws ArgumentError species(comm, 2)
    end

end
