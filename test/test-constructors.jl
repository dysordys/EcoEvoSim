using Test
using EcoEvoSim


numTests = 50


@testset "tests_of_constructors" begin

    @testset "testing_PopulationSize_constructors" begin
        for _ in 1:numTests
            T = Float64
            stageClasses = rand(1:10)
            popsizeVals = rand(Float64, stageClasses)
            popsize = Vector{T}(popsizeVals)
            ps = PopulationSize{T}(popsize)
            @test ps.popsize == popsize
        end
            # Test batch constructor from matrix
            mat = rand(Float64, 4, 2)
            pss = PopulationSize(mat)
            @test length(pss) == 4
            @test all(ps -> ps isa PopulationSize, pss)
            @test all(ps -> length(ps.popsize) == 2, pss)
        # Test invalid: non-positive StageClasses
        @test_throws ArgumentError PopulationSize{Float64}(Float64[])

        # Outer constructor: accept a regular Vector and convert to SVector
        for _ in 1:numTests
            T = Float64
            stageClasses = rand(1:10)
            popsizeVals = rand(Float64, stageClasses)
            ps = PopulationSize(popsizeVals)
            @test ps.popsize == Vector{T}(popsizeVals)
        end
    end


    @testset "testing_Phenotype_constructor" begin
        for _ in 1:numTests
            T = Float64
            traitDim = rand(1:10)
            traitVals = rand(Float64, traitDim)
            trait = Vector{T}(traitVals)
            ph = Phenotype{T}(trait)
            @test ph.trait == trait
        end
        # Test invalid: non-positive TraitDim
        @test_throws ArgumentError Phenotype{Float64}(Float64[])
            # Test batch constructor from matrix
            mat = rand(Float64, 5, 3)
            phs = Phenotype(mat)
            @test length(phs) == 5
            @test all(ph -> ph isa Phenotype, phs)
            @test all(ph -> length(ph.trait) == 3, phs)
    end
        # Test batch constructor from matrices
        popMat = rand(Float64, 3, 2)
        traitMat = rand(Float64, 3, 4)
        sps = Species(popMat, traitMat)
        @test length(sps) == 3
        @test all(sp -> sp isa Species, sps)
        @test all(sp -> length(sp.popsize.popsize) == 2, sps)
        @test all(sp -> length(sp.trait.trait) == 4, sps)

        # Test batch constructor from vectors of vectors
        popVecs = [rand(Float64, 2) for _ in 1:3]
        traitVecs = [rand(Float64, 4) for _ in 1:3]
        sps2 = Species(popVecs, traitVecs)
        @test length(sps2) == 3
        @test all(sp -> sp isa Species, sps2)
        @test all(sp -> length(sp.popsize.popsize) == 2, sps2)
        @test all(sp -> length(sp.trait.trait) == 4, sps2)


    @testset "testing_Species_constructor" begin
        for _ in 1:numTests
            T = Float64
            stageClasses = rand(1:5)
            traitDim = rand(1:5)
            popsize = PopulationSize{T}(Vector{T}(rand(Float64, stageClasses)))
            trait = Phenotype{T}(Vector{T}(rand(Float64, traitDim)))
            sp = Species{T}(popsize, trait)
            @test sp.popsize === popsize
            @test sp.trait === trait
        end
    end


    @testset "testing_Species_constructor_from_single_objects" begin
        for _ in 1:numTests
            T = Float64
            stageClasses = rand(1:5)
            traitDim = rand(1:5)
            popsize = PopulationSize{T}(Vector{T}(
                rand(Float64, stageClasses)))
            phenotype = Phenotype{T}(Vector{T}(
                rand(Float64, traitDim)))
            sp = Species(popsize, phenotype)
            @test sp.popsize === popsize
            @test sp.trait === phenotype
        end
    end


    @testset "testing_Community_constructor" begin
        for _ in 1:numTests
            T = Float64
            stageClasses = rand(1:3)
            traitDim = rand(1:3)
            auxClasses = rand(0:3)
            sp = Species{T}(
                PopulationSize{T}(Vector{T}(rand(Float64, stageClasses))),
                Phenotype{T}(Vector{T}(rand(Float64, traitDim))))
            aux = [PopulationSize{T}([rand(Float64)])
                for _ in 1:auxClasses]
            comm = Community{T, auxClasses}([sp], aux)
            @test length(comm.species) == 1
            @test length(comm.aux) == auxClasses
        end
    end

    @testset "testing_Community_constructor_with_explicit_time" begin
        for _ in 1:numTests
            T = Float64
            stageClasses = rand(1:3)
            traitDim = rand(1:3)
            auxClasses = rand(0:3)
            time = rand(Float64)
            sp = Species{T}(
                PopulationSize{T}(Vector{T}(rand(Float64, stageClasses))),
                Phenotype{T}(Vector{T}(rand(Float64, traitDim))))
            aux = [PopulationSize{T}([rand(Float64)])
                for _ in 1:auxClasses]
            comm = Community{T, auxClasses}(
                [sp], aux, time)
            @test length(comm.species) == 1
            @test length(comm.aux) == auxClasses
            @test comm.time == time
        end
    end


    @testset "testing_Community_constructor_with_no_type_info" begin
        for _ in 1:numTests
            T = Float64
            stageClasses = rand(1:3)
            traitDim = rand(1:3)
            auxClasses = rand(0:3)
            sp = Species{T}(
                PopulationSize{T}(Vector{T}(rand(Float64, stageClasses))),
                Phenotype{T}(Vector{T}(rand(Float64, traitDim))))
            aux = [PopulationSize{T}([rand(Float64)])
                for _ in 1:auxClasses]
            comm = Community([sp], aux)
            @test length(comm.species) == 1
            @test length(comm.aux) == auxClasses
            @test comm.time == zero(T)
        end
    end


    @testset "testing_Community_constructor_from_values" begin
        for _ in 1:numTests
            T = Float64
            numSpecies = rand(1:5)
            popvals = rand(T, numSpecies)
            traitvals = rand(T, numSpecies) # must match length of popvals
            auxvals = rand(T, rand(0:3))

            comm = Community(popvals, traitvals, auxvals)
            @test length(comm.species) == numSpecies
            for (i, sp) in enumerate(comm.species)
                @test sp.popsize.popsize == [popvals[i]]
                @test sp.trait.trait == [traitvals[i]]
            end
            if isempty(auxvals)
                @test length(comm.aux) == 0
            else
                @test length(comm.aux) == 1  # Single structured auxiliary variable
                @test comm.aux[1].popsize == auxvals
            end
        end

        # Invalid: pop and trait must have same length
        @test_throws ArgumentError Community(rand(Float64, 3), rand(Float64, 4), Float64[])
    end

    @testset "testing_Community_constructor_from_values_without_aux" begin
        for _ in 1:numTests
            T = Float64
            numSpecies = rand(1:5)
            popvals = rand(T, numSpecies)
            traitvals = rand(T, numSpecies) # must match length of popvals

            comm = Community(popvals, traitvals)
            @test length(comm.species) == numSpecies
            @test length(comm.aux) == 0
            for (i, sp) in enumerate(comm.species)
                @test sp.popsize.popsize == [popvals[i]]
                @test sp.trait.trait == [traitvals[i]]
            end
        end

        # Verify equivalence with the version that explicitly provides empty aux
        numSpecies = 3
        popvals = [0.5, 1.0, 0.2]
        traitvals = [-0.1, 0.5, 0.3]
        comm_no_aux = Community(popvals, traitvals)
        comm_with_empty_aux = Community(popvals, traitvals, Float64[])
        @test length(comm_no_aux.species) == length(comm_with_empty_aux.species)
        @test length(comm_no_aux.aux) == length(comm_with_empty_aux.aux) == 0
        for i in 1:numSpecies
            @test comm_no_aux.species[i].popsize.popsize == comm_with_empty_aux.species[i].popsize.popsize
            @test comm_no_aux.species[i].trait.trait == comm_with_empty_aux.species[i].trait.trait
        end

        # Invalid: pop and trait must have same length
        @test_throws ArgumentError Community(rand(Float64, 3), rand(Float64, 4))
    end

    @testset "testing_Community_constructor_from_popMat_and_traitVec" begin
        for _ in 1:numTests
            T = Float64
            numSpecies = rand(1:5)
            stageClasses = rand(1:5)
            popMat = rand(T, numSpecies, stageClasses)
            traitVec = rand(T, numSpecies)

            comm = Community(popMat, traitVec)
            @test length(comm.species) == numSpecies
            @test length(comm.aux) == 0
            for (i, sp) in enumerate(comm.species)
                @test sp.popsize.popsize == vec(popMat[i, :])
                @test sp.trait.trait == [traitVec[i]]
            end
        end

        # Invalid: popMat and traitVec must have matching number of species
        @test_throws ArgumentError Community(rand(Float64, 3, 2), rand(Float64, 4))
    end

    @testset "testing_Community_constructor_from_popMat_traitVec_and_aux" begin
        for _ in 1:numTests
            T = Float64
            numSpecies = rand(1:5)
            stageClasses = rand(1:5)
            popMat = rand(T, numSpecies, stageClasses)
            traitVec = rand(T, numSpecies)
            auxvals = rand(T, rand(0:3))

            comm = Community(popMat, traitVec, auxvals)
            @test length(comm.species) == numSpecies
            if isempty(auxvals)
                @test length(comm.aux) == 0
            else
                @test length(comm.aux) == 1  # Single structured auxiliary variable
            end
            for (i, sp) in enumerate(comm.species)
                @test sp.popsize.popsize == vec(popMat[i, :])
                @test sp.trait.trait == [traitVec[i]]
            end
            if !isempty(auxvals)
                @test comm.aux[1].popsize == auxvals
            end
        end

        # Invalid: popMat and traitVec must have matching number of species
        @test_throws ArgumentError Community(rand(Float64, 3, 2), rand(Float64, 4), Float64[])
    end

    @testset "testing_Community_constructor_from_popVec_and_traitMat" begin
        for _ in 1:numTests
            T = Float64
            numSpecies = rand(1:5)
            traitDim = rand(1:5)
            popVec = rand(T, numSpecies)
            traitMat = rand(T, numSpecies, traitDim)

            comm = Community(popVec, traitMat)
            @test length(comm.species) == numSpecies
            @test length(comm.aux) == 0
            for (i, sp) in enumerate(comm.species)
                @test sp.popsize.popsize == [popVec[i]]
                @test sp.trait.trait == vec(traitMat[i, :])
            end
        end

        # Invalid: popVec and traitMat must have matching number of species
        @test_throws ArgumentError Community(rand(Float64, 3), rand(Float64, 4, 2))
    end

    @testset "testing_Community_constructor_from_popVec_traitMat_and_aux" begin
        for _ in 1:numTests
            T = Float64
            numSpecies = rand(1:5)
            traitDim = rand(1:5)
            popVec = rand(T, numSpecies)
            traitMat = rand(T, numSpecies, traitDim)
            auxvals = rand(T, rand(0:3))

            comm = Community(popVec, traitMat, auxvals)
            @test length(comm.species) == numSpecies
            if isempty(auxvals)
                @test length(comm.aux) == 0
            else
                @test length(comm.aux) == 1  # Single structured auxiliary variable
            end
            for (i, sp) in enumerate(comm.species)
                @test sp.popsize.popsize == [popVec[i]]
                @test sp.trait.trait == vec(traitMat[i, :])
            end
            if !isempty(auxvals)
                @test comm.aux[1].popsize == auxvals
            end
        end

        # Invalid: popVec and traitMat must have matching number of species
        @test_throws ArgumentError Community(rand(Float64, 3), rand(Float64, 4, 2), Float64[])
    end

    @testset "testing_Community_constructor_from_popMat_and_traitMat" begin
        for _ in 1:numTests
            T = Float64
            numSpecies = rand(1:5)
            stageClasses = rand(1:5)
            traitDim = rand(1:5)
            popMat = rand(T, numSpecies, stageClasses)
            traitMat = rand(T, numSpecies, traitDim)

            comm = Community(popMat, traitMat)
            @test length(comm.species) == numSpecies
            @test length(comm.aux) == 0
            for (i, sp) in enumerate(comm.species)
                @test sp.popsize.popsize == vec(popMat[i, :])
                @test sp.trait.trait == vec(traitMat[i, :])
            end
        end

        # Invalid: popMat and traitMat must have matching number of species
        @test_throws ArgumentError Community(rand(Float64, 3, 2), rand(Float64, 4, 2))
    end

    @testset "testing_Community_constructor_from_popMat_traitMat_and_aux" begin
        for _ in 1:numTests
            T = Float64
            numSpecies = rand(1:5)
            stageClasses = rand(1:5)
            traitDim = rand(1:5)
            popMat = rand(T, numSpecies, stageClasses)
            traitMat = rand(T, numSpecies, traitDim)
            auxvals = rand(T, rand(0:3))

            comm = Community(popMat, traitMat, auxvals)
            @test length(comm.species) == numSpecies
            if isempty(auxvals)
                @test length(comm.aux) == 0
            else
                @test length(comm.aux) == 1  # Single structured auxiliary variable
            end
            for (i, sp) in enumerate(comm.species)
                @test sp.popsize.popsize == vec(popMat[i, :])
                @test sp.trait.trait == vec(traitMat[i, :])
            end
            if !isempty(auxvals)
                @test comm.aux[1].popsize == auxvals
            end
        end

        # Invalid: popMat and traitMat must have matching number of species
        @test_throws ArgumentError Community(rand(Float64, 3, 2), rand(Float64, 4, 2), Float64[])
    end


    @testset "testing_EvoHistory_constructor" begin
        for _ in 1:numTests
            T = Float64
            stageClasses = rand(1:3)
            traitDim = rand(1:3)
            auxClasses = rand(0:3)
            sp = Species{T}(
                PopulationSize{T}(Vector{T}(rand(Float64, stageClasses))),
                Phenotype{T}(Vector{T}(rand(Float64, traitDim))))
            aux = [PopulationSize{T}([rand(Float64)])
                for _ in 1:auxClasses]
            comm = Community{T, auxClasses}([sp], aux)
            hist = EvoHistory(comm)
            @test length(hist.history) == 1
            @test hist.history[1] === comm
        end
    end


    @testset "testing_EvoHistory_constructor_from_vector" begin
        for _ in 1:numTests
            T = Float64
            stageClasses = rand(1:3)
            traitDim = rand(1:3)
            auxClasses = rand(0:3)
            numComms = rand(2:5)
            comms = Community{T, auxClasses}[]
            for _ in 1:numComms
                sp = Species{T}(
                    PopulationSize{T}(Vector{T}(rand(Float64, stageClasses))),
                    Phenotype{T}(Vector{T}(rand(Float64, traitDim))))
                aux = [PopulationSize{T}([rand(Float64)])
                    for _ in 1:auxClasses]
                comm = Community{T, auxClasses}([sp], aux)
                push!(comms, comm)
            end
            hist = EvoHistory(comms)
            @test length(hist.history) == numComms
            @test all(hist.history[i] === comms[i] for i in 1:numComms)
        end
    end


    @testset "testing_emptyCommunity_and_emptyEvoHistory" begin
        ec = emptyCommunity()
        @test ec isa Community{Float64, 0}
        @test numSpecies(ec) == 0
        @test length(ec.aux) == 0
        @test ec.time == 0.0
        @test speciesList(ec) == Species{Float64}[]

        ec32 = emptyCommunity(Float32)
        @test ec32 isa Community{Float32, 0}
        @test numSpecies(ec32) == 0
        @test ec32.time == 0f0

        eh = emptyEvoHistory()
        @test eh isa EvoHistory{Float64, 0}
        @test length(eh.history) == 1
        @test eh.history[1] isa Community{Float64, 0}
        @test numSpecies(eh.history[1]) == 0

        eh32 = emptyEvoHistory(Float32)
        @test eh32 isa EvoHistory{Float32, 0}
        @test length(eh32.history) == 1
        @test eh32.history[1] isa Community{Float32, 0}
    end

end
