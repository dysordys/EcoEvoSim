using Test
using EcoEvoSim


numTests = 50


@testset "tests_of_basic_utils" begin

    @testset "testing_addSpecies_single_species" begin
        for _ in 1:numTests
            T = Float64
            auxClasses = rand(0:3)
            numSpeciesInit = rand(1:5)
            sps = Species{T}[]
            for _ in 1:numSpeciesInit
                ps = PopulationSize(rand(T, rand(1:3)))
                ph = Phenotype(rand(T, rand(1:3)))
                push!(sps, Species(ps, ph))
            end
            aux = [PopulationSize(rand(T)) for _ in 1:auxClasses]
            comm = Community(sps, aux)

            # Create a new species to add
            new_ps = PopulationSize(rand(T, rand(1:3)))
            new_ph = Phenotype(rand(T, rand(1:3)))
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
            auxClasses = rand(0:3)
            numSpeciesInit = rand(1:5)
            sps = Species{T}[]
            for _ in 1:numSpeciesInit
                ps = PopulationSize(rand(T, rand(1:3)))
                ph = Phenotype(rand(T, rand(1:3)))
                push!(sps, Species(ps, ph))
            end
            aux = [PopulationSize(rand(T)) for _ in 1:auxClasses]
            comm = Community(sps, aux)

            # Create multiple new species to add
            numNew = rand(1:3)
            new_sps = Species{T}[]
            for _ in 1:numNew
                ps = PopulationSize(rand(T, rand(1:3)))
                ph = Phenotype(rand(T, rand(1:3)))
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
            auxClasses = rand(1:3)
            numSpecies = rand(1:5)
            sps = Species{T}[]
            for _ in 1:numSpecies
                ps = PopulationSize(rand(T, rand(1:3)))
                ph = Phenotype(rand(T, rand(1:3)))
                push!(sps, Species(ps, ph))
            end
            aux = [PopulationSize(rand(T)) for _ in 1:auxClasses]
            time = rand(T)
            comm = Community(sps, aux, time)

            # Create a new species to add
            new_ps = PopulationSize(rand(T, rand(1:3)))
            new_ph = Phenotype(rand(T, rand(1:3)))
            new_sp = Species(new_ps, new_ph)

            comm_new = addSpecies(comm, new_sp)

            # Check that auxiliaries and time are preserved
            @test auxs(comm_new) == aux
            @test comm_new.time == time
        end
    end


    @testset "testing_speciesBelowThreshold" begin
        for _ in 1:numTests
            T = Float64
            auxClasses = rand(0:3)
            numSp = rand(3:6)

            # Create community with known population sizes
            sps = Species{T}[]
            pop_totals = Float64[]
            for _ in 1:numSp
                # Create species with single stage class for simplicity
                pop = rand(T)
                push!(pop_totals, pop)
                ps = PopulationSize([pop])
                ph = Phenotype(rand(T, rand(1:3)))
                push!(sps, Species(ps, ph))
            end
            aux = [PopulationSize(rand(T)) for _ in 1:auxClasses]
            comm = Community(sps, aux)

            # Test with threshold that captures some species
            sorted_pops = sort(pop_totals)
            threshold = sorted_pops[div(numSp, 2) + 1]  # Middle value
            indices = speciesBelowThreshold(comm, threshold)

            # Verify correct species identified
            for i in 1:numSp
                if pop_totals[i] < threshold
                    @test i in indices
                else
                    @test i ∉ indices
                end
            end
        end
    end


    @testset "testing_speciesBelowThreshold_edge_cases" begin
        T = Float64

        # Test empty community
        empty_comm = Community(Species{T}[], PopulationSize{T}[])
        @test speciesBelowThreshold(empty_comm, 0.5) == Int[]

        # Test threshold of zero (no species below)
        ps1 = PopulationSize([0.5])
        ph1 = Phenotype([0.1])
        comm = Community([Species(ps1, ph1)], PopulationSize{T}[])
        @test isempty(speciesBelowThreshold(comm, 0.0))

        # Test negative threshold throws error
        @test_throws ArgumentError speciesBelowThreshold(comm, -0.1)

        # Test all species below threshold
        ps1 = PopulationSize([0.1])
        ps2 = PopulationSize([0.2])
        ph = Phenotype([0.5])
        comm = Community([Species(ps1, ph), Species(ps2, ph)], PopulationSize{T}[])
        @test speciesBelowThreshold(comm, 1.0) == [1, 2]
    end


    @testset "testing_removeExtinct" begin
        for _ in 1:numTests
            T = Float64
            auxClasses = rand(0:3)
            numSp = rand(3:6)

            # Create community with known population sizes
            sps = Species{T}[]
            pop_totals = Float64[]
            for _ in 1:numSp
                # Single stage class for simplicity
                pop = rand(T)
                push!(pop_totals, pop)
                ps = PopulationSize([pop])
                ph = Phenotype(rand(T, rand(1:3)))
                push!(sps, Species(ps, ph))
            end
            aux = [PopulationSize(rand(T)) for _ in 1:auxClasses]
            time = rand(T)
            comm = Community(sps, aux, time)

            # Test removal with threshold
            sorted_pops = sort(pop_totals)
            threshold = sorted_pops[div(numSp, 2) + 1]  # Middle value
            comm_pruned = removeExtinct(comm, threshold)

            # Count expected survivors
            expected_survivors = sum(pop_totals .>= threshold)
            @test numSpecies(comm_pruned) == expected_survivors

            # Verify all remaining species are above threshold
            for i in 1:numSpecies(comm_pruned)
                pop_vec = popsizes(comm_pruned, i)
                @test sum(pop_vec) >= threshold
            end

            # Verify auxiliaries and time preserved
            @test auxs(comm_pruned) == aux
            @test comm_pruned.time == time
        end
    end


    @testset "testing_removeExtinct_edge_cases" begin
        T = Float64

        # Test empty community
        empty_comm = Community(Species{T}[], PopulationSize{T}[])
        @test removeExtinct(empty_comm, 0.5) === empty_comm

        # Test no species removed (all above threshold)
        ps1 = PopulationSize([1.0])
        ps2 = PopulationSize([2.0])
        ph = Phenotype([0.5])
        aux = [PopulationSize([0.1])]
        comm = Community([Species(ps1, ph), Species(ps2, ph)], aux)
        comm_pruned = removeExtinct(comm, 0.5)
        @test numSpecies(comm_pruned) == 2
        @test comm_pruned === comm  # Should return same object

        # Test all species removed
        ps1 = PopulationSize([0.1])
        ps2 = PopulationSize([0.2])
        comm = Community([Species(ps1, ph), Species(ps2, ph)], aux)
        comm_pruned = removeExtinct(comm, 1.0)
        @test numSpecies(comm_pruned) == 0
        @test auxs(comm_pruned) == aux  # Auxiliaries preserved
    end


    @testset "testing_orderByTrait" begin
        for _ in 1:numTests
            T = Float64
            auxClasses = rand(0:3)
            numSp = rand(2:5)
            traitDim = rand(2:4)

            # Create community with known trait values
            sps = Species{T}[]
            trait_matrix = rand(T, numSp, traitDim)
            for i in 1:numSp
                ps = PopulationSize(rand(T, rand(1:3)))
                ph = Phenotype(Vector{T}(trait_matrix[i, :]))
                push!(sps, Species(ps, ph))
            end
            aux = [PopulationSize(rand(T)) for _ in 1:auxClasses]
            time = rand(T)
            comm = Community(sps, aux, time)

            # Test ordering by different trait dimensions
            for n in 1:traitDim
                comm_ordered = orderByTrait(comm, n)

                # Verify correct number of species
                @test numSpecies(comm_ordered) == numSp

                # Verify ordering
                for i in 1:(numSp-1)
                    trait_i = traits(comm_ordered, i)
                    trait_next = traits(comm_ordered, i+1)
                    @test trait_i[n] <= trait_next[n]
                end

                # Verify auxiliaries and time preserved
                @test auxs(comm_ordered) == aux
                @test comm_ordered.time == time
            end

            # Test default (first trait dimension)
            comm_default = orderByTrait(comm)
            comm_first = orderByTrait(comm, 1)
            for i in 1:numSp
                @test traits(comm_default, i) == traits(comm_first, i)
            end
        end
    end


    @testset "testing_orderByTrait_edge_cases" begin
        T = Float64

        # Test empty community
        empty_comm = Community(Species{T}[], PopulationSize{T}[])
        @test orderByTrait(empty_comm) === empty_comm
        @test orderByTrait(empty_comm, 1) === empty_comm

        # Test single species
        ps = PopulationSize([1.0])
        ph = Phenotype([0.5, 0.3])
        comm = Community([Species(ps, ph)], PopulationSize{T}[])
        @test numSpecies(orderByTrait(comm)) == 1
        @test numSpecies(orderByTrait(comm, 2)) == 1

        # Test invalid trait dimension
        ps1 = PopulationSize([1.0])
        ps2 = PopulationSize([2.0])
        ph = Phenotype([0.5])
        comm = Community([Species(ps1, ph), Species(ps2, ph)], PopulationSize{T}[])
        @test_throws ArgumentError orderByTrait(comm, 0)
        @test_throws ArgumentError orderByTrait(comm, 2)

        # Test species already ordered
        ph1 = Phenotype([0.1, 0.5])
        ph2 = Phenotype([0.2, 0.3])
        ph3 = Phenotype([0.3, 0.1])
        comm = Community([Species(ps1, ph1), Species(ps2, ph2), Species(ps1, ph3)],
                          PopulationSize{T}[])
        comm_ordered = orderByTrait(comm, 1)
        @test traits(comm_ordered, 1)[1] ≈ 0.1
        @test traits(comm_ordered, 2)[1] ≈ 0.2
        @test traits(comm_ordered, 3)[1] ≈ 0.3
    end

end
