using Test
using EcoEvoSim
using Random


numTests = 50


@testset "tests_of_basic_utilities" begin

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


    @testset "testing_selectTraitDim" begin
        for _ in 1:numTests
            T = Float64
            auxClasses = rand(0:3)
            numSp = rand(2:5)
            traitDim = rand(2:4)

            # Create community with multidimensional traits
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

            # Test selecting each trait dimension
            for dimIndex in 1:traitDim
                comm_selected = selectTraitDim(comm, dimIndex)

                # Check that the trait space is now 1-dimensional
                @test traitSpaceDim(comm_selected) == 1

                # Check that each species has the correct trait value
                for i in 1:numSp
                    selected_trait = traits(comm_selected, i)[1]
                    original_trait = trait_matrix[i, dimIndex]
                    @test selected_trait ≈ original_trait
                end

                # Check that population sizes are preserved
                for i in 1:numSp
                    @test popsizes(comm_selected, i) == popsizes(comm, i)
                end

                # Check that auxiliaries and time are preserved
                @test auxs(comm_selected) == aux
                @test comm_selected.time == time

                # Check number of species is preserved
                @test numSpecies(comm_selected) == numSp
            end
        end
    end


    @testset "testing_selectTraitDim_1D_traits" begin
        # Test that selectTraitDim works correctly with 1D traits
        T = Float64
        for _ in 1:numTests
            numSp = rand(2:5)
            sps = Species{T}[]
            trait_vals = rand(T, numSp)
            for i in 1:numSp
                ps = PopulationSize(rand(T, rand(1:3)))
                ph = Phenotype(trait_vals[i])
                push!(sps, Species(ps, ph))
            end
            comm = Community(sps, PopulationSize{T}[])

            # Selecting the first (only) dimension should preserve the community structure
            comm_selected = selectTraitDim(comm, 1)

            @test traitSpaceDim(comm_selected) == 1
            @test numSpecies(comm_selected) == numSp
            for i in 1:numSp
                @test traits(comm_selected, i)[1] ≈ trait_vals[i]
                @test popsizes(comm_selected, i) == popsizes(comm, i)
            end
        end
    end


    @testset "testing_selectTraitDim_edge_cases" begin
        T = Float64

        # Test with empty community
        empty_comm = Community(Species{T}[], PopulationSize{T}[])
        @test_throws ArgumentError selectTraitDim(empty_comm, 1)

        # Test with invalid dimension index (too small)
        ps = PopulationSize([1.0])
        ph = Phenotype([0.5, 0.3, 0.8])
        comm = Community([Species(ps, ph)], PopulationSize{T}[])
        @test_throws ArgumentError selectTraitDim(comm, 0)

        # Test with invalid dimension index (too large)
        @test_throws ArgumentError selectTraitDim(comm, 4)

        # Test with valid dimension indices
        comm1 = selectTraitDim(comm, 1)
        @test traits(comm1, 1)[1] ≈ 0.5

        comm2 = selectTraitDim(comm, 2)
        @test traits(comm2, 1)[1] ≈ 0.3

        comm3 = selectTraitDim(comm, 3)
        @test traits(comm3, 1)[1] ≈ 0.8
    end


    @testset "testing_selectTraitDim_with_multiple_stage_classes" begin
        # Test that selectTraitDim preserves multiple stage classes correctly
        T = Float64
        for _ in 1:numTests
            numSp = rand(2:4)
            numStageClasses = rand(2:3)
            traitDim = rand(2:4)

            sps = Species{T}[]
            for i in 1:numSp
                ps = PopulationSize(rand(T, numStageClasses))
                ph = Phenotype(rand(T, traitDim))
                push!(sps, Species(ps, ph))
            end
            comm = Community(sps, PopulationSize{T}[])

            # Select each trait dimension and verify stage classes are preserved
            for dimIndex in 1:traitDim
                comm_selected = selectTraitDim(comm, dimIndex)

                for i in 1:numSp
                    original_popsize = popsizes(comm, i)
                    selected_popsize = popsizes(comm_selected, i)
                    @test length(selected_popsize) == numStageClasses
                    @test selected_popsize ≈ original_popsize
                end
            end
        end
    end


    @testset "testing_selectTraitDim_multiple_dimensions" begin
        T = Float64
        for _ in 1:numTests
            numSp = rand(2:5)
            traitDim = rand(3:5)

            # Create community with multidimensional traits
            sps = Species{T}[]
            trait_matrix = rand(T, numSp, traitDim)
            for i in 1:numSp
                ps = PopulationSize(rand(T, rand(1:3)))
                ph = Phenotype(Vector{T}(trait_matrix[i, :]))
                push!(sps, Species(ps, ph))
            end
            aux = [PopulationSize(rand(T)) for _ in 1:rand(0:2)]
            time = rand(T)
            comm = Community(sps, aux, time)

            # Test selecting multiple dimensions
            numDimsToSelect = rand(2:traitDim-1)
            dimIndices = sort(randperm(traitDim)[1:numDimsToSelect])

            comm_selected = selectTraitDim(comm, dimIndices)

            # Check that the trait space has correct dimensionality
            @test traitSpaceDim(comm_selected) == numDimsToSelect

            # Check that each species has the correct trait values
            for i in 1:numSp
                selected_traits = traits(comm_selected, i)
                @test length(selected_traits) == numDimsToSelect
                for (j, dimIdx) in enumerate(dimIndices)
                    @test selected_traits[j] ≈ trait_matrix[i, dimIdx]
                end
            end

            # Check that population sizes are preserved
            for i in 1:numSp
                @test popsizes(comm_selected, i) == popsizes(comm, i)
            end

            # Check that auxiliaries and time are preserved
            @test auxs(comm_selected) == aux
            @test comm_selected.time == time

            # Check number of species is preserved
            @test numSpecies(comm_selected) == numSp
        end
    end


    @testset "testing_selectTraitDim_multiple_dimensions_reordering" begin
        # Test that dimension reordering works correctly
        T = Float64
        ps = PopulationSize([1.0])
        ph = Phenotype([0.1, 0.2, 0.3, 0.4, 0.5])
        comm = Community([Species(ps, ph)], PopulationSize{T}[])

        # Select dimensions in reverse order
        comm_rev = selectTraitDim(comm, [5, 4, 3, 2, 1])
        @test traitSpaceDim(comm_rev) == 5
        @test traits(comm_rev, 1) ≈ [0.5, 0.4, 0.3, 0.2, 0.1]

        # Select non-contiguous dimensions
        comm_noncontig = selectTraitDim(comm, [1, 3, 5])
        @test traitSpaceDim(comm_noncontig) == 3
        @test traits(comm_noncontig, 1) ≈ [0.1, 0.3, 0.5]

        # Select same dimensions in different order
        comm_order1 = selectTraitDim(comm, [2, 4])
        comm_order2 = selectTraitDim(comm, [4, 2])
        @test traits(comm_order1, 1) ≈ [0.2, 0.4]
        @test traits(comm_order2, 1) ≈ [0.4, 0.2]
    end


    @testset "testing_selectTraitDim_multiple_dimensions_edge_cases" begin
        T = Float64

        # Test with empty community
        empty_comm = Community(Species{T}[], PopulationSize{T}[])
        @test_throws ArgumentError selectTraitDim(empty_comm, [1, 2])

        # Test with empty dimension vector
        ps = PopulationSize([1.0])
        ph = Phenotype([0.5, 0.3, 0.8])
        comm = Community([Species(ps, ph)], PopulationSize{T}[])
        @test_throws ArgumentError selectTraitDim(comm, Int[])

        # Test with invalid dimension index (too small)
        @test_throws ArgumentError selectTraitDim(comm, [0, 1])
        @test_throws ArgumentError selectTraitDim(comm, [1, 0, 2])

        # Test with invalid dimension index (too large)
        @test_throws ArgumentError selectTraitDim(comm, [1, 4])
        @test_throws ArgumentError selectTraitDim(comm, [1, 2, 3, 4])

        # Test with duplicate indices
        @test_throws ArgumentError selectTraitDim(comm, [1, 2, 1])
        @test_throws ArgumentError selectTraitDim(comm, [2, 2])

        # Test selecting all dimensions
        comm_all = selectTraitDim(comm, [1, 2, 3])
        @test traitSpaceDim(comm_all) == 3
        @test traits(comm_all, 1) ≈ [0.5, 0.3, 0.8]

        # Test selecting single dimension using vector notation
        comm_single = selectTraitDim(comm, [2])
        @test traitSpaceDim(comm_single) == 1
        @test traits(comm_single, 1)[1] ≈ 0.3
    end


    @testset "testing_historyToTable_basic_functionality" begin
        # Create a simple history with single stage class and trait dimension
        T = Float64
        comm1 = Community([1.0, 2.0], [0.1, 0.2], Float64[])

        sp1 = Species(PopulationSize([1.5]), Phenotype([0.15]))
        sp2 = Species(PopulationSize([2.5]), Phenotype([0.25]))
        sp3 = Species(PopulationSize([3.0]), Phenotype([0.35]))
        comm2 = Community([sp1, sp2, sp3], PopulationSize{T}[], 10.0)

        sp4 = Species(PopulationSize([1.8]), Phenotype([0.18]))
        sp5 = Species(PopulationSize([2.8]), Phenotype([0.28]))
        comm3 = Community([sp4, sp5], PopulationSize{T}[], 20.0)
        history = EvoHistory([comm1, comm2, comm3])

        table = historyToTable(history)

        # Check basic structure (long format: one row per species)
        @test table isa AbstractDict
        @test haskey(table, "mutNo")
        @test haskey(table, "time")
        @test haskey(table, "species")
        @test length(table["mutNo"]) == 7  # 2 + 3 + 2 species across events

        # Check first community (2 species)
        @test table["mutNo"][1:2] == [0, 0]
        @test table["time"][1:2] == [0.0, 0.0]
        @test table["species"][1:2] == [1, 2]
        @test table["popsize_1"][1:2] == [1.0, 2.0]
        @test table["trait_1"][1:2] == [0.1, 0.2]

        # Check second community (3 species)
        @test table["mutNo"][3:5] == [1, 1, 1]
        @test table["time"][3:5] == [10.0, 10.0, 10.0]
        @test table["species"][3:5] == [1, 2, 3]
        @test table["popsize_1"][3:5] == [1.5, 2.5, 3.0]
        @test table["trait_1"][3:5] == [0.15, 0.25, 0.35]

        # Check third community (2 species)
        @test table["mutNo"][6:7] == [2, 2]
        @test table["time"][6:7] == [20.0, 20.0]
        @test table["species"][6:7] == [1, 2]
        @test table["popsize_1"][6:7] == [1.8, 2.8]
        @test table["trait_1"][6:7] == [0.18, 0.28]
    end


    @testset "testing_historyToTable_with_stage_classes" begin
        # Create history with multiple stage classes
        T = Float64
        ps1 = PopulationSize([1.0, 1.5])
        ps2 = PopulationSize([2.0, 2.5])
        ph1 = Phenotype([0.1])
        ph2 = Phenotype([0.2])
        comm1 = Community([Species(ps1, ph1), Species(ps2, ph2)], PopulationSize{T}[])

        ps3 = PopulationSize([1.2, 1.6])
        ps4 = PopulationSize([2.2, 2.6])
        comm2 = Community([Species(ps3, ph1), Species(ps4, ph2)], PopulationSize{T}[], 5.0)

        history = EvoHistory([comm1, comm2])
        table = historyToTable(history)

        # Check that both stage classes are present (4 total rows: 2 species × 2 events)
        @test haskey(table, "popsize_1")
        @test haskey(table, "popsize_2")
        @test length(table["mutNo"]) == 4

        # First community, species 1
        @test table["popsize_1"][1] == 1.0
        @test table["popsize_2"][1] == 1.5
        # First community, species 2
        @test table["popsize_1"][2] == 2.0
        @test table["popsize_2"][2] == 2.5
        # Second community, species 1
        @test table["popsize_1"][3] == 1.2
        @test table["popsize_2"][3] == 1.6
        # Second community, species 2
        @test table["popsize_1"][4] == 2.2
        @test table["popsize_2"][4] == 2.6
    end


    @testset "testing_historyToTable_with_multidimensional_traits" begin
        # Create history with 2D trait space
        T = Float64
        ph1 = Phenotype([0.1, 0.5])
        ph2 = Phenotype([0.2, 0.6])
        ps = PopulationSize([1.0])
        comm1 = Community([Species(ps, ph1), Species(ps, ph2)], PopulationSize{T}[])

        ph3 = Phenotype([0.15, 0.55])
        comm2 = Community([Species(ps, ph3)], PopulationSize{T}[], 3.0)

        history = EvoHistory([comm1, comm2])
        table = historyToTable(history)

        # Check both trait dimensions (3 rows total: 2 species + 1 species)
        @test haskey(table, "trait_1")
        @test haskey(table, "trait_2")
        @test length(table["mutNo"]) == 3

        # First community, species 1
        @test table["trait_1"][1] == 0.1
        @test table["trait_2"][1] == 0.5
        # First community, species 2
        @test table["trait_1"][2] == 0.2
        @test table["trait_2"][2] == 0.6
        # Second community, species 1
        @test table["trait_1"][3] == 0.15
        @test table["trait_2"][3] == 0.55
    end


    @testset "testing_historyToTable_with_auxiliary_variables" begin
        # Create history with auxiliary variables
        T = Float64
        ps = PopulationSize([1.0])
        ph = Phenotype([0.1])
        aux1 = [PopulationSize([10.0, 20.0])]
        comm1 = Community([Species(ps, ph)], aux1)

        aux2 = [PopulationSize([15.0, 25.0])]
        comm2 = Community([Species(ps, ph)], aux2, 2.0)

        history = EvoHistory([comm1, comm2])
        table = historyToTable(history)

        # Check auxiliary variable columns (2 rows total: 1 species × 2 events)
        @test haskey(table, "aux_1_1")
        @test haskey(table, "aux_1_2")
        @test length(table["mutNo"]) == 2
        # Aux variables are repeated for each species row
        @test table["aux_1_1"] == [10.0, 15.0]
        @test table["aux_1_2"] == [20.0, 25.0]
    end


    @testset "testing_historyToTable_empty_history_error" begin
        # Test that empty history throws error
        history = EvoHistory(Community{Float64, 0}[])
        @test_throws ArgumentError historyToTable(history)
    end


    @testset "testing_historyToTable_all_extinct" begin
        # Test with a step where all species go extinct
        T = Float64
        comm1 = Community([1.0, 2.0], [0.1, 0.2], Float64[])
        comm2 = Community(Species{T}[], PopulationSize{T}[], 5.0)  # All extinct

        sp1 = Species(PopulationSize([1.5]), Phenotype([0.15]))
        comm3 = Community([sp1], PopulationSize{T}[], 10.0)  # One survivor

        history = EvoHistory([comm1, comm2, comm3])
        table = historyToTable(history)

        # Only 3 rows total: 2 from comm1, 0 from comm2, 1 from comm3
        @test length(table["mutNo"]) == 3
        @test haskey(table, "popsize_1")

        # First community: 2 species
        @test table["mutNo"][1:2] == [0, 0]
        @test table["species"][1:2] == [1, 2]
        @test table["popsize_1"][1:2] == [1.0, 2.0]

        # Second community contributes no rows (all extinct)
        # Third community: 1 species
        @test table["mutNo"][3] == 2
        @test table["species"][3] == 1
        @test table["popsize_1"][3] == 1.5
    end


    @testset "testing_filterHistory_predicate_function" begin
        for _ in 1:numTests
            T = Float64
            n = rand(5:20)
            communities = Community{T, 0}[]
            for i in 1:n
                sp = Species(PopulationSize(T(i)), Phenotype(T(i)))
                push!(communities, Community([sp], PopulationSize{T}[], T(i)))
            end
            history = EvoHistory(communities)

            # Test filtering for first half
            filtered = filterHistory(history, i -> i <= div(n, 2))
            @test length(filtered) == div(n, 2)
            for i in 1:div(n, 2)
                @test historyList(filtered, i).time == T(i)
            end

            # Test filtering for indices matching a pattern
            filtered_odd = filterHistory(history, i -> isodd(i))
            @test length(filtered_odd) == ceil(Int, n / 2)
            for (idx, i) in enumerate(1:2:n)
                @test historyList(filtered_odd, idx).time == T(i)
            end

            # Test filtering for last k snapshots
            k = rand(1:div(n, 2))
            filtered_last = filterHistory(history, i -> i > n - k)
            @test length(filtered_last) == k
            for (idx, i) in enumerate(n-k+1:n)
                @test historyList(filtered_last, idx).time == T(i)
            end
        end
    end


    @testset "testing_filterHistory_index_vector" begin
        for _ in 1:numTests
            T = Float64
            n = rand(10:20)
            communities = Community{T, 0}[]
            for i in 1:n
                sp = Species(PopulationSize(T(i)), Phenotype(T(i)))
                push!(communities, Community([sp], PopulationSize{T}[], T(i)))
            end
            history = EvoHistory(communities)

            # Test selecting specific indices
            indices = sort(unique(rand(1:n, rand(3:min(n, 8)))))
            filtered = filterHistory(history, indices)
            @test length(filtered) == length(indices)
            for (idx, i) in enumerate(indices)
                @test historyList(filtered, idx).time == T(i)
            end

            # Test selecting a single index
            single_idx = rand(1:n)
            filtered_single = filterHistory(history, [single_idx])
            @test length(filtered_single) == 1
            @test historyList(filtered_single, 1).time == T(single_idx)

            # Test that order is preserved
            selected_indices = [n, 1, div(n, 2), n - 1]
            filtered_ordered = filterHistory(history, selected_indices)
            @test length(filtered_ordered) == length(selected_indices)
            for (idx, i) in enumerate(selected_indices)
                @test historyList(filtered_ordered, idx).time == T(i)
            end
        end
    end


    @testset "testing_filterHistory_range" begin
        for _ in 1:numTests
            T = Float64
            n = rand(10:20)
            communities = Community{T, 0}[]
            for i in 1:n
                sp = Species(PopulationSize(T(i)), Phenotype(T(i)))
                push!(communities, Community([sp], PopulationSize{T}[], T(i)))
            end
            history = EvoHistory(communities)

            # Test selecting a contiguous range
            start = rand(1:div(n, 2))
            stop = rand(start:n)
            filtered = filterHistory(history, start:stop)
            @test length(filtered) == stop - start + 1
            for (idx, i) in enumerate(start:stop)
                @test historyList(filtered, idx).time == T(i)
            end

            # Test selecting all (1:n)
            filtered_all = filterHistory(history, 1:n)
            @test length(filtered_all) == n
            for i in 1:n
                @test historyList(filtered_all, i).time == T(i)
            end

            # Test single-element range
            single = rand(1:n)
            filtered_single = filterHistory(history, single:single)
            @test length(filtered_single) == 1
            @test historyList(filtered_single, 1).time == T(single)
        end
    end


    @testset "testing_filterHistory_stride" begin
        for _ in 1:numTests
            T = Float64
            n = rand(10:30)
            communities = Community{T, 0}[]
            for i in 1:n
                sp = Species(PopulationSize(T(i)), Phenotype(T(i)))
                push!(communities, Community([sp], PopulationSize{T}[], T(i)))
            end
            history = EvoHistory(communities)

            # Test stride of 1 (should return all)
            filtered_stride1 = filterHistory(history, 1)
            @test length(filtered_stride1) == n
            for i in 1:n
                @test historyList(filtered_stride1, i).time == T(i)
            end

            # Test stride of 2 (every other)
            filtered_stride2 = filterHistory(history, 2)
            expected_length = ceil(Int, n / 2)
            @test length(filtered_stride2) == expected_length
            for (idx, i) in enumerate(1:2:n)
                @test historyList(filtered_stride2, idx).time == T(i)
            end

            # Test stride larger than history
            large_stride = n + 5
            filtered_large_stride = filterHistory(history, large_stride)
            @test length(filtered_large_stride) == 1
            @test historyList(filtered_large_stride, 1).time == T(1)

            # Test with various strides
            for stride in 3:min(n, 10)
                filtered = filterHistory(history, stride)
                expected = collect(1:stride:n)
                @test length(filtered) == length(expected)
                for (idx, i) in enumerate(expected)
                    @test historyList(filtered, idx).time == T(i)
                end
            end
        end
    end


    @testset "testing_filterHistory_error_handling" begin
        T = Float64
        n = rand(5:10)
        communities = Community{T, 0}[]
        for i in 1:n
            sp = Species(PopulationSize(T(i)), Phenotype(T(i)))
            push!(communities, Community([sp], PopulationSize{T}[], T(i)))
        end
        history = EvoHistory(communities)

        # Test out-of-bounds indices with vector
        @test_throws ArgumentError filterHistory(history, [0])
        @test_throws ArgumentError filterHistory(history, [n + 1])
        @test_throws ArgumentError filterHistory(history, [1, n + 1])

        # Test out-of-bounds range
        @test_throws ArgumentError filterHistory(history, 0:5)
        @test_throws ArgumentError filterHistory(history, 1:n+1)
        @test_throws ArgumentError filterHistory(history, -1:5)

        # Test invalid stride
        @test_throws ArgumentError filterHistory(history, 0)
        @test_throws ArgumentError filterHistory(history, -1)

        # Test predicate that selects nothing
        @test_throws ArgumentError filterHistory(history, i -> false)
    end


    @testset "testing_filterHistory_preserves_community_content" begin
        for _ in 1:numTests
            T = Float64
            n = rand(5:15)
            # Create communities with varying numbers of species
            communities = Community{T, 0}[]
            for i in 1:n
                num_sp = rand(1:5)
                sps = [Species(PopulationSize(T(i + j)), Phenotype(T(i - j)))
                       for j in 1:num_sp]
                push!(communities, Community(sps, PopulationSize{T}[], T(i)))
            end
            history = EvoHistory(communities)

            # Filter and check community content is unchanged
            indices = [1, div(n, 2), n]
            filtered = filterHistory(history, indices)

            for (idx, orig_idx) in enumerate(indices)
                orig_comm = historyList(history, orig_idx)
                filtered_comm = historyList(filtered, idx)
                @test numSpecies(filtered_comm) == numSpecies(orig_comm)
                @test speciesList(filtered_comm) == speciesList(orig_comm)
                @test filtered_comm.time == orig_comm.time
            end
        end
    end


    @testset "testing_lastCommunity" begin
        for _ in 1:numTests
            n = rand(2:10)
            communities = [Community([Species(float(i), rand(1))],
                                     PopulationSize{Float64}[], float(i))
                           for i in 1:n]
            history = EvoHistory(communities)

            last_comm = lastCommunity(history)

            # Must equal the final entry
            @test last_comm.time == communities[end].time
            @test speciesList(last_comm) == speciesList(communities[end])

            # Single-element history works too
            single = EvoHistory(communities[1:1])
            @test lastCommunity(single).time == communities[1].time

            # end-indexing works
            @test lastCommunity(history).time == history[end].time
        end
    end

end
