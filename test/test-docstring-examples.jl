using EcoEvoSim
using Test
using Random
using OrdinaryDiffEq
using SteadyStateDiffEq
using Plots
using Distributions

@testset "Docstring Examples" begin

    @testset "Basic Types" begin
        # PopulationSize examples
        @test_nowarn PopulationSize(10.5)
        @test_nowarn PopulationSize([5.0, 15.0, 8.0])

        # Phenotype examples
        @test_nowarn Phenotype(0.5)
        @test_nowarn Phenotype([0.5, -0.2, 1.0])

        # Species examples
        @test_nowarn Species(10.0, 0.5)
        @test_nowarn Species([5.0, 15.0], [0.5, -0.2])

        # Community examples
        @test_nowarn Community([1.0, 2.0, 3.0], [0.1, 0.2, 0.3], Float64[])
        species_vec = [Species(1.0, 0.1), Species(2.0, 0.2)]
        aux_vec = [PopulationSize(100.0)]
        @test_nowarn Community(species_vec, aux_vec, 0.0)

        # EvoHistory examples
        initial_community = Community([1.0, 2.0], [0.1, 0.2], Float64[])
        @test_nowarn EvoHistory(initial_community)
    end

    @testset "Basic Constructors" begin
        # emptyCommunity examples
        @test_nowarn emptyCommunity()
        @test_nowarn emptyCommunity(Float32)

        # emptyEvoHistory examples
        @test_nowarn emptyEvoHistory()
        @test_nowarn emptyEvoHistory(Float32)
    end

    @testset "Basic Selectors" begin
        comm = Community([1.0, 2.0, 3.0], [0.1, 0.2, 0.3], Float64[])

        # numSpecies examples
        n = numSpecies(comm)
        @test n == 3

        # speciesIndices examples
        indices = speciesIndices(comm)
        @test indices == 1:3

        # speciesList examples
        @test_nowarn speciesList(comm)
        @test_nowarn speciesList(comm, 1)
        @test_nowarn speciesList(comm, [1, 3])

        # popsizes examples
        @test_nowarn popsizes(comm)
        @test_nowarn popsizes(comm, 1)
        @test_nowarn popsizes(comm, [1, 3])

        # traitSpaceDim examples
        comm2 = Community([1.0, 2.0], [0.1, 0.2], Float64[])
        dim = traitSpaceDim(comm2)
        @test dim == 1

        # traits examples
        @test_nowarn traits(comm)
        @test_nowarn traits(comm, 1)
        @test_nowarn traits(comm, [1, 3])

        # auxs examples
        species = [Species(1.0, 0.1), Species(2.0, 0.2)]
        aux = [PopulationSize(100.0)]
        comm_with_aux = Community(species, aux, 0.0)
        @test_nowarn auxs(comm_with_aux)
        @test_nowarn auxs(comm_with_aux, 1)

        # randomSpecies examples
        Random.seed!(42)  # For reproducibility
        @test_nowarn randomSpecies(comm)
        @test_nowarn randomSpecies(comm, 2)

        # weightedRandomSpecies examples
        comm3 = Community([1.0, 10.0, 5.0], [0.1, 0.2, 0.3], Float64[])
        Random.seed!(42)
        @test_nowarn weightedRandomSpecies(comm3)

        # speciesBelowThreshold examples
        comm4 = Community([0.001, 2.0, 0.005], [0.1, 0.2, 0.3], Float64[])
        extinct_indices = speciesBelowThreshold(comm4, 0.01)
        @test extinct_indices == [1, 3]

        # numStages examples
        comm_stages = Community([1.0 2.0; 3.0 4.0], [0.1, 0.2], Float64[])
        nstages = numStages(comm_stages)
        @test nstages == 2

        # historyList examples
        hist_comms = [Community([1.0 * i, 2.0 * i], [0.1, 0.2], Float64[]) for i in 1:5]
        test_hist = EvoHistory(hist_comms)
        all_communities = historyList(test_hist)
        @test all_communities isa Vector
        @test length(all_communities) == 5
        single_snap = historyList(test_hist, 3)
        @test single_snap isa Community
        subset_snaps = historyList(test_hist, [1, 3, 5])
        @test length(subset_snaps) == 3
    end

    @testset "Basic Utils" begin
        # addSpecies examples
        comm = Community([1.0, 2.0], [0.1, 0.2], Float64[])
        new_sp = Species(1.0, 0.5)
        comm2 = addSpecies(comm, new_sp)
        @test numSpecies(comm2) == 3

        # changePopsizes examples
        comm = Community([1.0, 2.0, 3.0], [0.1, 0.2, 0.3], Float64[])
        comm2 = changePopsizes(comm, [2.0, 3.0, 4.0])
        @test popsizes(comm2, 1) == [2.0]

        # changeTraits examples
        comm = Community([1.0, 2.0, 3.0], [0.1, 0.2, 0.3], Float64[])
        comm2 = changeTraits(comm, [0.15, 0.25, 0.35])
        @test traits(comm2, 1) == [0.15]

        # removeSpecies examples
        comm = Community([1.0, 2.0, 3.0], [0.1, 0.2, 0.3], Float64[])
        comm2 = removeSpecies(comm, 2)
        @test numSpecies(comm2) == 2
        comm3 = removeSpecies(comm, [1, 3])
        @test numSpecies(comm3) == 1

        # removeExtinct examples
        comm = Community([0.001, 2.0, 0.005, 3.0], [0.1, 0.2, 0.3, 0.4], Float64[])
        comm2 = removeExtinct(comm, 0.01)
        @test numSpecies(comm2) == 2

        # orderByTrait examples
        comm = Community([1.0, 2.0, 3.0], [0.3, 0.1, 0.2], Float64[])
        sorted_comm = orderByTrait(comm)
        @test traits(sorted_comm, 1) == [0.1]
        @test traits(sorted_comm, 3) == [0.3]

        # orderByTrait with dimension argument
        comm_2d = Community([1.0, 2.0], [0.5 0.9; 0.3 0.7])
        sorted_2d = orderByTrait(comm_2d, 2)  # Sort by second trait dimension
        @test traits(sorted_2d, 1)[2] < traits(sorted_2d, 2)[2]

        # length(EvoHistory) examples
        simple_hist = EvoHistory([Community([1.0], [0.1], Float64[]),
                                   Community([1.1], [0.1], Float64[])])
        @test length(simple_hist) == 2

        # filterHistory examples
        filt_comms = [Community([Float64(i)], [0.1], Float64[]) for i in 1:10]
        filt_hist = EvoHistory(filt_comms)

        filtered_fn = filterHistory(filt_hist, i -> isodd(i))
        @test length(filtered_fn) == 5

        filtered_idx = filterHistory(filt_hist, [1, 5, 10])
        @test length(filtered_idx) == 3

        filtered_range = filterHistory(filt_hist, 3:7)
        @test length(filtered_range) == 5

        filtered_step = filterHistory(filt_hist, 2)  # every other snapshot
        @test length(filtered_step) == 5

        # changePopsizes matrix variant (multiple stage classes)
        comm_stages = Community([1.0 2.0; 3.0 4.0], [0.1, 0.2])
        comm_stages2 = changePopsizes(comm_stages, [1.5 2.5; 3.5 4.5])
        @test popsizes(comm_stages2, 1) == [1.5, 2.5]
        @test popsizes(comm_stages2, 2) == [3.5, 4.5]

        # changeTraits multidimensional variant
        comm_multi_trait = Community([1.0, 2.0], [0.1 0.5; 0.2 0.6])
        comm_mt2 = changeTraits(comm_multi_trait, [0.15 0.55; 0.25 0.65])
        @test traits(comm_mt2, 1) == [0.15, 0.55]
        @test traits(comm_mt2, 2) == [0.25, 0.65]

        # selectTraitDim single dimension
        comm_2dtrait = Community([1.0, 2.0], [0.3 0.5; 0.1 0.8])
        comm_1dtrait = selectTraitDim(comm_2dtrait, 1)
        @test traitSpaceDim(comm_1dtrait) == 1
        @test traits(comm_1dtrait, 1) == [0.3]
        comm_2nd = selectTraitDim(comm_2dtrait, 2)
        @test traits(comm_2nd, 1) == [0.5]

        # selectTraitDim multiple dimensions
        comm_4dtrait = Community([1.0, 2.0], [0.1 0.2 0.3 0.4; 0.5 0.6 0.7 0.8])
        comm_2dtrait2 = selectTraitDim(comm_4dtrait, [1, 3])
        @test traitSpaceDim(comm_2dtrait2) == 2
        @test traits(comm_2dtrait2, 1) == [0.1, 0.3]
        @test traits(comm_2dtrait2, 2) == [0.5, 0.7]
    end

    @testset "Eco-Evo Configuration" begin
        # IntegrationParams examples
        @test_nowarn IntegrationParams(maxTime = 50.0)
        @test_nowarn IntegrationParams(maxTime = 100.0, algorithm = Tsit5(),
                                       abstol = 1e-10, reltol = 1e-8)
        @test_nowarn IntegrationParams(maxTime = 100.0, algorithm = DynamicSS())

        # EcoEvoConfig examples
        growthFn(z) = 1.0 - sum(z.^2)
        kernelFn(z_i, z_j) = -exp(-sum((z_i .- z_j).^2) / 0.15^2)
        mutationGen = generateMutant(invaderPopsize=0.001, variance=0.01)

        config = EcoEvoConfig(
            ecoDyn = lotkaVolterra(growthFn, kernelFn),
            mutationGenerator = mutationGen,
            integrationParams = IntegrationParams(maxTime = 50.0),
            extThreshold = 0.003
        )
        @test config.extThreshold == 0.003
    end

    @testset "Eco-Evo Functions" begin
        # Set up configuration
        growthFn(z) = 1.0 - sum(z.^2)
        kernelFn(z_i, z_j) = -exp(-sum((z_i .- z_j).^2) / 0.15^2)

        config = EcoEvoConfig(
            ecoDyn = lotkaVolterra(growthFn, kernelFn),
            mutationGenerator = generateMutant(invaderPopsize=0.001, variance=0.01),
            integrationParams = IntegrationParams(maxTime = 10.0, abstol=1e-6, reltol=1e-4),
            extThreshold = 0.003
        )

        comm = Community([1.0, 1.0], [0.1, 0.2], Float64[])

        # unpackCommunity examples
        @test_nowarn unpackCommunity(comm)
        u = unpackCommunity(comm)
        @test length(u) == 2

        # packCommunity examples
        @test_nowarn packCommunity(u, comm, 0.0)

        # ecoDyn examples
        # Note: ecoDyn may produce solver warnings for some parameter combinations
        # We just verify it doesn't throw an error
        result = ecoDyn(comm, config)
        @test result isa Community

        # generateMutant examples
        Random.seed!(42)
        @test_nowarn generateMutant(; invaderPopsize=0.001, variance=0.01)
        covMat = [0.01;;]  # 1x1 matrix
        @test_nowarn generateMutant(; invaderPopsize=0.001, covMat=covMat)

        # generateMutantWeighted examples
        Random.seed!(42)
        @test_nowarn generateMutantWeighted(; invaderPopsize=0.001, variance=0.01)
        @test_nowarn generateMutantWeighted(; invaderPopsize=0.001, covMat=covMat)

        # singleEvoStep examples
        Random.seed!(42)
        result = singleEvoStep(comm, config)
        @test result isa Community

        # evolve! examples
        Random.seed!(42)
        history = EvoHistory(comm)
        EcoEvoSim.evolve!(history, config, 5, showProgress=false)
        @test length(history.history) == 6  # Initial + 5 steps

        # Alternative evolve! form
        Random.seed!(42)
        comm_new = Community([1.0, 1.0], [0.1, 0.2], Float64[])
        history2 = EcoEvoSim.evolve!(comm_new, config, 5, showProgress=false)
        @test length(history2.history) == 6

        # generateMutantSpatial examples (1 species, 2 patches)
        Random.seed!(42)
        spatial_comm = Community([1.0 1.0], [0.0])
        gen_spatial = generateMutantSpatial(invaderPopsize=0.001, variance=0.01^2)
        mutant_spatial = gen_spatial(spatial_comm)
        @test numSpecies(mutant_spatial) == 2

        # generateMutantSpatialWeighted examples (2 species, 2 patches each)
        Random.seed!(42)
        spatial_comm2 = Community([1.0 1.0; 10.0 10.0], [0.0, 0.3])
        gen_spatial_w = generateMutantSpatialWeighted(invaderPopsize=0.001, variance=0.01^2)
        mutant_spatial2 = gen_spatial_w(spatial_comm2)
        @test numSpecies(mutant_spatial2) == 3
    end

    @testset "Models" begin
        # lotkaVolterra examples (from docstring)
        growthFn(z) = 1.0 - sum(z.^2)
        kernelFn(z_i, z_j) = -exp(-sum((z_i .- z_j).^2) / 0.15^2)

        community = Community([1.0, 1.0, 1.0], [-0.2, 0.0, 0.3], Float64[])

        config = EcoEvoConfig(
            ecoDyn = lotkaVolterra(growthFn, kernelFn),
            mutationGenerator = generateMutant(invaderPopsize=0.001, variance=0.01),
            integrationParams = IntegrationParams(maxTime = 10.0, abstol=1e-6, reltol=1e-4),
            extThreshold = 0.003
        )

        result = ecoDyn(community, config)
        @test result isa Community

        # unstructuredModel examples
        r_fn(z) = 1 - sum(z.^2)
        α_fn(zi, zj) = exp(-sum((zi .- zj).^2) / 0.04)

        ecology_unstr = unstructuredModel() do i, n, z, nSpecies
            n[i] * (r_fn(z[i]) - sum(α_fn(z[i], z[j]) * n[j] for j in 1:nSpecies))
        end
        @test ecology_unstr isa Function

        ecology_unstr_pre = unstructuredModel(
            precompute = (z, nSpecies) -> (
                b = [r_fn(z[i]) for i in 1:nSpecies],
                A = [α_fn(z[i], z[j]) for i in 1:nSpecies, j in 1:nSpecies]
            )
        ) do i, n, z, nSpecies, pre
            n[i] * (pre.b[i] - sum(pre.A[i, j] * n[j] for j in 1:nSpecies))
        end
        @test ecology_unstr_pre isa Function

        # structuredModel examples (requires Distributions)
        d_val = 1.0; mu_val = 0.1; alpha_val = 1.0
        y_val = [d_val / 2, -d_val / 2]

        ecology_str = structuredModel() do i, j, N, z, R, nSpecies, nPatches
            growth = pdf(Normal(0, 1), z[i][1] - y_val[j])
            dd = alpha_val * sum(N[k, j] for k in 1:nSpecies)
            (growth - dd - mu_val) * N[i, j] + mu_val * sum(N[i, k] for k in 1:nPatches if k != j)
        end
        @test ecology_str isa Function

        ecology_str_pre = structuredModel(
            precompute = (z, nSpecies, nPatches) ->
                [pdf(Normal(0, 1), z[i][1] - y_val[j]) for i in 1:nSpecies, j in 1:nPatches]
        ) do i, j, N, z, R, nSpecies, nPatches, pre
            dd = alpha_val * sum(N[k, j] for k in 1:nSpecies)
            (pre[i, j] - dd - mu_val) * N[i, j] + mu_val * sum(N[i, k] for k in 1:nPatches if k != j)
        end
        @test ecology_str_pre isa Function
    end

    @testset "Visualization" begin
        # niceTickInterval examples
        @test niceTickInterval(1500) == 200
        @test niceTickInterval(95000) == 20000  # Actual result is 20000, not 10000

        # plotEvo examples (just verify it doesn't error, don't display)
        growthFn(z) = 1.0 - sum(z.^2)
        kernelFn(z_i, z_j) = -exp(-sum((z_i .- z_j).^2) / 0.15^2)

        config = EcoEvoConfig(
            ecoDyn = lotkaVolterra(growthFn, kernelFn),
            mutationGenerator = generateMutant(invaderPopsize=0.001, variance=0.01),
            integrationParams = IntegrationParams(maxTime = 10.0, abstol=1e-6, reltol=1e-4),
            extThreshold = 0.003
        )

        Random.seed!(42)
        comm = Community([1.0, 1.0], [0.1, 0.2], Float64[])
        history = evolve(comm, config, 3, showProgress=false)

        @test_nowarn plotEvo(history)
    end

    @testset "History to Table" begin
        # historyToTable examples
        growthFn(z) = 1.0 - sum(z.^2)
        kernelFn(z_i, z_j) = -exp(-sum((z_i .- z_j).^2) / 0.15^2)

        config = EcoEvoConfig(
            ecoDyn = lotkaVolterra(growthFn, kernelFn),
            mutationGenerator = generateMutant(invaderPopsize=0.001, variance=0.01),
            integrationParams = IntegrationParams(maxTime = 10.0, abstol=1e-6, reltol=1e-4),
            extThreshold = 0.003
        )

        Random.seed!(42)
        comm = Community([1.0, 1.0], [0.1, 0.2], Float64[])
        history = evolve(comm, config, 3, showProgress=false)

        table = historyToTable(history)
        @test haskey(table, "mutNo")
        @test haskey(table, "time")
        @test haskey(table, "species")
        @test haskey(table, "popsize_1")
        @test haskey(table, "trait_1")
    end

end
