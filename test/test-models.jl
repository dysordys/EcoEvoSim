using Test
using EcoEvoSim
using LinearAlgebra


@testset "tests_of_models" begin

    @testset "testing_lotkaVolterra_basic_functionality" begin
        # Define simple growth and interaction functions
        growthFn = traits -> 1.0 .- traits.^2
        kernelFn = (z_i, z_j) -> -exp(-((z_i - z_j) / 0.15)^2)

        # Create a community
        comm = Community([1.0, 1.0, 1.0], [-0.2, 0.0, 0.3], Float64[])

        # Create the dynamics function
        ecoDynFn = lotkaVolterra(growthFn, kernelFn)(comm)

        # Check that it returns a function
        @test ecoDynFn isa Function

        # Test that the dynamics function has the correct signature
        u = [1.0, 1.0, 1.0]
        p = nothing
        t = 0.0
        dudt = ecoDynFn(u, p, t)
        @test dudt isa Vector{Float64}
        @test length(dudt) == 3

        # Test that dynamics makes sense: with positive growth rates and competition,
        # rates should be finite
        @test all(isfinite.(dudt))
    end


    @testset "testing_lotkaVolterra_with_single_species" begin
        # Single species with positive growth
        growthFn = z -> 0.5
        kernelFn = (z_i, z_j) -> -1.0

        comm = Community([1.0], [0.0], Float64[])
        ecoDynFn = lotkaVolterra(growthFn, kernelFn)(comm)

        # At low density, should grow approximately at intrinsic rate
        u = [0.1]
        dudt = ecoDynFn(u, nothing, 0.0)
        # dN/dt = N * (b - a*N) = 0.1 * (0.5 - 1.0*0.1) = 0.04
        @test dudt[1] ≈ 0.04 rtol=1e-10

        # At equilibrium (N* = b/a = 0.5), should have zero growth
        u = [0.5]
        dudt = ecoDynFn(u, nothing, 0.0)
        @test dudt[1] ≈ 0.0 atol=1e-10
    end


    @testset "testing_lotkaVolterra_two_species_competition" begin
        # Two species with identical traits (strong competition)
        growthFn = z -> 1.0
        kernelFn = (z_i, z_j) -> -1.0  # Equal competition

        comm = Community([0.5, 0.5], [0.0, 0.0], Float64[])
        ecoDynFn = lotkaVolterra(growthFn, kernelFn)(comm)

        u = [0.5, 0.5]
        dudt = ecoDynFn(u, nothing, 0.0)

        # Both should have same dynamics: N * (1 - 0.5 - 0.5) = 0
        @test dudt[1] ≈ 0.0 atol=1e-10
        @test dudt[2] ≈ 0.0 atol=1e-10
    end


    @testset "testing_lotkaVolterra_trait_dependent_growth" begin
        # Growth rate depends on trait value
        growthFn = z -> 1.0 - z^2
        kernelFn = (z_i, z_j) -> -0.5  # Weak competition

        comm = Community([1.0, 1.0], [0.0, 0.5], Float64[])
        ecoDynFn = lotkaVolterra(growthFn, kernelFn)(comm)

        u = [0.1, 0.1]
        dudt = ecoDynFn(u, nothing, 0.0)

        # Species 1 (trait 0.0): b = 1.0, rate = 0.1 * (1.0 - 0.5*0.1 - 0.5*0.1) = 0.09
        # Species 2 (trait 0.5): b = 0.75, rate = 0.1 * (0.75 - 0.5*0.1 - 0.5*0.1) = 0.065
        @test dudt[1] ≈ 0.09 rtol=1e-10
        @test dudt[2] ≈ 0.065 rtol=1e-10

        # Species with trait 0 should have higher growth rate
        @test dudt[1] > dudt[2]
    end


    @testset "testing_lotkaVolterra_trait_dependent_competition" begin
        # Competition depends on trait distance (niche differentiation)
        growthFn = z -> 1.0
        competitionWidth = 0.2
        kernelFn = (z_i, z_j) -> -exp(-((z_i - z_j) / competitionWidth)^2)

        # Two species with different traits
        comm = Community([1.0, 1.0], [-0.3, 0.3], Float64[])
        ecoDynFn = lotkaVolterra(growthFn, kernelFn)(comm)

        u = [0.5, 0.5]
        dudt = ecoDynFn(u, nothing, 0.0)

        # Competition coefficient between species:
        # a_12 = exp(-((0.3 - (-0.3))/0.2)^2) = exp(-(0.6/0.2)^2) = exp(-9) ≈ 1.234e-4
        # This is very small, so species experience weak interspecific competition
        a_12 = exp(-((0.3 - (-0.3)) / 0.2)^2)

        # Intraspecific competition is 1.0 (same trait)
        # Rate = N * (1 - 1*N_same - a_12*N_other)
        expected_rate = 0.5 * (1.0 - 1.0*0.5 - a_12*0.5)
        @test dudt[1] ≈ expected_rate rtol=1e-6
        @test dudt[2] ≈ expected_rate rtol=1e-6
    end


    @testset "testing_lotkaVolterra_integration_with_EcoEvoConfig" begin
        # Test full integration with EcoEvoConfig
        growthFn = z -> 1.0 - z^2
        kernelFn = (z_i, z_j) -> -exp(-((z_i - z_j) / 0.15)^2)

        comm = Community([1.0, 1.0, 1.0], [-0.2, 0.0, 0.3], Float64[])
        ecoDynFactory = lotkaVolterra(growthFn, kernelFn)

        # Create a valid EcoEvoConfig
        config = EcoEvoConfig(
            ecoDyn = ecoDynFactory,
            mutationGenerator = (x) -> x .+ 0.01,
            integrationParams = IntegrationParams(maxTime = 10.0),
            invaderPopsize = 0.001,
            extThreshold = 0.003
        )

        @test config.ecoDyn === ecoDynFactory
        @test config.integrationParams.maxTime == 10.0

        # Run dynamics
        finalComm = ecoDyn(comm, config)
        @test finalComm isa Community
        @test numSpecies(finalComm) <= numSpecies(comm)  # May have extinctions
        @test finalComm.time > comm.time
    end


    @testset "testing_lotkaVolterra_with_keyword_constructor" begin
        # Test using keyword constructor for EcoEvoConfig
        growthFn = z -> 1.0 - z^2
        kernelFn = (z_i, z_j) -> -exp(-((z_i - z_j) / 0.15)^2)

        comm = Community([1.0, 1.0], [0.0, 0.5], Float64[])
        ecoDynFactory = lotkaVolterra(growthFn, kernelFn)

        # Use keyword constructor
        config = EcoEvoConfig(
            ecoDyn = ecoDynFactory,
            mutationGenerator = (x) -> x .+ 0.01,
            integrationParams = IntegrationParams(maxTime = 10.0),
            invaderPopsize = 0.001,
            extThreshold = 0.003
        )

        @test config.ecoDyn === ecoDynFactory
        @test config.invaderPopsize == 0.001
        @test config.extThreshold == 0.003
    end


    @testset "testing_lotkaVolterra_different_numeric_types" begin
        # Test with Float32
        growthFn = z -> Float32(1.0 - z^2)
        kernelFn = (z_i, z_j) -> -Float32(exp(-((z_i - z_j) / 0.15)^2))

        comm = Community(Float32[1.0, 1.0], Float32[0.0, 0.5], Float32[])
        ecoDynFn = lotkaVolterra(growthFn, kernelFn)(comm)

        u = Float32[0.5, 0.5]
        dudt = ecoDynFn(u, nothing, 0.0)
        @test eltype(dudt) == Float32
        @test all(isfinite.(dudt))
    end


    @testset "testing_lotkaVolterra_zero_population" begin
        # Test behavior with zero populations
        growthFn = z -> 1.0
        kernelFn = (z_i, z_j) -> -1.0

        comm = Community([1.0, 1.0], [0.0, 0.5], Float64[])
        ecoDynFn = lotkaVolterra(growthFn, kernelFn)(comm)

        # One species at zero
        u = [0.0, 0.5]
        dudt = ecoDynFn(u, nothing, 0.0)
        @test dudt[1] == 0.0  # Zero population stays zero
        @test dudt[2] != 0.0  # Other species has dynamics
    end

end
