using Test
using EcoEvoSim
using LinearAlgebra


@testset "tests_of_models" begin

    @testset "testing_lotkaVolterra_basic_functionality" begin
        # Define simple growth and interaction functions
        growthFn = traits -> 1.0 - sum(traits.^2)
        kernelFn = (z_i, z_j) -> -exp(-sum((z_i .- z_j).^2) / 0.15^2)

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
        growthFn = z -> 1.0 - sum(z.^2)
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
        kernelFn = (z_i, z_j) -> -exp(-sum((z_i .- z_j).^2) / competitionWidth^2)

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
        growthFn = z -> 1.0 - sum(z.^2)
        kernelFn = (z_i, z_j) -> -exp(-sum((z_i .- z_j).^2) / 0.15^2)

        comm = Community([1.0, 1.0, 1.0], [-0.2, 0.0, 0.3], Float64[])
        ecoDynFactory = lotkaVolterra(growthFn, kernelFn)

        # Create a valid EcoEvoConfig
        config = EcoEvoConfig(
            ecoDyn = ecoDynFactory,
            mutationGenerator = (x) -> x .+ 0.01,
            integrationParams = IntegrationParams(maxTime = 10.0),
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
        growthFn = z -> 1.0 - sum(z.^2)
        kernelFn = (z_i, z_j) -> -exp(-sum((z_i .- z_j).^2) / 0.15^2)

        comm = Community([1.0, 1.0], [0.0, 0.5], Float64[])
        ecoDynFactory = lotkaVolterra(growthFn, kernelFn)

        # Use keyword constructor
        config = EcoEvoConfig(
            ecoDyn = ecoDynFactory,
            mutationGenerator = (x) -> x .+ 0.01,
            integrationParams = IntegrationParams(maxTime = 10.0),
            extThreshold = 0.003
        )

        @test config.ecoDyn === ecoDynFactory
        @test config.extThreshold == 0.003
    end


    @testset "testing_lotkaVolterra_different_numeric_types" begin
        # Test with Float32
        growthFn = z -> Float32(1.0 - sum(z.^2))
        kernelFn = (z_i, z_j) -> -Float32(exp(-sum((z_i .- z_j).^2) / 0.15^2))

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


@testset "tests_of_unstructuredModel" begin

    @testset "basic_functionality" begin
        r(z) = 1.0 - sum(z.^2)
        α(zi, zj) = exp(-sum((zi .- zj).^2) / 0.04)

        ecology = unstructuredModel() do i, n, z, nSpecies
            n[i] * (r(z[i]) - sum(α(z[i], z[j]) * n[j] for j in 1:nSpecies))
        end

        comm = Community([1.0, 1.0, 1.0], [-0.2, 0.0, 0.3], Float64[])
        ode_fn = ecology(comm)

        @test ode_fn isa Function

        u = [1.0, 1.0, 1.0]
        dudt = ode_fn(u, nothing, 0.0)
        @test dudt isa Vector{Float64}
        @test length(dudt) == 3
        @test all(isfinite.(dudt))
    end


    @testset "single_species_logistic" begin
        # dn/dt = n * r * (1 - n/K)
        r_val = 1.0; K = 10.0

        ecology = unstructuredModel() do i, n, z, nSpecies
            n[i] * r_val * (1.0 - n[i] / K)
        end

        comm = Community([5.0], [0.0], Float64[])
        ode_fn = ecology(comm)

        u = [5.0]
        dudt = ode_fn(u, nothing, 0.0)
        @test dudt[1] ≈ 5.0 * 1.0 * (1.0 - 5.0/10.0) rtol=1e-10

        # At carrying capacity, growth should be zero
        u_eq = [K]
        dudt_eq = ode_fn(u_eq, nothing, 0.0)
        @test dudt_eq[1] ≈ 0.0 atol=1e-10
    end


    @testset "matches_lotkaVolterra" begin
        # unstructuredModel should match lotkaVolterra for the same equations
        growthFn = z -> 1.0 - sum(z.^2)
        kernelFn = (zi, zj) -> -exp(-sum((zi .- zj).^2) / 0.15^2)

        lv_factory = lotkaVolterra(growthFn, kernelFn)

        ecology = unstructuredModel() do i, n, z, nSpecies
            n[i] * (growthFn(z[i]) + sum(kernelFn(z[i], z[j]) * n[j] for j in 1:nSpecies))
        end

        comm = Community([1.0, 1.0, 1.0], [-0.2, 0.0, 0.3], Float64[])
        u = unpackCommunity(comm)

        dudt_lv = lv_factory(comm)(u, nothing, 0.0)
        dudt_fn = ecology(comm)(u, nothing, 0.0)

        @test dudt_fn ≈ dudt_lv rtol=1e-10
    end


    @testset "zero_population" begin
        ecology = unstructuredModel() do i, n, z, nSpecies
            n[i] * (1.0 - n[i])
        end

        comm = Community([0.0, 0.5], [0.0, 0.0], Float64[])
        ode_fn = ecology(comm)

        dudt = ode_fn([0.0, 0.5], nothing, 0.0)
        @test dudt[1] == 0.0
        @test dudt[2] ≈ 0.5 * 0.5 rtol=1e-10
    end


    @testset "integration_with_EcoEvoConfig" begin
        ecology = unstructuredModel() do i, n, z, nSpecies
            n[i] * (1.0 - sum(n[j] for j in 1:nSpecies))
        end

        comm = Community([0.1, 0.2], [0.0, 0.5], Float64[])
        config = EcoEvoConfig(
            ecoDyn = ecology,
            mutationGenerator = (x) -> x .+ 0.01,
            integrationParams = IntegrationParams(maxTime = 10.0),
            extThreshold = 0.003
        )

        finalComm = ecoDyn(comm, config)
        @test finalComm isa Community
        @test finalComm.time > comm.time
    end


    @testset "precompute_matches_no_precompute" begin
        growthFn = z -> 1.0 - sum(z.^2)
        kernelFn = (zi, zj) -> -exp(-sum((zi .- zj).^2) / 0.15^2)

        # Without precompute
        ecology_plain = unstructuredModel() do i, n, z, nSpecies
            n[i] * (growthFn(z[i]) + sum(kernelFn(z[i], z[j]) * n[j] for j in 1:nSpecies))
        end

        # With precompute
        ecology_pre = unstructuredModel(
            precompute = (z, nSpecies) -> (
                b = [growthFn(z[i]) for i in 1:nSpecies],
                A = [kernelFn(z[i], z[j]) for i in 1:nSpecies, j in 1:nSpecies]
            )
        ) do i, n, z, nSpecies, pre
            n[i] * (pre.b[i] + sum(pre.A[i, j] * n[j] for j in 1:nSpecies))
        end

        comm = Community([1.0, 1.0, 1.0], [-0.2, 0.0, 0.3], Float64[])
        u = unpackCommunity(comm)

        dudt_plain = ecology_plain(comm)(u, nothing, 0.0)
        dudt_pre = ecology_pre(comm)(u, nothing, 0.0)

        @test dudt_pre ≈ dudt_plain rtol=1e-10
    end


    @testset "precompute_matches_lotkaVolterra" begin
        growthFn = z -> 1.0 - sum(z.^2)
        kernelFn = (zi, zj) -> -exp(-sum((zi .- zj).^2) / 0.15^2)

        lv_factory = lotkaVolterra(growthFn, kernelFn)

        ecology_pre = unstructuredModel(
            precompute = (z, nSpecies) -> (
                b = [growthFn(z[i]) for i in 1:nSpecies],
                A = [kernelFn(z[i], z[j]) for i in 1:nSpecies, j in 1:nSpecies]
            )
        ) do i, n, z, nSpecies, pre
            n[i] * (pre.b[i] + sum(pre.A[i, j] * n[j] for j in 1:nSpecies))
        end

        comm = Community([1.0, 1.0, 1.0], [-0.2, 0.0, 0.3], Float64[])
        u = unpackCommunity(comm)

        dudt_lv = lv_factory(comm)(u, nothing, 0.0)
        dudt_pre = ecology_pre(comm)(u, nothing, 0.0)

        @test dudt_pre ≈ dudt_lv rtol=1e-10
    end


    @testset "with_auxDynamics" begin
        eta = 1.0; K = 10.0

        ecology = unstructuredModel(
            auxDynamics = (R, n, z, nSpecies) ->
                [eta * (K - R[1]) - sum(n[i] * R[1] for i in 1:nSpecies)]
        ) do i, n, z, nSpecies
            n[i] * (R_unused = 0.0; 1.0 - n[i])  # simple logistic, ignore R here
        end

        # Need R accessible in the do block — rewrite using the full signature
        ecology = unstructuredModel(
            auxDynamics = (R, n, z, nSpecies) ->
                [eta * (K - R[1]) - sum(n[i] * R[1] for i in 1:nSpecies)]
        ) do i, n, z, nSpecies
            n[i] * (1.0 - n[i])
        end

        comm = Community([0.5, 0.3], [0.0, 1.0], [5.0])
        ode_fn = ecology(comm)

        u = unpackCommunity(comm)
        dudt = ode_fn(u, nothing, 0.0)

        # 2 species + 1 aux
        @test length(dudt) == 3

        # Species dynamics: n[i]*(1 - n[i])
        @test dudt[1] ≈ 0.5 * (1.0 - 0.5) rtol=1e-10
        @test dudt[2] ≈ 0.3 * (1.0 - 0.3) rtol=1e-10

        # Aux dynamics: eta*(K - R) - (n1 + n2)*R = 1*(10-5) - (0.5+0.3)*5 = 5 - 4 = 1
        @test dudt[3] ≈ 1.0 * (10.0 - 5.0) - (0.5 + 0.3) * 5.0 rtol=1e-10
    end


    @testset "precompute_with_auxDynamics" begin
        eta = 1.0; K = 10.0

        ecology = unstructuredModel(
            auxDynamics = (R, n, z, nSpecies, pre) ->
                [eta * (K - R[1]) - sum(pre.a[i] * n[i] * R[1] for i in 1:nSpecies)],
            precompute = (z, nSpecies) -> (
                a = [1.0 + z[i][1]^2 for i in 1:nSpecies],
            )
        ) do i, n, z, nSpecies, pre
            pre.a[i] * n[i] * (1.0 - n[i])
        end

        comm = Community([0.5, 0.3], [0.0, 1.0], [5.0])
        ode_fn = ecology(comm)

        u = unpackCommunity(comm)
        dudt = ode_fn(u, nothing, 0.0)

        @test length(dudt) == 3

        # pre.a = [1 + 0^2, 1 + 1^2] = [1.0, 2.0]
        # Species: pre.a[i]*n[i]*(1-n[i])
        @test dudt[1] ≈ 1.0 * 0.5 * 0.5 rtol=1e-10       # a=1.0
        @test dudt[2] ≈ 2.0 * 0.3 * 0.7 rtol=1e-10       # a=2.0

        # Aux: eta*(K-R) - sum(a[i]*n[i]*R) = 1*(10-5) - (1*0.5*5 + 2*0.3*5) = 5 - 5.5 = -0.5
        @test dudt[3] ≈ 1.0 * (10.0 - 5.0) - (1.0*0.5*5.0 + 2.0*0.3*5.0) rtol=1e-10
    end

end


@testset "tests_of_structuredModel" begin

    @testset "basic_two_patch" begin
        mu = 0.1; alpha_val = 1.0
        y = [0.5, -0.5]

        ecology = structuredModel() do i, j, N, z, R, nSpecies, nPatches
            growth = exp(-0.5 * (z[i][1] - y[j])^2)
            dd = alpha_val * sum(N[k, j] for k in 1:nSpecies)
            (growth - dd - mu) * N[i, j] + mu * sum(N[i, k] for k in 1:nPatches if k != j)
        end

        comm = Community([1.0 1.0;], [-0.2])
        ode_fn = ecology(comm)

        @test ode_fn isa Function

        u = unpackCommunity(comm)
        dudt = ode_fn(u, nothing, 0.0)
        @test length(dudt) == 2
        @test all(isfinite.(dudt))
    end


    @testset "two_species_two_patches" begin
        mu = 0.0  # no migration for simpler test

        ecology = structuredModel() do i, j, N, z, R, nSpecies, nPatches
            (1.0 - N[i, j]) * N[i, j]  # logistic per patch, independent
        end

        # Two species, each with [patch1 patch2] populations
        comm = Community([0.5 0.5; 0.2 0.8], [0.0, 1.0])
        ode_fn = ecology(comm)

        u = unpackCommunity(comm)
        dudt = ode_fn(u, nothing, 0.0)
        @test length(dudt) == 4

        # Species 1 patch 1: 0.5*(1-0.5) = 0.25
        @test dudt[1] ≈ 0.25 rtol=1e-10
        # Species 1 patch 2: 0.5*(1-0.5) = 0.25
        @test dudt[2] ≈ 0.25 rtol=1e-10
        # Species 2 patch 1: 0.2*(1-0.2) = 0.16
        @test dudt[3] ≈ 0.16 rtol=1e-10
        # Species 2 patch 2: 0.8*(1-0.8) = 0.16
        @test dudt[4] ≈ 0.16 rtol=1e-10
    end


    @testset "migration_conservation" begin
        # With only migration (no growth), total population should be conserved
        mu = 0.5

        ecology = structuredModel() do i, j, N, z, R, nSpecies, nPatches
            -mu * N[i, j] + mu * sum(N[i, k] for k in 1:nPatches if k != j) / (nPatches - 1)
        end

        comm = Community([3.0 1.0;], [0.0])
        ode_fn = ecology(comm)

        u = unpackCommunity(comm)
        dudt = ode_fn(u, nothing, 0.0)

        # Total derivative should sum to zero (conservation)
        @test sum(dudt) ≈ 0.0 atol=1e-10
    end


    @testset "with_auxDynamics" begin
        eta = 1.0; chi = 1.0; gamma = 1.0; beta = 1.0

        ecology = structuredModel(
            auxDynamics = (R, N, z, nSpecies, nPatches) ->
                [R[k] * (eta - chi * R[k]) - gamma * sum(N[i, k] for i in 1:nSpecies) * R[k]
                 for k in 1:nPatches]
        ) do i, j, N, z, R, nSpecies, nPatches
            (beta * R[j] - 1.0) * N[i, j]
        end

        comm = Community([1.0 1.0;], [0.0], [2.0, 3.0])
        ode_fn = ecology(comm)

        u = unpackCommunity(comm)
        dudt = ode_fn(u, nothing, 0.0)

        # 2 species state entries + 2 aux entries
        @test length(dudt) == 4

        # Species dynamics: (beta*R[j] - 1)*N[i,j]
        @test dudt[1] ≈ (1.0 * 2.0 - 1.0) * 1.0 rtol=1e-10  # patch 1
        @test dudt[2] ≈ (1.0 * 3.0 - 1.0) * 1.0 rtol=1e-10  # patch 2

        # Resource dynamics: R*(eta - chi*R) - gamma*total_consumers*R
        # Patch 1: R=2, consumers=1 → 2*(1-2) - 1*1*2 = -4
        @test dudt[3] ≈ 2.0 * (1.0 - 1.0*2.0) - 1.0*1.0*2.0 rtol=1e-10
        # Patch 2: R=3, consumers=1 → 3*(1-3) - 1*1*3 = -9
        @test dudt[4] ≈ 3.0 * (1.0 - 1.0*3.0) - 1.0*1.0*3.0 rtol=1e-10
    end


    @testset "no_auxDynamics_with_aux_vars" begin
        # Aux vars present but no auxDynamics — aux portion of du is untouched
        ecology = structuredModel() do i, j, N, z, R, nSpecies, nPatches
            (1.0 - N[i, j]) * N[i, j]
        end

        comm = Community([0.5 0.5;], [0.0], [1.0, 1.0])
        ode_fn = ecology(comm)

        u = unpackCommunity(comm)
        dudt = ode_fn(u, nothing, 0.0)
        @test length(dudt) == 4
        # Species dynamics still work
        @test dudt[1] ≈ 0.25 rtol=1e-10
        @test dudt[2] ≈ 0.25 rtol=1e-10
    end


    @testset "integration_with_EcoEvoConfig" begin
        mu = 0.1

        ecology = structuredModel() do i, j, N, z, R, nSpecies, nPatches
            growth = exp(-0.5 * (z[i][1] - 0.0)^2)
            (growth - mu) * N[i, j] + mu * sum(N[i, k] for k in 1:nPatches if k != j)
        end

        comm = Community([1.0 1.0;], [-0.2])
        config = EcoEvoConfig(
            ecoDyn = ecology,
            mutationGenerator =
                generateMutantSpatial(invaderPopsize=0.001, variance=0.003^2),
            integrationParams = IntegrationParams(maxTime = 10.0),
            extThreshold = 0.003
        )

        finalComm = ecoDyn(comm, config)
        @test finalComm isa Community
        @test finalComm.time > comm.time
    end


    @testset "precompute_matches_no_precompute" begin
        mu = 0.1; alpha_val = 1.0
        y = [0.5, -0.5]

        # Without precompute
        ecology_plain = structuredModel() do i, j, N, z, R, nSpecies, nPatches
            growth = exp(-0.5 * (z[i][1] - y[j])^2)
            dd = alpha_val * sum(N[k, j] for k in 1:nSpecies)
            (growth - dd - mu) * N[i, j] + mu * sum(N[i, k] for k in 1:nPatches if k != j)
        end

        # With precompute
        ecology_pre = structuredModel(
            precompute = (z, nSpecies, nPatches) ->
                [exp(-0.5 * (z[i][1] - y[j])^2) for i in 1:nSpecies, j in 1:nPatches]
        ) do i, j, N, z, R, nSpecies, nPatches, pre
            dd = alpha_val * sum(N[k, j] for k in 1:nSpecies)
            (pre[i, j] - dd - mu) * N[i, j] + mu * sum(N[i, k] for k in 1:nPatches if k != j)
        end

        comm = Community([1.0 1.0; 0.5 0.5], [-0.2, 0.3])
        u = unpackCommunity(comm)

        dudt_plain = ecology_plain(comm)(u, nothing, 0.0)
        dudt_pre = ecology_pre(comm)(u, nothing, 0.0)

        @test dudt_pre ≈ dudt_plain rtol=1e-10
    end


    @testset "precompute_with_auxDynamics" begin
        eta = 1.0; chi = 1.0; gamma = 1.0; beta = 1.0

        ecology = structuredModel(
            auxDynamics = (R, N, z, nSpecies, nPatches, pre) ->
                [R[k] * (eta - chi * R[k]) - gamma * sum(N[i, k] for i in 1:nSpecies) * R[k]
                 for k in 1:nPatches],
            precompute = (z, nSpecies, nPatches) ->
                (mort = [0.5 * z[i][1]^2 for i in 1:nSpecies],)
        ) do i, j, N, z, R, nSpecies, nPatches, pre
            (beta * R[j] - pre.mort[i]) * N[i, j]
        end

        comm = Community([1.0 1.0;], [0.0], [2.0, 3.0])
        ode_fn = ecology(comm)

        u = unpackCommunity(comm)
        dudt = ode_fn(u, nothing, 0.0)

        @test length(dudt) == 4
        # Species: (beta*R[j] - 0.5*z^2)*N = (1*R[j] - 0)*1
        @test dudt[1] ≈ (1.0 * 2.0 - 0.0) * 1.0 rtol=1e-10
        @test dudt[2] ≈ (1.0 * 3.0 - 0.0) * 1.0 rtol=1e-10
        # Aux still works
        @test dudt[3] ≈ 2.0 * (1.0 - 1.0*2.0) - 1.0*1.0*2.0 rtol=1e-10
        @test dudt[4] ≈ 3.0 * (1.0 - 1.0*3.0) - 1.0*1.0*3.0 rtol=1e-10
    end

end
