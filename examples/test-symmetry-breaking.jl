using EcoEvoSim
using Statistics
using Random

# Test if the asymmetry is random (sometimes dim1, sometimes dim2) or systematic

println("Testing symmetry breaking across multiple replicates...")
println("=" ^ 70)

results = []

for replicate in 1:10
    Random.seed!(replicate * 1000)  # Different seed for each replicate

    community = Community([Species(1.0, [0.0, 0.0])], PopulationSize{Float64}[], 0.0)

    growthFn = z -> 1.0 - sum(z.^2)
    kernelFn = (zi, zj) -> -exp(-sum((zi .- zj).^2) / 0.25^2)

    config = EcoEvoConfig(
        ecoDyn = lotkaVolterra(growthFn, kernelFn),
        mutationGenerator = (c, cfg) -> generateMutant(c, cfg, 0.002^2),
        integrationParams = IntegrationParams(maxTime = 1.0e8),
        invaderPopsize = 0.001,
        extThreshold = 0.003
    )

    history = evolve!(community, config, 100, showProgress=false)

    # Calculate standard deviations
    allTraits1 = [traits(comm, i)[1] for comm in history.history for i in 1:numSpecies(comm)]
    allTraits2 = [traits(comm, i)[2] for comm in history.history for i in 1:numSpecies(comm)]

    std1 = std(allTraits1)
    std2 = std(allTraits2)

    ratio = std1 / std2

    push!(results, (seed=replicate*1000, std1=std1, std2=std2, ratio=ratio))

    println("Replicate $replicate: Std1 = $(round(std1, digits=4)), Std2 = $(round(std2, digits=4)), Ratio = $(round(ratio, digits=2))")
end

println("=" ^ 70)
println("\nSummary:")
ratios = [r.ratio for r in results]
println("Mean ratio (Std1/Std2): $(round(mean(ratios), digits=2))")
println("If mean ratio ≈ 1: asymmetry is random (sometimes dim1, sometimes dim2)")
println("If mean ratio >> 1: dimension 1 always evolves more")
println("If mean ratio << 1: dimension 2 always evolves more")
println("\nIndividual ratios: $(round.(ratios, digits=2))")
