using Test
using EcoEvoSim.EcoEvoCore
using .TestUtils: generateCommunity, communityGen

@testset "addSpecies! / removeSpecies!" begin
    st = generateCommunity(traitDim=2, stageCls=1, nAux=2, nSpecies=2)

    @test nSpecies(st) == 2
    @test totalBiomass(st) ≈ 15.0  # 10 + 5 (default numbers from utils)

    newSp = makeSpecies([1.0, 0.7], [2.0])
    addSpecies!(st, newSp)

    @test nSpecies(st) == 3
    @test totalBiomass(st) ≈ 17.0
    @test traits(st.species[3]) == @SVector [1.0, 0.7]

    removed = removeSpecies!(st, 2)
    @test removed === st.species[2]   # the removed object is returned
    @test nSpecies(st) == 2
    @test totalBiomass(st) ≈ 12.0   # 10 + 2 after removal
end
