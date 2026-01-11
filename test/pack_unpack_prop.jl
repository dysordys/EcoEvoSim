using Test
using .EcoEvoCore
using .TestUtils: communityGen, sampleCommunity
using PropCheck
using PropCheck: @propcheck


function prop_pack_unpack_roundtrip(; ntrials=100)
    for _ in 1:ntrials
        # Draw a random Community (sampled without PropCheck shrinking)
        st = sampleCommunity()

        u, p = EcoEvoCore.packState(st)
        st2 = EcoEvoCore.unpackState(u, st)

        if !(st.time == st2.time &&
             all(abs.(EcoEvoCore.traits(st2) .- EcoEvoCore.traits(st)) .< 1e-12) &&
             EcoEvoCore.totalBiomass(st2) ≈ EcoEvoCore.totalBiomass(st))
            return false
        end
    end
    true
end

@testset "PropCheck round‑trip for pack/unpack" begin
    # Run the property with the default number of trials (≈100)
    @test prop_pack_unpack_roundtrip()
end
