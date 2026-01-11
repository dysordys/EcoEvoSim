# examples/usage_example.jl
# Small demo that exercises the primary EcoEvoCore APIs.

using EcoEvoSim.EcoEvoCore

# Create two species (trait dim = 1, stage classes = 2)
sp1 = makeSpecies([1.0], [10.0, 5.0])
sp2 = makeSpecies([2.0], [3.0, 7.0])

# Build a SystemState (explicit time to avoid makeSystemState default-time corner case)
aux = Vector{Population{Float64,0}}()  # typed empty aux vector
st = makeSystemState([sp1, sp2], aux, 0.0)

println("Initial state:")
println(st)
println("Total biomass: ", totalBiomass(st))

# Pack/unpack roundtrip
u, p = packState(st)
println("Packed state vector u = ", u)
println("Traits matrix p.traits = \n", p.traits)

st2 = unpackState(u, st)
println("Unpacked state:")
println(st2)

# Compare trait matrices and biomass to ensure pack/unpack preserved data
trmat1 = hcat([traits(s) for s in st.species]...)
trmat2 = hcat([traits(s) for s in st2.species]...)
@assert all(abs.(trmat1 .- trmat2) .< 1e-12)
@assert totalBiomass(st2) ≈ totalBiomass(st)

# Mutate the packed vector and apply with updateState!
newu = copy(u)
newu[1] -= 1.0  # reduce first species first stage
updateState!(st, newu)
println("After updateState!: ")
println(st)

# Add and remove species
addSpecies!(st, makeSpecies([0.5], [1.0, 1.0]))
println("After addSpecies!: ")
println(st)

removeSpecies!(st, 2)  # remove the previously second species
println("After removeSpecies!(idx=2): ")
println(st)

# Create a small EvoHistory and push the current state
h = makeEvoHistory(Float64, Val{1}(), Val{2}(), Val{0}())
push!(h.states, st)
push!(h.times, st.time)
println("EvoHistory entries: ", length(h.states))

println("Demo finished successfully.")
