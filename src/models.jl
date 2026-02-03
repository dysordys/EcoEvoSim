# Pre-defined ecological models


"""
    lotkaVolterra(growthFn, kernelFn)

Factory function for creating a multispecies Lotka-Volterra competition model.

Returns a function that takes a `Community` and produces the ecological dynamics
function (`ecoDyn`) for that community.

# Arguments
- `growthFn`: Function that takes a trait vector (for one species) and returns
  a scalar intrinsic growth rate. Signature: `Vector{T} -> T`
- `kernelFn`: Function that takes two trait vectors and returns the
  interaction coefficient. Signature: `(Vector{T}, Vector{T}) -> T`

# Returns
A function with signature `Community -> Function` that takes a community and
returns the dynamics function for the ODE solver.

# Model
The Lotka-Volterra model has the form:
```
dN_i/dt = N_i * (b_i - sum_j a_ij * N_j)
```
where:
- `N_i` is the population size of species i
- `b_i` is the intrinsic growth rate (from `growthFn`)
- `a_ij` is the interaction coefficient between species i and j (from `kernelFn`)

# Example
```julia
using EcoEvoSim

# Define growth rate as a function of trait
growthFn = traits -> 1.0 - sum(traits.^2)

# Define competition based on trait distance with Gaussian kernel
kernelFn = (z_i, z_j) -> -exp(-sum((z_i .- z_j).^2) / 0.15^2)

# Create a community
community = Community([1.0, 1.0, 1.0], [-0.2, 0.0, 0.3], Float64[])

# Create configuration with desired settings
config = EcoEvoConfig(
    ecoDyn = lotkaVolterra(growthFn, kernelFn),
    mutationGenerator = (x) -> x .+ 0.01,
    integrationParams = IntegrationParams(maxTime = 50.0),
    invaderPopsize = 0.001,
    extThreshold = 0.003
)

# Run dynamics
finalCommunity = ecoDyn(community, config)
```
"""


function lotkaVolterra(growthFn, kernelFn)
    # Return a function that creates the dynamics function for a given community
    function createDynamics(community::Community{T, AuxClasses}) where {T<:Real, AuxClasses}
        # Extract traits from the community
        nSpecies = numSpecies(community)
        spTraits = [traits(community, i) for i in 1:nSpecies]

        # Compute intrinsic growth rates using the provided function
        b = [growthFn(spTraits[i]) for i in 1:nSpecies]

        # Precompute interaction coefficient matrix using the provided function
        a = zeros(T, nSpecies, nSpecies)
        for i in 1:nSpecies
            for j in 1:nSpecies
                a[i, j] = kernelFn(spTraits[i], spTraits[j])
            end
        end

        # Return the dynamics function with proper signature for ODEProblem
        # Note: We can safely close over a and b because this function is recreated
        # for each new community in ecoDyn (which calls config.ecoDyn(community))
        return (u, p, t) -> u .* (b .+ (a * u))
    end

    return createDynamics
end


# User-friendly ecological model specification via real Julia syntax

"""
    @unstructuredModel block

Create an unstructured ecological model (no population structure) from dynamics equations.

The macro captures equations specifying how populations change over time. Variables are
indexed by species: `n[i]` for population of species i, `z[i]` for trait of species i.
Parameters should be passed via closure.

# Arguments
- `block`: A `begin...end` block (or single expression) containing one assignment:
  - `dn[i] = <expression>`: The per-capita growth rate expression for species i

The expression can reference:
- `n[i]` and `n[j]`: Population sizes of species (only `n` is the state variable)
- `z[i]` and `z[j]`: Trait vectors of species i and j
- `nSpecies`: Number of species (automatically available)
- Loop variables like `j` in comprehensions
- External parameters (captured in closure)

# Returns
A function with signature `Community -> (u, p, t) -> du` suitable for `EcoEvoConfig.ecoDyn`.

# Example
```julia
# Define intrinsic growth rate and interaction kernel as closures
r(z) = 1 - sum(z.^2) / 0.25
alpha(z_i, z_j) = exp(-sum((z_i .- z_j).^2) / 0.04)

# Create model
ecology = @unstructuredModel begin
    dn[i] = n[i] * (r(z[i]) - sum(alpha(z[i], z[j]) * n[j] for j in 1:nSpecies))
end

# Use in config
config = EcoEvoConfig(
    ecoDyn = ecology,
    mutationGenerator = ...,
    integrationParams = ...,
    extThreshold = ...
)
```
"""
macro unstructuredModel(block::Expr)
    # Parse the block to extract the equation
    # We expect either a single assignment or a begin...end block with one assignment

    eqs = Dict{Symbol, Expr}()

    if block.head == :block
        # Multiple expressions; find the assignment
        for expr in block.args
            if expr isa LineNumberNode
                continue
            elseif expr isa Expr && expr.head == :(=)
                # It's an assignment: dn[i] = rhs
                lhs = expr.args[1]
                rhs = expr.args[2]

                # Extract variable name from lhs
                if lhs isa Expr && lhs.head == :ref
                    # dn[i]
                    varname = lhs.args[1]::Symbol
                    eqs[varname] = rhs
                else
                    error("Expected equation of form dn[i] = ..., got $(lhs)")
                end
            else
                error("Expected assignment in @unstructuredModel, got $(expr.head)")
            end
        end
    elseif block.head == :(=)
        # Single assignment
        lhs = block.args[1]
        rhs = block.args[2]

        if lhs isa Expr && lhs.head == :ref
            varname = lhs.args[1]::Symbol
            eqs[varname] = rhs
        else
            error("Expected equation of form dn[i] = ..., got $(lhs)")
        end
    else
        error("@unstructuredModel expects assignment(s), got $(block.head)")
    end

    # Should have exactly one equation (dn)
    length(eqs) == 1 || error("@unstructuredModel expects exactly one equation (dn[i] = ...)")

    dn_expr = eqs[:dn]

    # Generate code that creates the factory function
    # Key optimization: Generate the actual ODE function code at macro time,
    # not at runtime, so we don't pay eval() costs in the integration loop

    # Build the ODE function body: a loop that unrolls to compute du[i] for each species
    # by substituting i into the dn_expr
    ode_body = Expr(:block)
    push!(ode_body.args, :(n = @view u[1:nSpecies]))
    push!(ode_body.args, :(du = similar(u)))

    # Generate code for each species index
    # We unroll the loop at macro time to avoid needing i as a variable
    loop = Expr(:for, :(i = 1:nSpecies), Expr(:block))
    push!(loop.args[2].args, :(du[i] = $(dn_expr)))
    push!(ode_body.args, loop)
    push!(ode_body.args, :(return du))

    # Generate the factory function that captures z and creates the ODE function
    # Use esc() so the generated code runs in the caller's module context
    # where parameters like r, alpha, etc. are defined
    esc(quote
        function _factory(community::Community{T, AuxClasses}) where {T<:Real, AuxClasses}
            nSpecies = numSpecies(community)
            z = [traits(community, i) for i in 1:nSpecies]

            function _ode_fn(u::Vector{T}, p, t) where T
                $(ode_body)
            end

            return _ode_fn
        end
        _factory
    end)
end
