"""
The supertype of all problems (including interdiction).
"""
abstract type AbstractProblem end

"""
The supertype of all pricing problems.
"""
abstract type PricingProblem <: AbstractProblem end

## Interface
"""
    num_items(prob::AbstractProblem) -> Int

Get the total number of items selectable by the followers (including tolled and toll-free items).

All items are indexed by the range `1:num_items(prob)`.
"""
function num_items end

"""
    tolled(prob::PricingProblem) -> BitSet

Get the set of tolled items.
"""
function tolled end

"""
    toll_free(prob::PricingProblem) -> BitSet

Get the set of toll-free items.
"""
function toll_free end

"""
    base_costs(prob::PricingProblem)::Vector{Float64}

Get the base costs of all selectable items.
They are the actual costs of the toll-free items.
The actual cost of a tolled item is the sum of its base cost and its toll price.
"""
function base_costs end

"""
    generate(type::Type{<:PricingProblem}, num_items::Int; kwargs...)

Generate a random pricing problem of the given `type` with `num_items` selectable items.
"""
function generate end

## Write to file
function Base.write(filename::AbstractString, prob::PricingProblem)
    open(filename, "w") do f
        JSON.print(f, Dict(:problem => prob))
    end
end

## Read from file
function _json_problem_value(::Type{T}, value) where {T}
    if T <: BitSet
        return BitSet(Int.(value))
    elseif T <: AbstractVector
        return [_json_problem_value(eltype(T), item) for item in value]
    elseif T <: Tuple
        return T(ntuple(i -> _json_problem_value(fieldtype(T, i), value[i]), fieldcount(T)))
    else
        return convert(T, value)
    end
end

function Base.read(filename::AbstractString, P::Type{<:PricingProblem})
    str = read(filename, String)
    data = JSON.parse(str)["problem"]
    values = (_json_problem_value(fieldtype(P, field), data[string(field)]) for field in fieldnames(P))
    return P(values...)
end
