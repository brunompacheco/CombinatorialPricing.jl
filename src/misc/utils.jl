expand_lift_type(::Type{VariableRef}) = AffExpr
expand_lift_type(::Type{<:Number}) = Float64

function expand_t(t::DenseAxisArray{T}, prob::AbstractProblem) where T
    L = expand_lift_type(T)
    n, i1 = num_items(prob), collect(tolled(prob))
    tfull = Vector{L}(zeros(n))
    tfull[i1] .= t.data
    return tfull
end

expand_t(t::AbstractVector, ::AbstractProblem) = t

make_ct(t, prob::PricingProblem) = base_costs(prob) .+ expand_t(t, prob)

convert_x_to_array(x::Vector{Bool}) = findall(x)
convert_x_to_array(x::Vector{Int}) = findall(==(1), x)
convert_x_to_array(x::Vector{Float64}) = findall(>(0.5), x)

convert_x_to_set(x) = BitSet(convert_x_to_array(x))

convert_set_to_x(s::BitSet, n::Int) = [i in s for i in 1:n]

union_sets!(target, sets) = reduce(union!, sets; init=target)
union_sets(sets) = union_sets!(BitSet(), sets)

function some(predicate, itr)
    for a in itr
        predicate(a) && return a
    end
    return nothing
end
