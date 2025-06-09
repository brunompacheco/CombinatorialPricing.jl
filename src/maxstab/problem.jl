struct MaxStableSetPricing <: PricingProblem
    vertices::Int
    edges::Vector{Tuple{Int,Int}}
    base_values::Vector{Float64}
    tolled::BitSet
end

# Interface
num_items(prob::MaxStableSetPricing) = prob.vertices
tolled(prob::MaxStableSetPricing) = prob.tolled
toll_free(prob::MaxStableSetPricing) = setdiff(BitSet(1:num_items(prob)), tolled(prob))
base_costs(prob::MaxStableSetPricing) = -prob.base_values

base_values(prob::MaxStableSetPricing) = prob.base_values

# Build graph (the result is cached, do not modify directly on the returned graph)
@memoize function graph(prob::MaxStableSetPricing)
    g = SimpleGraph(prob.vertices)
    for e in prob.edges
        add_edge!(g, e[1], e[2])
    end
    return g
end

follower_A(prob::MaxStableSetPricing) = incidence_matrix(graph(prob))
follower_b(prob::MaxStableSetPricing) = ones(ne(graph(prob)))

toll_bounds(prob::MaxStableSetPricing) = prob.base_values

# Pretty print
function Base.show(io::IO, prob::MaxStableSetPricing)
    i1 = length(prob.tolled)
    nv = prob.vertices
    ne = length(prob.edges)
    print(io, "$(typeof(prob)) with $nv vertices ($i1 tolled), $ne edges")
end
