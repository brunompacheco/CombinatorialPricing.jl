struct KnapsackPricing <: PricingProblem
    weights::Vector{Int}
    capacity::Int
    base_values::Vector{Float64}
    tolled::BitSet
end

# Interface
num_items(prob::KnapsackPricing) = length(prob.weights)
tolled(prob::KnapsackPricing) = prob.tolled
toll_free(prob::KnapsackPricing) = setdiff(BitSet(1:num_items(prob)), tolled(prob))
base_costs(prob::KnapsackPricing) = -prob.base_values

base_values(prob::KnapsackPricing) = prob.base_values
weights(prob::KnapsackPricing) = prob.weights
capacity(prob::KnapsackPricing) = prob.capacity
density(prob::KnapsackPricing) = prob.base_values ./ prob.weights

follower_A(prob::KnapsackPricing) = weights(prob)'
follower_b(prob::KnapsackPricing) = capacity(prob)

toll_bounds(prob::KnapsackPricing) = prob.base_values

# Pretty print
function Base.show(io::IO, prob::KnapsackPricing)
    i1 = length(prob.tolled)
    ni = num_items(prob)
    print(io, "$(typeof(prob)) with $ni items ($i1 tolled)")
end
