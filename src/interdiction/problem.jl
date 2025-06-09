struct KnapsackInterdiction <: AbstractProblem
    lower_cap::Int
    upper_cap::Int
    lower_weights::Vector{Int}
    upper_weights::Vector{Int}
    profits::Vector{Int}
end

# Interface
num_items(prob::KnapsackInterdiction) = length(prob.lower_weights)

lower_cap(prob::KnapsackInterdiction) = prob.lower_cap
upper_cap(prob::KnapsackInterdiction) = prob.upper_cap

lower_weights(prob::KnapsackInterdiction) = prob.lower_weights
upper_weights(prob::KnapsackInterdiction) = prob.upper_weights

profits(prob::KnapsackInterdiction) = prob.profits
lower_density(prob::KnapsackInterdiction) = profits(prob) ./ lower_weights(prob)
upper_density(prob::KnapsackInterdiction) = profits(prob) ./ upper_weights(prob)

follower_A(prob::KnapsackInterdiction) = lower_weights(prob)'
follower_b(prob::KnapsackInterdiction) = lower_cap(prob)
follower_c(prob::KnapsackInterdiction) = -profits(prob)

function ct(t, prob::KnapsackInterdiction)
    p = profits(prob)
    t_full = expand_t(t, prob)
    return -p .* (1 .- t_full)
end

# Read from file
function Base.read(filename::AbstractString, ::Type{KnapsackInterdiction})
    lines = readlines(filename)
    CL = parse(Int, lines[2])
    CU = parse(Int, lines[3])
    wL = parse.(Int, split(lines[4]))
    wU = parse.(Int, split(lines[5]))
    p = parse.(Int, split(lines[6]))
    return KnapsackInterdiction(CL, CU, wL, wU, p)
end

# Pretty print
function Base.show(io::IO, prob::KnapsackInterdiction)
    ni = num_items(prob)
    CL = lower_cap(prob)
    CU = upper_cap(prob)
    print(io, "$(typeof(prob)) with (n = $ni, C_L = $CL, C_U = $CU)")
end
