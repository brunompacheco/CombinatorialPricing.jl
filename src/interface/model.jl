"""
    toll_bounds(prob::PricingProblem) -> Vector{Float64}

Get the upper bounds of the toll prices. The entries corresponding to the toll-free items are ignored.
"""
function toll_bounds end

"""
    follower_A(prob::PricingProblem) -> Matrix{Float64}

Get of the follower's constraint matrix. This matrix represents the coefficients of the items in the follower's problem, which is typically a knapsack or similar combinatorial problem.
"""
function follower_A end

"""
    follower_b(prob::PricingProblem) -> Vector{Float64}

The right-hand side of the follower's constraints. This is typically the capacity of the knapsack or similar combinatorial problem.
"""
function follower_b end
