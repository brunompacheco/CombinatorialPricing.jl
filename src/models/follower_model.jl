function add_primal!(model::Model, prob::AbstractProblem)
    x = model[:x]
    @constraint(model, primal, follower_A(prob) .≤ follower_b(prob))
    return
end

function follower_model(prob::AbstractProblem; silent=false, threads=nothing)
    model = blank_model(; silent, threads)
    model[:prob] = prob

    n = num_items(prob)

    @variable(model, x[i=1:n], Bin)

    # Follower objective
    @objective(model, Min, follower_c(prob)' * x)

    # Stable set constraints
    add_primal!(model, prob)
    
    return model
end

function set_toll!(model::Model, prob::AbstractProblem, toll)
    x = model[:x]
    @objective(model, Min, ct(toll, prob)' * x)
end
