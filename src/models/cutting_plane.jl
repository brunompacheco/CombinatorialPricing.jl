function add_cutting_plane_callback!(cutgen::Function, model; threads=nothing, reset_follower=false)
    model[:follower] = follower_model(model[:prob]; silent=true, threads)
    model[:vf_x] = Set{Vector{Int}}()
    model[:callback_time] = 0.
    model[:subproblem_time] = 0.
    model[:subproblem_times] = Float64[]
    model[:callback_calls] = 0

    provider = haskey(model, :heur_provider) ? model[:heur_provider] : nothing

    # If enabled, reset the warmstart solution of the subproblem
    reset_fn = () -> nothing
    if reset_follower
        # Only applied to Gurobi optimizer
        # reset_fn = () -> GRBreset(unsafe_backend(model[:follower]), 0)
    end

    function lazy_callback(cb_data)
        status = callback_node_status(cb_data, model)
        if status != MOI.CALLBACK_NODE_STATUS_INTEGER
            return
        end

        model[:callback_time] += @elapsed begin
            subproblem_time = @elapsed begin
                # Solve the follower's problem with the given toll
                reset_fn()
                x̂, t̂, follower_obj, current_obj = solve_callback_solution(model, cb_data)
            end
            model[:callback_calls] += 1
            model[:subproblem_time] += subproblem_time
            push!(model[:subproblem_times], subproblem_time)
            not_optimal = follower_obj < current_obj - 1e-4
            
            @debug "Callback $(round.((follower_obj, current_obj), digits=1)) $(not_optimal ? " " : "*") $(round(Int, subproblem_time * 1000)) ms"
    
            # Add a constraint if the current_obj is not optimal
            if not_optimal
                cutgen(cb_data, x̂)
                push!(model[:vf_x], x̂)
            end
    
            # Log the current solution as heuristic
            isnothing(provider) || push!(provider, x̂, t̂)
        end
    end

    MOI.set(model, MOI.LazyConstraintCallback(), lazy_callback)
end

function solve_callback_solution(model, cb_data)
    t, g, follower, prob = model[:t], model[:g], model[:follower], model[:prob]

    t̂ = callback_value.(cb_data, t)
    set_toll!(follower, prob, t̂)
    optimize!(follower)

    follower_obj = objective_value(follower)
    current_obj = callback_value(cb_data, g)

    x̂ = round.(Int, value.(follower[:x]))

    return x̂, t̂, follower_obj, current_obj
end
