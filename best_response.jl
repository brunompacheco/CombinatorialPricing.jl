using JuMP, Gurobi, CombinatorialPricing, Printf, Random

Random.seed!(42)

# Import a problem from a file
file = "./problems/knapsack/expset-2/kpp2-n30-01.json"
prob = read(file, KnapsackPricing)


""" BUILD INITIAL DIAGRAM """
numpairs = 100

# Sample some solutions
sampler = MaximalKnapsackSampler(prob)
samples = rand(sampler, numpairs)

# Extract random pairs from the samples
pairs = random_pair.(samples)

# Create a selection diagram
sd = sdgraph_from_pairs(prob, pairs)

# Add empty arcs from last layer to terminal node (necessary for separation problem)
l = length(sd.layers) - 1
for u in sd.layers[l]
    empty_arc = SDArc((l, u), sink_node(sd), DPAction())
    push!(sd.arcs, empty_arc)
end

# println("| Iter |   Bound   | Follower Obj | Violation |")
# println("|------|-----------|--------------|-----------|")

for iter in 1:100
    """ FIND CANDIDATE SOLUTION """
    # Create baseline SD model
    model_sd = base_model(prob, silent=true)
    CombinatorialPricing.add_sdgraph_dual!(model_sd, sd)

    # Solve the relaxation
    set_silent(model_sd)
    optimize!(model_sd)

    # candidate solution
    x_hat = value.(model_sd[:x])
    t_hat = value.(model_sd[:t])

    """ FIND NEW CUT THROUGH BEST RESPONSE """
    model_f = follower_model(prob, silent=true)

    set_toll!(model_f, prob, value.(t_hat))

    inp = Dict(model_f[:x][i] => value(x_hat[i]) for i in 1:num_items(prob))
    follower_obj_hat = value(i -> inp[i], objective_function(model_f))

    # solve the follower's problem with fixed leader assignment (t_hat)
    set_silent(model_f)
    optimize!(model_f)

    best_response_obj = objective_value(model_f)

    ### STOPPING CRITERION ###
    if abs(best_response_obj - follower_obj_hat) < 1e-4
        println("Stopping criterion met")
        break
    end
    ### STOPPING CRITERION ###

    # follower's best response
    x_br = value.(model_f[:x])

    # compute deepest arc
    x_br_set = CombinatorialPricing.convert_x_to_set(x_br)

    # worst-case = root-to-terminal arc
    arc_best_response = SDArc(source_node(sd), sink_node(sd), x_br_set)

    added = false
    for l in 2:-1:1
        layer = sd.layers[l]
        for s in layer
            (s.selected ⊆ x_br_set) || continue
            arc_label = setdiff(x_br_set, s.selected)
            arc_best_response = SDArc((l, s), sink_node(sd), arc_label)
            added = true
            break
        end
        added && break
    end

    ct = base_costs(prob) .+ CombinatorialPricing.expand_t(t_hat, prob)
    cut_violation = -value(model_sd[:y][arc_best_response.src]) + sum(ct[i] for i in arc_best_response.action)

    # log
    if iter == 1
        @printf("| %4s | %9s | %12s | %10s |\n", "Iter", "Bound", "Follower Obj", "Violation")
        @printf("|%s|%s|%s|%s|\n", "-"^6, "-"^11, "-"^14, "-"^12)
    end
    @printf("| %4d | % 9.3f | % 12.3f | % 10.3f |\n", iter, objective_value(model_sd), -best_response_obj, -cut_violation)

    push!(sd.arcs, arc_best_response)
end
