using JuMP, CombinatorialPricing, Printf, Random

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

    """ FIND NEW CUT THROUGH SEPARATION PROBLEM """
    
    # problem parameters

    S = zeros(Int, num_items(prob), nv(sd))  # state matrix
    for (j, v) in enumerate(collect(vertices(sd)))
        for i in v[2].selected
            S[i,j] = 1
        end
    end

    y = zeros(nv(sd)) # dual variables vector

    for (j, v) in enumerate(vertices(sd))
        y[j] = -value(model_sd[:y][v])
    end

    # build separation problem
    model_sep = CombinatorialPricing.blank_model(silent=true)

    @variable(model_sep, u[1:nv(sd)], Bin)
    @variable(model_sep, ell[1:num_items(prob)], Bin)

    @constraint(model_sep, prob.weights' * (S * u .+ ell) <= prob.capacity)
    @constraint(model_sep, S * u .+ ell .<= 1)
    @constraint(model_sep, sum(u) == 1)
    @constraint(model_sep, u[end] == 0)

    ct = base_costs(prob) .+ CombinatorialPricing.expand_t(t_hat, prob)
    @objective(model_sep, Min, y' * u + ct' * ell)

    set_silent(model_sep)
    optimize!(model_sep)

    ### STOPPING CRITERION ###
    if objective_value(model_sep) >= -1e-4
        println("Stopping criterion met")
        break
    end
    ### STOPPING CRITERION ###

    # build arc
    u_star = collect(vertices(sd))[findfirst(round.(Bool, value.(u)))]
    arc_label = Set(findall(round.(Bool, value.(ell))))
    arc_separation = SDArc(u_star, sink_node(sd), arc_label)

    # solution quality
    m = Model()
    @variable(m, x[i=1:num_items(prob)], Bin)
    c = base_costs(prob)
    t = CombinatorialPricing.expand_t(value.(t_hat), prob)
    f_obj_func = (c + t)' * x
    function get_x_index(var)
        m = match(r"\[(\d+)\]", name(var))
        return parse(Int, m.captures[1])
    end
    response_obj = value(x_i -> get_x_index(x_i) in arc_separation.src[2].selected ∪ arc_separation.action, f_obj_func)

    # violation
    cut_violation = -value(model_sd[:y][arc_separation.src]) + sum(ct[i] for i in arc_separation.action)

    # log
    if iter == 1
        @printf("| %4s | %9s | %12s | %10s |\n", "Iter", "Bound", "Follower Obj", "Violation")
        @printf("|%s|%s|%s|%s|\n", "-"^6, "-"^11, "-"^14, "-"^12)
    end
    @printf("| %4d | % 9.3f | % 12.3f | % 10.3f |\n", iter, objective_value(model_sd), -response_obj, -cut_violation)

    push!(sd.arcs, arc_separation)
end
