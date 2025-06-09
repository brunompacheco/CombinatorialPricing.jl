using JuMP, CombinatorialPricing, Printf, Random, CSV, DataFrames, Glob, SCIP

const numpairs = 1000
const maxiter = 25

Random.seed!(42)

function find_best_response_arc(prob, t_hat, y_hat, sd)
    model_f = follower_model(prob, silent=true)
    set_toll!(model_f, prob, value.(t_hat))

    # solve the follower's problem with fixed leader assignment (t_hat)
    set_silent(model_f)
    optimize!(model_f)

    best_response_obj = objective_value(model_f)

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

    x_sep = CombinatorialPricing.convert_set_to_x(arc_best_response.action, num_items(prob))
    cut_violation = -y_hat[arc_best_response.src] + CombinatorialPricing.ct(t_hat, prob)' * x_sep

    return arc_best_response, cut_violation, best_response_obj
end

function find_separation_arc(prob, t_hat, y_hat, sd)
    # problem parameters
    S = zeros(Int, num_items(prob), nv(sd))  # state matrix
    for (j, v) in enumerate(collect(vertices(sd)))
        for i in v[2].selected
            S[i,j] = 1
        end
    end

    y = zeros(nv(sd)) # dual variables vector

    for (j, v) in enumerate(vertices(sd))
        y[j] = -y_hat[v]
    end

    # build separation problem
    # model_sep = CombinatorialPricing.blank_model(silent=true)
    model_sep = Model(SCIP.Optimizer)
    set_optimizer_attribute(model_sep, MOI.Silent(), true)

    @variable(model_sep, u[1:nv(sd)], Bin)
    @variable(model_sep, ell[1:num_items(prob)], Bin)
    @variable(model_sep, x[1:num_items(prob)], Bin)

    @constraint(model_sep, x .== S * u .+ ell)

    add_primal!(model_sep, prob)

    @constraint(model_sep, sum(u) == 1)
    @constraint(model_sep, u[end] == 0)  # not allowed to select terminal node

    @objective(model_sep, Min, y' * u + CombinatorialPricing.ct(t_hat, prob)' * ell)

    set_silent(model_sep)
    optimize!(model_sep)

    # build arc
    u_star = collect(vertices(sd))[findfirst(round.(Bool, value.(u)))]
    arc_label = Set(findall(round.(Bool, value.(ell))))
    arc_separation = SDArc(u_star, sink_node(sd), arc_label)

    # solution quality
    m = Model()
    @variable(m, x[i=1:num_items(prob)], Bin)
    f_obj_func = CombinatorialPricing.ct(t_hat, prob)' * x
    function get_x_index(var)
        m = match(r"\[(\d+)\]", name(var))
        return parse(Int, m.captures[1])
    end
    response_obj = value(x_i -> get_x_index(x_i) in arc_separation.src[2].selected ∪ arc_separation.action, f_obj_func)

    # violation
    x_sep = CombinatorialPricing.convert_set_to_x(arc_separation.action, num_items(prob))
    cut_violation = -y_hat[arc_separation.src] + CombinatorialPricing.ct(t_hat, prob)' * x_sep

    return arc_separation, cut_violation, response_obj
end

function run_test(file, method = :best_response)
    Random.seed!(42)

    # load problem
    if occursin("knapsack", file)
        prob = read(file, KnapsackPricing)
        sampler = MaximalKnapsackSampler(prob)
    elseif occursin("maxstab", file)
        prob = read(file, MaxStableSetPricing)
        sampler = MaximalStableSetSampler(prob)
    elseif occursin("mincover", file)
        prob = read(file, MinSetCoverPricing)
        sampler = MinimalSetCoverSampler(prob)
    elseif occursin("interdiction-flatten", file)
        prob = read(file, KnapsackInterdiction)
        sampler = MaximalKnapsackInterdictionSampler(prob)
    else
        error("Unknown problem type in file: $file")
    end

    """ BUILD INITIAL DIAGRAM """
    # Sample some solutions
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

    # Prepare to collect log data
    log_data = Vector{NamedTuple{(:iter, :dual_bound, :follower_obj, :violation), Tuple{Int, Float64, Float64, Float64}}}()

    for iter in 1:maxiter
        """ FIND CANDIDATE SOLUTION """
        # Create baseline SD model
        model_sd = base_model(prob, silent=true)
        CombinatorialPricing.add_sdgraph_dual!(model_sd, sd)

        # Solve the relaxation
        set_silent(model_sd)
        optimize!(model_sd)

        if ~is_solved_and_feasible(model_sd)
            println("Instance $file is not feasible. Stopping...")
            break
        end

        # candidate solution
        x_hat = value.(model_sd[:x])
        t_hat = value.(model_sd[:t])

        """ FIND NEW CUT """
        if method == :best_response
            new_arc, cut_violation, response_obj = find_best_response_arc(prob, t_hat, value.(model_sd[:y]), sd)
        elseif method == :separation
            new_arc, cut_violation, response_obj = find_separation_arc(prob, t_hat, value.(model_sd[:y]), sd)
        else
            error("Unknown method: $method")
        end

        ### STOPPING CRITERION ###
        m = Model()
        @variable(m, x[i=1:num_items(prob)], Bin)
        f_obj_func = CombinatorialPricing.ct(t_hat, prob)' * x
        function get_x_index(var)
            m = match(r"\[(\d+)\]", name(var))
            return parse(Int, m.captures[1])
        end
        inp = Dict(i => value(x_hat[i]) for i in 1:num_items(prob))
        follower_obj_hat = value(x_i -> inp[get_x_index(x_i)], f_obj_func)

        if abs(response_obj - follower_obj_hat) < 1e-4
            # println("Stopping criterion met")
            break
        end
        ### STOPPING CRITERION ###

        Base.GC.gc()

        # # log
        # if iter == 1
        #     @printf("| %4s | %9s | %12s | %10s |\n", "Iter", "Bound", "Follower Obj", "Violation")
        #     @printf("|%s|%s|%s|%s|\n", "-"^6, "-"^11, "-"^14, "-"^12)
        # end
        # @printf("| %4d | % 9.3f | % 12.3f | % 10.3f |\n", iter, objective_value(model_sd), -response_obj, -cut_violation)

        # Save log data
        push!(log_data, (iter=iter, dual_bound=objective_value(model_sd), follower_obj=-response_obj, violation=-cut_violation))

        push!(sd.arcs, new_arc)
    end

    return DataFrame(log_data)
end

results_dir = "results/"
for arg in ARGS
    if endswith(arg, ".json")
        println("Solving instance $arg")

        instance = splitext(basename(arg))[1]

        for method in [:best_response, :separation]
            # println("Method: $method")

            df_fpath = results_dir * instance * "_" * string(method) * ".csv"
            if isfile(df_fpath)
                # println("Results already exist. Skipping...")
            else
                df = run_test(arg, method)
                CSV.write(df_fpath, df)
            end
        end
    end
end
