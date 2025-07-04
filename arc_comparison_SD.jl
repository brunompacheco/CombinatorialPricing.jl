using JuMP, CombinatorialPricing, Printf, Random, CSV, DataFrames, Glob, SCIP

const numpairs = 1000
const maxiter = 25

Random.seed!(42)

function solve_best_response(prob, t_hat, sd)
    model_f = follower_model(prob, silent=true)
    set_toll!(model_f, prob, value.(t_hat))

    # solve the follower's problem with fixed leader assignment (t_hat)
    set_silent(model_f)
    optimize!(model_f)

    # follower's best response
    x_br = value.(model_f[:x])

    # compute deepest arc
    x_br_set = CombinatorialPricing.convert_x_to_set(x_br)

    # worst-case = root-to-terminal arc
    u_star = source_node(sd)
    ell_star = x_br_set

    added = false
    for l in 2:-1:1
        layer = sd.layers[l]
        for s in layer
            (s.selected ⊆ x_br_set) || continue
            ell_star = setdiff(x_br_set, s.selected)
            u_star = (l, s)
            added = true
            break
        end
        added && break
    end

    return u_star, ell_star, model_f
end

function solve_separation(prob, t_hat, y_hat, sd)
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

    u_star = collect(vertices(sd))[findfirst(round.(Bool, value.(u)))]
    # ell_star = round.(Bool, value.(ell))
    ell_star = CombinatorialPricing.convert_x_to_set(round.(Int, value.(ell)))
    
    return u_star, ell_star, model_sep
end

function solve_separation_extended(prob, t_hat, y_hat, sd; relaxed=false)
    # problem parameters
    θ = Dict()
    for v in vertices(sd)
        θ[v] = -y_hat[v]
    end

    labels = Dict(a => CombinatorialPricing.convert_set_to_x(a.action, num_items(prob)) for a in sd.arcs)
    regular_nodes = [u for u in vertices(sd) if u != sink_node(sd) && u != source_node(sd)]

    # build separation problem
    model_sep2 = Model(SCIP.Optimizer)

    if relaxed
        @variable(model_sep2, x[1:num_items(prob)])
    else
        @variable(model_sep2, x[1:num_items(prob)], Bin)
    end
    @variable(model_sep2, υ[sd.arcs] .>= 0)
    @variable(model_sep2, ι[1:num_items(prob)], Bin)

    selected_arc_label = sum(labels[a] .* υ[a] for a in sd.arcs if a.dst == sink_node(sd))
    @constraint(model_sep2, CombinatorialPricing.follower_A(prob) * (x - selected_arc_label + ι) <= CombinatorialPricing.follower_b(prob))
    @constraint(model_sep2, x - selected_arc_label + ι .<= 1)

    # SD primal
    @constraint(model_sep2, sum(υ[a] for a in sd.arcs if a.src == source_node(sd)) == 1)
    @constraint(model_sep2, sum(υ[a] for a in sd.arcs if a.dst == sink_node(sd)) == 1)
    @constraint(model_sep2, [u in regular_nodes], sum(υ[a] for a in sd.arcs if a.src == u) == sum(υ[a] for a in sd.arcs if a.dst == u))
    @constraint(model_sep2, x .== sum(labels[a] .* υ[a] for a in sd.arcs))

    @objective(model_sep2, Min, sum(θ[a.src] .* υ[a] for a in sd.arcs if a.dst == sink_node(sd)) + CombinatorialPricing.ct(t_hat, prob)' * ι)

    set_silent(model_sep2)
    optimize!(model_sep2)

    selected_target_arcs = [a for a in sd.arcs if (a.dst == sink_node(sd)) && (value(υ[a]) > 1e-3)]
    @assert relaxed || (length(selected_target_arcs) == 1) "Expected exactly one target arc, but found $(length(selected_target_arcs))"
    u_star = first(selected_target_arcs).src
    ι_star = CombinatorialPricing.convert_x_to_set(round.(Int, value.(ι)))
    # ι_star = Set(findall(round.(Bool, value.(ι))))

    return u_star, ι_star, model_sep2
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
    if method == :best_response
        log_data = Vector{NamedTuple{(:iter, :dual_bound, :solve_time), Tuple{Int, Float64, Float64}}}()
    elseif method == :separation
        log_data = Vector{NamedTuple{(:iter, :dual_bound, :solve_time, :solve_time_extended, :solve_time_extended_relaxed, :violation, :violation_extended, :violation_extended_relaxed), Tuple{Int, Float64, Float64, Float64, Float64, Float64, Float64, Float64}}}()
    else
        error("Unknown method: $method")
    end

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
        y_hat = value.(model_sd[:y])

        """ FIND NEW CUT """
        if method == :best_response
            u_star, ell_star, model_cut = solve_best_response(prob, t_hat, sd)
        elseif method == :separation
            u_star, ell_star, model_cut = solve_separation(prob, t_hat, y_hat, sd)
            u_star_extended, ell_star_extended, model_cut_extended = solve_separation_extended(prob, t_hat, y_hat, sd, relaxed=false)
            u_star_relaxed, ell_star_relaxed, model_cut_relaxed = solve_separation_extended(prob, t_hat, y_hat, sd, relaxed=true)
        end

        new_arc = SDArc(u_star, sink_node(sd), ell_star)

        # solution quality
        m = Model()
        @variable(m, x[i=1:num_items(prob)], Bin)
        f_obj_func = CombinatorialPricing.ct(t_hat, prob)' * x
        function get_x_index(var)
            m = match(r"\[(\d+)\]", name(var))
            return parse(Int, m.captures[1])
        end
        response_obj = value(x_i -> get_x_index(x_i) in new_arc.src[2].selected ∪ new_arc.action, f_obj_func)

        ### STOPPING CRITERION ###
        m = Model()
        @variable(m, x[i=1:num_items(prob)], Bin)
        f_obj_func = CombinatorialPricing.ct(t_hat, prob)' * x
        inp = Dict(i => value(x_hat[i]) for i in 1:num_items(prob))
        follower_obj_hat = value(x_i -> inp[get_x_index(x_i)], f_obj_func)

        if abs(response_obj - follower_obj_hat) < 1e-4
            # println("Stopping criterion met")
            break
        end
        ### STOPPING CRITERION ###

        Base.GC.gc()

        # log
        if iter == 1
            @printf("| %4s | %9s | %12s |\n", "Iter", "Bound", "Follower Obj")
            @printf("|%s|%s|%s|\n", "-"^6, "-"^11, "-"^14)
        end
        @printf("| %4d | % 9.3f | % 12.3f |\n", iter, objective_value(model_sd), -response_obj)

        # Save log data
        if method == :best_response
            push!(log_data, (
                iter=iter,
                dual_bound=objective_value(model_sd),
                solve_time=solve_time(model_sd),
            ))
        elseif method == :separation
            push!(log_data, (
                iter=iter,
                dual_bound=objective_value(model_sd),
                solve_time=solve_time(model_cut),
                solve_time_extended=solve_time(model_cut_extended),
                solve_time_extended_relaxed=solve_time(model_cut_relaxed),
                violation=objective_value(model_cut),
                violation_extended=objective_value(model_cut_extended),
                violation_extended_relaxed=objective_value(model_cut_relaxed),
            ))
        end

        push!(sd.arcs, new_arc)
    end

    return DataFrame(log_data)
end

results_dir = "results_SD/"
for arg in ARGS
    if endswith(arg, ".json") || endswith(arg, ".ki")
        println("Solving instance $arg")

        instance = splitext(basename(arg))[1]

        for method in [:best_response, :separation]
            df_fpath = results_dir * instance * "_" * string(method) * ".csv"
            if isfile(df_fpath) && method != :separation
                println("Results already exist for $method. Skipping...")
            else
                df = run_test(arg, method)
                CSV.write(df_fpath, df)
            end
        end
    end
end
