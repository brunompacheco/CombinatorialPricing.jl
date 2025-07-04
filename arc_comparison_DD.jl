using JuMP, CombinatorialPricing, Printf, Random, CSV, DataFrames, Glob, SCIP

const ddwidth = 50
const grouping = 2
const maxiter = 25

Random.seed!(42)

function solve_best_response(prob, t_hat, dd)
    model_f = follower_model(prob, silent=true)
    set_toll!(model_f, prob, value.(t_hat))

    # solve the follower's problem with fixed leader assignment (t_hat)
    set_silent(model_f)
    optimize!(model_f)

    # follower's best response
    x_br = value.(model_f[:x])

    # compute path over DD nodes that implements x_br 
    x_br_set = CombinatorialPricing.convert_x_to_set(x_br)
    path = structured_path(dd, x_br_set)

    # new arcs
    new_arcs = [a for a in path if ~(a in dd.arcs)]

    return new_arcs
end

function solve_separation_uv(prob, t_hat, θ_hat, dd)
    # problem parameters
    I = zeros(Int, num_items(prob), nv(dd))  # indicator of unassigned items
    θ = zeros(nv(dd))  # dual variables vector
    s = zeros(nv(dd), length(CombinatorialPricing.follower_b(prob)))  # state of the vertices
    for (j, (i, s_v)) in enumerate(collect(vertices(dd)))
        I[(i*grouping+1):end, j] .= 1
        θ[j] = θ_hat[(i, s_v)]
        # TODO: this only works for knapsack! How to handle the state of the other problems?
        s[j, :] .= s_v.remaining
    end

    # build separation problem
    # model_sep = CombinatorialPricing.blank_model(silent=true)
    model_sep = Model(SCIP.Optimizer)
    set_optimizer_attribute(model_sep, MOI.Silent(), true)

    @variable(model_sep, u[1:nv(dd)], Bin)
    @variable(model_sep, v[1:nv(dd)], Bin)
    @variable(model_sep, l[1:num_items(prob)], Bin)

    @constraint(model_sep, l .<= I * u - I * v)
    @constraint(model_sep, s' * v .== s' * u .+ CombinatorialPricing.follower_A(prob) * l)
    @constraint(model_sep, sum((I * u - I * v)[i] for i in 1:num_items(prob)) .>= 1)

    @constraint(model_sep, sum(u) == 1)
    @constraint(model_sep, sum(v) == 1)
    @constraint(model_sep, u[end] == 0)  # not allowed to select terminal node
    @constraint(model_sep, v[end] == 0)  # not allowed to select terminal node

    @objective(model_sep, Min, θ' * u - θ' * v + CombinatorialPricing.ct(t_hat, prob)' * l)

    set_silent(model_sep)
    optimize!(model_sep)

    if -objective_value(model_sep) < 1e-5
        return DPArc[]  # no new arc found
    else
        new_arc = DPArc(
            first(collect(vertices(dd))[value.(u) .≈ 1]),
            first(collect(vertices(dd))[value.(v) .≈ 1]),
            CombinatorialPricing.convert_x_to_set(value.(l)),
        )

        # println("New u-v arc found: $new_arc with value $(objective_value(model_sep))")
        # new_arc_slack = θ_hat[new_arc.src] - θ_hat[new_arc.dst] + CombinatorialPricing.ct(t_hat, prob)' * value.(l)
        # println("New arc slack: $new_arc_slack")
        # u_set = CombinatorialPricing.convert_x_to_set(value.(u))
        # println("u = $(u_set)")
        # println("theta_hat u = $(θ_hat[new_arc.src])")
        # println("theta u = $(θ' * value.(u))")
        # println("theta u = $(θ[first(u_set)])")
        # v_set = CombinatorialPricing.convert_x_to_set(value.(v))
        # println("v = $(v_set)")
        # println("theta_hat v = $(θ_hat[new_arc.dst])")
        # println("theta v = $(θ' * value.(v))")
        # println("theta v = $(θ[first(v_set)])")
        
        return [new_arc]
    end
end

function solve_separation_ut(prob, t_hat, θ_hat, dd)
    # problem parameters
    I = zeros(Int, num_items(prob), nv(dd))  # indicator of unassigned items
    θ = zeros(nv(dd))  # dual variables vector
    s = zeros(nv(dd), length(CombinatorialPricing.follower_b(prob)))  # state of the vertices
    for (j, (i, s_v)) in enumerate(collect(vertices(dd)))
        I[(i*grouping+1):end, j] .= 1
        θ[j] = θ_hat[(i, s_v)]
        # TODO: this only works for knapsack! How to handle the state of the other problems?
        s[j, :] .= s_v.remaining
    end

    # build separation problem
    # model_sep = CombinatorialPricing.blank_model(silent=true)
    model_sep = Model(SCIP.Optimizer)
    set_optimizer_attribute(model_sep, MOI.Silent(), true)

    @variable(model_sep, u[1:nv(dd)], Bin)
    @variable(model_sep, l[1:num_items(prob)], Bin)

    @constraint(model_sep, l .<= I * u)
    @constraint(model_sep, s' * u .+ CombinatorialPricing.follower_A(prob) * l .<= CombinatorialPricing.follower_b(prob))

    @constraint(model_sep, sum(u) == 1)
    @constraint(model_sep, u[end] == 0)  # not allowed to select terminal node

    @objective(model_sep, Min, θ' * u + CombinatorialPricing.ct(t_hat, prob)' * l)

    set_silent(model_sep)
    optimize!(model_sep)

    if -objective_value(model_sep) < 1e-5
        return DPArc[]  # no new arc found
    else
        new_arc = DPArc(
            first(collect(vertices(dd))[value.(u) .≈ 1]),
            sink_node(dd),
            CombinatorialPricing.convert_x_to_set(value.(l)),
        )
        
        # println("New u-t arc found: $new_arc with value $(objective_value(model_sep))")
        # new_arc_slack = θ_hat[new_arc.src] - θ_hat[new_arc.dst] + CombinatorialPricing.ct(t_hat, prob)' * value.(l)
        # println("New arc slack: $new_arc_slack")
        # u_set = CombinatorialPricing.convert_x_to_set(value.(u))
        # println("u = $(u_set)")
        # println("theta_hat u = $(θ_hat[new_arc.src])")
        # println("theta u = $(θ' * value.(u))")
        # println("theta u = $(θ[first(u_set)])")
        
        return [new_arc]
    end
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
    samples = rand(sampler, ddwidth)

    # Create a parition of items into layers
    partition = Partition(prob, DefaultPartitioning(grouping))
    # WARNING: the default partitioning is important for the current implementation of the
    # separation problem (construction of I matrix)!!!

    # Create an empty decision diagram
    dd = DPGraph(prob, partition)

    # Populate the DD with nodes using the above samples
    populate_nodes!(dd, BitSet.(samples))

    # Prepare to collect log data
    log_data = Vector{NamedTuple{(:iter, :dual_bound, :solve_time, :violation, :cut_time), Tuple{Int, Float64, Float64, Float64, Float64}}}()

    for iter in 1:maxiter
        """ FIND CANDIDATE SOLUTION """
        # Create baseline DD model
        model_dd = base_model(dd.prob)
        CombinatorialPricing.add_dpgraph_dual!(model_dd, dd)

        # Solve the relaxation
        set_silent(model_dd)
        optimize!(model_dd)

        if ~is_solved_and_feasible(model_dd)
            println("Instance $file is not feasible. Stopping...")
            break
        end

        # candidate solution
        x_hat = value.(model_dd[:x])
        t_hat = value.(model_dd[:t])
        θ_hat = -value.(model_dd[:y])

        """ FIND NEW CUT """
        start_time = time()
        if method == :best_response
            new_arcs = solve_best_response(prob, t_hat, dd)
        elseif method == :separation
            new_uv_arcs = solve_separation_uv(prob, t_hat, θ_hat, dd)
            new_ut_arcs = solve_separation_ut(prob, t_hat, θ_hat, dd)
            new_arcs = vcat(new_uv_arcs, new_ut_arcs)
        end
        cut_time = time() - start_time

        violation = 0.0
        for a in new_arcs
            # add new arcs
            push!(dd.arcs, a)

            # maximal violation over all arcs added
            ell = CombinatorialPricing.convert_set_to_x(a.action, num_items(prob))
            new_arc_slack = θ_hat[a.src] - θ_hat[a.dst] + CombinatorialPricing.ct(t_hat, prob)' * ell
            println("New arc $a with slack $new_arc_slack")
            violation = max(violation, -new_arc_slack)
        end

        Base.GC.gc()

        # log
        if iter == 1
            @printf("| %4s | %9s | %12s |\n", "Iter", "Bound", "Violation")
            @printf("|%s|%s|%s|\n", "-"^6, "-"^11, "-"^14)
        end
        @printf("| %4d | % 9.3f | % 12.3f |\n", iter, objective_value(model_dd), violation)

        # Save log data
        push!(log_data, (
            iter=iter,
            dual_bound=objective_value(model_dd),
            solve_time=solve_time(model_dd),
            violation=violation,
            cut_time=cut_time,
        ))

        ### STOPPING CRITERION ###
        if violation < 1e-5
            println("Stopping criterion met")
            break
        end
        ### STOPPING CRITERION ###
    end

    return DataFrame(log_data)
end

results_dir = "results_DD/"
for arg in ARGS
    if endswith(arg, ".json") || endswith(arg, ".ki")
        println("Solving instance $arg")

        instance = splitext(basename(arg))[1]

        for method in [:best_response, :separation]
            df_fpath = results_dir * instance * "_" * string(method) * ".csv"
            # if isfile(df_fpath) && method != :separation
            if isfile(df_fpath)
                println("Results already exist for $method. Skipping...")
            else
                println("Running $method method")
                df = run_test(arg, method)
                CSV.write(df_fpath, df)
            end
        end
    end
end
