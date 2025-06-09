module CombinatorialPricing

using Printf, UnicodeGraphics
using DataStructures, SparseArrays
using JSON, Unmarshal, UnPack, Memoization
using JuMP
using Graphs
using Random, Distributions, StatsBase
using Parameters: @with_kw
using GraphPlot, Colors
using DataFrames

using JuMP: JuMP.Containers.DenseAxisArray

include("interface/problem.jl")
include("interface/model.jl")
include("interface/dp_model.jl")
include("interface/sampler.jl")

include("misc/lazyenv.jl")
include("misc/blank_model.jl")
include("misc/print_utils.jl")
include("misc/matrices.jl")
include("misc/utils.jl")
include("misc/random_pair.jl")

include("dp_model/partition.jl")
include("dp_model/components.jl")
include("dp_model/graph.jl")
include("dp_model/type_utils.jl")
include("dp_model/unstructured_path.jl")
include("dp_model/structured_path.jl")
include("dp_model/populate.jl")
include("dp_model/stats.jl")
include("dp_model/unique.jl")
include("dp_model/selection_diagram.jl")

include("models/base_model.jl")
include("models/base_follower.jl")
include("models/follower_model.jl")
include("models/trivial_solutions.jl")
include("models/heuristic.jl")
include("models/cutting_plane.jl")
include("models/value_function.jl")
include("models/dpgraph_model.jl")
include("models/colgen_state.jl")
include("models/colgen_model.jl")
include("models/sdgraph_model.jl")

include("samplers/bilevel_feasible.jl")

include("knapsack/problem.jl")
include("knapsack/probgen.jl")
include("knapsack/dp_model.jl")
include("knapsack/maximal_sampler.jl")

include("maxstab/problem.jl")
include("maxstab/probgen.jl")
include("maxstab/dp_model.jl")
include("maxstab/maximal_sampler.jl")
include("maxstab/analysis.jl")

include("mincover/problem.jl")
include("mincover/probgen.jl")
include("mincover/dp_model.jl")
include("mincover/minimal_sampler.jl")

include("interdiction/problem.jl")
include("interdiction/dp_model.jl")
include("interdiction/maximal_sampler.jl")
include("interdiction/base_model.jl")

export AbstractProblem, PricingProblem, num_items, tolled, toll_free, base_costs, generate
export toll_bounds, add_primal!
export DPState, source_state, sink_state, transition, is_valid_transition
export SolutionSampler, problem

export Partition
export AbstractPartitioning, DefaultPartitioning, RandomPartitioning, TolledFirstPartitioning
export DPNode, DPAction, DPArc, action
export DPGraph, source_node, sink_node, nl, layer, partition
export node_stats, count_paths, arc_distribution
export unique_paths, unique_arcs, count_unique_paths, unique_graph

export unstructured_path
export structured_path
export populate_nodes!, populate_arcs!

export SDState, SDArc, SDGraph
export sdgraph_from_pairs, random_pair

export base_model
export follower_model, set_toll!
export toll_free_solution, null_toll_solution, toll_free_cost, null_toll_cost, difference_estimate
export value_function_model, add_value_function_constraint!
export dpgraph_model
export colgen_model
export sdgraph_model

export BilevelFeasibleSampler

export KnapsackPricing
export base_values, weights, capacity, density
export MaximalKnapsackSampler

export MaxStableSetPricing
export graph, base_values
export MaximalStableSetSampler
export plot_solution

export MinSetCoverPricing
export num_elements, restrict
export MinimalSetCoverSampler

export KnapsackInterdiction
export lower_cap, upper_cap, lower_weights, upper_weights, profits, upper_density, lower_density
export MaximalKnapsackInterdictionSampler

# Exports from JuMP
export optimize!, value, objective_value

# Exports from Graphs
export nv, ne, vertices, edges, src, dst

end # module CombinatorialPricing
