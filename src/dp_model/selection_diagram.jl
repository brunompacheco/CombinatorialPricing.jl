# State
struct SDState <: DPState{AbstractProblem}
    selected::BitSet
end

Base.show(io::IO, s::SDState) = print(io, Tuple(s.selected))
Base.hash(s::SDState, h::UInt) = hash(s.selected, h)
Base.:(==)(s::SDState, t::SDState) = s.selected == t.selected

const SDArc = DPArc{DPState}

# Graph
struct SDGraph{P<:AbstractProblem} <: AbstractGraph{SDArc}
    prob::P
    layers::Vector{Vector{SDState}}  # Not include the source layer (layer 0)
    arcs::Vector{SDArc}
end

source_state(::SDGraph) = SDState(BitSet())
sink_state(::SDGraph) = SDState(BitSet())

source_node(g::SDGraph) = (0, source_state(g))
sink_node(g::SDGraph) = (nl(g), sink_state(g))

Graphs.nv(g::SDGraph) = 1 + sum(length.(g.layers))      # 1 extra node for the source
Graphs.vertices(g::SDGraph) = Iterators.flatten(vcat((source_node(g),), [((l, st) for st in layer) for (l, layer) in enumerate(g.layers)]))
Graphs.has_vertex(g::SDGraph, u) = u == source_node(g) || any(layer -> u in layer, g.layers)

Graphs.ne(g::SDGraph) = length(g.arcs)
Graphs.edges(g::SDGraph) = g.arcs
Graphs.edgetype(g::SDGraph) = SDArc
Graphs.has_edge(g::SDGraph, u, v) = any(e -> e == SDArc(u, v), g.arcs)

Graphs.inneighbors(g::SDGraph, u) = [e.src for e in g.arcs if e.dst == u]
Graphs.outneighbors(g::SDGraph, u) = [e.dst for e in g.arcs if e.src == u]
Graphs.is_directed(::Type{SDGraph})::Bool = true
Graphs.is_directed(g::SDGraph)::Bool = true

nl(g::SDGraph) = length(g.layers)                    # Not include the source layer (layer 0)
layer(g::SDGraph, l) = l == 0 ? [source_state(g)] : g.layers[l]
Base.getindex(g::SDGraph, l) = layer(g, l)

num_items(g::SDGraph) = num_items(g.prob)

Base.show(io::IO, g::SDGraph{P}) where {P} = print(io, "SDGraph{$P} with $(nl(g)) layers, $(nv(g)) nodes, and $(ne(g)) arcs")

Base.widen(::Type{SDArc}) = SDArc

# Construction
function sdgraph_from_pairs(prob::AbstractProblem, pairs)
    layer1 = SDState.(BitSet.(unique!(collect(Iterators.flatten(pairs)))))
    layer2 = SDState.(BitSet.(unique!(pairs)))
    layer3 = [SDState(BitSet())]
    layers = [layer1, layer2, layer3]

    arcs = SDArc[]
    g = SDGraph(prob, layers, arcs)
    source = source_node(g)
    for s in layer1
        push!(arcs, SDArc(source, (1, s), s.selected))
    end
    for s in layer1, t in layer2
        issubset(s.selected, t.selected) || continue
        push!(arcs, SDArc((1, s), (2, t), setdiff(t.selected, s.selected)))
    end

    return g
end

