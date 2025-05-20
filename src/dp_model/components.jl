const DPNode{S} = Tuple{Int,S} where S<:DPState{<:AbstractProblem}
const DPAction = BitSet

struct DPArc{S<:DPState} <: AbstractEdge{S}
    src::DPNode{S}
    dst::DPNode{S}
    action::DPAction
end

Graphs.src(a::DPArc) = a.src
Graphs.dst(a::DPArc) = a.dst
action(a::DPArc) = a.action

Base.reverse(a::DPArc) = typeof(a)(a.dst, a.src, a.action)

Base.show(io::IO, a::DPArc) = print(io, "DPArc $(src(a)) => $(dst(a)) $(tuple(action(a)...))")
Base.hash(a::DPArc, h::UInt) = hash(src(a), hash(dst(a), hash(action(a), h)))
Base.:(==)(a::DPArc, b::DPArc) = (a.src == b.src) && (a.dst == b.dst) && (a.action == b.action)

Base.merge(a::DPArc, b::DPArc) = DPArc(src(a), dst(b), action(a) ∪ action(b))
Base.length(a::DPArc) = dst(a)[1] - src(a)[1]
