module GraphPlotExt

using Colors
using CombinatorialPricing
using GraphPlot
using Graphs

function CombinatorialPricing.plot_solution(prob::CombinatorialPricing.MaxStableSetPricing, x_set::BitSet)
    g = CombinatorialPricing.graph(prob)
    n, i1 = Graphs.nv(g), CombinatorialPricing.tolled(prob)

    nodelabel = string.(1:n)
    nodesize = 0.25 / sqrt(n)
    nodestrokec = [i in x_set ? colorant"red" : nothing for i in 1:n]
    nodefillc = [i in i1 ? colorant"orange" : colorant"gray" for i in 1:n]

    return GraphPlot.gplot(g; nodelabel, nodesize, nodestrokec, nodefillc, nodestrokelw=1.0)
end

end
