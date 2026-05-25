using CombinatorialPricing
using Test

@testset "pricing problem JSON round trip" begin
    mktemp() do path, _
        prob = KnapsackPricing([2, 3, 5], 8, [4.0, 6.0, 9.0], BitSet([1, 3]))
        write(path, prob)

        restored = read(path, KnapsackPricing)

        @test weights(restored) == weights(prob)
        @test capacity(restored) == capacity(prob)
        @test base_values(restored) == base_values(prob)
        @test tolled(restored) == tolled(prob)
    end

    mktemp() do path, _
        prob = MaxStableSetPricing(4, [(1, 2), (2, 3), (3, 4)], [1.0, 2.0, 3.0, 4.0], BitSet([2]))
        write(path, prob)

        restored = read(path, MaxStableSetPricing)

        @test restored.vertices == prob.vertices
        @test restored.edges == prob.edges
        @test base_values(restored) == base_values(prob)
        @test tolled(restored) == tolled(prob)
    end

    mktemp() do path, _
        prob = MinSetCoverPricing(4, [BitSet([1, 2]), BitSet([2, 4])], [1.5, 2.5], BitSet([1]))
        write(path, prob)

        restored = read(path, MinSetCoverPricing)

        @test num_elements(restored) == num_elements(prob)
        @test restored.sets == prob.sets
        @test base_costs(restored) == base_costs(prob)
        @test tolled(restored) == tolled(prob)
    end
end
