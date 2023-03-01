using PeriodicGraphs, LinearAlgebra, Graphs, StaticArrays
using SnoopPrecompile

@static if VERSION < v"1.7-"
    Returns(a) = (_ -> a)
end

macro _precompile_g(N, str)
    g = gensym("graph")
    ret = quote
        $g = PeriodicGraph($str)
        $g == parse(PeriodicGraph, string($g))
        PeriodicGraph{$N}($str) == parse(PeriodicGraph{$N}, string($g))
        $g != PeriodicGraph{$N}()
        vert = PeriodicVertex(3, $(ntuple(i -> ifelse(i == 1, 1, 0), N)))
        edg = PeriodicEdge(1, vert)
        PeriodicEdge{$N}(1, vert) == PeriodicEdge(1, 3, $(ntuple(Returns(0), N)))
        PeriodicEdge(1, 3, zero(SVector{$N,Int})) == PeriodicEdge{$N}(1, 3, $(ntuple(Returns(0), N)))
        PeriodicEdge{$N}(1, 3, zero(SVector{$N,Int}))
        add_edge!($g, 2, vert)
        rem_edge!($g, 2, vert)
        truncated_graph($g)
        quotient_graph($g)
        coordination_sequence($g, 2, 5)
        first(Iterators.drop(edges($g), 2)) < edg
        edg in edges($g)
        edges($g) == edges($g)
        edges($g) < edges($g)
        find_edges($g, 2, 4)
        has_edge($g, 1, 3)
        $g[[1,4]]
        $g[[3,1,4,2]]
        $g[[true, true, false, true]]
        connected_components($g)
        add_vertex!($g)
        add_vertices!($g, 5)
        rem_vertex!($g, 10)
        rem_vertices!($g, [8,9])
        rem_vertices!($g, [5,6,7], true)
        rings($g, 3)
        strong_rings($g, 3)
        strong_erings($g, 3)
        RingAttributions($g, true, 3)[1][1]
        RingAttributions($g, 3)[1][1]
        RingAttributions($g, true)[1][1]
        RingAttributions($g)[1][1]
    end
    if N != 0
        append!(ret.args, (quote
            slice_graph($g, [$N])
            slice_graph($g, ($(fld1(N,2)),))
            make_supercell($g, $(ntuple(Returns(2), N)))
            swap_axes!($g, $([i for i in 1:N]))
            swap_axes!($g, $(SVector{N,Int}([i for i in 1:N])))
            offset_representatives!($g, [SVector{$N,Int}(i+j*(-1)^i for j in 1:$N) for i in 1:nv($g)])
        end).args)
    end
    # dimensionality detection and reduction
    for i in 0:(N-1)
        push!(ret.args, :(PeriodicGraph{$i}(slice_graph($g, $(collect(1:(N-i)))))))
    end
    # dimensionality increase
    for j in N:3
        push!(ret.args, :(PeriodicGraph{$j}($g)))
    end

    return ret
end

@precompile_setup begin
    @precompile_all_calls begin
        @_precompile_g 0 "0   1 2   1 3   1 4  2 3   3 4"
        @_precompile_g 1 "1   1 2  1   1 3  0   1 4  0   2 3  0   3 4  -1"
        @_precompile_g 2 "2   1 2  1 0   1 3  0 0   1 4  0 0   2 3  0 0   3 4  -1 0"
        @_precompile_g 3 "3   1 2  1 0 0   1 3  0 0 0   1 4  0 0 0   2 3  0 0 -1   3 4  -1 0 1"

        wno = PeriodicGraph("3 1 2 0 0 0 1 2 0 1 0 1 3 -1 0 0 1 3 0 0 0 1 4 0 0 0 1 5 0 0 0 1 6 0 0 0 2 3 0 -1 0 2 4 0 0 0 2 5 0 -1 -1 2 6 0 0 0 2 6 1 -1 -1 3 4 0 1 1 3 4 1 0 0 3 5 0 0 0 3 6 1 0 -1 4 5 0 -1 -1 4 5 0 0 -1 4 6 0 0 -1 5 6 0 0 0 5 6 1 0 0")
        RingAttributions(wno)
        RingAttributions(wno, true)

        PeriodicGraphs.IterativeGaussianEliminationNone([3,4])
        PeriodicGraphs.IterativeGaussianEliminationLength([2,7])
        PeriodicGraphs.IterativeGaussianEliminationDecomposition([1,6])
        PeriodicGraphs.gaussian_elimination!(PeriodicGraphs.IterativeGaussianEliminationNone(), Int[2])
        PeriodicGraphs.gaussian_elimination!(PeriodicGraphs.IterativeGaussianEliminationLength(), Int[2])
        PeriodicGraphs.gaussian_elimination!(PeriodicGraphs.IterativeGaussianEliminationDecomposition(), Int[2])
    end
end
