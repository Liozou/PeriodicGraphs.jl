using Test
using PeriodicGraphs

using Aqua
Aqua.test_all(PeriodicGraphs; ambiguities=(; exclude=[sort!])) # TODO: fix upsrteam at Graphs.jl

using Graphs
using PeriodicGraphs: vertex_permutation, LoopException, extended_gcd
using StaticArrays: SVector, SMatrix
using Combinatorics

function naive_periodicgraph(nv::Integer, t::AbstractVector{PeriodicEdge{N}}) where N
    sort!(t); unique!(t)
    g = PeriodicGraph{N}(nv)
    for e in t
        add_edge!(g, e)
    end
    return g
end

function is_well_formed(g::PeriodicGraph{N}) where N
    isempty(g.nlist) || sum(length, g.nlist) == 2*g.ne[] || return false
    for (i, l) in enumerate(g.nlist)
        if isempty(l)
            g.directedgestart[i] == 1 || return false
            continue
        end
        issorted(l) || return false
        allunique(l) || return false
        any(==(PeriodicVertex{N}(i)), l) && return false
        dedgestart = g.directedgestart[i]
        if dedgestart > length(l)
            dedgestart == length(l) + 1 || return false
            isdirectedge(PeriodicEdge{N}(i, l[end])) && return false
            continue
        end
        isdirectedge(PeriodicEdge{N}(i, l[dedgestart])) || return false
        if dedgestart != 1 && isdirectedge(PeriodicEdge{N}(i, l[dedgestart-1]))
            return false
        end
        for j in dedgestart:length(l)
            rev = PeriodicVertex{N}(i, .- l[j].ofs)
            any(==(rev), g.nlist[l[j].v]) || return false
        end
    end
    g.width[] == -1 || g.width[] > 0 || return false
    return true
end

@testset "PeriodicVertex" begin
    @test_throws MethodError PeriodicVertex2D() # Always require a vertex identifier
    v, ofs = PeriodicVertex2D(1)
    @test v == 1 && iszero(ofs)
    @test v == first(PeriodicVertex3D(1, (3,4,5)))
    @test ndims(PeriodicVertex2D(1)) == 2
    @test PeriodicVertex2D(1) != PeriodicVertex3D(1)
    @test PeriodicVertex{4}(5) == PeriodicVertex{4}(5)
    @test PeriodicVertex1D(1, SVector{1,Int}(3)) == PeriodicVertex1D(1, [3])
    @test PeriodicVertex(3, (1,2)) == PeriodicVertex2D(3, (1,2)) == PeriodicVertex(3, SVector{2,Int}(1,2))
    @test PeriodicVertex{0}(7) == PeriodicVertex{0}(7, ()) == PeriodicVertex(7, ())
    h = hash(PeriodicVertex(3, (1,2)))
    @test h == hash(PeriodicVertex(3, SVector{2,Int}(1,2)))
    @test hash(PeriodicVertex2D(3)) != h != hash(PeriodicVertex3D(3, (1,2,0))) != hash(PeriodicVertex(2, (1,2)))
    @test PeriodicVertex[(1, (-1,0))] == [PeriodicVertex2D(1, (-1,0))] == PeriodicVertex2D[(1, (-1,0))]
    @test_throws MethodError PeriodicVertex(2, [4]) # Performance trap if the dimension is unknown
    @test_throws DimensionMismatch PeriodicVertex2D(1, (2,))
    @test_throws DimensionMismatch PeriodicVertex3D(3, (3, 1, 4, 0))
    @test ndims(PeriodicVertex(0, SVector{2,Int}(1,4))) == 2 # Inference from SVector
    @test ndims(PeriodicVertex{0}(1)) == 0 # No need of a second argument here
    @test PeriodicVertex{0}(1) < PeriodicVertex{0}(2)
    @test PeriodicVertex(2, (3,)) >  PeriodicVertex(2, (1,))
    @test string(PeriodicVertex3D(12, (1,0,-1))) == "PeriodicVertex3D(12, (1,0,-1))"
    @test string(PeriodicVertex2D[(1, (0,0)), (1, (-1,0))]) == "PeriodicVertex2D[(1, (0,0)), (1, (-1,0))]"
    @test string(PeriodicVertex1D(2, (1,))) == "PeriodicVertex1D(2, (1,))"
    @test string(PeriodicVertex{0}(12)) == "PeriodicVertex{0}(12)"
    @test string(PeriodicVertex{0}[1, 4]) == "PeriodicVertex{0}[1, 4]"
    @test last(PeriodicVertex(1, (2,1))) == SVector{2,Int}(2,1)
end

@testset "PeriodicEdge" begin
    @test_throws LoopException PeriodicEdge2D(1, 1, (0, 0))
    @test_throws LoopException PeriodicEdge(2, PeriodicVertex3D(2))
    @test_throws LoopException PeriodicEdge{0}(3, 3, Int[])
    @test_throws LoopException PeriodicEdge(1, 1, ())
    try
        PeriodicEdge(1, PeriodicVertex2D(1))
        @test false
    catch e
        @test e isa LoopException
        io = IOBuffer()
        showerror(io, e)
        @test occursin("LoopException", String(take!(io)))
    end
    # There should be no default zero offset to prevent programmers from forgetting the offset
    # The only exception is for PeriodicVertex and PeriodicEdge with explicit dimension 0
    @test_throws MethodError PeriodicEdge1D(1, 2)
    @test_throws MethodError PeriodicEdge(3, 6)
    @test PeriodicEdge2D(1, 2, (0,0)) != PeriodicEdge3D(1, 2, (0,0,0))
    @test PeriodicEdge(2, PeriodicVertex2D(1, (1,3))) == PeriodicEdge2D(2, 1, (1,3))
    @test PeriodicEdge(1, 1, SVector{1,Int}(2)) == PeriodicEdge1D(1, 1, SVector{1,Int}(2))
    @test PeriodicEdge(2, 1, ()) == PeriodicEdge{0}((2, 1, SVector{0,Int}())) == PeriodicEdge{0}(2, 1, ())
    @test PeriodicEdge((1, 2, (1,0,0))) == PeriodicEdge3D((1, 2, (1,0,0)))
    @test hash(PeriodicEdge((1, 2, (1,0,0)))) == hash(PeriodicEdge3D((1, 2, (1,0,0))))
    @test hash(PeriodicEdge((1, 2, (-1,)))) != hash(PeriodicEdge((2, 1, (1,))))
    @test directedge(PeriodicEdge((2, 1, (1,)))) == directedge(PeriodicEdge((1, 2, (-1,)))) == PeriodicEdge((1, 2, (-1,)))
    @test directedge(PeriodicVertex1D(3, (1,)), PeriodicVertex1D(2, (0,))) == PeriodicEdge(2, 3, (1,))
    @test directedge(1, 1, (1, 0)) == directedge(1, 1, (-1, 0)) == PeriodicEdge(1, 1, (1, 0))
    @test directedge(3, 1, ()) == PeriodicEdge{0}(1, 3, ()) == PeriodicEdge{0}(1, 3)
    @test convert(PeriodicEdge, (3, PeriodicVertex{1}(2, (1,)))) == convert(PeriodicEdge1D, (3, PeriodicVertex{1}(2, (1,))))
    @test PeriodicEdge1D[(1,2,SVector{1,Int}(3))] == PeriodicEdge[(1,2,(3,))] == [PeriodicEdge(1, 2, (3,))]
    @test PeriodicEdge{0}[(2,3)] == PeriodicEdge{0}[(2, PeriodicVertex{0}(3))] == [PeriodicEdge{0}((2, 3))]
    @test string(PeriodicEdge3D[(1, 2, (1,0,0))]) == "PeriodicEdge3D[(1, 2, (1,0,0))]"
    @test string(PeriodicEdge(3,4, (0,0))) == "PeriodicEdge2D(3, 4, (0,0))"
    @test string(PeriodicEdge{1}(1, 1, (2,))) == "PeriodicEdge1D(1, 1, (2,))"
    @test string(PeriodicEdge1D[(3, PeriodicVertex1D(4)), (2, 2, SVector{1,Int}(-1))]) == "PeriodicEdge1D[(3, 4, (0,)), (2, 2, (-1,))]"
    @test string(PeriodicEdge{0}[(1, 2), (2, 3)]) == "PeriodicEdge{0}[(1, 2), (2, 3)]"
    @test reverse(reverse(PeriodicEdge(1, 1, (2,3,1)))) == PeriodicEdge3D(1, 1, (2,3,1))
    @test reverse(PeriodicEdge(1, 2, (1,-2))) == PeriodicEdge(2, 1, (-1,2))
    @test reverse(PeriodicEdge(3, 3, (1,0))) == PeriodicEdge(3, 3, (-1,0))
    @test src(PeriodicEdge(1, 2, (-1, 0))) == PeriodicVertex2D(1)
    @test first(PeriodicEdge(4, 3, ())) == 4
    @test dst(PeriodicEdge(2, 2, (1,))) == last(PeriodicEdge(1, 2, (1,))) == PeriodicVertex1D(2, (1,))
    _s, (_d, _ofs) = PeriodicEdge(12, 4, (-1,-2))
    @test _s == 12 && _d == 4 && _ofs == [-1,-2]
    @test PeriodicEdge2D(1, 1, (2,3)) < PeriodicEdge2D(1, 1, (3,0)) < PeriodicEdge2D(1, 2, (3,0)) < PeriodicEdge2D(2, 2, (3,0))
    @test cmp(PeriodicEdge{0}(2, 3, ()), PeriodicEdge(1, 3, ())) == 1
    @test cmp(PeriodicEdge2D(1, 1, (0,1)), PeriodicEdge(1, 3, (0,0))) == -1
    @test cmp(PeriodicEdge2D(7, 7, (1,2)), PeriodicEdge(7, 7, (2,0))) == -1
    @test cmp(PeriodicEdge1D(1, 2, (0,)), PeriodicEdge(1, 2, (0,))) == 0
    @test ndims(PeriodicEdge(1, 2, ())) == 0
    @test ndims(PeriodicEdge(1, 1, (2,0,0,0))) == 4
end

@testset "Empty graphs" begin
    @test_throws MethodError PeriodicGraph() # Missing dimension
    g = PeriodicGraph3D()
    @test ne(g) == nv(g) == 0
    @test ndims(g) == 3
    @test !rem_edge!(g, PeriodicEdge3D(1, 2, (0,0,0)))
    @test_throws ArgumentError rem_vertex!(g, 0)
    @test_throws ArgumentError rem_vertex!(g, 1)
    @test isempty(rem_vertices!(g, Int[]))
    @test isempty(edges(g))
    @test isempty(vertices(g))
    gg = vertex_permutation(g, [])
    @test ne(gg) == nv(gg) == 0 && isempty(edges(gg))
    @test g == gg
    gg, vlist = induced_subgraph(g, Int[])
    @test isempty(vlist)
    @test ne(gg) == nv(gg) == 0 && isempty(edges(gg))
    @test g == gg
    gg = offset_representatives!(g, Int[])
    @test ne(gg) == nv(gg) == 0 && isempty(edges(gg))
    @test ne(g) == nv(g) == 0 && isempty(edges(g))
    @test g == gg
    @test_throws DimensionMismatch swap_axes!(g, [2,1])
    gg = swap_axes!(g, Int[1,2,3])
    @test ne(gg) == nv(gg) == 0 && isempty(edges(gg))
    @test ne(g) == nv(g) == 0 && isempty(edges(g))
    @test g == gg
    @test isempty(connected_components(g))
    @test edgetype(g) == PeriodicEdge3D
    @test eltype(g) == PeriodicVertex3D
    @test edges(g) isa PeriodicGraphs.PeriodicEdgeIter{3}
    @test_throws BoundsError coordination_sequence(g, 1, 0)
    @test_throws BoundsError coordination_sequence(g, 1, 0) # once width is set
    @test !has_vertex(g, 0)
    @test !has_vertex(g, 1)
end

@testset "Basic graph construction and modification" begin
    @test_throws MethodError PeriodicGraph(2) # missing dimension
    g = PeriodicGraph3D()
    @test is_well_formed(g)

    @testset "add_vertex!, rem_vertex!" begin
        @test add_vertex!(g)
        @test is_well_formed(g)
        @test ne(g) == 0 && nv(g) == 1
        @test g == PeriodicGraph3D(1)
        @test add_vertex!(g)
        @test ne(g) == 0 && nv(g) == 2
        @test g == PeriodicGraph3D(2)
        @test_throws ArgumentError rem_vertex!(g, 0)
        @test_throws ArgumentError rem_vertex!(g, 3)
        @test rem_vertex!(g, 1)
        @test is_well_formed(g)
        @test ne(g) == 0 && nv(g) == 1
        @test g == PeriodicGraph3D(1)
        @test rem_vertex!(g, 1)
        @test ne(g) == 0 && nv(g) == 0
        @test g == PeriodicGraph3D(0)
    end

    @testset "add_vertices!, rem_vertices!" begin
        @test isempty(edges(g))
        @test isempty(vertices(g))
        @test_throws ArgumentError rem_vertices!(g, Int[1])
        @test isempty(rem_vertices!(g, Int[]))
        @test add_vertices!(g, 0) == 0
        @test nv(g) == ne(g) == 0 && g == PeriodicGraph3D()
        @test add_vertices!(g, 1) == 1
        @test ne(g) == 0 && nv(g) == 1 && g == PeriodicGraph3D(1)
        @test only(vertices(g)) == 1
        @test rem_vertices!(g, [1]) == Int[]
        @test add_vertices!(g, 5) == 5
        @test is_well_formed(g)
        @test ne(g) == 0 && nv(g) == 5 && g == PeriodicGraph3D(5)
        @test vertices(g) == 1:5
        @test issetequal(rem_vertices!(g, [2,3]), [1,4,5])
        @test ne(g) == 0 && nv(g) == 3 && g == PeriodicGraph3D(3)
        @test vertices(g) == 1:3
        @test rem_vertex!(g, 2)
        @test is_well_formed(g)
        @test ne(g) == 0 && nv(g) == 2 && g == PeriodicGraph3D(2)
        @test vertices(g) == [1,2]
        @test rem_vertices!(g, [1,2]) == Int[]
        @test add_vertices!(g, 5) == 5
        @test rem_vertices!(g, [2,4], true) == [1,3,5]
        @test is_well_formed(g)
        @test g == PeriodicGraph3D(3)
        @test rem_vertices!(g, [1,2,3], true) == Int[]
        @test g == PeriodicGraph3D()
    end

    @testset "Single edge" begin
        @test isempty(edges(g))
        @test isempty(vertices(g))
        @test !rem_edge!(g, 1, PeriodicVertex(2, (0,0,0)))
        @test g == PeriodicGraph3D()
        @test !add_edge!(g, PeriodicEdge(1, 1, (-1,0,0)))
        @test g == PeriodicGraph3D()
        @test !has_edge(g, 1, 1)
        @test !has_edge(g, PeriodicEdge3D(1, 2, (0,0,0)))
        @test add_vertex!(g)
        @test !has_edge(g, 1, 1)
        @test !has_edge(g, PeriodicEdge3D(1, 1, (0,1,1)))
        @test add_edge!(g, 1, PeriodicVertex(1, (-1,0,0)))
        @test is_well_formed(g)
        @test has_edge(g, 1, 1)
        @test !add_edge!(g, PeriodicEdge(1, 1, (-1,0,0)))
        @test has_edge(g, 1, 1)
        @test has_edge(g, PeriodicEdge3D(1, 1, (1,0,0)))
        @test has_edge(g, 1, PeriodicVertex(1, (-1,0,0)))
        @test g == PeriodicGraph3D(PeriodicEdge3D[(1, 1, (1,0,0))])
        @test g != PeriodicGraph3D(1)
        @test rem_edge!(g, PeriodicEdge(1, 1, (-1,0,0)))
        @test is_well_formed(g)
        @test g == PeriodicGraph3D(1)
        @test !has_edge(g, 1, 1)
        @test !has_edge(g, 1, PeriodicVertex(1, (1,0,0)))
        @test !has_edge(g, PeriodicEdge3D(1, 1, (-1,0,0)))
        @test add_edge!(g, 1, PeriodicVertex(1, (2,-1,0)))
        @test !add_edge!(g, PeriodicEdge(1, 1, (-2,1,0)))
        @test g == PeriodicGraph3D(PeriodicEdge3D[(1, 1, (-2,1,0))])
        @test has_edge(g, 1, 1)
        @test has_edge(g, PeriodicEdge3D(1, 1, (2,-1,0)))
        @test !has_edge(g, PeriodicEdge3D(1, 1, (1,0,0)))
        @test rem_edge!(g, 1, PeriodicVertex(1, (2,-1,0)))
        @test g == PeriodicGraph3D(1)
        @test add_edge!(g, PeriodicEdge(1, 1, (1,0,0)))
        @test rem_vertex!(g, 1)
        @test g == PeriodicGraph3D(0)
        @test add_vertices!(g, 2) == 2
        @test !has_edge(g, 1, 2)
        @test g == PeriodicGraph3D(2)
        @test add_edge!(g, PeriodicEdge(1, 2, (0,0,0)))
        @test has_edge(g, 1, 2)
        @test has_edge(g, PeriodicEdge(1, 2, (0,0,0)))
        @test !has_edge(g, PeriodicEdge(1, 2, (1,0,0)))
        @test g == PeriodicGraph3D(PeriodicEdge3D[(1, 2, (0,0,0))])
        @test g != PeriodicGraph3D(2)
        @test rem_edge!(g, PeriodicEdge(1, 2, (0,0,0)))
        @test !has_edge(g, 1, 2)
        @test g == PeriodicGraph3D(2)
        @test add_edge!(g, 1, PeriodicVertex(2, (-1,0,0)))
        @test g == PeriodicGraph3D(PeriodicEdge3D[(1, 2, (-1,0,0))])
        @test g == PeriodicGraph3D(PeriodicEdge3D[(2, 1, (1,0,0))])
        @test g != PeriodicGraph3D(PeriodicEdge3D[(1, 2, (1,0,0))])
        @test has_edge(g, 1, 2)
        @test has_edge(g, 1, PeriodicVertex(2, (-1,0,0)))
        @test has_edge(g, PeriodicEdge(2, 1, (1,0,0)))
        @test rem_vertex!(g, 1)
        @test g == PeriodicGraph3D(1)
        @test !has_edge(g, 1, 2)
        @test !has_edge(g, 1, 1)
        @test add_vertex!(g)
        @test nv(g) == 2
        @test add_edge!(g, PeriodicEdge(1, 2, (0,1,1)))
        @test ne(g) == 1
        @test rem_vertices!(g, [2]) == [1]
        @test ne(g) == 0
        @test nv(g) == 1
    end

    @testset "Multiple edges" begin
        edgs1 = PeriodicEdge2D[(1, 1, (0,1)), (1, 1, (-1,0)), (1, 1, (1,0))]
        g1 = PeriodicGraph(edgs1)
        @test nv(g1) == 1
        @test ne(g1) == 2
        @test is_well_formed(g1)
        @test g1 == naive_periodicgraph(1, edgs1) == PeriodicGraph(1, edgs1)
        _edgs1 = collect(edges(g1))
        @test g1 == PeriodicGraph(_edgs1)

        edgs2 = PeriodicEdge2D[(3, 1, (0,0)), (2, 3, (0,0)), (1, 1, (2,0)), (1, 3, (0,1))]
        g2 = PeriodicGraph(edgs2)
        @test nv(g2) == 3
        @test ne(g2) == 4
        @test is_well_formed(g2)
        @test g2 == naive_periodicgraph(3, edgs2) == PeriodicGraph(3, edgs2)
        _edgs2 = collect(edges(g2))
        @test g2 == PeriodicGraph(_edgs2)
        sort!(_edgs2)
        @test g2 == PeriodicGraphs.from_edges(_edgs2) == PeriodicGraphs.from_edges(3, _edgs2)

        edgs3 = PeriodicEdge2D[(1, 1, (0,1)), (1, 1, (1,0)), (1, 1, (2,0)), (1, 2, (0,1)),
                               (1, 3, (-1,0)), (1, 3, (0,1)), (2, 3, (0,0)), (3, 4, (0,0))]
        g3 = PeriodicGraph(edgs3)
        @test nv(g3) == 4
        @test ne(g3) == 8
        @test is_well_formed(g3)
        @test g3 == naive_periodicgraph(4, edgs3) == PeriodicGraph(4, edgs3)
        _edgs3 = collect(edges(g3))
        @test g3 == PeriodicGraph(_edgs3)
        sort!(_edgs3)
        @test g3 == PeriodicGraphs.from_edges(_edgs3) == PeriodicGraphs.from_edges(4, _edgs3)

        edgs4 = PeriodicEdge3D[(2, 3, (-1,0,0)), (2, 3, (0,1,0)), (4, 4, (1,-1,0)),
                               (1, 3, (0,0,1)), (3, 2, (0,0,0))]
        g4 = PeriodicGraph(edgs4)
        @test nv(g4) == 4
        @test ne(g4) == 5
        @test is_well_formed(g4)
        @test g4 == naive_periodicgraph(4, edgs4) == PeriodicGraph(4, edgs4)
        _edgs4 = collect(edges(g4))
        @test g4 == PeriodicGraph(_edgs4)
        sort!(_edgs4)
        @test g4 == PeriodicGraphs.from_edges(_edgs4) == PeriodicGraphs.from_edges(4, _edgs4)
    end
end

@testset "String graph construction" begin
    @test_throws MethodError length(PeriodicGraphs.KeyString{Int}("3 1 2 0 0 0"))
    @test isempty(PeriodicGraphs.KeyString{Int}("  "))
    @test_throws ArgumentError iterate(PeriodicGraphs.KeyString{Int}("  "), 1)
    @test isempty(PeriodicGraphs.KeyString{Int}(""))
    @test PeriodicGraph("0 1 2 1 3
    ") == PeriodicGraph(PeriodicEdge{0}[(1, 2, ()), (1, 3, ())])
    @test PeriodicGraph("3") == PeriodicGraph3D("3") == PeriodicGraph3D()
    @test ne(PeriodicGraph("3  1 1  1 0 0  1 1 -1 0 0  1 1 3 0 0")) == 2
    @test nv(PeriodicGraph("2  11 12 0 0")) == 12
    @test collect(PeriodicGraphs.KeyString{Int}("3 2 1 0 1 -1")) == [3, 2, 1, 0, 1, -1]
    @test_throws ArgumentError PeriodicGraph2D("")
    @test_throws ArgumentError PeriodicGraph("")
    @test_throws ArgumentError PeriodicGraph3D("a")
    @test_throws ArgumentError PeriodicGraph2D("2 a")
    @test_throws DimensionMismatch PeriodicGraph2D("4 1 2 0 0 0 1")
    @test_throws DimensionMismatch PeriodicGraph3D("2  1 2  0 0 0   1 3  0 1 0")
    @test_throws ArgumentError PeriodicGraph("1 2 1")
    @test_throws ArgumentError PeriodicGraph{4}("4 1")
    @test_throws ArgumentError PeriodicGraph("2 1 1 0 1 0")
    @test_throws LoopException PeriodicGraph("3 1 1 0 0 0")
    g = PeriodicGraph(" 2 1 2 0 0 2 3 1 0  ")
    sg = string(g)
    @test sg == "2 1 2 0 0 2 3 1 0"
    @test parse(PeriodicGraph, sg) == g
    @test parse(PeriodicGraph2D, sg) == g
    @test PeriodicGraph(sg) == g
    @test string(PeriodicGraph("1    23 4  5    10  9 43")) == "1 4 23 -5 9 10 -43"
    @test string(PeriodicGraph{5}()) == "5"
end

@testset "Basic graph representation and functions" begin
    @test zero(PeriodicGraph{6}) == PeriodicGraph{6}()
    io = IOBuffer()
    show(io, PeriodicGraph("2 1 2 0 0 2 3 1 0"))
    g1 = eval(Meta.parse(String(take!(io))))
    @test g1 isa PeriodicGraph2D
    @test string(g1) == "2 1 2 0 0 2 3 1 0"
    g2 = PeriodicGraph(g1)
    g3 = deepcopy(g1)
    @test g1 == g2 == g3
    @test hash(g1) == hash(g2) == hash(g3)
    rem_vertex!(g1, 1)
    @test g1 != g2 == g3
    @test hash(g1) != hash(g2) == hash(g3)
    @test add_edge!(g2, PeriodicEdge(1, 2, (-1,0)))
    @test add_edge!(g2, 1, PeriodicVertex(1, (1,1)))
    @test add_edge!(g2, PeriodicEdge(1, 3, (0,0)))
    @test find_edges(g2, 1, 1) == PeriodicVertex2D[(1, (-1,-1)), (1, (1,1))]
    @test find_edges(g2, 1, 2) == PeriodicVertex2D[(2, (-1,0)), (2, (0,0))]
    @test find_edges(g2, 3, 2) == PeriodicVertex2D[(2, (-1,0))]
    @test inneighbors(g2, 2) == outneighbors(g2, 2)
    @test has_vertex(g1, 1)
end

@testset "Complex graphs modifications" begin
    edgs = PeriodicEdge3D[(1, 1, (1,0,0)),  (1, 2, (0,0,0)), (1, 3, (0,0,1)),
                          (1, 3, (0,1,-1)), (2, 3, (0,0,0)), (2, 4, (0,-1,0))]
    g = PeriodicGraph3D(edgs)
    @test nv(g) == 4
    @test ne(g) == length(edges(g)) == length(edgs)
    @test collect(edges(g)) == edgs
    insert!(edgs, 5, PeriodicEdge(1, 4, (-1,-1,0)))
    @test add_edge!(g, edgs[5])
    @test ne(g) == length(edges(g)) == length(edgs)
    @test collect(edges(g)) == edgs
    @test rem_edge!(g, edgs[4])
    deleteat!(edgs, 4)
    @test ne(g) == length(edges(g)) == length(edgs)
    @test collect(edges(g)) == edgs
    @test g[[1,2,3,4]] == g
    gg = g[[1,4,2,3]]
    @test gg != g
    @test gg[[1,3,4,2]] == g
    gg = g[[3,1]]
    @test nv(gg) == 2
    @test collect(edges(gg)) == PeriodicEdge3D[(1, 2, (0,0,-1)), (2, 2, (1,0,0))]
    @test gg == g[[true, false, true, false]][[2,1]]
    @test_throws ArgumentError g[[3,3]]
    @test_throws ArgumentError offset_representatives!(gg, [2])
    gcopy = PeriodicGraph(g)
    @test offset_representatives!(g, [0,0,0,0]) == gcopy == g
    @test is_well_formed(g)
    ggg = offset_representatives!(gg, [(1,-1,1),0])
    @test ggg == gg
    @test collect(edges(gg)) == PeriodicEdge3D[(1, 2, (1,-1,0)), (2, 2, (1,0,0))]
    gcopy = PeriodicGraph3D(gg)
    offset_representatives!(gg, [(2,1,-3), (2,1,-3)])
    @test gg == gcopy
    @test_throws DimensionMismatch swap_axes!(g, [1,3])
    ggg = swap_axes!(gg, (1,2,3))
    @test is_well_formed(ggg)
    @test gcopy == gg == ggg
    ggg = swap_axes!(gg, (2,1,3))
    @test gcopy != gg == ggg
    @test collect(edges(gg)) == PeriodicEdge3D[(1, 2, (-1,1,0)), (2, 2, (0,1,0))]
    @test edges(gg) < edges(gcopy) # lexicographical ordering on the list of edges
end

const sqc7399 = PeriodicGraph3D("3 1 6 0 -1 -1 1 10 -1 -1 0 1 10 -1 0 0 2 5 0 0 0 2 7 -1 0 0 2 7 0 0 0 3 4 0 1 0 3 7 -1 1 0 3 7 0 0 -1 4 10 -1 0 1 4 10 0 -1 0 5 10 -1 0 0 5 10 0 0 0 6 7 0 0 0 6 7 0 1 0 7 7 0 0 1 7 7 0 1 0 7 7 1 -1 -1 7 7 1 0 0 7 8 0 0 0 7 8 0 0 1 7 10 -1 0 1 7 10 0 -1 0 7 10 0 0 0 7 10 0 0 1 8 9 0 -1 -1 9 10 0 0 0 9 10 0 0 1 10 10 0 0 1 10 10 0 1 0 10 10 1 -1 -1 10 10 1 0 0");
const sqc10673 = PeriodicGraph3D("3 1 3 0 0 0 1 25 -1 -1 -1 1 26 -1 -1 -1 2 5 0 0 0 2 23 -1 0 0 2 24 -1 0 0 3 4 0 0 -1 3 6 0 -1 0 3 7 0 0 0 4 5 0 0 0 4 16 -1 0 0 4 19 0 0 1 5 6 0 0 0 5 8 0 0 0 6 16 0 1 0 6 19 -1 0 0 7 9 0 0 0 7 10 0 0 0 8 11 0 0 0 8 12 0 0 0 9 13 0 0 0 9 17 0 0 0 10 14 0 0 0 10 18 0 0 0 11 13 0 0 0 11 20 0 0 1 12 14 0 0 0 12 15 0 1 0 13 14 0 0 0 13 21 0 0 0 14 22 0 0 0 15 16 0 0 0 15 26 0 -1 0 16 17 0 0 0 17 23 0 0 0 18 19 0 0 0 18 24 0 0 0 19 20 0 0 0 20 25 0 0 -1 21 22 0 0 0 21 23 0 0 0 21 25 0 0 0 22 24 0 0 0 22 26 0 0 0")

@testset "Neighbors" begin
    g = PeriodicGraph3D(1)
    @test add_vertices!(g, 2) == 2
    @test add_edge!(g, PeriodicEdge3D(2, 3, (0, 0, 1)))
    @test is_well_formed(g)
    @test only(neighbors(g, 2)) == PeriodicVertex3D(3, (0, 0, 1))
    @test only(neighbors(g, 3)) == PeriodicVertex3D(2, (0, 0, -1))
    @test !has_edge(g, 2, 2)
    @test add_edge!(g, PeriodicEdge3D(2, 2, (1, 0, 0)))
    @test only(neighbors(g, 3)) == PeriodicVertex3D(2, (0, 0, -1))
    @test length(neighbors(g, 2)) == 3
    @test has_edge(g, 2, 2)
    @test neighborhood(g, 1, 1) == pushfirst!(collect(neighbors(g, 1)), PeriodicVertex3D(1, (0, 0, 0)))
    @test coordination_sequence(g, 1, 3) == zeros(Int, 3)
    @test coordination_sequence(g, 3, 6) == [1, 2, 4, 4, 4, 4]

    g = PeriodicGraph{0}(4)
    @test add_edge!(g, PeriodicEdge{0}(1, 2, ()))
    @test add_edge!(g, PeriodicEdge(1, 3, ()))
    @test add_edge!(g, 2, PeriodicVertex(4, ()))
    @test coordination_sequence(g, 1, 10) == [2, 1, 0, 0, 0, 0, 0, 0, 0, 0,]

    g = PeriodicGraph1D(3)
    @test add_edge!(g, PeriodicEdge1D(1, 1, (1,)))
    @test add_edge!(g, PeriodicEdge(1, 2, (0,)))
    @test add_edge!(g, 2, PeriodicVertex(3, (0,)))
    @test add_edge!(g, 3, PeriodicVertex1D(3, (1,)))
    @test coordination_sequence(g, 1, 5) == coordination_sequence(g, 3, 5) == [3, 5, 6, 6, 6]
    @test coordination_sequence(g, 2, 5) == [2, 4, 6, 6, 6]

    g = PeriodicGraph{4}(1)
    @test add_edge!(g, PeriodicEdge{4}(1, 1, (1,0,0,0)))
    @test add_edge!(g, PeriodicEdge(1, 1, (0,1,0,0)))
    @test add_edge!(g, 1, PeriodicVertex(1, (0,0,1,0)))
    @test add_edge!(g, 1, PeriodicVertex{4}(1, (0,0,0,1)))
    @test coordination_sequence(g, 1, 30) == [8i*(i^2+2)/3 for i in 1:30] # sequence A008412 of the OEIS

    @static if VERSION < v"1.7.0-"
        @test PeriodicGraphs.savedhashes3 isa Matrix{Int16}
    else
        @test PeriodicGraphs.savedhashes3 isa Array{Int16,3}
    end
    @test coordination_sequence(sqc7399, 1, 10) == [3, 35, 198, 501, 923, 1465, 2127, 2909, 3811, 4833]
    @test coordination_sequence(sqc10673, 3, 10) == [4, 10, 22, 40, 64, 84, 130, 162, 234, 270]
end

@testset "Graph reduction" begin
    g = PeriodicGraph3D("3  1 2  0 0 0  1 3  0 1 0   2 3  1 0 -1
                            4 4  0 0 1  4 5  0 -1 0  5 6  1 0 0")
    c = truncated_graph(g)
    p = quotient_graph(g)
    @test nv(g) == nv(p) == nv(c)
    @test issubset(edges(c), edges(p))
    for e in edges(p)
        @test has_edge(g, src(e), dst(e))
    end
    @test connected_components(p) == connected_components(g)

    @test quotient_graph(slice_graph(sqc10673, SVector{3,Int}(1,2,3))) == truncated_graph(sqc10673)
    @test quotient_graph(slice_graph(sqc7399, [1,2,3])) == truncated_graph(sqc7399)

    g1 = slice_graph(g, (2,))
    @test g1 == PeriodicGraph2D("2   1 2  0 0   2 3  1 -1
                                     4 4  0 1   5 6  1 0")
    g13 = slice_graph(g, (1,3))
    @test g13 == PeriodicGraph1D(6, PeriodicEdge1D[(1, 2, (0,)), (1, 3, (1,)), (4, 5, (-1,))])

    g23 = slice_graph(g, [2,3])
    @test g23 == PeriodicGraph3D("3 1 2  0 0 0   5 6  1 0 0")
end

@testset "Periodic vertex hash" begin
    ## dim = 0
    @test hash_position(SVector{0,Int}()) == 0
    ## dim = 1
    @test hash_position(SVector{1,Int}(0)) == 0
    @test hash_position(SVector{1,Int}(-1)) == 1
    @test hash_position(SVector{1,Int}(1)) == 2
    @test hash_position(SVector{1,Int}(-2)) == 3
    @test hash_position(SVector{1,Int}(2)) == 4
    @test hash_position(SVector{1,Int}(3)) == 6
    ## dim = 2
    n::Int = 4
    seen = falses((2n-1)^2)
    for i in -n+1:n-1, j in -n+1:n-1
        x = hash_position(SVector{2,Int}(i,j))+1
        @test !seen[x]
        seen[x] = true
        @test reverse_hash_position(x-1, Val(2)) == SVector{2,Int}(i,j)
    end
    @test all(seen)
    append!(seen, falses((2n+1)^2-(2n-1)^2))
    for i in (-n,n), j in -n:n
        x = hash_position(SVector{2,Int}(i,j))+1
        @test !seen[x]
        seen[x] = true
        @test reverse_hash_position(x-1, Val(2)) == SVector{2,Int}(i,j)
    end
    for i in -n+1:n-1, j in (-n,n)
        x = hash_position(SVector{2,Int}(i,j))+1
        @test !seen[x]
        seen[x] = true
        @test reverse_hash_position(x-1, Val(2)) == SVector{2,Int}(i,j)
    end
    @test all(seen)
    ## dim = 3
    seen = falses((2n-1)^3)
    for i in -n+1:n-1, j in -n+1:n-1, k in -n+1:n-1
        x = hash_position(SVector{3,Int}(i,j,k))+1
        @test !seen[x]
        seen[x] = true
        @test reverse_hash_position(x-1, Val(3)) == SVector{3,Int}(i,j,k)
    end
    @test all(seen)
    append!(seen, falses((2n+1)^3-(2n-1)^3))
    for i in (-n,n), j in -n:n, k in -n:n
        x = hash_position(SVector{3,Int}(i,j,k))+1
        @test !seen[x]
        seen[x] = true
        @test reverse_hash_position(x-1, Val(3)) == SVector{3,Int}(i,j,k)
    end
    for i in -n+1:n-1, j in (-n,n), k in -n:n
        x = hash_position(SVector{3,Int}(i,j,k))+1
        @test !seen[x]
        seen[x] = true
        @test reverse_hash_position(x-1, Val(3)) == SVector{3,Int}(i,j,k)
    end
    for i in -n+1:n-1, j in -n+1:n-1, k in (-n,n)
        x = hash_position(SVector{3,Int}(i,j,k))+1
        @test !seen[x]
        seen[x] = true
        @test reverse_hash_position(x-1, Val(3)) == SVector{3,Int}(i,j,k)
    end

    @test all(seen)

    for x in (SVector{0,Int}(), SVector{1,Int}(47),
              SVector{2,Int}(9,3), SVector{2,Int}(13,-34),
              SVector{3,Int}(8,7,6), SVector{3,Int}(4,-8,-1), SVector{3,Int}(13,-10,24))
        @test reverse_hash_position(hash_position(x), Val(length(x))) == x
        y = PeriodicVertex(7, x)
        @test reverse_hash_position(hash_position(y, 13), 13, Val(length(x))) == y
    end

    g = PeriodicGraph3D([PeriodicEdge3D(1, 1, (0, 0, 1)),
                         PeriodicEdge3D(1, 1, (1, 1, 0))])
    @test neighborhood(g, 1, 1) == PeriodicVertex3D[(1, (0, 0, 0)),
            (1, (-1, -1, 0)), (1, (0, 0, -1)), (1, (0, 0, 1)), (1, (1, 1, 0))]
    @test neighborhood(g, 1, 2) == PeriodicVertex3D[(1, (0, 0, 0)), (1, (-1, -1, 0)),
            (1, (0, 0, -1)), (1, (0, 0, 1)), (1, (1, 1, 0)), (1, (-2, -2, 0)),
            (1, (-1, -1, -1)), (1, (-1, -1, 1)), (1, (0, 0, -2)), (1, (1, 1, -1)),
            (1, (0, 0, 2)), (1, (1, 1, 1)), (1, (2, 2, 0))]
end

@testset "Vertex and edge iteration" begin
    g::PeriodicGraph3D = PeriodicGraph3D([PeriodicEdge3D(2, 3, (-1, 0, 0)),
                                          PeriodicEdge3D(2, 3, (0, 1, 0)),
                                          PeriodicEdge3D(4, 4, (-1, 1, 0)),
                                          PeriodicEdge3D(1, 3, (0, 0, 1)),
                                          PeriodicEdge3D(3, 2, (0, 0, 0))])
    expected_edges = PeriodicEdge3D[(1, 3, (0, 0, 1)), (2, 3, (-1, 0, 0)),
                                    (2, 3, (0, 0, 0)), (2, 3, (0, 1, 0)),
                                    (4, 4, (1, -1, 0))]
    @test collect(edges(g)) == expected_edges
    @test PeriodicEdge3D(4, 4, (1, -1, 0)) ∈ edges(g)
    add_vertex!(g)
    @test collect(edges(g)) == expected_edges
    @test length(edges(g)) == ne(g)

    expected_neighbors = PeriodicVertex3D[(1, (7, -9, 1)), (2, (7, -10, 2)),
                                          (2, (7, -9, 2)), (2, (8, -9, 2))]
    ofslist = neighbors(g, PeriodicVertex(3, (7, -9, 2)))
    @test ofslist == collect(ofslist) == ofslist[begin:end] == expected_neighbors
    @test keys(ofslist) == 1:4
    @test neighbors(g, PeriodicVertex3D(2)) == neighbors(g, 2)
    @test isempty(outneighbors(g, PeriodicVertex(5, (1, 0, 1))))
    @test only(inneighbors(g, PeriodicVertex(1, (4, 5, 6)))) == PeriodicVertex(3, (4, 5, 7))
    @test indegree(g, 2) == indegree(g, PeriodicVertex3D(2)) == outdegree(g, PeriodicVertex3D(2)) == degree(g, PeriodicVertex3D(2))
    @test string(ofslist) == string(expected_neighbors)
end

@testset "Vertex removal" begin
    g::PeriodicGraph3D = PeriodicGraph3D([PeriodicEdge3D(1, 3, (0, 0, 1)),
                                          PeriodicEdge3D(2, 3, (0, 1, 0)),
                                          PeriodicEdge3D(2, 3, (0, -1, 0)),
                                          PeriodicEdge3D(3, 2, (0, 0, 0)),
                                          PeriodicEdge3D(4, 4, (1, 1, 0))])
    gg = copy(g)
    @test rem_vertices!(g, Int[]) == [1, 2, 3, 4]
    @test is_well_formed(g)
    @test g == gg
    @test rem_vertices!(g, [2, 4], false) == [1, 3]
    @test is_well_formed(g)
    @test nv(g) == 2
    @test ne(g) == 1
    @test neighbors(g, 1) == PeriodicVertex3D[(2, (0, 0, 1))]
    @test neighbors(g, 2) == PeriodicVertex3D[(1, (0, 0, -1))]
    @test rem_vertex!(gg, 2)
    @test rem_vertex!(gg, 2)
    @test g == gg

    g1 = PeriodicGraph("2
            1 1  1 0
            1 2  0 0
            1 2  0 1
            1 3  0 0
            2 3  1 0
            2 4  0 1
            3 4  0 0
            2 5  0 0
            4 5  1 1")
    g2 = copy(g1)
    @test g2 == parse(PeriodicGraph, string(g1))
    @test rem_vertices!(g1, [2,4], false) == [1, 5, 3]
    @test g1 == PeriodicGraph("2 1 1 1 0 1 3 0 0")
    @test rem_vertices!(g2, [2,4], true) == [1, 3, 5]
    gref = PeriodicGraph("2 1 1 1 0 1 2 0 0")
    @test add_vertex!(gref)
    @test g2 == gref
end

function equivalent_dict(d1, d2)
    keys(d1) == keys(d2) || return false
    for k in keys(d1)
        s1 = Set(Set(x) for x in d1[k])
        s2 = Set(Set(x) for x in d2[k])
        s1 == s2 || return false
    end
    return true
end

@testset "Dimensionality" begin
    @test dimensionality(PeriodicGraph{0}(5)) == Dict(0 => [collect(1:5)])
    g::PeriodicGraph3D = PeriodicGraph3D([PeriodicEdge3D(1, 3, (0, 0, 1)),
                                          PeriodicEdge3D(2, 3, (0, 1, 0)),
                                          PeriodicEdge3D(2, 3, (0, -1, 0)),
                                          PeriodicEdge3D(4, 4, (1, 1, 0))])
    @test equivalent_dict(dimensionality(g), Dict(1 => connected_components(g)))
    @test add_edge!(g, PeriodicEdge(3, 2, (0,0,0)))
    @test equivalent_dict(dimensionality(g), Dict(1 => connected_components(g)))
    @test add_edge!(g, PeriodicEdge(2, 2, (0,0,1)))
    @test equivalent_dict(dimensionality(g), Dict(1 => [[4]], 2 => [[1,2,3]]))
    @test rem_edge!(g, PeriodicEdge(1, 3, (0,0,1)))
    @test equivalent_dict(dimensionality(g), Dict(0 => [[1]], 1 => [[4]], 2 => [[2,3]]))
    @test rem_edge!(g, 3, PeriodicVertex(2, (0,0,0)))
    @test equivalent_dict(dimensionality(g), Dict(0 => [[1]], 1 => [[4]], 2 => [[2,3]]))
    @test rem_edge!(g, PeriodicEdge(2, 3, (0,1,0)))
    @test equivalent_dict(dimensionality(g), Dict(0 => [[1]], 1 => [[2,3],[4]]))
    @test add_edge!(g, 2, PeriodicVertex3D(3, (1,0,-1)))
    @test equivalent_dict(dimensionality(g), Dict(0 => [[1]], 1 => [[4]], 2 => [[2,3]]))
    @test add_edge!(g, 2, PeriodicVertex(3, (2,0,-1)))
    @test equivalent_dict(dimensionality(g), Dict(0 => [[1]], 1 => [[4]], 3 => [[2,3]]))
end

@testset "Dimension change" begin
    g = PeriodicGraph("3
            1 2 0 0 0
            2 3 1 0 0
            3 4 -1 0 0
            4 1 0 0 1
            4 4 0 0 1")
    @test equivalent_dict(dimensionality(g), Dict(1 => [[1,2,3,4]]))
    gg = PeriodicGraph1D(g)
    @test string(gg) == "1 1 2 0 1 4 -1 2 3 0 3 4 0 4 4 1"
    @test PeriodicGraph2D(gg) == PeriodicGraph2D(g) ==
        PeriodicGraph("2 1 2 0 0 1 4 -1 0 2 3 0 0 3 4 0 0 4 4 1 0")
    @test_throws DimensionMismatch PeriodicGraph{0}(g)

    s = [14,28,-17,-34,12]
    d, λ = extended_gcd(s)
    @test d == 1
    @test reduce(+, s .* λ) == d
end

@testset "ConstMiniBitSet" begin
    CMBitSet = PeriodicGraphs.ConstMiniBitSet
    @test CMBitSet() == CMBitSet{UInt64}()
    @test only(CMBitSet(52)) == 52
    @test CMBitSet([3, 4]) == PeriodicGraphs.push(CMBitSet{UInt64}(4), 3)
    for (m,T) in zip((7, 15, 31, 63, 127), (UInt8, UInt16, UInt32, UInt64, UInt128))
        c = CMBitSet{T}([2,3,m])
        @test 3 ∈ c
        @test length(c) == 3
        @test 4 ∈ symdiff(c, CMBitSet{T}(4))
        @test 4 == only(setdiff(CMBitSet{T}([m, 4]), c))
        @test intersect(c, CMBitSet{T}([m,3,m])) == CMBitSet{T}([3,m])
        @test isempty(CMBitSet{T}())
        @test collect(c) == [2, 3, m]
        @test minimum(c) == 2
        @test maximum(c) == m
        @test occursin("[2, 3, $m]", string(c))
    end
end

@testset "Iterative gaussian elimination" begin
    gausslengths = PeriodicGraphs.IterativeGaussianEliminationLength([3, 4, 7])
    @test !PeriodicGraphs.gaussian_elimination!(gausslengths, [3, 6])
    @test !PeriodicGraphs.gaussian_elimination!(gausslengths, [4, 7])
    @test !PeriodicGraphs.gaussian_elimination!(gausslengths, [1, 3, 5])
    @test length(gausslengths.track) == length(gausslengths.rings) == 4
    @test first(PeriodicGraphs.gaussian_elimination(gausslengths, [1, 4, 5, 7]))
    @test PeriodicGraphs.gaussian_elimination!(gausslengths, [1, 4, 5, 7])
    @test length(gausslengths.track) == length(gausslengths.rings) == 4
    @test !PeriodicGraphs.gaussian_elimination!(gausslengths, [3, 5, 6, 7])
    notindependent, info = PeriodicGraphs.gaussian_elimination(gausslengths, [2, 3])
    @test !notindependent
    @test !PeriodicGraphs.gaussian_elimination!(gausslengths, [2, 3], notindependent, info)
    @test !PeriodicGraphs.gaussian_elimination!(gausslengths, [1, 2, 3])
    @test length(gausslengths.track) == length(gausslengths.rings) == 7

    gaussnone = PeriodicGraphs.IterativeGaussianElimination([3, 4, 7])
    @test !PeriodicGraphs.gaussian_elimination!(gaussnone, [3, 6])
    @test !PeriodicGraphs.gaussian_elimination!(gaussnone, [4, 7])
    @test !PeriodicGraphs.gaussian_elimination!(gaussnone, [1, 3, 5])
    @test length(gaussnone.rings) == 4
    @test first(PeriodicGraphs.gaussian_elimination(gaussnone, [1, 4, 5, 7]))
    @test PeriodicGraphs.gaussian_elimination!(gaussnone, [1, 4, 5, 7])
    @test !PeriodicGraphs.gaussian_elimination!(gaussnone, [3, 5, 6, 7])
    @test !PeriodicGraphs.gaussian_elimination!(gaussnone, [2, 3])
    @test !first(PeriodicGraphs.gaussian_elimination(gaussnone, [1, 2, 3]))
    @test !PeriodicGraphs.gaussian_elimination!(gaussnone, [1, 2, 3])
    @test PeriodicGraphs.gaussian_elimination!(gaussnone, [1, 2, 3])
    @test gaussnone.shortcuts == Int32[4, 6, 1, 2, 5, 3, 7]
    length(gaussnone.rings) == 7

    gausstrack = PeriodicGraphs.IterativeGaussianEliminationDecomposition([3, 4, 7])
    @test !PeriodicGraphs.gaussian_elimination!(gausstrack, [3, 6])
    @test !PeriodicGraphs.gaussian_elimination!(gausstrack, [4, 7])
    @test !PeriodicGraphs.gaussian_elimination!(gausstrack, [1, 3, 5])
    @test length(gausstrack.rings) == 4
    notindependent, info = PeriodicGraphs.gaussian_elimination(gausstrack, [1, 4, 5, 7])
    @test notindependent
    @test PeriodicGraphs.gaussian_elimination!(gausstrack, [1, 4, 5, 7], notindependent, info)
    @test PeriodicGraphs.retrieve_track!(gausstrack) == Int32[5, 4, 1]
    @test !PeriodicGraphs.gaussian_elimination!(gausstrack, [3, 5, 6, 7])
    @test !PeriodicGraphs.gaussian_elimination!(gausstrack, [2, 3])
    @test !first(PeriodicGraphs.gaussian_elimination(gausstrack, [1, 2, 3]))
    @test !PeriodicGraphs.gaussian_elimination!(gausstrack, [1, 2, 3])
    @test PeriodicGraphs.gaussian_elimination!(gausstrack, [1, 2, 3])
    @test PeriodicGraphs.retrieve_track!(gausstrack) == Int32[9, 8]
    @test gausstrack.shortcuts == Int32[4, 7, 1, 2, 6, 3, 8]
    @test length(gausstrack.rings) == 9

    @test PeriodicGraphs.IterativeGaussianElimination().track == PeriodicGraphs.IterativeGaussianElimination{Nothing}().track == nothing
end

const cpi = PeriodicGraph("2 1 2 0 0 1 3 0 0 1 4 0 0 1 5 0 0 2 3 0 -1 2 5 1 -1 3 4 1 0 4 5 0 -1");
const nab = PeriodicGraph("3 1 2 0 -1 0 1 3 0 -1 -1 1 4 -1 0 0 1 5 -1 0 -1 2 3 0 0 -1 2 4 -1 0 0 2 5 -1 0 -1 3 4 0 0 0 3 5 0 0 0 4 5 0 0 0");
const cai = PeriodicGraph("3 1 1 1 0 0 1 2 0 0 0 1 3 0 0 0 1 4 0 0 0 1 5 0 0 0 2 2 1 0 0 2 5 0 -1 0 2 6 0 0 0 3 3 1 0 0 3 4 0 0 1 3 6 0 0 1 4 4 1 0 0 4 6 0 1 0 5 5 1 0 0 5 6 0 1 1 6 6 1 0 0");
const sny = PeriodicGraph("3 1 1 1 0 0 1 2 0 0 0 1 2 1 0 0 1 3 0 0 0 1 4 0 0 0 2 2 1 0 0 2 5 0 0 0 2 6 0 0 0 3 3 1 0 0 3 4 0 0 0 3 6 0 0 1 3 6 1 0 1 4 4 1 0 0 4 5 0 1 0 4 5 1 1 0 5 5 1 0 0 5 6 0 0 0 6 6 1 0 0");
const wno = PeriodicGraph("3 1 2 0 0 0 1 2 0 1 0 1 3 -1 0 0 1 3 0 0 0 1 4 0 0 0 1 5 0 0 0 1 6 0 0 0 2 3 0 -1 0 2 4 0 0 0 2 5 0 -1 -1 2 6 0 0 0 2 6 1 -1 -1 3 4 0 1 1 3 4 1 0 0 3 5 0 0 0 3 6 1 0 -1 4 5 0 -1 -1 4 5 0 0 -1 4 6 0 0 -1 5 6 0 0 0 5 6 1 0 0");
const txt = PeriodicGraph("3 1 2 -2 1 1 1 2 -1 1 1 1 2 0 0 0 1 3 -2 1 1 1 3 -1 1 0 1 3 0 0 0 1 4 -2 1 1 1 4 -1 0 1 1 4 0 0 0 1 5 -2 1 1 1 5 -1 0 0 1 5 0 0 0 1 6 -1 1 0 1 6 -1 1 1 1 6 0 0 0 1 7 -1 0 1 1 7 -1 1 1 1 7 0 0 0 1 8 -1 0 0 1 8 -1 1 0 1 8 0 0 0 1 9 -1 0 0 1 9 -1 0 1 1 9 0 0 0 1 10 -1 1 1 1 10 0 0 0 1 10 0 1 0 1 11 -1 1 1 1 11 0 0 0 1 11 0 0 1 1 12 -1 1 0 1 12 0 0 0 1 12 0 1 0 1 13 -1 0 1 1 13 0 0 0 1 13 0 0 1");
const afi = PeriodicGraph("3 1 8 0 1 0 1 12 0 0 0 1 14 0 0 0 1 21 0 0 0 2 9 0 0 0 2 12 0 0 0 2 13 -1 0 0 2 22 0 0 0 3 9 0 0 0 3 11 0 0 0 3 14 0 0 0 3 19 0 0 0 4 11 0 0 0 4 13 0 0 0 4 14 0 0 0 4 18 0 0 0 5 11 0 0 -1 5 16 0 0 0 5 18 0 0 0 5 19 0 0 0 6 8 0 0 0 6 11 0 0 0 6 12 0 -1 0 6 16 0 0 0 7 8 0 0 0 7 9 1 0 0 7 13 0 0 0 7 17 0 0 0 8 15 0 0 1 9 10 0 0 1 10 17 -1 0 0 10 19 0 0 0 10 22 0 0 0 12 23 0 0 1 13 24 0 0 1 14 20 0 0 1 15 16 0 0 0 15 17 0 0 0 15 21 0 -1 0 16 23 0 -1 0 17 24 0 0 0 18 20 0 0 0 18 24 0 0 0 19 20 0 0 0 20 21 0 0 0 21 23 0 0 0 22 23 0 0 0 22 24 -1 0 0");
const lta = PeriodicGraph("3 1 3 0 0 0 1 7 0 -1 0 1 9 0 0 -1 1 21 -1 0 -1 2 4 0 0 0 2 8 0 -1 0 2 9 0 0 0 2 21 -1 0 0 3 5 0 0 0 3 10 0 0 -1 3 22 -1 0 -1 4 6 0 0 0 4 10 0 0 0 4 22 -1 0 0 5 6 0 0 0 5 11 0 0 0 5 23 -1 0 0 6 12 0 0 0 6 24 -1 0 0 7 8 0 0 0 7 11 0 0 0 7 23 -1 0 0 8 12 0 0 0 8 24 -1 0 0 9 10 0 0 0 9 16 0 -1 0 10 13 0 0 0 11 12 0 0 0 11 14 0 0 0 12 15 0 0 0 13 14 0 0 1 13 15 0 0 0 13 17 0 0 0 14 16 0 0 -1 14 18 0 0 0 15 16 0 0 0 15 19 0 0 0 16 20 0 0 0 17 18 0 0 1 17 19 0 0 0 17 22 0 0 0 18 20 0 0 -1 18 23 0 0 0 19 20 0 0 0 19 24 0 0 0 20 21 0 1 0 21 22 0 0 0 23 24 0 0 0");
const mtn = PeriodicGraph("3 1 2 0 0 0 1 8 0 0 -1 1 11 0 -1 0 1 29 -1 0 0 2 7 0 0 0 2 12 0 -1 0 2 30 -1 0 0 3 5 0 0 0 3 7 0 0 0 3 9 0 0 0 3 31 -1 0 0 4 6 0 0 0 4 8 0 0 0 4 9 0 0 0 4 32 -1 0 0 5 10 0 0 -1 5 11 0 0 0 5 33 -1 0 0 6 10 0 0 0 6 12 0 0 0 6 34 -1 0 0 7 13 0 0 -1 7 15 0 0 0 8 13 0 0 0 8 16 0 0 0 9 14 0 0 0 9 17 0 0 0 10 13 0 0 0 10 18 0 0 0 11 14 0 0 0 11 20 0 0 -1 12 14 0 0 0 12 19 0 0 0 13 21 0 0 0 14 22 0 0 0 15 17 0 0 0 15 19 0 -1 0 15 23 0 0 0 16 17 0 0 0 16 20 0 -1 0 16 24 0 0 0 17 25 0 0 0 18 19 0 0 0 18 20 0 0 0 18 26 0 0 0 19 27 0 0 0 20 28 0 0 0 21 23 0 0 1 21 24 0 0 0 21 26 0 0 0 22 25 0 0 0 22 27 0 0 0 22 28 0 0 -1 23 29 0 0 0 23 33 0 -1 0 24 30 0 0 0 24 34 0 -1 0 25 29 0 0 0 25 30 0 0 0 26 31 0 0 0 26 32 0 0 0 27 31 0 0 0 27 33 0 0 0 28 32 0 0 0 28 34 0 0 0 29 32 0 0 -1 30 31 0 0 0 33 34 0 0 -1");

@testset "Particular examples" begin
    @test is_well_formed(cpi)
    @test nv(cpi) == 5 && ne(cpi) == 8
    @test degree(cpi) == [4; fill(3, 4)]
    @test is_well_formed(nab)
    @test ne(nab) == 2*nv(nab) == 10
    @test all(==(4), degree(nab))
    @test is_well_formed(cai)
    @test nv(cai) == 6 && ne(cai) == 16
    @test degree(cai) == [6; fill(5, 4); 6]
    @test is_well_formed(sny)
    @test nv(sny) == 6 && ne(sny) == 18
    @test is_well_formed(wno)
    @test nv(wno) == 6 && ne(wno) == 21
    @test is_well_formed(txt)
    @test nv(txt) == 13 && ne(txt) == 36
    @test is_well_formed(afi)
    @test ne(afi) == 2*nv(afi) == 48
    @test all(==(4), degree(afi))
    @test is_well_formed(lta)
    @test ne(lta) == 2*nv(lta) == 48
    @test all(==(4), degree(lta))
    @test is_well_formed(mtn)
    @test ne(mtn) == 2*nv(mtn) == 68
    @test all(==(4), degree(mtn))
end

@testset "Ring statistics" begin
    dag, vertexnums = PeriodicGraphs.arcs_list(cpi, 1, 10, nothing, Bool[0;1;1;0;0])
    @test string(dag[1]) == "JunctionNode{$Int}([])"
    rea = PeriodicGraphs.RingsEndingAt(dag, 12)
    @test_throws MethodError length(rea)
    @test isempty(rea)
    dag, vertexnums = PeriodicGraphs.arcs_list(cpi, 1, 10)
    rea = PeriodicGraphs.RingsEndingAt(dag, 12)
    @test only(rea) == [5, 1, 4, 10, 12]
    @test eltype(rea) == Vector{Int}

    rings_nab = RingAttributions(nab)
    lens_nab = [sort!(length.(x)) for x in rings_nab]
    @test lens_nab == [[3; 4; 8; 8; fill(9, 14)], [3; 4; 8; 8; fill(9, 14)], [3; 3; fill(9, 16)], [3; 4; 8; 8; fill(9, 14)], [3; 4; 8; 8; fill(9, 14)]]
    for r in permutations(1:5)
        nabr = nab[r]
        rings_nab_r = RingAttributions(nab[r])
        @test length(rings_nab_r.rings) == 12
        @test sort!(length.(rings_nab_r.rings)) == [3; 3; 4; 8; fill(9, 8)]
        for (j, ringaroundi) in zip(r, rings_nab_r)
            @test sort!(length.(ringaroundi)) == lens_nab[j]
        end
    end
    nab_rotated = nab[[3,1,2,4,5]]
    rings_nab_rotated = RingAttributions(nab_rotated)
    @test sort!(RingAttributions(nab, true).rings) == sort(rings_nab.rings)
    @test sort!(RingAttributions(nab_rotated, true).rings) == sort(rings_nab_rotated.rings)

    rings_cai = RingAttributions(cai)
    cai_rotated = cai[[3,1,2,4,5,6]]
    rings_cai_rotated = RingAttributions(cai_rotated)
    @test length(rings_cai.rings) == length(rings_cai_rotated.rings) == 14
    @test length(RingAttributions(cai_rotated, true).rings) == length(RingAttributions(cai, true).rings)

    rings_sny_strong = RingAttributions(sny, true, 2) # test for orig_1
    @test length(rings_sny_strong.rings) == 14
    @test all(==(8), length.(rings_sny_strong.attrs))
    @test sort!(length.(rings_sny_strong.rings)) == [fill(3, 8); fill(4, 6)]
    @test sort!(length.(RingAttributions(sny, true).rings)) == [fill(3, 8); fill(4, 6); fill(12, 20)]

    for r in permutations(1:6)
        rings_wno_r = RingAttributions(wno[r], false, 5)
        @test length(rings_wno_r.rings) == 35
        @test all(==(29), length.(rings_wno_r.attrs))
        @test sort!(length.(rings_wno_r.rings)) == [fill(3, 12); 4; 4; 4; fill(5, 6); 6; 6; fill(7, 12)]
    end

    rings_txt = RingAttributions(txt, 4)
    _r1 = [2; 1; 3:13]
    rings_txt_r1 = RingAttributions(txt[_r1], false, 4)
    _r2 = [1; 8; 3:7; 2; 9:13]
    rings_txt_r2 = RingAttributions(txt[_r2], 4)
    _r3 = [1; 2; 8; 4:7; 3; 9:13]
    rings_txt_r3 = RingAttributions(txt[_r3], false, 4)
    @test length(rings_txt.rings) == length(rings_txt_r1.rings) == length(rings_txt_r2.rings) == length(rings_txt_r3.rings) == 1566
    @test length.(rings_txt) == [5208; fill(434, 12)]
    for i in 1:13
        @test length(rings_txt[i]) == length(rings_txt_r1[_r1[i]]) == length(rings_txt_r2[_r2[i]]) == length(rings_txt_r3[_r3[i]])
    end

    @test length(rings(afi, 7)[1]) == 72

    rings_lta = RingAttributions(lta, false, 9)
    @test length(rings_lta.rings) == 113
    @test all(==(59), length.(rings_lta.attrs))
    rings_lta_strong = RingAttributions(lta, true, 9)
    @test length(rings_lta_strong.rings) == 29
    for i in 1:nv(lta)
        rs = collect(rings_lta_strong[i])
        @test length(rs) == 6
        sort!(rs; by=length)
        @test length(rs[1]) == length(rs[2]) == length(rs[3]) == 4
        @test length(rs[4]) == length(rs[5]) == 6
        @test length(rs[6]) == 8
    end

    rings_mtn = RingAttributions(mtn)
    @test length(rings_mtn.rings) == 64
    @test sort!(length.(rings_mtn.attrs)) == [fill(12, 24); fill(15, 8); 18; 18]
    @test sort!(length.(rings_mtn[34])) == [fill(5,5); 6; fill(10, 6)]
    @test sort!(length.(rings_mtn[10])) == [fill(5,6); fill(10, 9)]
    @test sort!(length.(rings_mtn[13])) == [fill(5,6); fill(10, 12)]
    @test eltype(rings_mtn) == PeriodicGraphs.RingIncluding{3}
    str_rings_mtn = sprint(io -> show(io, MIME("text/plain"), rings_mtn))
    @test occursin("rings per node: [12", str_rings_mtn)
    @test eltype(rings_mtn[1]) == PeriodicGraphs.OffsetVertexIterator{3}
    @test isempty(rings(mtn, 0)[1])
    @test isempty(strong_rings(mtn, 0)[1])

    house = PeriodicGraph{0}("0  1 2  2 3  3 4  4 1
                                 1 5  2 6  3 7  4 8
                                 5 6  6 7  7 8  8 5
                                 5 9  6 9  7 9  8 9   ");
    @test length(rings(house, 5)[1]) == 14
    @test length(strong_rings(house, 5)[1]) == 9

    @test PeriodicGraphs.symdiff_cycles([3,4], [3,5,6]) == PeriodicGraphs.symdiff_cycles([4,5,6,7,8], [7,8]) == [4,5,6]
    @test PeriodicGraphs.intersect_cycles([3,4], [3,5,6]) == PeriodicGraphs.intersect_cycles([2,3,5], [1,3]) == [3]
    @test PeriodicGraphs.intersect_cycles([3,4,5], [4,5,6]) == PeriodicGraphs.intersect_cycles([4,5,7], [2,4,5,6]) == [4,5]
    @test PeriodicGraphs.union_cycles([3,4], [4,5,6]) == PeriodicGraphs.union_cycles([4,5,6], [3,4]) == [3,4,5,6]
    @test PeriodicGraphs.union_cycles([4,6], [3,5,6]) == PeriodicGraphs.union_cycles([3,5,6], [3,4,6]) == [3,4,5,6]

    # keep track of the limitations
    @test_throws ErrorException rings(lta, 63)
    very_high_degree = PeriodicGraph{0}(128)
    for i in 2:128
        add_edge!(very_high_degree, 1, PeriodicVertex{0}(i))
    end
    @test_throws ErrorException rings(very_high_degree, 1)
end


module PeriodicSymmetries
    # heavily inspired from PeriodicGraphEmbeddings.jl
    using PeriodicGraphs, Graphs, StaticArrays

    struct PeriodicSymmetry3D{T} <: PeriodicGraphs.AbstractSymmetry
        vmap::SubArray{PeriodicVertex3D,1,Matrix{PeriodicVertex3D},Tuple{Base.Slice{Base.OneTo{Int}},Int},true}
        rotation::SMatrix{3,3,Int,9}
        translation::SVector{3,T}
    end
    Base.getindex(symm::PeriodicSymmetry3D, i::Integer) = symm.vmap[i].v
    function Base.getindex(symm::PeriodicSymmetry3D, x::PeriodicVertex3D)
        dst = symm.vmap[x.v]
        _ofs = muladd(symm.rotation, x.ofs, dst.ofs)
        PeriodicVertex3D(dst.v, _ofs)
    end
    function Base.isequal(x::PeriodicSymmetry3D{T}, y::PeriodicSymmetry3D{T}) where T
        x.vmap == y.vmap && x.rotation == y.rotation && x.translation == y.translation
    end

    struct SymmetryGroup3D{T} <: PeriodicGraphs.AbstractSymmetryGroup{PeriodicSymmetry3D{T}}
        vmaps::Matrix{PeriodicVertex3D}
        rotations::Vector{SMatrix{3,3,Int,9}}
        translations::Vector{SVector{3,T}}
        uniquemap::Vector{Int}
        uniques::Vector{Int}
        hasmirror::Bool
    end

    (s::SymmetryGroup3D)(i::Integer) = s.uniques[s.uniquemap[i]]
    Base.unique(s::SymmetryGroup3D) = s.uniques
    function Base.getindex(s::SymmetryGroup3D{T}, i::Integer) where {T}
        PeriodicSymmetry3D{T}((@view s.vmaps[:,i]), s.rotations[i], s.translations[i])
    end
    Base.iterate(s::SymmetryGroup3D, state=1) = state > length(s) ? nothing : (s[state], state+1)
    Base.eltype(::Type{SymmetryGroup3D{T}}) where {T} = PeriodicSymmetry3D{T}
    Base.length(s::SymmetryGroup3D) = length(s.rotations)
    function Base.one(s::SymmetryGroup3D{T}) where T
        n = length(s.uniquemap)
        PeriodicSymmetry3D{T}((@view reshape(collect(PeriodicVertex3D.(Base.OneTo(n))), n, 1)[:,1]),
                            one(SMatrix{3,3,Int,9}), zero(SVector{3,T}))
    end

end

const symmetries_sny = PeriodicSymmetries.SymmetryGroup3D{Rational{Int16}}(
    PeriodicVertex3D[(2, (0,-1,-1)) (2, (1,-1,-1)) (1, (-1,0,0)) (2, (0,-1,-1)) (1, (0,0,0)) (1, (-1,0,0)) (2, (1,-1,-1)); (1, (0,-1,-1)) (1, (0,-1,-1)) (2, (0,0,0)) (1, (0,-1,-1)) (2, (0,0,0)) (2, (0,0,0)) (1, (0,-1,-1)); (6, (0,-1,-1)) (6, (1,-1,-1)) (3, (-1,0,0)) (5, (0,-1,-1)) (4, (0,0,0)) (4, (-1,0,0)) (5, (1,-1,-1)); (5, (0,-1,-1)) (5, (1,-1,-1)) (4, (-1,0,0)) (6, (0,-1,-1)) (3, (0,0,0)) (3, (-1,0,0)) (6, (1,-1,-1)); (4, (0,-1,-1)) (4, (0,-1,-1)) (5, (0,0,0)) (3, (0,-1,-1)) (6, (0,0,0)) (6, (0,0,0)) (3, (0,-1,-1)); (3, (0,-1,-1)) (3, (0,-1,-1)) (6, (0,0,0)) (4, (0,-1,-1)) (5, (0,0,0)) (5, (0,0,0)) (4, (0,-1,-1))],
    SMatrix{3, 3, Int, 9}[[-1 0 0; 0 -1 0; 0 0 -1], [1 0 0; 0 -1 0; 0 0 -1], [-1 0 0; 0 1 0; 0 0 1], [-1 0 0; 0 0 -1; 0 -1 0], [1 0 0; 0 0 1; 0 1 0], [-1 0 0; 0 0 1; 0 1 0], [1 0 0; 0 0 -1; 0 -1 0]],
    SVector{3, Rational{Int16}}[[1//2, 0//1, 0//1], [1//2, 0//1, 0//1], [0//1, 0//1, 0//1], [1//2, 0//1, 0//1], [0//1, 0//1, 0//1], [0//1, 0//1, 0//1], [1//2, 0//1, 0//1]],
    [1, 1, 2, 2, 2, 2],
    [1, 3],
    true
)

const symmetries_afi = PeriodicSymmetries.SymmetryGroup3D{Rational{Int64}}(
    PeriodicVertex3D[(16, (-1,-1,-1)) (6, (-1,-1,0)) (21, (0,0,-1)) (9, (-1,-1,-1)) (24, (0,0,1)) (13, (0,0,-1)) (10, (-1,-1,1)); (17, (-1,-1,-1)) (7, (-1,-1,0)) (22, (0,0,-1)) (12, (-1,-1,-1)) (15, (0,0,1)) (8, (0,0,-1)) (23, (-1,-1,1)); (18, (-1,-1,-1)) (4, (-1,-1,0)) (19, (0,0,-1)) (14, (-1,-1,-1)) (5, (0,0,1)) (11, (0,0,-1)) (20, (-1,-1,1)); (19, (-1,-1,-1)) (3, (-1,-1,0)) (18, (0,0,-1)) (11, (-1,-1,-1)) (20, (0,0,1)) (14, (0,0,-1)) (5, (-1,-1,1)); (14, (-1,-1,-1)) (20, (-1,-1,0)) (11, (0,0,-1)) (18, (-1,-1,0)) (3, (0,0,0)) (19, (0,0,0)) (4, (-1,-1,0)); (21, (-1,-1,-1)) (1, (-1,-1,0)) (16, (0,0,-1)) (13, (-1,-1,-1)) (10, (0,0,1)) (9, (0,0,-1)) (24, (-1,-1,1)); (22, (-1,-1,-1)) (2, (-1,-1,0)) (17, (0,0,-1)) (8, (-1,-1,-1)) (23, (0,0,1)) (12, (0,0,-1)) (15, (-1,-1,1)); (23, (-1,-1,-1)) (12, (-1,-1,0)) (15, (0,0,-1)) (7, (-1,-1,-1)) (22, (0,0,1)) (2, (0,0,-1)) (17, (-1,-1,1)); (24, (-1,-1,-1)) (13, (-1,-1,0)) (10, (0,0,-1)) (1, (-1,-1,-1)) (16, (0,0,1)) (6, (0,0,-1)) (21, (-1,-1,1)); (13, (-1,-1,-1)) (24, (-1,-1,0)) (9, (0,0,-1)) (21, (-1,-1,0)) (6, (0,0,0)) (16, (0,0,0)) (1, (-1,-1,0)); (20, (-1,-1,-1)) (14, (-1,-1,0)) (5, (0,0,-1)) (4, (-1,-1,-1)) (19, (0,0,1)) (3, (0,0,-1)) (18, (-1,-1,1)); (15, (-1,-1,-1)) (8, (-1,-1,0)) (23, (0,0,-1)) (2, (-1,-1,-1)) (17, (0,0,1)) (7, (0,0,-1)) (22, (-1,-1,1)); (10, (-1,-1,-1)) (9, (-1,-1,0)) (24, (0,0,-1)) (6, (-1,-1,-1)) (21, (0,0,1)) (1, (0,0,-1)) (16, (-1,-1,1)); (5, (-1,-1,-1)) (11, (-1,-1,0)) (20, (0,0,-1)) (3, (-1,-1,-1)) (18, (0,0,1)) (4, (0,0,-1)) (19, (-1,-1,1)); (12, (-1,-1,-1)) (23, (-1,-1,0)) (8, (0,0,-1)) (17, (-1,-1,0)) (2, (0,0,0)) (22, (0,0,0)) (7, (-1,-1,0)); (1, (-1,-1,-1)) (21, (-1,-1,0)) (6, (0,0,-1)) (24, (-1,-1,0)) (9, (0,0,0)) (10, (0,0,0)) (13, (-1,-1,0)); (2, (-1,-1,-1)) (22, (-1,-1,0)) (7, (0,0,-1)) (15, (-1,-1,0)) (12, (0,0,0)) (23, (0,0,0)) (8, (-1,-1,0)); (3, (-1,-1,-1)) (19, (-1,-1,0)) (4, (0,0,-1)) (5, (-1,-1,0)) (14, (0,0,0)) (20, (0,0,0)) (11, (-1,-1,0)); (4, (-1,-1,-1)) (18, (-1,-1,0)) (3, (0,0,-1)) (20, (-1,-1,0)) (11, (0,0,0)) (5, (0,0,0)) (14, (-1,-1,0)); (11, (-1,-1,-1)) (5, (-1,-1,0)) (14, (0,0,-1)) (19, (-1,-1,0)) (4, (0,0,0)) (18, (0,0,0)) (3, (-1,-1,0)); (6, (-1,-1,-1)) (16, (-1,-1,0)) (1, (0,0,-1)) (10, (-1,-1,0)) (13, (0,0,0)) (24, (0,0,0)) (9, (-1,-1,0)); (7, (-1,-1,-1)) (17, (-1,-1,0)) (2, (0,0,-1)) (23, (-1,-1,0)) (8, (0,0,0)) (15, (0,0,0)) (12, (-1,-1,0)); (8, (-1,-1,-1)) (15, (-1,-1,0)) (12, (0,0,-1)) (22, (-1,-1,0)) (7, (0,0,0)) (17, (0,0,0)) (2, (-1,-1,0)); (9, (-1,-1,-1)) (10, (-1,-1,0)) (13, (0,0,-1)) (16, (-1,-1,0)) (1, (0,0,0)) (21, (0,0,0)) (6, (-1,-1,0))],
    SMatrix{3, 3, Int64, 9}[[-1 0 0; 0 -1 0; 0 0 -1], [-1 0 0; 0 -1 0; 0 0 1], [1 0 0; 0 1 0; 0 0 -1], [0 -1 0; -1 0 0; 0 0 -1], [0 1 0; 1 0 0; 0 0 1], [0 1 0; 1 0 0; 0 0 -1], [0 -1 0; -1 0 0; 0 0 1]],
    SVector{3, Rational{Int64}}[[0//1, 0//1, 0//1], [0//1, 0//1, 0//1], [0//1, 0//1, 0//1], [0//1, 0//1, 1//2], [0//1, 0//1, 1//2], [0//1, 0//1, 1//2], [0//1, 0//1, 1//2]],
    [1, 2, 3, 3, 3, 1, 2, 2, 1, 1, 3, 2, 1, 3, 2, 1, 2, 3, 3, 3, 1, 2, 2, 1],
    [1, 2, 3],
    true
)

const symmetries_lta = PeriodicSymmetries.SymmetryGroup3D{Rational{Int32}}(
    PeriodicVertex3D[(4, (0,0,0)) (18, (0,0,0)) (15, (0,0,0)) (3, (0,0,0)) (2, (0,0,0)) (14, (0,0,0)) (19, (0,0,0)) (4, (0,0,0)) (1, (0,0,0)) (19, (0,0,0)) (14, (0,0,0)) (2, (0,0,0)) (3, (0,0,0)) (15, (0,0,0)) (18, (0,0,0)) (11, (0,0,0)) (24, (0,0,0)) (10, (0,0,0)) (21, (0,0,0)) (12, (0,0,0)) (23, (0,0,0)) (9, (0,0,0)) (22, (0,0,0)) (24, (0,0,0)) (11, (0,0,0)) (22, (0,0,0)) (9, (0,0,0)) (23, (0,0,0)) (12, (0,0,0)) (21, (0,0,0)) (10, (0,0,0)) (16, (0,0,0)) (17, (0,0,0)) (8, (0,0,0)) (5, (0,0,0)) (20, (0,0,0)) (13, (0,0,0)) (7, (0,0,0)) (6, (0,0,0)) (17, (0,0,0)) (16, (0,0,0)) (6, (0,0,0)) (7, (0,0,0)) (13, (0,0,0)) (20, (0,0,0)) (5, (0,0,0)) (8, (0,0,0)); (3, (0,0,0)) (19, (0,0,0)) (14, (0,0,0)) (4, (0,0,0)) (1, (0,0,0)) (15, (0,0,0)) (18, (0,0,0)) (3, (0,0,0)) (2, (0,0,0)) (18, (0,0,0)) (15, (0,0,0)) (1, (0,0,0)) (4, (0,0,0)) (14, (0,0,0)) (19, (0,0,0)) (23, (0,0,0)) (12, (0,0,0)) (22, (0,0,0)) (9, (0,0,0)) (24, (0,0,0)) (11, (0,0,0)) (21, (0,0,0)) (10, (0,0,0)) (12, (0,0,0)) (23, (0,0,0)) (10, (0,0,0)) (21, (0,0,0)) (11, (0,0,0)) (24, (0,0,0)) (9, (0,0,0)) (22, (0,0,0)) (13, (0,1,0)) (20, (0,-1,0)) (6, (0,1,0)) (7, (0,-1,0)) (17, (0,1,0)) (16, (0,-1,0)) (5, (0,1,0)) (8, (0,-1,0)) (20, (0,-1,0)) (13, (0,1,0)) (8, (0,-1,0)) (5, (0,1,0)) (16, (0,-1,0)) (17, (0,1,0)) (7, (0,-1,0)) (6, (0,1,0)); (2, (0,0,0)) (14, (0,0,0)) (19, (0,0,0)) (1, (0,0,0)) (4, (0,0,0)) (18, (0,0,0)) (15, (0,0,0)) (2, (0,0,0)) (3, (0,0,0)) (15, (0,0,0)) (18, (0,0,0)) (4, (0,0,0)) (1, (0,0,0)) (19, (0,0,0)) (14, (0,0,0)) (12, (0,0,0)) (23, (0,0,0)) (9, (0,0,0)) (22, (0,0,0)) (11, (0,0,0)) (24, (0,0,0)) (10, (0,0,0)) (21, (0,0,0)) (23, (0,0,0)) (12, (0,0,0)) (21, (0,0,0)) (10, (0,0,0)) (24, (0,0,0)) (11, (0,0,0)) (22, (0,0,0)) (9, (0,0,0)) (20, (0,0,0)) (13, (0,0,0)) (7, (0,0,0)) (6, (0,0,0)) (16, (0,0,0)) (17, (0,0,0)) (8, (0,0,0)) (5, (0,0,0)) (13, (0,0,0)) (20, (0,0,0)) (5, (0,0,0)) (8, (0,0,0)) (17, (0,0,0)) (16, (0,0,0)) (6, (0,0,0)) (7, (0,0,0)); (1, (0,0,0)) (15, (0,0,0)) (18, (0,0,0)) (2, (0,0,0)) (3, (0,0,0)) (19, (0,0,0)) (14, (0,0,0)) (1, (0,0,0)) (4, (0,0,0)) (14, (0,0,0)) (19, (0,0,0)) (3, (0,0,0)) (2, (0,0,0)) (18, (0,0,0)) (15, (0,0,0)) (24, (0,0,0)) (11, (0,0,0)) (21, (0,0,0)) (10, (0,0,0)) (23, (0,0,0)) (12, (0,0,0)) (22, (0,0,0)) (9, (0,0,0)) (11, (0,0,0)) (24, (0,0,0)) (9, (0,0,0)) (22, (0,0,0)) (12, (0,0,0)) (23, (0,0,0)) (10, (0,0,0)) (21, (0,0,0)) (17, (0,1,0)) (16, (0,-1,0)) (5, (0,1,0)) (8, (0,-1,0)) (13, (0,1,0)) (20, (0,-1,0)) (6, (0,1,0)) (7, (0,-1,0)) (16, (0,-1,0)) (17, (0,1,0)) (7, (0,-1,0)) (6, (0,1,0)) (20, (0,-1,0)) (13, (0,1,0)) (8, (0,-1,0)) (5, (0,1,0)); (8, (0,-1,0)) (11, (0,0,0)) (24, (0,0,0)) (7, (0,-1,0)) (6, (0,0,0)) (23, (0,0,0)) (12, (0,0,0)) (8, (0,-1,0)) (5, (0,0,0)) (12, (0,0,0)) (23, (0,0,0)) (6, (0,0,0)) (7, (0,-1,0)) (24, (0,0,0)) (11, (0,0,0)) (15, (0,0,0)) (18, (0,0,0)) (16, (0,-1,0)) (17, (0,0,0)) (14, (0,0,0)) (19, (0,0,0)) (13, (0,0,0)) (20, (0,-1,0)) (18, (0,0,0)) (15, (0,0,0)) (20, (0,-1,0)) (13, (0,0,0)) (19, (0,0,0)) (14, (0,0,0)) (17, (0,0,0)) (16, (0,-1,0)) (21, (0,1,0)) (10, (0,0,0)) (1, (0,1,0)) (4, (0,0,0)) (9, (0,1,0)) (22, (0,0,0)) (2, (0,1,0)) (3, (0,0,0)) (10, (0,0,0)) (21, (0,1,0)) (3, (0,0,0)) (2, (0,1,0)) (22, (0,0,0)) (9, (0,1,0)) (4, (0,0,0)) (1, (0,1,0)); (7, (0,-1,0)) (12, (0,0,0)) (23, (0,0,0)) (8, (0,-1,0)) (5, (0,0,0)) (24, (0,0,0)) (11, (0,0,0)) (7, (0,-1,0)) (6, (0,0,0)) (11, (0,0,0)) (24, (0,0,0)) (5, (0,0,0)) (8, (0,-1,0)) (23, (0,0,0)) (12, (0,0,0)) (19, (0,0,0)) (14, (0,0,0)) (20, (0,-1,0)) (13, (0,0,0)) (18, (0,0,0)) (15, (0,0,0)) (17, (0,0,0)) (16, (0,-1,0)) (14, (0,0,0)) (19, (0,0,0)) (16, (0,-1,0)) (17, (0,0,0)) (15, (0,0,0)) (18, (0,0,0)) (13, (0,0,0)) (20, (0,-1,0)) (22, (0,1,0)) (9, (0,0,0)) (3, (0,1,0)) (2, (0,0,0)) (10, (0,1,0)) (21, (0,0,0)) (4, (0,1,0)) (1, (0,0,0)) (9, (0,0,0)) (22, (0,1,0)) (1, (0,0,0)) (4, (0,1,0)) (21, (0,0,0)) (10, (0,1,0)) (2, (0,0,0)) (3, (0,1,0)); (6, (0,-1,0)) (23, (-1,0,0)) (12, (1,0,0)) (5, (0,-1,0)) (8, (0,0,0)) (11, (1,0,0)) (24, (-1,0,0)) (6, (0,-1,0)) (7, (0,0,0)) (24, (-1,0,0)) (11, (1,0,0)) (8, (0,0,0)) (5, (0,-1,0)) (12, (1,0,0)) (23, (-1,0,0)) (14, (0,0,1)) (19, (0,0,-1)) (13, (0,-1,0)) (20, (0,0,0)) (15, (0,0,-1)) (18, (0,0,1)) (16, (0,0,0)) (17, (0,-1,0)) (19, (0,0,-1)) (14, (0,0,1)) (17, (0,-1,0)) (16, (0,0,0)) (18, (0,0,1)) (15, (0,0,-1)) (20, (0,0,0)) (13, (0,-1,0)) (9, (1,1,0)) (22, (-1,0,0)) (2, (0,1,-1)) (3, (0,0,1)) (21, (-1,1,0)) (10, (1,0,0)) (1, (0,1,1)) (4, (0,0,-1)) (22, (-1,0,0)) (9, (1,1,0)) (4, (0,0,-1)) (1, (0,1,1)) (10, (1,0,0)) (21, (-1,1,0)) (3, (0,0,1)) (2, (0,1,-1)); (5, (0,-1,0)) (24, (-1,0,0)) (11, (1,0,0)) (6, (0,-1,0)) (7, (0,0,0)) (12, (1,0,0)) (23, (-1,0,0)) (5, (0,-1,0)) (8, (0,0,0)) (23, (-1,0,0)) (12, (1,0,0)) (7, (0,0,0)) (6, (0,-1,0)) (11, (1,0,0)) (24, (-1,0,0)) (18, (0,0,1)) (15, (0,0,-1)) (17, (0,-1,0)) (16, (0,0,0)) (19, (0,0,-1)) (14, (0,0,1)) (20, (0,0,0)) (13, (0,-1,0)) (15, (0,0,-1)) (18, (0,0,1)) (13, (0,-1,0)) (20, (0,0,0)) (14, (0,0,1)) (19, (0,0,-1)) (16, (0,0,0)) (17, (0,-1,0)) (10, (1,1,0)) (21, (-1,0,0)) (4, (0,1,-1)) (1, (0,0,1)) (22, (-1,1,0)) (9, (1,0,0)) (3, (0,1,1)) (2, (0,0,-1)) (21, (-1,0,0)) (10, (1,1,0)) (2, (0,0,-1)) (3, (0,1,1)) (9, (1,0,0)) (22, (-1,1,0)) (1, (0,0,1)) (4, (0,1,-1)); (22, (-1,0,-1)) (20, (0,0,0)) (13, (0,0,-1)) (22, (-1,0,0)) (9, (0,0,-1)) (13, (0,0,0)) (20, (0,0,-1)) (10, (0,0,-1)) (21, (-1,0,0)) (17, (0,0,-1)) (16, (0,0,0)) (21, (-1,0,-1)) (10, (0,0,0)) (16, (0,0,-1)) (17, (0,0,0)) (7, (1,0,0)) (6, (0,0,0)) (3, (1,0,1)) (2, (0,0,0)) (6, (1,0,0)) (7, (0,0,0)) (2, (1,0,0)) (3, (0,0,1)) (8, (0,0,0)) (5, (1,0,0)) (4, (0,0,0)) (1, (1,0,1)) (5, (0,0,0)) (8, (1,0,0)) (1, (0,0,1)) (4, (1,0,0)) (14, (0,1,1)) (19, (0,-1,0)) (12, (0,1,0)) (23, (-1,-1,0)) (19, (0,1,0)) (14, (0,-1,1)) (23, (-1,1,0)) (12, (0,-1,0)) (18, (0,-1,1)) (15, (0,1,0)) (24, (-1,-1,0)) (11, (0,1,0)) (15, (0,-1,0)) (18, (0,1,1)) (11, (0,-1,0)) (24, (-1,1,0)); (21, (-1,0,-1)) (16, (0,0,0)) (17, (0,0,-1)) (21, (-1,0,0)) (10, (0,0,-1)) (17, (0,0,0)) (16, (0,0,-1)) (9, (0,0,-1)) (22, (-1,0,0)) (13, (0,0,-1)) (20, (0,0,0)) (22, (-1,0,-1)) (9, (0,0,0)) (20, (0,0,-1)) (13, (0,0,0)) (8, (1,0,0)) (5, (0,0,0)) (1, (1,0,1)) (4, (0,0,0)) (5, (1,0,0)) (8, (0,0,0)) (4, (1,0,0)) (1, (0,0,1)) (7, (0,0,0)) (6, (1,0,0)) (2, (0,0,0)) (3, (1,0,1)) (6, (0,0,0)) (7, (1,0,0)) (3, (0,0,1)) (2, (1,0,0)) (18, (0,1,1)) (15, (0,-1,0)) (11, (0,1,0)) (24, (-1,-1,0)) (15, (0,1,0)) (18, (0,-1,1)) (24, (-1,1,0)) (11, (0,-1,0)) (14, (0,-1,1)) (19, (0,1,0)) (23, (-1,-1,0)) (12, (0,1,0)) (19, (0,-1,0)) (14, (0,1,1)) (12, (0,-1,0)) (23, (-1,1,0)); (24, (-1,-1,0)) (7, (0,0,0)) (6, (1,0,0)) (23, (-1,-1,0)) (12, (0,0,0)) (5, (1,0,0)) (8, (0,0,0)) (12, (0,-1,0)) (23, (-1,0,0)) (6, (0,0,0)) (7, (1,0,0)) (24, (-1,0,0)) (11, (0,-1,0)) (8, (1,0,0)) (5, (0,0,0)) (16, (0,0,0)) (17, (0,0,-1)) (14, (0,-1,1)) (19, (0,0,0)) (13, (0,0,-1)) (20, (0,0,0)) (15, (0,0,0)) (18, (0,-1,1)) (20, (0,0,-1)) (13, (0,0,0)) (19, (0,-1,0)) (14, (0,0,1)) (17, (0,0,0)) (16, (0,0,-1)) (18, (0,0,1)) (15, (0,-1,0)) (1, (1,1,1)) (4, (0,0,0)) (9, (0,1,-1)) (22, (-1,0,0)) (2, (0,1,0)) (3, (1,0,1)) (21, (-1,1,0)) (10, (0,0,-1)) (3, (0,0,1)) (2, (1,1,0)) (22, (-1,0,-1)) (9, (0,1,0)) (4, (1,0,0)) (1, (0,1,1)) (10, (0,0,0)) (21, (-1,1,-1)); (23, (-1,-1,0)) (8, (0,0,0)) (5, (1,0,0)) (24, (-1,-1,0)) (11, (0,0,0)) (6, (1,0,0)) (7, (0,0,0)) (11, (0,-1,0)) (24, (-1,0,0)) (5, (0,0,0)) (8, (1,0,0)) (23, (-1,0,0)) (12, (0,-1,0)) (7, (1,0,0)) (6, (0,0,0)) (20, (0,0,0)) (13, (0,0,-1)) (18, (0,-1,1)) (15, (0,0,0)) (17, (0,0,-1)) (16, (0,0,0)) (19, (0,0,0)) (14, (0,-1,1)) (16, (0,0,-1)) (17, (0,0,0)) (15, (0,-1,0)) (18, (0,0,1)) (13, (0,0,0)) (20, (0,0,-1)) (14, (0,0,1)) (19, (0,-1,0)) (3, (1,1,1)) (2, (0,0,0)) (10, (0,1,-1)) (21, (-1,0,0)) (4, (0,1,0)) (1, (1,0,1)) (22, (-1,1,0)) (9, (0,0,-1)) (1, (0,0,1)) (4, (1,1,0)) (21, (-1,0,-1)) (10, (0,1,0)) (2, (1,0,0)) (3, (0,1,1)) (9, (0,0,0)) (22, (-1,1,-1)); (20, (-1,-1,-1)) (9, (0,1,0)) (22, (0,0,-1)) (20, (-1,-1,0)) (13, (0,0,-1)) (22, (0,0,0)) (9, (0,1,-1)) (16, (0,-1,-1)) (17, (-1,0,0)) (10, (0,0,-1)) (21, (0,1,0)) (17, (-1,0,-1)) (16, (0,-1,0)) (21, (0,1,-1)) (10, (0,0,0)) (2, (1,1,0)) (3, (0,0,0)) (7, (1,-1,1)) (6, (0,0,0)) (3, (1,0,0)) (2, (0,1,0)) (6, (1,0,0)) (7, (0,-1,1)) (1, (0,1,0)) (4, (1,0,0)) (8, (0,-1,0)) (5, (1,0,1)) (4, (0,0,0)) (1, (1,1,0)) (5, (0,0,1)) (8, (1,-1,0)) (23, (0,1,1)) (12, (0,-1,0)) (14, (0,1,0)) (19, (-1,-1,0)) (12, (0,1,0)) (23, (0,-1,1)) (19, (-1,1,0)) (14, (0,-1,0)) (11, (0,-1,1)) (24, (0,1,0)) (18, (-1,-1,0)) (15, (0,1,0)) (24, (0,-1,0)) (11, (0,1,1)) (15, (0,-1,0)) (18, (-1,1,0)); (19, (-1,-1,0)) (1, (0,1,0)) (4, (1,0,0)) (18, (-1,-1,0)) (15, (0,0,0)) (3, (1,0,0)) (2, (0,1,0)) (15, (0,-1,0)) (18, (-1,0,0)) (4, (0,0,0)) (1, (1,1,0)) (19, (-1,0,0)) (14, (0,-1,0)) (2, (1,1,0)) (3, (0,0,0)) (9, (0,1,0)) (22, (0,0,-1)) (11, (0,-1,1)) (24, (0,0,0)) (10, (0,0,-1)) (21, (0,1,0)) (12, (0,0,0)) (23, (0,-1,1)) (21, (0,1,-1)) (10, (0,0,0)) (24, (0,-1,0)) (11, (0,0,1)) (22, (0,0,0)) (9, (0,1,-1)) (23, (0,0,1)) (12, (0,-1,0)) (7, (1,0,1)) (6, (0,0,0)) (16, (0,0,-1)) (17, (-1,0,0)) (8, (0,0,0)) (5, (1,0,1)) (20, (-1,0,0)) (13, (0,0,-1)) (5, (0,0,1)) (8, (1,0,0)) (17, (-1,0,-1)) (16, (0,0,0)) (6, (1,0,0)) (7, (0,0,1)) (13, (0,0,0)) (20, (-1,0,-1)); (18, (-1,-1,0)) (2, (0,1,0)) (3, (1,0,0)) (19, (-1,-1,0)) (14, (0,0,0)) (4, (1,0,0)) (1, (0,1,0)) (14, (0,-1,0)) (19, (-1,0,0)) (3, (0,0,0)) (2, (1,1,0)) (18, (-1,0,0)) (15, (0,-1,0)) (1, (1,1,0)) (4, (0,0,0)) (21, (0,1,0)) (10, (0,0,-1)) (23, (0,-1,1)) (12, (0,0,0)) (22, (0,0,-1)) (9, (0,1,0)) (24, (0,0,0)) (11, (0,-1,1)) (9, (0,1,-1)) (22, (0,0,0)) (12, (0,-1,0)) (23, (0,0,1)) (10, (0,0,0)) (21, (0,1,-1)) (11, (0,0,1)) (24, (0,-1,0)) (5, (1,1,1)) (8, (0,-1,0)) (13, (0,1,-1)) (20, (-1,-1,0)) (6, (0,1,0)) (7, (1,-1,1)) (17, (-1,1,0)) (16, (0,-1,-1)) (7, (0,-1,1)) (6, (1,1,0)) (20, (-1,-1,-1)) (13, (0,1,0)) (8, (1,-1,0)) (5, (0,1,1)) (16, (0,-1,0)) (17, (-1,1,-1)); (17, (-1,-1,-1)) (21, (-1,1,0)) (10, (1,0,-1)) (17, (-1,-1,0)) (16, (0,0,-1)) (10, (1,0,0)) (21, (-1,1,-1)) (13, (0,-1,-1)) (20, (-1,0,0)) (22, (-1,0,-1)) (9, (1,1,0)) (20, (-1,0,-1)) (13, (0,-1,0)) (9, (1,1,-1)) (22, (-1,0,0)) (1, (1,1,1)) (4, (0,0,-1)) (5, (1,-1,1)) (8, (0,0,0)) (4, (1,0,-1)) (1, (0,1,1)) (8, (1,0,0)) (5, (0,-1,1)) (2, (0,1,-1)) (3, (1,0,1)) (6, (0,-1,0)) (7, (1,0,1)) (3, (0,0,1)) (2, (1,1,-1)) (7, (0,0,1)) (6, (1,-1,0)) (11, (1,1,1)) (24, (-1,-1,0)) (15, (0,1,-1)) (18, (-1,-1,1)) (24, (-1,1,0)) (11, (1,-1,1)) (18, (-1,1,1)) (15, (0,-1,-1)) (23, (-1,-1,1)) (12, (1,1,0)) (19, (-1,-1,-1)) (14, (0,1,1)) (12, (1,-1,0)) (23, (-1,1,1)) (14, (0,-1,1)) (19, (-1,1,-1)); (16, (-1,-1,-1)) (10, (0,1,0)) (21, (0,0,-1)) (16, (-1,-1,0)) (17, (0,0,-1)) (21, (0,0,0)) (10, (0,1,-1)) (20, (0,-1,-1)) (13, (-1,0,0)) (9, (0,0,-1)) (22, (0,1,0)) (13, (-1,0,-1)) (20, (0,-1,0)) (22, (0,1,-1)) (9, (0,0,0)) (4, (1,1,0)) (1, (0,0,0)) (8, (1,-1,1)) (5, (0,0,0)) (1, (1,0,0)) (4, (0,1,0)) (5, (1,0,0)) (8, (0,-1,1)) (3, (0,1,0)) (2, (1,0,0)) (7, (0,-1,0)) (6, (1,0,1)) (2, (0,0,0)) (3, (1,1,0)) (6, (0,0,1)) (7, (1,-1,0)) (24, (0,1,1)) (11, (0,-1,0)) (18, (0,1,0)) (15, (-1,-1,0)) (11, (0,1,0)) (24, (0,-1,1)) (15, (-1,1,0)) (18, (0,-1,0)) (12, (0,-1,1)) (23, (0,1,0)) (14, (-1,-1,0)) (19, (0,1,0)) (23, (0,-1,0)) (12, (0,1,1)) (19, (0,-1,0)) (14, (-1,1,0)); (15, (-1,-1,0)) (3, (0,1,0)) (2, (1,0,0)) (14, (-1,-1,0)) (19, (0,0,0)) (1, (1,0,0)) (4, (0,1,0)) (19, (0,-1,0)) (14, (-1,0,0)) (2, (0,0,0)) (3, (1,1,0)) (15, (-1,0,0)) (18, (0,-1,0)) (4, (1,1,0)) (1, (0,0,0)) (10, (0,1,0)) (21, (0,0,-1)) (12, (0,-1,1)) (23, (0,0,0)) (9, (0,0,-1)) (22, (0,1,0)) (11, (0,0,0)) (24, (0,-1,1)) (22, (0,1,-1)) (9, (0,0,0)) (23, (0,-1,0)) (12, (0,0,1)) (21, (0,0,0)) (10, (0,1,-1)) (24, (0,0,1)) (11, (0,-1,0)) (8, (1,0,1)) (5, (0,0,0)) (20, (0,0,-1)) (13, (-1,0,0)) (7, (0,0,0)) (6, (1,0,1)) (16, (-1,0,0)) (17, (0,0,-1)) (6, (0,0,1)) (7, (1,0,0)) (13, (-1,0,-1)) (20, (0,0,0)) (5, (1,0,0)) (8, (0,0,1)) (17, (0,0,0)) (16, (-1,0,-1)); (14, (-1,-1,0)) (4, (0,1,0)) (1, (1,0,0)) (15, (-1,-1,0)) (18, (0,0,0)) (2, (1,0,0)) (3, (0,1,0)) (18, (0,-1,0)) (15, (-1,0,0)) (1, (0,0,0)) (4, (1,1,0)) (14, (-1,0,0)) (19, (0,-1,0)) (3, (1,1,0)) (2, (0,0,0)) (22, (0,1,0)) (9, (0,0,-1)) (24, (0,-1,1)) (11, (0,0,0)) (21, (0,0,-1)) (10, (0,1,0)) (23, (0,0,0)) (12, (0,-1,1)) (10, (0,1,-1)) (21, (0,0,0)) (11, (0,-1,0)) (24, (0,0,1)) (9, (0,0,0)) (22, (0,1,-1)) (12, (0,0,1)) (23, (0,-1,0)) (6, (1,1,1)) (7, (0,-1,0)) (17, (0,1,-1)) (16, (-1,-1,0)) (5, (0,1,0)) (8, (1,-1,1)) (13, (-1,1,0)) (20, (0,-1,-1)) (8, (0,-1,1)) (5, (1,1,0)) (16, (-1,-1,-1)) (17, (0,1,0)) (7, (1,-1,0)) (6, (0,1,1)) (20, (0,-1,0)) (13, (-1,1,-1)); (13, (-1,-1,-1)) (22, (-1,1,0)) (9, (1,0,-1)) (13, (-1,-1,0)) (20, (0,0,-1)) (9, (1,0,0)) (22, (-1,1,-1)) (17, (0,-1,-1)) (16, (-1,0,0)) (21, (-1,0,-1)) (10, (1,1,0)) (16, (-1,0,-1)) (17, (0,-1,0)) (10, (1,1,-1)) (21, (-1,0,0)) (3, (1,1,1)) (2, (0,0,-1)) (6, (1,-1,1)) (7, (0,0,0)) (2, (1,0,-1)) (3, (0,1,1)) (7, (1,0,0)) (6, (0,-1,1)) (4, (0,1,-1)) (1, (1,0,1)) (5, (0,-1,0)) (8, (1,0,1)) (1, (0,0,1)) (4, (1,1,-1)) (8, (0,0,1)) (5, (1,-1,0)) (12, (1,1,1)) (23, (-1,-1,0)) (19, (0,1,-1)) (14, (-1,-1,1)) (23, (-1,1,0)) (12, (1,-1,1)) (14, (-1,1,1)) (19, (0,-1,-1)) (24, (-1,-1,1)) (11, (1,1,0)) (15, (-1,-1,-1)) (18, (0,1,1)) (11, (1,-1,0)) (24, (-1,1,1)) (18, (0,-1,1)) (15, (-1,1,-1)); (10, (-1,0,-1)) (17, (0,1,0)) (16, (0,-1,-1)) (10, (-1,0,0)) (21, (0,0,-1)) (16, (0,-1,0)) (17, (0,1,-1)) (22, (0,0,-1)) (9, (-1,0,0)) (20, (0,-1,-1)) (13, (0,1,0)) (9, (-1,0,-1)) (22, (0,0,0)) (13, (0,1,-1)) (20, (0,-1,0)) (5, (1,1,0)) (8, (0,-1,0)) (4, (1,0,1)) (1, (0,0,0)) (8, (1,-1,0)) (5, (0,1,0)) (1, (1,0,0)) (4, (0,0,1)) (6, (0,1,0)) (7, (1,-1,0)) (3, (0,0,0)) (2, (1,0,1)) (7, (0,-1,0)) (6, (1,1,0)) (2, (0,0,1)) (3, (1,0,0)) (15, (0,1,1)) (18, (0,-1,0)) (24, (0,1,0)) (11, (-1,-1,0)) (18, (0,1,0)) (15, (0,-1,1)) (11, (-1,1,0)) (24, (0,-1,0)) (19, (0,-1,1)) (14, (0,1,0)) (12, (-1,-1,0)) (23, (0,1,0)) (14, (0,-1,0)) (19, (0,1,1)) (23, (0,-1,0)) (12, (-1,1,0)); (9, (-1,0,-1)) (13, (0,1,0)) (20, (0,-1,-1)) (9, (-1,0,0)) (22, (0,0,-1)) (20, (0,-1,0)) (13, (0,1,-1)) (21, (0,0,-1)) (10, (-1,0,0)) (16, (0,-1,-1)) (17, (0,1,0)) (10, (-1,0,-1)) (21, (0,0,0)) (17, (0,1,-1)) (16, (0,-1,0)) (6, (1,1,0)) (7, (0,-1,0)) (2, (1,0,1)) (3, (0,0,0)) (7, (1,-1,0)) (6, (0,1,0)) (3, (1,0,0)) (2, (0,0,1)) (5, (0,1,0)) (8, (1,-1,0)) (1, (0,0,0)) (4, (1,0,1)) (8, (0,-1,0)) (5, (1,1,0)) (4, (0,0,1)) (1, (1,0,0)) (19, (0,1,1)) (14, (0,-1,0)) (23, (0,1,0)) (12, (-1,-1,0)) (14, (0,1,0)) (19, (0,-1,1)) (12, (-1,1,0)) (23, (0,-1,0)) (15, (0,-1,1)) (18, (0,1,0)) (11, (-1,-1,0)) (24, (0,1,0)) (18, (0,-1,0)) (15, (0,1,1)) (24, (0,-1,0)) (11, (-1,1,0)); (12, (-1,-1,0)) (5, (0,1,0)) (8, (1,-1,0)) (11, (-1,-1,0)) (24, (0,0,0)) (7, (1,-1,0)) (6, (0,1,0)) (24, (0,-1,0)) (11, (-1,0,0)) (8, (0,-1,0)) (5, (1,1,0)) (12, (-1,0,0)) (23, (0,-1,0)) (6, (1,1,0)) (7, (0,-1,0)) (13, (0,1,0)) (20, (0,-1,-1)) (15, (0,-1,1)) (18, (0,0,0)) (16, (0,-1,-1)) (17, (0,1,0)) (14, (0,0,0)) (19, (0,-1,1)) (17, (0,1,-1)) (16, (0,-1,0)) (18, (0,-1,0)) (15, (0,0,1)) (20, (0,-1,0)) (13, (0,1,-1)) (19, (0,0,1)) (14, (0,-1,0)) (2, (1,1,1)) (3, (0,0,0)) (21, (0,1,-1)) (10, (-1,0,0)) (1, (0,1,0)) (4, (1,0,1)) (9, (-1,1,0)) (22, (0,0,-1)) (4, (0,0,1)) (1, (1,1,0)) (10, (-1,0,-1)) (21, (0,1,0)) (3, (1,0,0)) (2, (0,1,1)) (22, (0,0,0)) (9, (-1,1,-1)); (11, (-1,-1,0)) (6, (0,1,0)) (7, (1,-1,0)) (12, (-1,-1,0)) (23, (0,0,0)) (8, (1,-1,0)) (5, (0,1,0)) (23, (0,-1,0)) (12, (-1,0,0)) (7, (0,-1,0)) (6, (1,1,0)) (11, (-1,0,0)) (24, (0,-1,0)) (5, (1,1,0)) (8, (0,-1,0)) (17, (0,1,0)) (16, (0,-1,-1)) (19, (0,-1,1)) (14, (0,0,0)) (20, (0,-1,-1)) (13, (0,1,0)) (18, (0,0,0)) (15, (0,-1,1)) (13, (0,1,-1)) (20, (0,-1,0)) (14, (0,-1,0)) (19, (0,0,1)) (16, (0,-1,0)) (17, (0,1,-1)) (15, (0,0,1)) (18, (0,-1,0)) (4, (1,1,1)) (1, (0,0,0)) (22, (0,1,-1)) (9, (-1,0,0)) (3, (0,1,0)) (2, (1,0,1)) (10, (-1,1,0)) (21, (0,0,-1)) (2, (0,0,1)) (3, (1,1,0)) (9, (-1,0,-1)) (22, (0,1,0)) (1, (1,0,0)) (4, (0,1,1)) (21, (0,0,0)) (10, (-1,1,-1))],
    SMatrix{3, 3, Int, 9}[[-1 0 0; 0 -1 0; 0 0 -1], [0 -1 0; 1 0 0; 0 0 1], [0 1 0; -1 0 0; 0 0 -1], [-1 0 0; 0 -1 0; 0 0 1], [1 0 0; 0 1 0; 0 0 -1], [0 1 0; -1 0 0; 0 0 1], [0 -1 0; 1 0 0; 0 0 -1], [1 0 0; 0 -1 0; 0 0 -1], [-1 0 0; 0 1 0; 0 0 1], [0 -1 0; -1 0 0; 0 0 -1], [0 1 0; 1 0 0; 0 0 1], [-1 0 0; 0 1 0; 0 0 -1], [1 0 0; 0 -1 0; 0 0 1], [0 1 0; 1 0 0; 0 0 -1], [0 -1 0; -1 0 0; 0 0 1], [0 0 1; 1 0 0; 0 1 0], [0 0 -1; -1 0 0; 0 -1 0], [0 0 1; 0 -1 0; 1 0 0], [0 0 -1; 0 1 0; -1 0 0], [0 0 1; -1 0 0; 0 -1 0], [0 0 -1; 1 0 0; 0 1 0], [0 0 1; 0 1 0; -1 0 0], [0 0 -1; 0 -1 0; 1 0 0], [0 0 -1; 1 0 0; 0 -1 0], [0 0 1; -1 0 0; 0 1 0], [0 0 -1; 0 -1 0; -1 0 0], [0 0 1; 0 1 0; 1 0 0], [0 0 -1; -1 0 0; 0 1 0], [0 0 1; 1 0 0; 0 -1 0], [0 0 -1; 0 1 0; 1 0 0], [0 0 1; 0 -1 0; -1 0 0], [0 1 0; 0 0 1; 1 0 0], [0 -1 0; 0 0 -1; -1 0 0], [1 0 0; 0 0 1; 0 -1 0], [-1 0 0; 0 0 -1; 0 1 0], [0 -1 0; 0 0 1; -1 0 0], [0 1 0; 0 0 -1; 1 0 0], [-1 0 0; 0 0 1; 0 1 0], [1 0 0; 0 0 -1; 0 -1 0], [0 -1 0; 0 0 -1; 1 0 0], [0 1 0; 0 0 1; -1 0 0], [-1 0 0; 0 0 -1; 0 -1 0], [1 0 0; 0 0 1; 0 1 0], [0 1 0; 0 0 -1; -1 0 0], [0 -1 0; 0 0 1; 1 0 0], [1 0 0; 0 0 -1; 0 1 0], [-1 0 0; 0 0 1; 0 -1 0]],
    SVector{3, Rational{Int32}}[[0//1, 1//4, 3//4], [5//8, 5//8, 0//1], [3//8, 5//8, 3//4], [0//1, 1//4, 0//1], [0//1, 0//1, 3//4], [3//8, 5//8, 0//1], [5//8, 5//8, 3//4], [0//1, 1//4, 3//4], [0//1, 0//1, 0//1], [5//8, 5//8, 3//4], [3//8, 5//8, 0//1], [0//1, 0//1, 3//4], [0//1, 1//4, 0//1], [3//8, 5//8, 3//4], [5//8, 5//8, 0//1], [1//8, 5//8, 1//4], [7//8, 5//8, 1//2], [1//8, 1//4, 7//8], [7//8, 0//1, 7//8], [1//8, 5//8, 1//2], [7//8, 5//8, 1//4], [1//8, 0//1, 7//8], [7//8, 1//4, 7//8], [7//8, 5//8, 1//2], [1//8, 5//8, 1//4], [7//8, 1//4, 7//8], [1//8, 0//1, 7//8], [7//8, 5//8, 1//4], [1//8, 5//8, 1//2], [7//8, 0//1, 7//8], [1//8, 1//4, 7//8], [3//8, 3//4, 7//8], [5//8, 1//2, 7//8], [0//1, 3//4, 1//2], [0//1, 1//2, 1//4], [5//8, 3//4, 7//8], [3//8, 1//2, 7//8], [0//1, 3//4, 1//4], [0//1, 1//2, 1//2], [5//8, 1//2, 7//8], [3//8, 3//4, 7//8], [0//1, 1//2, 1//2], [0//1, 3//4, 1//4], [3//8, 1//2, 7//8], [5//8, 3//4, 7//8], [0//1, 1//2, 1//4], [0//1, 3//4, 1//2]],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1],
    true
)

# utility to compare RingIncluding
function canonicalize_ri(ri::PeriodicGraphs.RingIncluding{D}) where D
    ret = [PeriodicVertex{D}[] for _ in 1:length(ri)]
    for i in 1:length(ri)
        cycle = collect(ri[i])
        fst = argmin(cycle)
        lenc = length(cycle)
        endreverse = cycle[mod1(fst-1, lenc)] > cycle[mod1(fst+1, lenc)]
        if fst != 1 || !endreverse
            reverse!(cycle, 1, fst - endreverse)
            reverse!(cycle, fst - endreverse + 1, lenc)
            endreverse && reverse!(cycle)
        end
        ret[i] = cycle
    end
    return sort!(ret)
end

@testset "Symmetry and rings" begin
    ras_sny2 = RingAttributions(sny, true, 2)
    rasym_sny2 = RingAttributions(sny, 2, symmetries_sny) # test orig_1 with symmetries
    @test length(ras_sny2) == length(rasym_sny2)
    for (r1, r2) in zip(ras_sny2, rasym_sny2)
        @test canonicalize_ri(r1) == canonicalize_ri(r2)
    end

    # AbstractArray interface
    _ri = first(ras_sny2)
    @test _ri == ras_sny2[1] == ras_sny2[begin]
    @test last(ras_sny2) == ras_sny2[length(ras_sny2)] == ras_sny2[end]

    @test _ri[1] == first(_ri) == _ri[begin]
    @test _ri[length(_ri)] == last(_ri) == _ri[end]

    ras_sny6 = RingAttributions(sny, true, 6)
    rasym_sny6 = RingAttributions(sny, symmetries_sny)
    @test length(ras_sny6) == length(rasym_sny6)
    for (r1, r2) in zip(ras_sny6, rasym_sny6)
        @test canonicalize_ri(r1) == canonicalize_ri(r2)
    end

    ras_afi = RingAttributions(afi, true, 10)
    rasym_afi = RingAttributions(afi, true, symmetries_afi)
    @test length(ras_afi) == length(rasym_afi)
    for (r1, r2) in zip(ras_afi, rasym_afi)
        @test canonicalize_ri(r1) == canonicalize_ri(r2)
    end

    ras_lta = RingAttributions(lta, true, 10)
    rasym_lta = RingAttributions(lta, true, 10, symmetries_lta)
    @test length(ras_lta) == length(rasym_lta)
    for (r1, r2) in zip(ras_lta, rasym_lta)
        @test canonicalize_ri(r1) == canonicalize_ri(r2)
    end

    withid = PeriodicGraphs.IncludingIdentity(symmetries_lta)
    @test eltype(withid) == PeriodicSymmetries.PeriodicSymmetry3D{Rational{Int32}}
    @test length(withid) == 48
    id, id2 = withid
    @test isequal(id, withid[1]) && isequal(id, one(symmetries_lta))
    @test isequal(id2, withid[2]) && isequal(id2, first(symmetries_lta))

    trivsymmgroup = PeriodicGraphs.IncludingIdentity(NoSymmetryGroup(lta))
    @test unique(trivsymmgroup) == 1:nv(lta)
    @test trivsymmgroup(13) == 13
    @test length(trivsymmgroup) == 1
    @test PeriodicGraphs.IncludingIdentity(trivsymmgroup) === trivsymmgroup
    @test_throws ErrorException one(trivsymmgroup)
    triv = only(trivsymmgroup)
    @test eltype(trivsymmgroup) == IdentitySymmetry
    @test triv[5] == 5
    @test triv[PeriodicVertex3D(14)] == PeriodicVertex3D(14)
    @test triv([1, 2, 4]) == [1, 2, 4]

    weaks, symmweaks = rings(lta, symmetries_lta)
    @test length(symmweaks) == length(symmetries_lta)
    for str in weaks
        id = symmweaks(str)
        for symm in symmweaks
            image = symm(str)
            @test image == first(PeriodicGraphs.normalize_cycle!(copy(image), nv(lta), Val(3)))
            imagevertices = [reverse_hash_position(x, lta) for x in image]
            @test image == hash_position.(first(PeriodicGraphs.normalize_cycle!(imagevertices)), nv(lta))
            @test image ∈ weaks
            @test symmweaks(image) == id
        end
    end

    strs, symmstrs = strong_rings(lta, symmetries_lta)
    @test length(symmstrs) == length(symmetries_lta)
    for str in strs
        id = symmstrs(str)
        for symm in symmstrs
            image = symm(str)
            @test image == first(PeriodicGraphs.normalize_cycle!(copy(image), nv(lta), Val(3)))
            imagevertices = [reverse_hash_position(x, lta) for x in image]
            @test image == hash_position.(first(PeriodicGraphs.normalize_cycle!(imagevertices)), nv(lta))
            @test image ∈ strs
            @test symmstrs(image) == id
        end
    end

    _strs, _symmstrs, estrs, kp = strong_erings(lta, symmetries_lta)
    @test _strs == strs
    @test unique(symmstrs) == unique(_symmstrs)
    @test length.(estrs) == length.(strs)
    for (str, estr) in Iterators.take(zip(strs, estrs), 5)
        @test issorted(estr)
        for e in estr
            _e = kp[e]
            @test get!(kp, _e) == kp[_e] == e
        end
        lenstr = length(str)
        _edgestr = [minmax(str[i], str[mod1(i+1,lenstr)]) for i in 1:lenstr]
        @test issetequal(_edgestr, minmax(hash_position.(kp[e], nv(lta))...) for e in estr)
    end
end

@testset "Simple symmetries" begin
    map = [2,1,4,3,5]
    symm = SimpleSymmetry(map)
    @test symm[3] == 4
    @test symm[5] == 5
    dict = Dict(:a => 1, :b => 2, :c => 3, :d => 4, :e => 5, :f => 1, :g => 2, :h => 3, :i => 4, :j => 5)
    ssymm = SimpleSymmetry(map, dict)
    @test ssymm[:c] == ssymm[:h] == ssymm[3] == 4
    @test ssymm[:e] == ssymm[:j] == ssymm[5] == 5
end
