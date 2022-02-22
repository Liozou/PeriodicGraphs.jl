using Test
using Graphs
using PeriodicGraphs
import PeriodicGraphs: ofs, vertex_permutation, dimensionality, LoopException, extended_gcd
import StaticArrays: SVector

# using Aqua
# Aqua.test_all(PeriodicGraphs)


@testset "PeriodicVertex" begin
    @test_throws MethodError PeriodicVertex2D() # Always require a vertex identifier
    @test iszero(ofs(PeriodicVertex2D(1)))
    @test ndims(PeriodicVertex2D(1)) == 2
    @test PeriodicVertex2D(1) != PeriodicVertex3D(1)
    @test PeriodicVertex{4}(5) == PeriodicVertex{4}(5)
    @test PeriodicVertex1D(1, SVector{1,Int}(3)) == PeriodicVertex1D(1, [3])
    @test PeriodicVertex(3, (1,2)) == PeriodicVertex2D(3, (1,2)) == PeriodicVertex(3, SVector{2,Int}(1,2))
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
    @static if VERSION < v"1.6.0-DEV"
        @test string(PeriodicVertex3D(12, (1,0,-1))) == "PeriodicVertex{3}(12, (1,0,-1))"
        @test string(PeriodicVertex2D[(1, (0,0)), (1, (-1,0))]) == "PeriodicVertex{2}[(1, (0,0)), (1, (-1,0))]"
        @test string(PeriodicVertex1D(2, (1,))) == "PeriodicVertex{1}(2, (1,))"
    else
        @test string(PeriodicVertex3D(12, (1,0,-1))) == "PeriodicVertex3D(12, (1,0,-1))"
        @test string(PeriodicVertex2D[(1, (0,0)), (1, (-1,0))]) == "PeriodicVertex2D[(1, (0,0)), (1, (-1,0))]"
        @test string(PeriodicVertex1D(2, (1,))) == "PeriodicVertex1D(2, (1,))"
    end
    @test ofs(PeriodicVertex(1, (2,1))) == SVector{2,Int}(2,1)
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
    @test_throws MethodError PeriodicEdge1D(1, 2)
    @test PeriodicEdge2D(1, 2, (0,0)) != PeriodicEdge3D(1, 2, (0,0,0))
    @test PeriodicEdge(2, PeriodicVertex2D(1, (1,3))) == PeriodicEdge2D(2, 1, (1,3))
    @test PeriodicEdge(1, 1, SVector{1,Int}(2)) == PeriodicEdge1D(1, 1, SVector{1,Int}(2))
    @test PeriodicEdge(2, 1, ()) == PeriodicEdge{0}((2, 1, SVector{0,Int}()))
    @test PeriodicEdge((1, 2, (1,0,0))) == PeriodicEdge3D((1, 2, (1,0,0)))
    @test hash(PeriodicEdge((1, 2, (1,0,0)))) == hash(PeriodicEdge3D((1, 2, (1,0,0))))
    @test hash(PeriodicEdge((1, 2, (1,)))) != hash(PeriodicEdge((2, 1, (-1,))))
    @test convert(PeriodicEdge, (3, PeriodicVertex{1}(2, (1,)))) == convert(PeriodicEdge1D, (3, PeriodicVertex{1}(2, (1,))))
    @test PeriodicEdge1D[(1,2,SVector{1,Int}(3))] == PeriodicEdge[(1,2,(3,))] == [PeriodicEdge(1, 2, (3,))]
    @static if VERSION < v"1.6.0-DEV"
        @test string(PeriodicEdge3D[(1, 2, (1,0,0))]) == "PeriodicEdge{3}[(1, 2, (1,0,0))]"
        @test string(PeriodicEdge(3,4, (0,0))) == "PeriodicEdge{2}(3, 4, (0,0))"
        @test string(PeriodicEdge{1}(1, 1, (2,))) == "PeriodicEdge{1}(1, 1, (2,))"
    else
        @test string(PeriodicEdge3D[(1, 2, (1,0,0))]) == "PeriodicEdge3D[(1, 2, (1,0,0))]"
        @test string(PeriodicEdge(3,4, (0,0))) == "PeriodicEdge2D(3, 4, (0,0))"
        @test string(PeriodicEdge{1}(1, 1, (2,))) == "PeriodicEdge1D(1, 1, (2,))"
    end
    @test reverse(reverse(PeriodicEdge(1, 1, (2,3,1)))) == PeriodicEdge3D(1, 1, (2,3,1))
    @test reverse(PeriodicEdge(1, 2, (1,-2))) == PeriodicEdge(2, 1, (-1,2))
    @test reverse(PeriodicEdge(3, 3, (1,0))) == PeriodicEdge(3, 3, (-1,0))
    @test src(PeriodicEdge(1, 2, (-1, 0))) == 1
    @test dst(PeriodicEdge(2, 2, (1,))) == 2
    @test ofs(PeriodicEdge(1, 1, (0,1,3))) == SVector{3,Int}(0,1,3)
    @test PeriodicEdge2D(1, 1, (2,3)) < PeriodicEdge2D(1, 1, (3,0)) < PeriodicEdge2D(1, 2, (3,0)) < PeriodicEdge2D(2, 2, (3,0))
    @test ndims(PeriodicEdge(1, 2, ())) == 0ndims(PeriodicEdge(1, 2, ())) == 0
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
            PeriodicGraphs.isindirectedge(PeriodicEdge{N}(i, l[end])) || return false
            continue
        end
        PeriodicGraphs.isindirectedge(PeriodicEdge{N}(i, l[dedgestart])) && return false
        if dedgestart != 1 && !PeriodicGraphs.isindirectedge(PeriodicEdge{N}(i, l[dedgestart-1]))
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
    @test PeriodicGraph("0 1 2 1 3") == PeriodicGraph(PeriodicEdge{0}[(1, 2, ()), (1, 3, ())])
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
    g = PeriodicGraph("2 1 2 0 0 2 3 1 0")
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
    @test_throws ArgumentError g[[3,3]]
    @test_throws ArgumentError offset_representatives!(gg, [2])
    gcopy = deepcopy(g)
    @test offset_representatives!(g, [0,0,0,0]) == gcopy == g
    @test is_well_formed(g)
    ggg = offset_representatives!(gg, [(1,-1,1),0])
    @test ggg == gg
    @test collect(edges(gg)) == PeriodicEdge3D[(1, 2, (1,-1,0)), (2, 2, (1,0,0))]
    gcopy = deepcopy(gg)
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
end

@testset "Graph reduction" begin
    g = PeriodicGraph3D("3  1 2  0 0 0  1 3  0 1 0   2 3  1 0 -1
                            4 4  0 0 1  4 5  0 -1 0  5 6  0 0 0")
    c = cellgraph(g)
    p = periodiccellgraph(g)
    @test nv(g) == nv(p) == nv(c)
    @test issubset(edges(c), edges(p))
    for e in edges(p)
        @test has_edge(g, src(e), dst(e))
    end
    @test connected_components(p) == connected_components(g)
end

@testset "Periodic vertex hash" begin
    ## dim = 0
    @test PeriodicGraphs.hash_position(SVector{0,Int}()) == 0
    ## dim = 1
    @test PeriodicGraphs.hash_position(SVector{1,Int}(0)) == 0
    @test PeriodicGraphs.hash_position(SVector{1,Int}(-1)) == 1
    @test PeriodicGraphs.hash_position(SVector{1,Int}(1)) == 2
    @test PeriodicGraphs.hash_position(SVector{1,Int}(-2)) == 3
    @test PeriodicGraphs.hash_position(SVector{1,Int}(2)) == 4
    @test PeriodicGraphs.hash_position(SVector{1,Int}(3)) == 6
    ## dim = 2
    n::Int = 4
    seen = falses((2n-1)^2)
    for i in -n+1:n-1, j in -n+1:n-1
        x = PeriodicGraphs.hash_position(SVector{2,Int}(i,j))+1
        @test !seen[x]
        seen[x] = true
    end
    @test all(seen)
    append!(seen, falses((2n+1)^2-(2n-1)^2))
    for i in (-n,n), j in -n:n
        x = PeriodicGraphs.hash_position(SVector{2,Int}(i,j))+1
        @test !seen[x]
        seen[x] = true
    end
    for i in -n+1:n-1, j in (-n,n)
        x = PeriodicGraphs.hash_position(SVector{2,Int}(i,j))+1
        @test !seen[x]
        seen[x] = true
    end
    @test all(seen)
    ## dim = 3
    seen = falses((2n-1)^3)
    for i in -n+1:n-1, j in -n+1:n-1, k in -n+1:n-1
        x = PeriodicGraphs.hash_position(SVector{3,Int}(i,j,k))+1
        @test !seen[x]
        seen[x] = true
    end
    @test all(seen)
    append!(seen, falses((2n+1)^3-(2n-1)^3))
    for i in (-n,n), j in -n:n, k in -n:n
        x = PeriodicGraphs.hash_position(SVector{3,Int}(i,j,k))+1
        @test !seen[x]
        seen[x] = true
    end
    for i in -n+1:n-1, j in (-n,n), k in -n:n
        x = PeriodicGraphs.hash_position(SVector{3,Int}(i,j,k))+1
        @test !seen[x]
        seen[x] = true
    end
    for i in -n+1:n-1, j in -n+1:n-1, k in (-n,n)
        x = PeriodicGraphs.hash_position(SVector{3,Int}(i,j,k))+1
        @test !seen[x]
        seen[x] = true
    end

    @test all(seen)
    g = PeriodicGraph3D([PeriodicEdge3D(1, 1, (0, 0, 1)),
                         PeriodicEdge3D(1, 1, (1, 1, 0))])
    @test neighborhood(g, 1, 1) == PeriodicVertex3D[(1, (0, 0, 0)),
            (1, (-1, -1, 0)), (1, (0, 0, -1)), (1, (0, 0, 1)), (1, (1, 1, 0))]
    @test neighborhood(g, 1, 2) == PeriodicVertex3D[(1, (0, 0, 0)), (1, (-1, -1, 0)),
            (1, (0, 0, -1)), (1, (0, 0, 1)), (1, (1, 1, 0)), (1, (-2, -2, 0)),
            (1, (-1, -1, -1)), (1, (-1, -1, 1)), (1, (0, 0, -2)), (1, (1, 1, -1)),
            (1, (0, 0, 2)), (1, (1, 1, 1)), (1, (2, 2, 0))]
end

@testset "Edge iteration" begin
    g::PeriodicGraph3D = PeriodicGraph3D([PeriodicEdge3D(2, 3, (-1, 0, 0)),
                                          PeriodicEdge3D(2, 3, (0, 1, 0)),
                                          PeriodicEdge3D(4, 4, (-1, 1, 0)),
                                          PeriodicEdge3D(1, 3, (0, 0, 1)),
                                          PeriodicEdge3D(3, 2, (0, 0, 0))])
    expected = PeriodicEdge3D[(1, 3, (0, 0, 1)), (2, 3, (-1, 0, 0)),
                              (2, 3, (0, 0, 0)), (2, 3, (0, 1, 0)),
                              (4, 4, (1, -1, 0))]
    @test collect(edges(g)) == expected
    @test PeriodicEdge3D(4, 4, (1, -1, 0)) ∈ edges(g)
    add_vertex!(g)
    @test collect(edges(g)) == expected
    @test length(edges(g)) == ne(g)
end

@testset "Vertex removal" begin
    g::PeriodicGraph3D = PeriodicGraph3D([PeriodicEdge3D(1, 3, (0, 0, 1)),
                                          PeriodicEdge3D(2, 3, (0, 1, 0)),
                                          PeriodicEdge3D(2, 3, (0, -1, 0)),
                                          PeriodicEdge3D(3, 2, (0, 0, 0)),
                                          PeriodicEdge3D(4, 4, (1, 1, 0))])
    gg = deepcopy(g)
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
    g2 = deepcopy(g1)
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
    gg = change_dimension(PeriodicGraph1D, g)
    @test string(gg) == "1 1 2 0 1 4 -1 2 3 0 3 4 0 4 4 1"
    @test change_dimension(PeriodicGraph2D, gg) == change_dimension(PeriodicGraph2D, g) ==
        PeriodicGraph("2 1 2 0 0 1 4 -1 0 2 3 0 0 3 4 0 0 4 4 1 0")
    @test_throws DimensionMismatch change_dimension(PeriodicGraph{0}, g)

    s = [14,28,-17,-34,12]
    d, λ = extended_gcd(s)
    @test d == 1
    @test reduce(+, s .* λ) == d
end
