module PeriodicGraphs

using LinearAlgebra
using StaticArrays
using Graphs

export PeriodicVertex, PeriodicEdge, PeriodicGraph,
       PeriodicVertex1D, PeriodicEdge1D, PeriodicGraph1D,
       PeriodicVertex2D, PeriodicEdge2D, PeriodicGraph2D,
       PeriodicVertex3D, PeriodicEdge3D, PeriodicGraph3D,
       coordination_sequence, cellgraph, periodiccellgraph,
       offset_representatives!, swap_axes!, find_edges,
       dimensionality, change_dimension

import Base: (==), isless, convert, show, showerror, eltype, iterate, zero,
             length, in, ndims, print, cmp
import Base.Order: Forward, Lt

"""
    PeriodicVertex{N}

Vertex type for an `N`-periodic graph.

A vertex is uniquely determined by the identifier of its representative in the
a fixed initial cell, and the offset of the cell containing the vertex compared
to the the initial cell.
Vertex identifiers start at 1.
"""
struct PeriodicVertex{N}
    v::Int
    ofs::SVector{N,Int}
end
PeriodicVertex{N}(n::Integer) where {N} = PeriodicVertex{N}(n, zero(SVector{N,Int}))
PeriodicVertex(v::Integer, ofs::Union{SVector{N},NTuple{N}}) where {N} = PeriodicVertex{N}(v, ofs)

function show(io::IO, x::PeriodicVertex{N}) where N
    if get(io, :typeinfo, Any) != PeriodicVertex{N}
        print(io, PeriodicVertex{N})
    end
    print(io, '(', x.v, ", (", join(x.ofs, ','))
    N == 1 && print(io, ',')
    print(io, ')', ')')
end
function convert(::Type{PeriodicVertex{N}}, (dst, offset)::Tuple{Integer,Any}) where N
    PeriodicVertex{N}(dst, offset)
end
function convert(::Type{PeriodicVertex}, (dst, offset)::Tuple{Integer,Union{SVector{N},NTuple{N}}}) where N
    PeriodicVertex{N}(dst, offset)
end
function cmp(x::PeriodicVertex{N}, y::PeriodicVertex{N}) where N
    c = cmp(x.v, y.v)
    iszero(c) || return c
    return cmp(x.ofs, y.ofs)
end
isless(x::PeriodicVertex{N}, y::PeriodicVertex{N}) where {N} = cmp(x, y) < 0
==(x::PeriodicVertex{N}, y::PeriodicVertex{M}) where {N,M} = N == M && iszero(cmp(x, y))

const PeriodicVertex1D = PeriodicVertex{1}
const PeriodicVertex2D = PeriodicVertex{2}
const PeriodicVertex3D = PeriodicVertex{3}


ZtoN(x::Signed) = -(x<0) + 2*abs(x)
function hash_position((x1,x2,x3)::SVector{3,<:Integer})
    x1 = ZtoN(x1); x2 = ZtoN(x2); x3 = ZtoN(x3)
    _b = x1 >= x2
    b1 = _b & (x1 >= x3)
    b2 = (!_b) & (x2 >= x3)
    b3 = (!b1) & (!b2)
    return b1*(x2*(x1 + 1) + x3 + x1^3) +
           b2*((x2 + 1)*(x1 + x2 + 1) + x3 + x2^3) +
           b3*((x3 + 1)*(2x3 + 1) + x3*(x1 + x3^2) + x2)
end

function hash_position((x1,x2)::SVector{2,<:Integer})
    x1 = ZtoN(x1); x2 = ZtoN(x2)
    b = x1 >= x2
    return b*(x2 + x1^2) + (!b)*(x1 + x2*(x2 + 1) + 1)
end

function hash_position((x,)::SVector{1,<:Integer})
    return ZtoN(x)
end

function hash_position(::SVector{0,<:Integer})
    return 0
end

"""
    hash_position(x::PeriodicVertex{N}, n::Integer) where {N}

Given a `PeriodicVertex{N}` and the number n of vertex identifiers in a graph, compute a
unique hash for the given vertex.

This hash function is a bijection between the set of all the vertices of the graph and
the set of positive integers. Its value is an integer between `1+n\\*(2d-1)^N`
(or `1` if `d == 0`) and `n*(2d+1)^N`, where `d = max.(abs.(x.ofs))`.
This means that when one unit cell A is further than another B (for the Manhattan
distance), all vertices in A will have a larger hash than all vertices in B.
"""
function hash_position(x::PeriodicVertex, n::Integer)
    return x.v + n*hash_position(x.ofs)
end

"""
    LoopException <: Exception

Error type for constructing an invalid `PeriodicEdge{N}` of the form (u, u, zeros(Int,N)).
Loops are not expected in the algorithms implemented in PeriodicGraphs.jl. If you
still want to construct them, use the `unsafe_edge{N}` constructor instead of
`PeriodicEdge{N}``.
"""
struct LoopException <: Exception
    src::Int
end
showerror(io::IO, e::LoopException) = print(io, "LoopException: a loop from vertex $(e.src) to itself in the same unit cell is a forbidden edges. Maybe the offset is wrong?")
@noinline __throw_loopexception(src) = throw(LoopException(src))

"""
    unsafe_edge{N}

Internal constructor for `PeriodicEdge{N}` that bypasses the loop check.
This function is not part of the official API, use with caution.
"""
struct unsafe_edge{N} end

"""
    PeriodicEdge{N} <: Graphs.SimpleGraphs.AbstractSimpleEdge{Int}

Edge type for an `N`-periodic graph.

An edge is uniquely determined by the vertex identifiers of its source and
destination, and the cell offset between the source vertex and the destination vertex.
"""
struct PeriodicEdge{N} <: Graphs.SimpleGraphs.AbstractSimpleEdge{Int}
    src::Int
    dst::PeriodicVertex{N}
    function (::Type{unsafe_edge{N}})(src, dst) where N
        return new{N}(src, dst)
    end
end
unsafe_edge{N}(src, dst, ofs) where {N} = unsafe_edge{N}(src, PeriodicVertex{N}(dst, ofs))
function PeriodicEdge{N}(src, dst::PeriodicVertex{N}) where N
    src == dst.v && all(iszero.(dst.ofs)) && __throw_loopexception(src)
    return unsafe_edge{N}(src, dst)
end
PeriodicEdge(src, dst::PeriodicVertex{N}) where {N} = PeriodicEdge{N}(src, dst)
function PeriodicEdge{N}(src, dst, offset) where N
    return PeriodicEdge{N}(src, PeriodicVertex{N}(dst, offset))
end
function PeriodicEdge(src, dst, offset::Union{SVector{N,T},NTuple{N,T}}) where {N,T}
    PeriodicEdge{N}(src, dst, offset)
end
function PeriodicEdge{N}((src, dst, offset)::Tuple{Any,Any,Union{SVector{N,T},NTuple{N,T}}}) where {N,T<:Integer}
    PeriodicEdge{N}(src, dst, offset)
end
function PeriodicEdge((src, dst, offset)::Tuple{Any,Any,Union{SVector{N,T},NTuple{N,T}}}) where {N,T<:Integer}
    PeriodicEdge{N}(src, dst, offset)
end

function convert(::Type{PeriodicEdge{N}}, (src, dst, offset)::Tuple{Any,Any,Any}) where {N}
    PeriodicEdge{N}(src, dst, offset)
end
function convert(::Type{PeriodicEdge}, (src, dst, offset)::Tuple{Any,Any,Union{SVector{N,T},NTuple{N,T}}}) where {N,T<:Integer}
    PeriodicEdge{N}(src, dst, offset)
end
function convert(::Union{Type{PeriodicEdge},Type{PeriodicEdge{N}}}, (src, v)::Tuple{Any,PeriodicVertex{N}}) where N
    PeriodicEdge{N}(src, v)
end

function show(io::IO, x::PeriodicEdge{N}) where N
    if get(io, :typeinfo, Any) != PeriodicEdge{N}
        print(io, PeriodicEdge{N})
    end
    print(io, '(', x.src, ", ", x.dst.v, ", (", join(x.dst.ofs, ','))
    N == 1 && print(io, ',')
    print(io, ')', ')')
end

const PeriodicEdge1D = PeriodicEdge{1}
const PeriodicEdge2D = PeriodicEdge{2}
const PeriodicEdge3D = PeriodicEdge{3}

function Graphs.reverse(e::PeriodicEdge{N}) where N
    return unsafe_edge{N}(e.dst.v, e.src, .-e.dst.ofs)
end
Graphs.src(e::PeriodicEdge) = e.src
Graphs.dst(e::PeriodicEdge) = e.dst.v

"""
    ofs(v::PeriodicVertex)
    ofs(e::PeriodicEdge)

The offset of a vertex or an edge.
"""
ofs(v::PeriodicVertex) = v.ofs
ofs(e::PeriodicEdge) = ofs(e.dst)
isindirectedge(e::PeriodicEdge) = e.src > e.dst.v || (e.src == e.dst.v && e.dst.ofs < zero(e.dst.ofs))
function cmp(x::PeriodicEdge{N}, y::PeriodicEdge{N}) where N
    c = cmp(x.src, y.src)
    iszero(c) || return c
    return cmp(x.dst, y.dst)
end
isless(x::PeriodicEdge{N}, y::PeriodicEdge{N}) where {N} = cmp(x, y) < 0
==(x::PeriodicEdge{N}, y::PeriodicEdge{M}) where {N,M} = N == M && iszero(cmp(x, y))


"""
    PeriodicGraph{N} <: AbstractGraph{Int}

Type representing an undirected `N`-periodic graph.
"""
struct PeriodicGraph{N} <: AbstractGraph{Int}
    ne::Base.RefValue{Int}
    nlist::Vector{Vector{PeriodicVertex{N}}} # For each vertex, the sorted list of its neighbors
    directedgestart::Vector{Int}
    width::Base.RefValue{Rational{Int}}
end

const PeriodicGraph1D = PeriodicGraph{1}
const PeriodicGraph2D = PeriodicGraph{2}
const PeriodicGraph3D = PeriodicGraph{3}

function PeriodicGraph{N}(nv::Integer, t, s) where N
    return PeriodicGraph{N}(Ref(nv), t, s, Ref(-1//1))
end

"""
    PeriodicGraph{N}(nv::Integer=0)

Construct a `PeriodicGraph{N}` with `nv` vertices and 0 edge.
"""
function PeriodicGraph{N}(n::Integer = 0) where N
    @assert n >= 0
    return PeriodicGraph{N}(0, [PeriodicVertex{N}[] for _ in 1:n], [1 for _ in 1:n])
end
"""
    PeriodicGraph([nv::Integer, ]edge_list::AbstractVector{PeriodicEdge{N}})
    PeriodicGraph{N}([nv::Integer, ]edge_list::AbstractVector{PeriodicEdge{N}})

Construct a `PeriodicGraph{N}` from a vector of edges.
If `nv` is unspecified, the number of vertices is the highest that is used
in an edge in `edge_list`.

### Implementation Notes
This constructor works the fastest when `edge_list` is sorted by the lexical ordering
and does not contain any duplicates.

## Examples
```jldoctest
julia> el = PeriodicEdge2D[(1, 1, (1,0)), (1, 3, (0,1)), (3, 1, (0,-1))];

julia> g = PeriodicGraph(el)
PeriodicGraph2D(3, PeriodicEdge2D[(1, 1, (1, 0)), (1, 3, (0, 1))])

julia> ne(g)
2
```
"""
function PeriodicGraph{N}(nv::Integer, t::AbstractVector{PeriodicEdge{N}}) where N
    sort!(t); unique!(t)
    ne = length(t)
    g = PeriodicGraph{N}(nv)
    for e in t
        add_edge!(g, e)
    end
    return g
end
PeriodicGraph(nv::Integer, t::AbstractVector{PeriodicEdge{N}}) where {N} = PeriodicGraph{N}(nv, t)
function PeriodicGraph{N}(t::AbstractVector{PeriodicEdge{N}}) where N
    @static if VERSION < v"1.6-"
        isempty(t) && return PeriodicGraph(0, t)
        return PeriodicGraph(maximum(max(e.src, e.dst.v) for e in t), t)
    else
        return PeriodicGraph(maximum(max(e.src, e.dst.v) for e in t; init=0), t)
    end
end
PeriodicGraph(t::AbstractVector{PeriodicEdge{N}}) where {N} = PeriodicGraph{N}(t)


struct KeyString{T,S<:AbstractString}
    x::S
    start::Base.RefValue{Int}
end
function KeyString{T}(x) where T
    KeyString{T,typeof(x)}(x, Ref(firstindex(x)))
end
iterate(::KeyString, ::Nothing) = nothing
function iterate(k::KeyString{T}, (next_char, idx)) where T
    start = idx
    state = nothing
    while isspace(next_char)
        state = iterate(k.x, idx)
        state isa Nothing && return nothing
        start = idx
        next_char, idx = state
    end
    stop = start
    tmp = idx
    while !isspace(next_char)
        state = iterate(k.x, idx)
        state isa Nothing && break
        stop = tmp
        tmp = idx
        next_char, idx = state
    end
    ret = tryparse(T, SubString(k.x, start:stop))
    if ret isa T
        return (ret, state)
    end
    throw(ArgumentError("Input string does not represent a graph"))
end
function iterate(k::KeyString)
    iterate(k, (' ', k.start[]))
end
function Base.popfirst!(k::KeyString)
    next = iterate(k)
    next isa Nothing && throw(ArgumentError("Input string does not represent a graph"))
    char, state = next
    k.start[] = state isa Nothing ? (ncodeunits(k.x) + 1) : last(state)
    return char
end
Base.isempty(k::KeyString) = isempty(SubString(k.x, k.start[]))
Base.IteratorSize(::Type{<:KeyString}) = Base.SizeUnknown()

"""
    PeriodicGraph(key::AbstractString)
    PeriodicGraph{N}(key::AbstractString)

Construct a `PeriodicGraph{N}` from a `key`, which is a string of whitespace-separated
values of the form
`"N src1 dst1 ofs1_1 ofs1_2 ... ofs1_N src2 dst2 ofs2_1 ofs2_2 ... ofs2_N  ...  srcm dstm ofsm_1 ofsm_2 ... ofsm_N"`
where `N` is the number of repeating dimensions of the graph, `m` is the number of edges
and for all `i` between `1` and `m`, the number of edges, the `i`-th edge is described as
- `srci`, the vertex identifier of the source vertex,
- `dsti`, the vertex identifier of the destination vertex and
- `(ofsi_1, ofsi_2, ..., ofsi_N)` the offset of the edge.

This compact representation of a graph can be obtained simply by `print`ing the graph
or with `string`.

## Examples
```jldoctest
julia> PeriodicGraph("2  1 2 0 0  2 1 1 0  1 1 0 -1")
PeriodicGraph2D(2, PeriodicEdge2D[(1, 1, (0,1)), (1, 2, (-1,0)), (1, 2, (0,0))])

julia> PeriodicGraph3D("3  1 1 0 0 1  1 1 0 1 0  1 1 1 0 0")
PeriodicGraph3D(1, PeriodicEdge3D[(1, 1, (0,0,1)), (1, 1, (0,1,0)), (1, 1, (1,0,0))])

julia> string(ans)
"3 1 1 0 0 1 1 1 0 1 0 1 1 1 0 0"

julia> string(PeriodicGraph(ans)) == ans
true
```
"""
function PeriodicGraph{N}(s::AbstractString) where N
    key = KeyString{Int}(s)
    M = popfirst!(key)
    M != N && throw(DimensionMismatch("Cannot construct a $N-dimensional graph from a $M-dimensional key"))
    edgs = PeriodicEdge{N}[]
    while !isempty(key)
        src = popfirst!(key)
        dst = popfirst!(key)
        ofs = SVector{N,Int}([popfirst!(key) for _ in 1:N])
        push!(edgs, PeriodicEdge{N}(src, dst, ofs))
    end
    return PeriodicGraph{N}(edgs)
end
function PeriodicGraph(s::AbstractString)
    N = first(KeyString{Int}(s))
    return PeriodicGraph{N}(s)
end
function PeriodicGraph{N}(g::PeriodicGraph{N}) where N
    PeriodicGraph{N}(g.ne[], [copy(x) for x in g.nlist], copy(g.directedgestart))
end
PeriodicGraph(g::PeriodicGraph{N}) where {N} = PeriodicGraph{N}(g)


function show(io::IO, g::PeriodicGraph{N}) where N
    if get(io, :typeinfo, Any) != PeriodicGraph{N}
        print(io, PeriodicGraph{N})
    end
    print(io, '(', nv(g), ',', ' ', collect(edges(g)), ')')
end
function print(io::IO, g::PeriodicGraph{N}) where N
    print(io, N)
    for e in edges(g)
        print(io, ' ', e.src, ' ', e.dst.v, ' ', join(e.dst.ofs, ' '))
    end
end

"""
    ==(g1::PeriodicGraph, g2::PeriodicGraph)

Test whether two `PeriodicGraph`s have the same representations.

!!! warning
    Two `PeriodicGraph`s may be isomorphic but still have different representations.
"""
function ==(g1::PeriodicGraph{N}, g2::PeriodicGraph{M}) where {N,M}
    return N == M && nv(g1) == nv(g2) && edges(g1) == edges(g2)
end

ndims(::PeriodicVertex{N}) where {N} = N
ndims(::PeriodicEdge{N}) where {N} = N
ndims(::PeriodicGraph{N}) where {N} = N

"""
    ndims(::PeriodicVertex{N}) where N
    ndims(::PeriodicEdge{N}) where N
    ndims(::PeriodicGraph{N}) where N

Return `N`, the number of repeating dimensions.

## Examples
```jldoctest
julia> ndims(PeriodicVertex{5}(2))
5

julia> ndims(PeriodicEdge3D(1, 3, (0,0,1)))
3

julia> ndims(PeriodicGraph2D(7))
2
```
"""
ndims

Graphs.ne(g::PeriodicGraph) = g.ne[]
Graphs.nv(g::PeriodicGraph) = length(g.nlist)
Graphs.vertices(g::PeriodicGraph) = Base.OneTo(nv(g))
Graphs.edges(g::PeriodicGraph{N}) where {N} = PeriodicEdgeIter{N}(g)
eltype(g::PeriodicGraph{N}) where {N} = PeriodicVertex{N}
Graphs.edgetype(::PeriodicGraph{N}) where {N} = PeriodicEdge{N}
function Graphs.has_edge(g::PeriodicGraph, s, d)
    ((s < 1) | (s > nv(g))) && return false
    #=@inbounds=# begin
        start = g.directedgestart[s]
        lo, hi = s > d ? (1, start-1) : (start, lastindex(g.nlist[s]))
        i = searchsortedfirst(g.nlist[s], d, lo, hi, Lt((x,y)->isless(x.v, y)))
        return i <= length(g.nlist[s]) && g.nlist[s][i].v == d
    end
end

"""
    find_edges(g::PeriodicGraph, s::Int, d::Int)

Return the set of PeriodicVertex `v` of graph `g` such that there is an edge
between a source vertex of identifier `s` and `v`, and the identifier of `v` is `d`.
"""
function find_edges(g::PeriodicGraph{N}, s::Int, d::Int) where N
    ((s < 1) | (s > nv(g))) && return false
    #=@inbounds=# begin
        start = g.directedgestart[s]
        lo, hi = s > d ? (1, start-1) : (start, lastindex(g.nlist[s]))
        rng = searchsorted(g.nlist[s], d, lo, hi, Lt((x,y)->isless(x isa Integer ? x : x.v, y isa Integer ? y : y.v)))
        if s == d
            rng = (2*first(rng) - last(rng) - 1):last(rng)
        end
        return g.nlist[s][rng]
    end
end
function Graphs.has_edge(g::PeriodicGraph, e::PeriodicEdge)
    s, d = e.src, e.dst
    ((s < 1) | (s > nv(g))) && return false
    #=@inbounds=# begin
        start = g.directedgestart[s]
        lo, hi = isindirectedge(e) ? (1, start-1) : (start, lastindex(g.nlist[s]))
        i = searchsortedfirst(g.nlist[s], d, lo, hi, Forward)
        return i <= length(g.nlist[s]) && g.nlist[s][i] == d
    end
end
Graphs.has_edge(g::PeriodicGraph, i, x::PeriodicVertex) = has_edge(g, PeriodicEdge(i, x))
Graphs.outneighbors(g::PeriodicGraph, v::Integer) = g.nlist[v]
Graphs.inneighbors(g::PeriodicGraph, v::Integer) = outneighbors(g, v)
zero(::Type{PeriodicGraph{N}}) where N = PeriodicGraph{N}(0)
Graphs.is_directed(::Type{<:PeriodicGraph}) = false
@static if isdefined(Graphs, :has_contiguous_vertices)
    @inline Graphs.has_contiguous_vertices(::Type{<:PeriodicGraph}) = true
end
Graphs.has_vertex(g::PeriodicGraph, v::Integer) = 1 <= v <= nv(g)
function Graphs.SimpleGraphs.add_vertices!(g::PeriodicGraph{N}, n::Integer) where N
    append!(g.nlist, [PeriodicVertex{N}[] for _ in 1:n])
    append!(g.directedgestart, [1 for _ in 1:n])
    g.width[] = -1
    return n
end
function Graphs.SimpleGraphs.add_vertex!(g::PeriodicGraph{N}) where N
    push!(g.nlist, PeriodicVertex{N}[])
    push!(g.directedgestart, 1)
    g.width[] = -1
    true
end

function _add_edge!(g::PeriodicGraph, e::PeriodicEdge, ::Val{check}) where check
    #=@inbounds=# begin
        s, dst = e.src, e.dst
        neigh = g.nlist[s]
        start = g.directedgestart[s]
        indirectedge = isindirectedge(e)
        lo, hi = indirectedge ? (1, start-1) : (start, lastindex(neigh))
        i = searchsortedfirst(neigh, dst, lo, hi, Forward)
        if check
            i <= length(neigh) && neigh[i] == dst && return false
        end
        g.directedgestart[s] += indirectedge
        insert!(neigh, i, dst)
        return true
    end
end
function Graphs.add_edge!(g::PeriodicGraph, e::PeriodicEdge)
    (src(e) < 1 || src(e) > nv(g) || dst(e) < 1 || dst(e) > nv(g)) && return false
    success = _add_edge!(g, e, Val(true)) && _add_edge!(g, reverse(e), Val(false))
    if success
        g.ne[] += 1
        g.width[] = -1
    end
    return success
end
Graphs.add_edge!(g::PeriodicGraph, i, x::PeriodicVertex) = add_edge!(g, PeriodicEdge(i, x))

function _rem_edge!(g::PeriodicGraph, e::PeriodicEdge, ::Val{check}) where check
    #=@inbounds=# begin
        s, dst = e.src, e.dst
        neigh = g.nlist[s]
        start = g.directedgestart[s]
        indirectedge = isindirectedge(e)
        lo, hi = indirectedge ? (1, start-1) : (start, lastindex(neigh))
        i = searchsortedfirst(neigh, dst, lo, hi, Forward)
        if check
            i <= length(neigh) && neigh[i] == dst || return false
        end
        g.directedgestart[s] -= indirectedge
        deleteat!(neigh, i)
        return true
    end
end
function Graphs.rem_edge!(g::PeriodicGraph, e::PeriodicEdge)
    (src(e) < 1 || src(e) > nv(g) || dst(e) < 1 || dst(e) > nv(g)) && return false
    success = _rem_edge!(g, e, Val(true)) && _rem_edge!(g, reverse(e), Val(false))
    if success
        g.ne[] -= 1
        g.width[] = -1
    end
    return success
end
Graphs.rem_edge!(g::PeriodicGraph, i, x::PeriodicVertex) = rem_edge!(g, PeriodicEdge(i, x))


function Graphs.SimpleGraphs.rem_vertices!(g::PeriodicGraph{N}, t::AbstractVector{<:Integer}, keep_order::Bool=false) where N
    isempty(t) && return collect(1:nv(g))
    sort!(t)
    (first(t) < 1 || last(t) > nv(g)) && throw(ArgumentError("Vertices to be removed must be in the range 1:nv(g)."))

    bt = falses(nv(g))
    bt[t] .= true

    vmap = Int[]
    rev_vmap = zeros(Int, nv(g))

    if keep_order
        append!(vmap, collect(1:nv(g)))
        deleteat!(vmap, bt)
        rev_vmap[vmap] .= 1:length(vmap)
        deleteat!(g.nlist, bt)
        deleteat!(g.directedgestart, bt)
    else
        i_next_vertex_to_del = 1
        next_vertex_to_del = t[i_next_vertex_to_del]
        t_end = length(t)
        g_end = length(g.nlist)
        i = 1
        while i <= g_end
            if i == next_vertex_to_del
                i_next_vertex_to_del += 1
                if i_next_vertex_to_del > t_end
                    next_vertex_to_del = nv(g) + 1
                else
                    next_vertex_to_del = t[i_next_vertex_to_del]
                end
                while t_end >= i_next_vertex_to_del && g_end == t[t_end]
                    t_end -= 1
                    g_end -= 1
                end
                if i < g_end
                    g.nlist[i], g.nlist[g_end] = g.nlist[g_end], g.nlist[i]
                    push!(vmap, g_end)
                    rev_vmap[g_end] = length(vmap)
                end
                g_end -= 1
            else
                push!(vmap, i)
                rev_vmap[i] = length(vmap)
            end
            i += 1
        end
        resize!(g.nlist, g_end)
        resize!(g.directedgestart, g_end)
    end

    counter_edges = 0

    for i in vertices(g)
        neighbors = g.nlist[i]
        remove_edges = falses(length(neighbors))
        startoffset = 1
        for (k, x) in enumerate(neighbors)
            if bt[x.v]
                remove_edges[k] = true
            else
                neigh = PeriodicVertex{N}(rev_vmap[x.v], x.ofs)
                neighbors[k] = neigh
                startoffset += isindirectedge(PeriodicEdge{N}(i, neigh))
            end
        end
        deleteat!(neighbors, remove_edges)
        sort!(neighbors)
        g.directedgestart[i] = startoffset
        counter_edges += startoffset - 1
    end

    g.ne[] = counter_edges
    g.width[] = -1

    return vmap
end

function Graphs.SimpleGraphs.rem_vertex!(g::PeriodicGraph, v::Integer)
    n = nv(g)
    return length(rem_vertices!(g, [v])) == n - 1
end


"""
    PeriodicEdgeIter{N} <: AbstractEdgeIter

Edge iterator type for undirected `N`-periodic graphs.

The iterator only yields edges in the form `(u, v, ofs)` with either `u < v` or
`u == v && ofs > zero(ofs)`.
This is possible because `PeriodicGraph`s are undirected, hence to each edge
`(u, v, ofs)` in the graph corresponds its reverse edge `(v, u, .-ofs)`. The iterator
thus yields each edge of the graph exactly once.
"""
struct PeriodicEdgeIter{N} <: AbstractEdgeIter
    g::PeriodicGraph{N}
end

eltype(::Type{PeriodicEdgeIter{N}}) where {N} = PeriodicEdge{N}
length(iter::PeriodicEdgeIter) = iter.g.ne[]

function iterate(iter::PeriodicEdgeIter{N}, (vertex, neigh)) where N
    nlists = iter.g.nlist
    n = length(nlists)
    #=@inbounds=# while vertex <= n
        if iszero(neigh)
            neigh = iter.g.directedgestart[vertex]
        end
        neighbors = nlists[vertex]
        if neigh > length(neighbors)
            vertex += 1
            neigh = 0
            continue
        end
        return (unsafe_edge{N}(vertex, neighbors[neigh]), (vertex, neigh+1))
    end
    return nothing
end

iterate(iter::PeriodicEdgeIter) = iterate(iter, (1, 0))

function in(edge, iter::PeriodicEdgeIter{N}) where N
    has_edge(iter.g, edge)
end

function cmp(it1::PeriodicEdgeIter{N}, it2::PeriodicEdgeIter{N}) where N
    n = length(it1)
    m = length(it2)
    n == m || return cmp(n,m)
    n == 0 && return 0
    cmpofs = 0
    e1::PeriodicEdge{N}, st1::Tuple{Int,Int} = iterate(it1)
    e2::PeriodicEdge{N}, st2::Tuple{Int,Int} = iterate(it2)
    for _ in 1:n-1
        c = cmp((e1.src, e1.dst.v), (e2.src, e2.dst.v))
        if iszero(c)
            if iszero(cmpofs)
                cmpofs = cmp(e1.dst.ofs, e2.dst.ofs)
            end
            e1, st1 = iterate(it1, st1)
            e2, st2 = iterate(it2, st2)
            continue
        end
        return c
    end
    c = cmp((e1.src, e1.dst.v), (e2.src, e2.dst.v))
    return iszero(c) ? (iszero(cmpofs) ? cmp(e1.dst.ofs, e2.dst.ofs) : cmpofs) : c
end

function isless(it1::PeriodicEdgeIter{N}, it2::PeriodicEdgeIter{N}) where N
    return cmp(it1, it2) < 0
end
function ==(it1::PeriodicEdgeIter{N}, it2::PeriodicEdgeIter{N}) where N
    return iszero(cmp(it1, it2))
end

"""
    vertex_permutation(g::PeriodicGraph, vlist)

Return the `PeriodicGraph` corresponding to `g` with its vertices identifiers
permuted according to `vlist`. `isperm(vlist)` must hold and will not be checked.

See also [`induced_subgraph`](@ref) for the more general case where `vlist` is not a
permutation.

!!! note
    The resulting graph is isomorphic to the initial one, only the representation has
    changed.
"""
function vertex_permutation(g::PeriodicGraph{N}, vlist) where N
    n = length(vlist)
    newvid = Vector{Int}(undef, n)
    for i in 1:n
        newvid[vlist[i]] = i
    end
    edges = Vector{Vector{PeriodicVertex{N}}}(undef, n)
    startoffsets = [1 for _ in 1:n]
    #=@inbounds=# for i in 1:n
        neighs = copy(g.nlist[vlist[i]])
        edges[i] = neighs
        for j in 1:length(neighs)
            dst = neighs[j]
            neigh = PeriodicVertex{N}(newvid[dst.v], dst.ofs)
            neighs[j] = neigh
            startoffsets[i] += isindirectedge(PeriodicEdge{N}(i, neigh))
        end
        sort!(neighs)
    end
    return PeriodicGraph{N}(Ref(g.ne[]), edges, startoffsets, Ref(g.width[]))
end

function Graphs.induced_subgraph(g::PeriodicGraph{N}, vlist::AbstractVector{U}) where {N, U<:Integer}
    allunique(vlist) || __throw_unique_vlist()
    n = length(vlist)
    n == nv(g) && return (vertex_permutation(g, vlist), vlist)
    newvid = zeros(Int, nv(g))
    for i in 1:n
        newvid[vlist[i]] = i
    end

    ne = 0
    edges = Vector{Vector{PeriodicVertex{N}}}(undef, n)
    startoffsets = [1 for _ in 1:n]
    for i in 1:n
        edges[i] = PeriodicVertex{N}[]
        startne = ne
        for dst in g.nlist[vlist[i]]
            v = newvid[dst.v]
            iszero(v) && continue
            neigh = PeriodicVertex{N}(v, dst.ofs)
            push!(edges[i], neigh)
            ne += isindirectedge(PeriodicEdge{N}(i, neigh))
        end
        startoffsets[i] = 1 + ne - startne
        sort!(edges[i])
    end
    return (PeriodicGraph{N}(ne, edges, startoffsets), vlist)
end
@noinline __throw_unique_vlist() = throw(ArgumentError("Vertices in subgraph list must be unique"))

"""
    offset_representatives!(g::PeriodicGraph, offsets)

In-place modifies graph `g` so that the `i`-th vertex of the new initial cell corresponds
to the `i`-th vertex in cell `offsets[i]` compared to the previous initial cell.

!!! note
    The resulting graph is isomorphic to the initial one, only the representation has
    changed.
"""
function offset_representatives!(g::PeriodicGraph{N}, offsets) where N
    n = nv(g)
    length(offsets) == n || __throw_invalid_offsets()
    for i in 1:n
        neighs = g.nlist[i]
        ofsi = offsets[i]
        startoffset = 1
        for j in 1:length(neighs)
            x = neighs[j]
            neigh = PeriodicVertex{N}(x.v, x.ofs .+ ofsi .- offsets[x.v])
            neighs[j] = neigh
            startoffset += isindirectedge(PeriodicEdge{N}(i, neigh))
        end
        sort!(neighs)
        g.directedgestart[i] = startoffset
        g.width[] = -1
    end
    g
end
@noinline __throw_invalid_offsets() = throw(ArgumentError("The size of offsets does not match the number of vertices"))

"""
    swap_axes!(g::PeriodicGraph, t)

In-place modifies graph `g` so that the new initial cell corresponds to the previous
one with its axes swapped according to the permutation `t`.

!!! note
    The resulting graph is isomorphic to the initial one, only the representation has
    changed.
"""
function swap_axes!(g::PeriodicGraph{N}, t) where N
    length(t) == N || __throw_invalid_axesswap()
    tt::Vector{Int} = t isa Vector ? t : collect(t)
    for i in vertices(g)
        neighs = g.nlist[i]
        for j in 1:length(neighs)
            x = neighs[j]
            neigh = PeriodicVertex{N}(x.v, x.ofs[tt])
            neighs[j] = neigh
        end
        sort!(neighs)
        # Note: g.width[] is unchanged
    end
    g
end
@noinline __throw_invalid_axesswap() = throw(DimensionMismatch("The number of axes must match the dimension of the graph"))


function Graphs.connected_components(g::PeriodicGraph)
    nvg = nv(g)
    label = zeros(Int, nvg)

    for u in vertices(g)
        label[u] != 0 && continue
        label[u] = u
        Q = Int[]
        push!(Q, u)
        #=@inbounds=# while !isempty(Q)
            src = popfirst!(Q)
            for dst in outneighbors(g, src)
                vertex = dst.v
                if label[vertex] == 0
                    push!(Q, vertex)
                    label[vertex] = u
                end
            end
        end
    end
    c, d = Graphs.components(label)
    return c
end

"""
    graph_width!(g::PeriodicGraph{N}) where N

Set the `width` internal field of the graph so that the for all `n` ∈ N\\*,
the `n`-th neighbor of any vertex `v` of the initial cell is in a cell
`(i_1, i_2, ..., i_N)` such that `max(abs.((i_1, i_2, ..., i_N))) ≤ 1 + fld((n - 1), width)`.

This function is meant for internal use and will be used whenever the `width` field
is required but unset. If you decide to modify the other internal fields of `g`,
it is probably a good idea to do `g.width[] = -1` so that this function gets
automatically called when needed, unless you are sure the width will not be
affected by your change.

It is never necessary to use this function or to touch the `width` field of the graph
if you only modify the graph using the official API (i.e. if you never directly touch
the fields).
"""
function graph_width!(g::PeriodicGraph{N}) where N
    distances = floyd_warshall_shortest_paths(cellgraph(g)).dists
    extremalpoints = NTuple{N,NTuple{2,Vector{Tuple{Int,Int}}}}([([],[]) for _ in 1:N])
    # a, x ∈ extremalpoints[i][j] where i ∈ ⟦1,N⟧ and j ∈ ⟦1,2⟧ means that
    # vertex x has a neighbor whose offset is a*(-1)^(j-1) along dimension i
    maxa = 1
    for e in edges(g)
        iszero(ofs(e)) && continue
        offset = ofs(e)
        for i in 1:N
            iszero(offset[i]) && continue
            j = signbit(offset[i]) + 1
            a = abs(offset[i])
            if a > maxa
                maxa = a
            end
            push!(extremalpoints[i][j], (a, src(e)))
            push!(extremalpoints[i][3-j], (a, dst(e)))
        end
    end

    width::Rational{Int} = Rational(nv(g)+1)
    for i in 1:N
        if any(isempty.(extremalpoints[i]))
            # TODO check this part
            if width > 1
                width = 1//1
            end
            continue
        end
        for (a1, x1) in extremalpoints[i][1], (a2, x2) in extremalpoints[i][2]
            dist = distances[x1, x2]
            if dist == typemax(Int)
                dist = 1
            end
            d = (dist + 1) // (a1 + a2)
            if d < width
                width = d
            end
        end
    end
    g.width[] = width == nv(g)+1 ? Rational(maxa) : width
end

function Graphs._neighborhood(g::Union{PeriodicGraph{0},PeriodicGraph{1},PeriodicGraph{2},PeriodicGraph{3}}, v::Integer, d::Real, distmx::AbstractMatrix{U}, neighborfn::Function) where U <: Real
    @assert typeof(neighborfn) === typeof(outneighbors)
    N = ndims(g)
    Q = Tuple{PeriodicVertex{N}, U}[]
    d < zero(U) && return Q
    start_vertex = PeriodicVertex{N}(v)
    push!(Q, (start_vertex, zero(U),) )
    n = nv(g)
    width = g.width[]
    if width == -1
        width = graph_width!(g)
    end
    seen_size = n*(2*(1 + fld(d-1, width)) + 1)^N
    seen = falses(seen_size)
    seen[hash_position(start_vertex, n)] = true
    #=@inbounds=# for (src, currdist) in Q
        currdist == d && continue # should be in Q but all its neighbours are too far
        for dst in outneighbors(g, src.v)
            dst = PeriodicVertex{N}(dst.v, dst.ofs .+ src.ofs)
            position = hash_position(dst, n)
            if !seen[position]
                seen[position] = true
                distance = currdist + distmx[src.v, dst.v]
                if distance <= d
                    push!(Q, (dst, distance))
                end
            end
        end
    end
    return Q
end

function Graphs._neighborhood(g::PeriodicGraph{N}, v::Integer, d::Real, distmx::AbstractMatrix{U}, neighborfn::Function) where {N,U <: Real}
    @assert typeof(neighborfn) === typeof(outneighbors)
    Q = Tuple{PeriodicVertex, U}[]
    d < zero(U) && return Q
    start_vertex = PeriodicVertex{N}(v)
    push!(Q, (start_vertex, zero(U),) )
    n = nv(g)
    seen = Set{PeriodicVertex{N}}()
    push!(seen, start_vertex)
    #=@inbounds=# for (src, currdist) in Q
        currdist == d && continue # should be in Q but all its neighbours are too far
        @simd for dst in outneighbors(g, src.v)
            dst = PeriodicVertex(dst.v, dst.ofs .+ src.ofs)
            if dst ∉ seen
                push!(seen, dst)
                distance = currdist + distmx[src.v, dst.v]
                if distance <= d
                    push!(Q, (dst, distance))
                end
            end
        end
    end
    return Q
end

"""
    coordination_sequence(g::PeriodicGraph, v::Integer, dmax)

Compute the list of numbers of `n`-th neighbors of vertex `v` in graph `g`, for
`1 ≤ n ≤ dmax`.
"""
function coordination_sequence(g::PeriodicGraph, v::Integer, dmax)
    Q = Graphs._neighborhood(g, v, dmax, weights(g), outneighbors)
    popfirst!(Q)
    ret = zeros(Int, dmax)
    for (_, d) in Q
        ret[d] += 1
    end
    return ret
end

"""
    cellgraph(g::PeriodicGraph)

Extract a simple graph from `g` by only keeping the edges that are strictly
within the initial cell.

See also [`periodiccellgraph`](@ref) to keep these edges.
"""
function cellgraph(g::PeriodicGraph)
    edgs = [Edge{Int}(src(x), dst(x)) for x in edges(g) if iszero(ofs(x))]
    ret = SimpleGraph(edgs)
    add_vertices!(ret, nv(g) - nv(ret))
    return ret
end

"""
    periodiccellgraph(g::PeriodicGraph)

Extract a simple graph from `g` by removing all indications of offset in the edges.
This means that edges that used to cross the boundaries of the initial cell now
bind the source vertex to the representative of the destination vertex that is in
the initial cell.

Note that these modified edges may turn into loops.

See also [`cellgraph`](@ref) to remove these edges.
"""
function periodiccellgraph(g::PeriodicGraph)
    ret = SimpleGraph([Edge{Int}(i, j.v) for i in vertices(g) for j in outneighbors(g, i)])
    add_vertices!(ret, nv(g) - nv(ret))
    return ret
end

include("dimensionality.jl")
include("precompile.jl")

end
