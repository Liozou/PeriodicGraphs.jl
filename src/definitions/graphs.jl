# PeriodicGraph definition and basic functions

export PeriodicGraph, PeriodicGraph1D, PeriodicGraph2D, PeriodicGraph3D

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

function PeriodicGraph{N}(ne::Integer, t, s) where N
    return PeriodicGraph{N}(Ref(ne), t, s, Ref(-1//1))
end

"""
    PeriodicGraph{N}(nv::Integer=0)

Construct a `PeriodicGraph{N}` with `nv` vertices and 0 edge.
"""
function PeriodicGraph{N}(n::Integer = 0) where N
    return PeriodicGraph{N}(0, [PeriodicVertex{N}[] for _ in 1:n], [1 for _ in 1:n])
end

function from_edges(nv, t::AbstractVector{PeriodicEdge{N}}) where {N}
    isempty(t) && return PeriodicGraph{N}(nv)

    directedgestart = ones(Int, nv)
    counterdirect = zeros(Int, nv)
    numdifferentindirect = zeros(Int, nv)
    last_e = first(t)
    numdifferentindirect[last_e.dst.v] = 1
    for e in t
        s = e.src
        d = e.dst.v
        directedgestart[d] += 1
        counterdirect[s] += 1
        numdifferentindirect[d] += ((s != last_e.src) | (d != last_e.dst.v))
        last_e = e
    end

    nlist = Vector{Vector{PeriodicVertex{N}}}(undef, nv)
    counterindirect = Vector{Vector{Int}}(undef, nv)
    for i in 1:nv
        nlist[i] = Vector{PeriodicVertex{N}}(undef, counterdirect[i] + directedgestart[i] - 1)
        counterdirect[i] = 0
        counterindirect[i] = zeros(Int, numdifferentindirect[i] + 1)
    end

    # First pass: set the direct edges
    idx_indirect = ones(Int, nv)
    last_e = first(t)
    for e in t
        s = e.src
        d = e.dst.v
        nlist[s][directedgestart[s] + counterdirect[s]] = e.dst
        counterdirect[s] += 1
        if ((s != last_e.src) | (d != last_e.dst.v))
            last_d = last_e.dst.v
            idx_indirect[last_d] += 1
            val = idx_indirect[last_d]
            counterindirect[last_d][val] += counterindirect[last_d][val-1]
        end
        counterindirect[d][idx_indirect[d]] += 1
        last_e = e
    end

    @simd for i in 1:nv; idx_indirect[i] = 1; end

    # Second pass: set the indirect edges
    last_e = first(t)
    current_idx = idx_indirect[last_e.dst.v]
    for e in t
        s = e.src
        d = e.dst.v
        if ((s != last_e.src) | (d != last_e.dst.v))
            idx_indirect[last_e.dst.v] += 1
            current_idx = idx_indirect[d]
        end
        nlist[d][counterindirect[d][current_idx]] = PeriodicVertex(s, .- e.dst.ofs)
        counterindirect[d][current_idx] -= 1
        last_e = e
    end

    return PeriodicGraph{N}(length(t), nlist, directedgestart)
end

function from_edges(t::AbstractVector{PeriodicEdge{N}}) where N
    from_edges(maximum(max(e.src, e.dst.v) for e in t; init=0), t)
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
PeriodicGraph2D(3, PeriodicEdge2D[(1, 1, (1,0)), (1, 3, (0,1))])

julia> ne(g)
2
```
"""
function PeriodicGraph{N}(nv::Integer, t::AbstractVector{PeriodicEdge{N}}) where N
    for (i, e) in enumerate(t)
        if !isdirectedge(e)
            t[i] = reverse(e)
        end
    end
    sort!(t); unique!(t)
    return from_edges(nv, t)
end

PeriodicGraph(nv::Integer, t::AbstractVector{PeriodicEdge{N}}) where {N} = PeriodicGraph{N}(nv, t)
function PeriodicGraph{N}(t::AbstractVector{PeriodicEdge{N}}) where N
    PeriodicGraph{N}(maximum(max(e.src, e.dst.v) for e in t; init=0), t)
end
PeriodicGraph(t::AbstractVector{PeriodicEdge{N}}) where {N} = PeriodicGraph{N}(t)

"""
    ==(g1::PeriodicGraph, g2::PeriodicGraph)

Test whether two `PeriodicGraph`s have the same representations.

!!! warning
    Two `PeriodicGraph`s may be isomorphic but still have different representations.
"""
function ==(g1::PeriodicGraph{N}, g2::PeriodicGraph{M}) where {N,M}
    return N == M && nv(g1) == nv(g2) && edges(g1) == edges(g2)
end

function hash(g::PeriodicGraph, h::UInt)
    ret = h
    for e in edges(g)
        ret = hash(e, ret)
    end
    return ret
end

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

function PeriodicGraph{N}(g::PeriodicGraph{N}) where N
    PeriodicGraph{N}(g.ne[], [copy(x) for x in g.nlist], copy(g.directedgestart))
end
PeriodicGraph(g::PeriodicGraph{N}) where {N} = PeriodicGraph{N}(g)
Base.copy(g::PeriodicGraph{N}) where {N} = PeriodicGraph{N}(g)
