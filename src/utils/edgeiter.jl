# Utility for iterating over the edges of a PeriodicGraph

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
