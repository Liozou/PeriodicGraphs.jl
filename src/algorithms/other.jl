# Other utilities specific to periodic graphs

export offset_representatives!, swap_axes!, cellgraph, periodiccellgraph

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
