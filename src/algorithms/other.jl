# Other utilities specific to periodic graphs

export offset_representatives!,
       swap_axes!,
       make_supercell,
       truncated_graph,
       quotient_graph,
       slice_graph

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
            startoffset += !isdirectedge(PeriodicEdge{N}(i, neigh))
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

struct MetaClock{T}
    maxs::T
    factors::Vector{Int}
    len::Int
end
function MetaClock(t)
    t isa AbstractVector || (t = collect(t))
    factors = Vector{Int}(undef, length(t))
    cumprod!(factors, t)
    len = pop!(factors)
    pushfirst!(factors, 1)
    return MetaClock(t, factors, len)
end
Base.length(x::MetaClock) = x.len
Base.eltype(::Type{MetaClock{T}}) where {T} = Vector{Int}
function Base.iterate(x::MetaClock, state::Vector{Int}=[-1; zeros(Int, length(x.maxs)-1)])
    maxs = x.maxs
    for i in 1:length(maxs)
        y = state[i] + 1
        if y == maxs[i]
            state[i] = 0
        else
            state[i] = y
            return (state, state)
        end
    end
    return nothing
end

"""
    make_supercell(g::PeriodicGraph, t)

Return a graph isomorphic to the input `g` whose its unit cell is a repetition of that of
`g`, each dimension `i` being repeated `t[i]` times.
It follows that the number of vertices of `make_supercell(g, t)` is `prod(t)*nv(g)`

`t` must be an interator over positive integers.
"""
function make_supercell(g::PeriodicGraph{N}, t::S) where {N,S}
    length(t) == N || __throw_invalid_axesswap()
    N == 0 && return g
    newedges = PeriodicEdge{N}[]
    __check_nonpositive_axe(minimum(t))
    n = nv(g)
    clock = MetaClock(t)
    for e in edges(g)
        (src, (dst, ofs)) = e
        for pos in clock
            factor = n*sum(prod, zip(clock.factors, pos); init=0)
            newofs, newpos = eachrow(reinterpret(reshape, Int, fldmod.(ofs .+ pos, t)))
            newfactor = n*sum(prod, zip(clock.factors, newpos); init=0)
            push!(newedges, PeriodicEdge{N}(factor+src, newfactor+dst, newofs))
        end
    end
    return PeriodicGraph{N}(length(clock)*n, newedges)
end
__check_nonpositive_axe(i) = i > 0 || throw(DomainError(i, "All supercell dimensions must be strictly positive"))

"""
    truncated_graph(g::PeriodicGraph)

Extract a simple graph from `g` by only keeping the edges that are strictly
within the initial cell.

See also [`quotient_graph`](@ref) to keep all of these edges and
[`slice_graph`](@ref) to keep only some of these edges.
"""
function truncated_graph(g::PeriodicGraph)
    edgs = [Edge{Int}(x.src, x.dst.v) for x in edges(g) if iszero(x.dst.ofs)]
    ret = SimpleGraph(edgs)
    add_vertices!(ret, nv(g) - nv(ret))
    return ret
end

"""
    quotient_graph(g::PeriodicGraph)

Extract a simple graph from `g` by removing all indications of offset in the edges.
This means that edges that used to cross the boundaries of the initial cell now
bind the source vertex to the representative of the destination vertex that is in
the initial cell.

Note that these modified edges may turn into loops.

See also [`truncated_graph`](@ref) to remove all of these edges and
[`slice_graph`](@ref) to keep only some of these edges.
"""
function quotient_graph(g::PeriodicGraph)
    ret = SimpleGraph([Edge{Int}(i, j.v) for i in vertices(g) for j in outneighbors(g, i)])
    add_vertices!(ret, nv(g) - nv(ret))
    return ret
end

"""
    slice_graph(g::PeriodicGraph{D}, remove::Union{SVector{N},NTuple{N}}) where {D,N}

Extract a `PeriodicGraph{D-N}` from `g` by removing all edges that have an offset `o` such
that `!iszero(o[remove])` and shrinking the resulting offsets. In other words, remove the
dimensions in `remove`.

To only remove the edges while keeping the same number of dimensions, use
[`slice_graph(g, collect(remove))`](@ref slice_graph(g::PeriodicGraph{D}, remove::Vector{<:Integer}) where D)

`remove` is assumed to be sorted and to contain unique elements.

!!! warning
    No verification of the previous assumption will be performed.
"""
function slice_graph(g::PeriodicGraph{D}, _remove::Union{SVector{N},NTuple{N}}) where {D,N}
    K = D - N
    remove = SVector{N,Int}(_remove)
    _map = Vector{Int}(undef, K)
    next_remove_j = 1
    next_remove = isempty(remove) ? 0 : first(remove)
    counter = 1
    for i in 1:D
        if i == next_remove
            next_remove_j += 1
            next_remove = next_remove_j > length(remove) ? 0 : remove[next_remove_j]
        else
            _map[counter] = i
            counter += 1
        end
    end
    map = SVector{K,Int}(_map)
    _edges = PeriodicEdge{K}[PeriodicEdge{K}(s, d, o[map]) for (s, (d, o)) in edges(g) if iszero(o[remove])]
    return PeriodicGraph{K}(nv(g), _edges)
end


"""
    slice_graph(g::PeriodicGraph{D}, remove::Vector{<:Integer}) where D

Extract a `PeriodicGraph{D}` from `g` by removing all edges that have an offset `o` such
that `!iszero(o[remove])`.

Contrarily to the [`slice_graph(g::PeriodicGraph{D}, remove::Union{SVector{N},NTuple{N}}) where {D,N}`](@ref)
method, the dimensions along which the edges are erased are still kept here.

`remove` is assumed to be sorted and to contain unique elements.

!!! warning
    No verification of the previous assumption will be performed.
"""
function slice_graph(g::PeriodicGraph{D}, remove::Vector{<:Integer}) where D
    _edges = PeriodicEdge{D}[e for e in edges(g) if iszero(last(last(e))[remove])]
    return PeriodicGraph{D}(nv(g), _edges)
end
