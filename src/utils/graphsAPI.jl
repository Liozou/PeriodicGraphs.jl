# Simple functions required by the AbstractGraph API

export find_edges

import Base.Order: Forward, Lt

Graphs.ne(g::PeriodicGraph) = g.ne[]
Graphs.nv(g::PeriodicGraph) = length(g.nlist)
Graphs.vertices(g::PeriodicGraph) = Base.OneTo(nv(g))
Graphs.edges(g::PeriodicGraph{N}) where {N} = PeriodicEdgeIter{N}(g)
# Note: the following is intentionally not eltype(::Type{PeriodicGraph{N}}) where N
eltype(::PeriodicGraph{N}) where {N} = PeriodicVertex{N}
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
        lo, hi = isdirectedge(e) ? (start, lastindex(g.nlist[s])) : (1, start-1)
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
    # Note: g.width[] is unchanged as long as there is no modification of the edges
    return n
end
function Graphs.SimpleGraphs.add_vertex!(g::PeriodicGraph{N}) where N
    push!(g.nlist, PeriodicVertex{N}[])
    push!(g.directedgestart, 1)
    # Note: g.width[] is unchanged as long as there is no modification of the edges
    true
end

function _add_edge!(g::PeriodicGraph, e::PeriodicEdge, ::Val{check}) where check
    #=@inbounds=# begin
        s, dst = e.src, e.dst
        neigh = g.nlist[s]
        start = g.directedgestart[s]
        _directedge = isdirectedge(e)
        lo, hi = _directedge ? (start, lastindex(neigh)) : (1, start-1)
        i = searchsortedfirst(neigh, dst, lo, hi, Forward)
        if check
            i <= length(neigh) && neigh[i] == dst && return false
        end
        g.directedgestart[s] += !_directedge
        insert!(neigh, i, dst)
        return true
    end
end
function Graphs.add_edge!(g::PeriodicGraph, e::PeriodicEdge)
    (e.src < 1 || e.src > nv(g) || e.dst.v < 1 || e.dst.v > nv(g)) && return false
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
        _directedge = isdirectedge(e)
        lo, hi = _directedge ? (start, lastindex(neigh)) : (1, start-1)
        i = searchsortedfirst(neigh, dst, lo, hi, Forward)
        if check
            i <= length(neigh) && neigh[i] == dst || return false
        end
        g.directedgestart[s] -= !_directedge
        deleteat!(neigh, i)
        return true
    end
end
function Graphs.rem_edge!(g::PeriodicGraph, e::PeriodicEdge)
    (e.src < 1 || e.src > nv(g) || e.dst.v < 1 || e.dst.v > nv(g)) && return false
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
                startoffset += !isdirectedge(PeriodicEdge{N}(i, neigh))
            end
        end
        deleteat!(neighbors, remove_edges)
        sort!(neighbors)
        g.directedgestart[i] = startoffset
        counter_edges += startoffset - 1
    end

    if g.ne[] != counter_edges
        g.ne[] = counter_edges
        g.width[] = -1
    end # otherwise, no modification to the width either

    return vmap
end

function Graphs.SimpleGraphs.rem_vertex!(g::PeriodicGraph, v::Integer)
    n = nv(g)
    return length(rem_vertices!(g, [v])) == n - 1
end


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
    return first(Graphs.components(label))
end


"""
    vertex_permutation(g::PeriodicGraph, vlist)

Return the `PeriodicGraph` corresponding to `g` with its vertices identifiers
permuted according to `vlist`. `isperm(vlist)` must hold and will not be checked.

See also `Graphs.induced_subgraph` for the more general case where `vlist` is not a
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
            startoffsets[i] += !isdirectedge(PeriodicEdge{N}(i, neigh))
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
            ne += !isdirectedge(PeriodicEdge{N}(i, neigh))
        end
        startoffsets[i] = 1 + ne - startne
        sort!(edges[i])
    end
    return (PeriodicGraph{N}(ne, edges, startoffsets), vlist)
end
@noinline __throw_unique_vlist() = throw(ArgumentError("Vertices in subgraph list must be unique"))

function Graphs.induced_subgraph(g::PeriodicGraph{N}, vlist::AbstractVector{Bool}) where N
    length(vlist) == nv(g) || throw(BoundsError(g, vlist))
    return induced_subgraph(g, findall(vlist))
end

"""
    OffsetVertexIterator{D} <: AbstractVector{PeriodicVertex{D}}
    OffsetVertexIterator(ofs::SVector{D,Int}, list::AbstractVector{PeriodicVertex{D}}) where D

Iterator type that yields the sequence of `PeriodicVertex` in `list`, each offset by the
input `ofs`.
"""
struct OffsetVertexIterator{D} <: AbstractVector{PeriodicVertex{D}}
    ofs::SVector{D,Int}
    nlist::Vector{PeriodicVertex{D}}
end

function Base.getindex(x::OffsetVertexIterator{D}, i::Int) where D
    neigh = x.nlist[i]
    return PeriodicVertex{D}(neigh.v, neigh.ofs .+ x.ofs)
end
Base.size(x::OffsetVertexIterator) = (length(x.nlist),)
Base.IndexStyle(::Type{OffsetVertexIterator{D}}) where {D} = Base.IndexLinear()

for (neigh, deg) in ((:neighbors, :degree),
                     (:inneighbors, :indegree),
                     (:outneighbors, :outdegree))
    @eval begin
        function (Graphs.$neigh)(g::PeriodicGraph{D}, u::PeriodicVertex{D}) where D
            OffsetVertexIterator{D}(u.ofs, g.nlist[u.v])
        end
        function (Graphs.$deg)(g::PeriodicGraph{D}, u::PeriodicVertex{D}) where D
            length(($neigh)(g, u))
        end
    end
end

function reverse_hash_position(hash::Integer, g::PeriodicGraph{D}) where D
    reverse_hash_position(hash, nv(g), Val(D))
end
