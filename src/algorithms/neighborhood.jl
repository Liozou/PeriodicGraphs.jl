# Extension of the Graphs._neighborhood function

export coordination_sequence

function _sortunique!(vec)
    sort!(vec)
    last_x = 0
    toremove = Int[]
    for (i, (x, _)) in enumerate(vec)
        if x == last_x
            push!(toremove, i-1)
        else
            last_x = x
        end
    end
    deleteat!(vec, toremove)
    vec
end

"""
    graph_width!(g::PeriodicGraph{N}) where N

Set the `width` internal field of the graph so that the for most `n` ∈ N\\*,
the `n`-th neighbor of any vertex `v` of the initial cell is in a cell
`(i_1, i_2, ..., i_N)` such that `max(abs.((i_1, i_2, ..., i_N))) ≤ 1 + fld((n - 1), width)`.
Return the new `width`.

This function is only a heuristic, it may produce a `width` that does not satisfy the
condition.

This function returns the current `width` if it is not equal to `-1` (internal value used
to mark an unset `width`). If you decide to modify the other internal fields of `g`, it is
probably a good idea to do `g.width[] = -1` so that this function gets automatically called
when needed, unless you are sure the width will not be affected by your change.
"""
function graph_width!(g::PeriodicGraph{N}) where N
    previous_width = g.width[]
    previous_width == -1 || return previous_width
    distances = floyd_warshall_shortest_paths(truncated_graph(g)).dists
    extremalpoints = SVector{N,NTuple{2,Vector{Tuple{Int,Int}}}}([(Tuple{Int,Int}[], Tuple{Int,Int}[]) for _ in 1:N])
    # x, a ∈ extremalpoints[i][j] where i ∈ ⟦1,N⟧ and j ∈ ⟦1,2⟧ means that
    # vertex x has a neighbor whose offset is a*(-1)^(j-1) along dimension i
    maxa = 0
    for e in edges(g)
        _, (_, offset) = e
        iszero(offset) && continue
        for i in 1:N
            iszero(offset[i]) && continue
            j = signbit(offset[i]) + 1
            a = abs(offset[i])
            maxa = max(maxa, a)
            push!(extremalpoints[i][j], (e.src, a))
            push!(extremalpoints[i][3-j], (e.dst.v, a))
        end
    end

    if maxa == 0# non-periodic graph with no edge crossing cells
        g.width[] = 1//0
        return 1//0
    end

    width::Rational{Int} = Rational(nv(g)+1)
    for i in 1:N
        extri_pre, extri_post = extremalpoints[i]
        _sortunique!(extri_pre)
        _sortunique!(extri_post)
        isempty(extri_pre) && continue
        for (x1, a1) in extri_pre, (x2, a2) in extri_post
            dist = distances[x1, x2]
            if dist == typemax(Int)
                dist = 1
            end
            d = (dist + 1) // (a1 + a2)
            width = min(width, d)
        end
    end
    g.width[] = width == nv(g)+1 ? Rational(maxa) : width
end

function Graphs._neighborhood(g::Union{PeriodicGraph{0},PeriodicGraph{1},PeriodicGraph{2},PeriodicGraph{3}}, v::Integer, d::Real, distmx::AbstractMatrix{U}, ::typeof(outneighbors)) where U <: Real
    N = ndims(g)
    Q = Tuple{PeriodicVertex{N}, U}[]
    d < zero(U) && return Q
    start_vertex = PeriodicVertex{N}(v)
    sizehint!(Q, floor(Int, 6N*(d*ne(g)/(N^2*(1+nv(g))))^N)) # heuristically determined
    push!(Q, (start_vertex, zero(U),) )
    n = nv(g)
    if n < 100 # heuristic constant
        width = graph_width!(g)
        seen_size = n*(2*(1 + fld(d-1, width)) + 1)^N
    else
        seen_size = n
    end
    seen = falses(seen_size)
    seen[v] = true
    @inbounds for (src, currdist) in Q
        currdist == d && continue # should be in Q but all its neighbours are too far
        for dst in outneighbors(g, src)
            position = hash_position(dst, n)
            if position > length(seen)
                append!(seen, falses(position - length(seen)))
            end
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

function Graphs._neighborhood(g::PeriodicGraph{N}, v::Integer, d::Real, distmx::AbstractMatrix{U}, ::typeof(outneighbors)) where {N,U <: Real}
    Q = Tuple{PeriodicVertex, U}[]
    d < zero(U) && return Q
    start_vertex = PeriodicVertex{N}(v)
    hintsize = floor(Int, 6N*(d*ne(g)/(N^2*(1+nv(g))))^N)  # heuristically determined
    sizehint!(Q, hintsize)
    push!(Q, (start_vertex, zero(U),) )
    seen = Set{PeriodicVertex{N}}()
    sizehint!(seen, hintsize)
    push!(seen, start_vertex)
    @inbounds for (src, currdist) in Q
        currdist == d && continue # should be in Q but all its neighbours are too far
        for dst in outneighbors(g, src.v)
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
