# Extension of the Graphs._neighborhood function

export coordination_sequence

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

function Graphs._neighborhood(g::Union{PeriodicGraph{0},PeriodicGraph{1},PeriodicGraph{2},PeriodicGraph{3}}, v::Integer, d::Real, distmx::AbstractMatrix{U}, ::typeof(outneighbors)) where U <: Real
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

function Graphs._neighborhood(g::PeriodicGraph{N}, v::Integer, d::Real, distmx::AbstractMatrix{U}, ::typeof(outneighbors)) where {N,U <: Real}
    Q = Tuple{PeriodicVertex, U}[]
    d < zero(U) && return Q
    start_vertex = PeriodicVertex{N}(v)
    push!(Q, (start_vertex, zero(U),) )
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
