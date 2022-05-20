export Tiling, tilingof


struct Tiling{D}
    rings::Vector{Vector{PeriodicVertex{D}}} # each sublist is the list of vertices
    erings::Vector{Vector{Int}} # each sublist is the list of edge indices in kp of a ring
    tiles::Vector{Vector{PeriodicVertex{D}}} # each sublist is a list of indices of rings
    tileedges::Vector{Vector{Int}} # each sublist is the list of edge indices in kp of a tile
    tiledict::Dict{Vector{Int},Int} # for each tile (represented as its unique edges), its index
    tilesofring::Vector{Tuple{PeriodicVertex{D},PeriodicVertex{D}}} # for each ring, the two tiles attached to it
    ringsofedge::Dict{PeriodicEdge{D},Vector{PeriodicVertex{D}}} # for each edge, the list of associated ring
    rgraph::PeriodicGraph{D} # The graph of rings: two rings are bonded if they share an edge
    kp::EdgeDict{D}

    function Tiling{D}(rings, erings, kp) where {D}
        n = length(erings)
        ringsofedge = Dict{PeriodicEdge{D},Vector{PeriodicVertex{D}}}()
        rgraph = PeriodicGraph{D}(n)
        for (i, ering) in enumerate(erings)
            for _e in ering
                u1, u2 = kp[_e]
                @assert (u1, u2) == minmax(u1, u2)
                e = PeriodicEdge{D}(u1.v, u2.v, u2.ofs .- u1.ofs)
                ringsofe = get!(ringsofedge, e, PeriodicVertex{D}[])
                for idx in ringsofe
                    newofs = u1.ofs - idx.ofs
                    if idx.v != i || !iszero(newofs)
                        add_edge!(rgraph, i, PeriodicVertex{D}(idx.v, newofs))
                    end
                end
                push!(ringsofe, PeriodicVertex{D}(i, u1.ofs))
            end
        end

        tiles = Vector{PeriodicVertex{D}}[]
        tilesofring = [(PeriodicVertex{D}(0),PeriodicVertex{D}(0)) for _ in 1:n]
        tiledict = Dict{Vector{Int},Int}()
        return new{D}(rings, erings, tiles, tiledict, tilesofring, ringsofedge, rgraph, kp)
    end
end

function countconsecutiveunique(l::Vector{T}) where T
    ret = Tuple{Int,T}[]
    count = 1
    lst = first(l)
    for x in Iterators.rest(l, 2)
        if x == lst
            count += 1
        else
            push!(ret, (count, lst))
            count = 1
            lst = x
        end
    end
    push!(ret, (count, lst))
    return ret
end

# const superscriptdigits = ('⁰', '¹', '²', '³', '⁴', '⁵', '⁶', '⁷', '⁸', '⁹')
# function tosuperscript(io::IO, x::Int)
#     ds = digits(x)
#     for d in Iterators.reverse(ds)
#         print(io, superscriptdigits[d+1])
#     end
#     nothing
# end

function Base.show(io::IO, ::MIME"text/plain", t::Tiling{D}) where D
    tiles = [sort!([length(t.rings[x.v]) for x in _t]) for _t in t.tiles]
    sort!(tiles)
    sort!(tiles; by=length)
    ctiles = countconsecutiveunique([countconsecutiveunique(t) for t in tiles])
    for (i, (c, tile)) in enumerate(ctiles)
        i == 1 || print(io, " + ")
        c == 1 || print(io, c)
        print(io, '[')
        for (j, (c2, x)) in enumerate(tile)
            j == 1 || print(io, ", ")
            print(io, x)
            if c2 != 1
                print(io, '^', c2)
            end
        end
        print(io, ']')
    end
    nothing
end



struct TilingAroundCycle{D}
    gauss::IterativeGaussianEliminationDecomposition
    encountered::Set{PeriodicVertex{D}}
    Q::Vector{Tuple{PeriodicVertex{D},Int}}
    restart::Base.RefValue{Int}
    tiles::Vector{Vector{PeriodicVertex{D}}}
end
function TilingAroundCycle{D}(i) where D
    gauss = IterativeGaussianEliminationDecomposition()
    encountered = Set{PeriodicVertex{D}}((PeriodicVertex{D}(i),))
    Q = [(PeriodicVertex{D}(i), 0)]
    tiles = Vector{PeriodicVertex{D}}[]
    return TilingAroundCycle{D}(gauss, encountered, Q, 1, tiles)
end

function cycle_at_pos(tiling::Tiling{D}, u::PeriodicVertex{D}) where D
    ring = tiling.rings[u.v]
    len = length(ring)
    ret = Vector{Int}(undef, len)
    convert_to_ering!(ret, ring, len, tiling.kp, u.ofs)
    return sort!(ret)
end

function explore_around_cycle!(tac::TilingAroundCycle{D}, tiling::Tiling{D}, untilfirstfound=false) where D
    Q = tac.Q
    restart = Q.restart[]
    maxdist = untilfirstfound ? typemax(Int) : isone(restart) ? 3 : last(Q[restart]) + 1
    newrestart = restart
    for (u, dist) in Iterators.rest(Q, restart)
        dist > maxdist && break
        newrestart += 1
        if gaussian_elimination!(tac.gauss, cycle_at_pos(tiling, u)) # a sum of previously encountered rings is empty
            track = retrieve_track!(tac.gauss)
            if last(track) == 1
                if untilfirstfound
                    maxdist = dist
                end
                push!(tac.tiles, sort!([first(Q[x]) for x in track]))
            end
        end
        for x in neighbors(tiling.rgraph, u)
            x ∈ tac.encountered && continue
            push!(tac.encountered, x)
            push!(Q, (x, dist+1))
        end
    end
    tac.restart[] = newrestart
    nothing
end

function unique_edges(tile::Vector{PeriodicVertex{D}}, tiling::Tiling{D}) where D
    doubleedges = VertexPair{D}[]
    for t in tile
        ring = tiling.rings[t.v]
        lst = PeriodicVertex{D}(last(ring).v, last(ring).ofs .+ t.ofs)
        for r in ring
            x = PeriodicVertex{D}(r.v, r.ofs .+ t.ofs)
            push!(doubleedges, minmax(lst, x))
            lst = x
        end
    end
    sort!(doubleedges)
    ret = doubleedges[1:2:end]
    bef = first(ret)
    for i in 2:length(ret)
        y = ret[i]
        bef ≥ y && return VertexPair{D}[] # error: one edge is included more than twice
        bef = y
    end
    return doubleedges
end

function canonical_tile!(tiling::Tiling{D}, tile::Vector{PeriodicVertex{D}}) where D
    normalized, ofs = normalize_cycle!(tile)
    uniqueedges = unique_edges(normalized, tiling)
    isempty(uniqueedges) && return PeriodicVertex{D}(0)
    n = length(tiling.tiles) + 1
    idx = get!(tiling.tiledict, uniqueedges, n)
    if idx == n
        push!(tiling.tiles, normalized)
        push!(tiling.tileedges, uniqueedges)
        for x in normalized
            fst, lst = tiling.tilesofring[x.v]
            new = PeriodicVertex{D}(idx, x.ofs + ofs)
            if iszero(fst.v)
                tiling.tilesofring[x.v] = (new, lst)
            elseif iszero(lst.v)
                tiling.tilesofring[x.v] = (fst, new)
            else
                if !(fst == new || lst == new)
                    @show fst
                    @show lst
                    @show new
                    @show x
                    @show normalized
                end
            end
        end
    elseif tiling.tiles[idx] != normalized
        # two tiles share the same edges but are made from different rings
        @error "Colliding tiles"
    end
    return PeriodicVertex{D}(idx, ofs)
end

function tilingof(g::PeriodicGraph{D}, depth=15, symmetries::AbstractSymmetryGroup=NoSymmetryGroup(g), dist::DistanceRecord=DistanceRecord(g,depth)) where D
    _rings, symms, erings, kp = strong_erings(g, depth, symmetries, dist)
    rings = Vector{PeriodicVertex{D}}[[reverse_hash_position(x, g) for x in r] for r in _rings]
    tiling = Tiling{D}(rings, erings, kp)
    exploration = [TilingAroundCycle{D}(i) for i in 1:n]
    for i in 1:length(erings)
        last(tiling.tilesofring[i]).v == 0 || continue
        ts = tiles_including_cycle(tiling, i)
        if length(ts) > 2
            println("Multiple tile candidates for ring ", i)
        else
            t1, t2 = ts
            canonical_tile!(tiling, t1)
            canonical_tile!(tiling, t2)
        end
    end
    return tiling
end
