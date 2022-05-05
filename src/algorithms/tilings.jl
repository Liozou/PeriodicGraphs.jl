struct IdxWithOfs{D}
    idx::Int
    ofs::SVector{D,Int}
end
IdxWithOfs{D}(idx) where {D} = IdxWithOfs{D}(idx, zero(SVector{D,Int}))

struct Tiling{D}
    rings::Vector{Vector{PeriodicVertex{D}}} # each sublist is the list of vertices
    erings::Vector{Vector{Int}} # each sublist is the list of known_pairs of edges
    tiles::Vector{Vector{IdxWithOfs{D}}} # each sublist is a list of indices of rings
    tilesofring::Vector{Tuple{IdxWithOfs{D},IdxWithOfs{D}}} # for each ring, the two tiles attached to it
    ringsofedge::Dict{PeriodicEdge{D},Vector{IdxWithOfs{D}}} # for each edge, the list of associated ring
    rgraph::PeriodicGraph{D} # The graph of rings: two rings are bonded if they share an edge
    prepared_known_pairs::Tuple{Vector{VertexPair{D}},Dict{VertexPair{D},Int}}

    function Tiling{D}(rings, erings, prepared_known_pairs) where {D}
        n = length(erings)
        tiles = Vector{IdxWithOfs{D}}[]
        tilesofring = [(IdxWithOfs{D}(0),IdxWithOfs{D}(0)) for _ in 1:n]
        ringsofedge = Dict{PeriodicEdge{D},Vector{IdxWithOfs{D}}}()
        rgraph = PeriodicGraph{D}(n)
        known_pairs = first(prepared_known_pairs)
        for (i, ering) in enumerate(erings)
            for _e in ering
                u1, u2 = known_pairs[_e]
                @assert (u1, u2) == minmax(u1, u2)
                e = PeriodicEdge{D}(u1.v, u2.v, u2.ofs .- u1.ofs)
                ringsofe = get!(ringsofedge, e, IdxWithOfs{D}[])
                for idx in ringsofe
                    newofs = idx.ofs .- u1.ofs # TODO: check
                    if idx.idx != i || !iszero(newofs)
                        add_edge!(rgraph, i, PeriodicVertex{D}(idx.idx, newofs))
                    end
                end
                push!(ringsofe, IdxWithOfs{D}(i, u1.ofs))
            end
        end

        return new{D}(rings, erings, tiles, tilesofring, ringsofedge, rgraph, prepared_known_pairs)
    end
end

function cycle_at_pos(tiling::Tiling{D}, u::PeriodicVertex{D}) where D
    ring = tiling.rings[u.v]
    len = length(ring)
    ret = Vector{Int}(undef, len)
    convert_to_ering!(ret, ring, len, tiling.prepared_known_pairs, u.ofs)
    return sort!(ret)
end

function tiles_including_cycle(tiling::Tiling{D}, i) where D
    gauss = IterativeGaussianElimination(tiling.erings[i])
    encountered = Set{PeriodicVertex{D}}((PeriodicVertex{D}(i),))
    Q = [(x, 1) for x in neighbors(tiling.rgraph, i)]
    maxdist = typemax(Int)
    tiles = Vector{PeriodicVertex{D}}[]
    for (u, dist) in Q
        dist > 40 && @show dist
        dist > maxdist && break
        if gaussian_elimination!(gauss, cycle_at_pos(tiling, u)) # a sum of previously encountered rings is empty
            track = gauss.lengths
            if i ∈ track
                if !isempty(tiles) # TODO: refine to actually account for tiles on the two sides of cycle i
                    maxdist = dist # TODO: check that this is correct / desired
                end
                push!(tiles, sort!([first(Q[x]) for x in track]))
            end
        end
        for x in neighbors(tiling.rgraph, u)
            x ∈ encountered && continue
            push!(encountered, x)
            push!(Q, (x, dist+1))
        end
    end
    @show maxdist
    @show tiles
end


function tiling(g::PeriodicGraph{D}, depth=15, symmetries::AbstractSymmetryGroup=NoSymmetryGroup(g), dist::DistanceRecord=DistanceRecord(g,depth)) where D
    _rings, symms, erings, prepared_known_pairs = strong_erings(g, depth, symmetries, dist)
    rings = Vector{PeriodicVertex{D}}[[reverse_hash_position(x, g) for x in r] for r in _rings]
    tiling = Tiling{D}(rings, erings, prepared_known_pairs)
    for i in 1:length(erings)
        tiling.tilesofring[i][2].idx == 0 || continue
        ts = tiles_including_cycle(tiling, i)
    end
end
