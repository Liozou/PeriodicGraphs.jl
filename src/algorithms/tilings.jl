struct IdxWithOfs{D}
    idx::Int
    ofs::SVector{D,Int}
end
IdxWithOfs{D}(idx) = IdxWithOfs{D}(idx, zero(SVector{D,Int}))

struct Tiling{D}
    rings::Vector{Vector{Int}} # each sublist is the list of hash_position of vertices
    erings::Vector{Vector{Int}} # each sublist is the list of known_pairs of edges
    tiles::Vector{Vector{IdxWithOfs{D}}} # each sublist is a list of indices of rings
    tilesofring::Vector{Tuple{IdxWithOfs{D},IdxWithOfs{D}}} # for each ring, the two tiles attached to it
    ringsofedge::Dict{PeriodicEdge{D},Vector{IdxWithOfs{D}}} # for each edge, the list of associated ring
    rgraph::PeriodicGraph{D} # The graph of rings: two rings are bonded if they share an edge

    function Tiling{D}(g, rings, erings, known_pairs)
        n = length(erings)
        tiles = Vector{IdxWithOfs{D}}[]
        tilesofring = [(IdxWithOfs{D}(0),IdxWithOfs{D}(0)) for _ in 1:n]
        ringsofedge = Dict{PeriodicEdge{D},Vector{IdxWithOfs{D}}}()
        rgraph = PeriodicGraph{D}(n)
        for (i, ering) in enumerate(erings)
            for _e in ering
                u1, u2 = known_pairs[_e]
                @assert u1, u2 == minmax(u1, u2)
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
        return new{D}(rings, erings, tiles, tilesofring, ringsofedge, rgraph)
    end
end

function ering_graph(t::Tiling{D}, known_pairs) where D
    g = PeriodicGraph{D}(length(t.erings))

end

function tiles_including_cycle(tiling::Tiling{D}, i) where D
    gauss = IterativeGaussianElimination(tiling.erings[i])
    encountered = Set{PeriodicVertex{D}}(PeriodicVertex{D}(i))
    Q = [(x, 1) for x in neighbors(tiling.rgraph, i)]
    maxdist = typemax(Int)
    tiles = Vector{Int}[]
    for (u, dist) in Q
        dist > maxdist && break
        if gaussian_elimination!(gauss, u) # a sum of previously encountered rings is empty
            track = gauss.lengths
            if i ∈ track
                if !isempty(tiles) # TODO: refine to actually account for tiles on the two sides of cycle i
                    maxdist = dist # TODO: check that this is correct / desired
                end
                push!(tiles, sort!(collect(track)))
            end
        end
        for x in neighbors(tiling.rgraph, u)
            x ∈ encountered && continue
            push!(encountered, x)
            push!(Q, (x, dist+1))
        end
    end
end


function tiling(g::PeriodicGraph{D}, depth=15, symmetries::AbstractSymmetryGroup=NoSymmetryGroup(g), dist::DistanceRecord=DistanceRecord(g,depth)) where D
    rings, symms, erings, prepared_known_pairs = strong_erings(g, depth, symmetries, dist)
    tiling = Tiling{D}(rings, erings)
    for (i, (ring, ering)) in enumerate(zip(rings, erings))
        tiling.tilesofring[i][2].v == 0 || continue
        ts = tiles_including_cycle(tiling, ering)
    end
end