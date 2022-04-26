macro enforce(expr) # strong @assert
    msg = string(expr)
    return :($(esc(expr)) ? $(nothing) : throw(AssertionError($msg)))
end

function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing

    @enforce Base.precompile(Tuple{Type{PeriodicGraph},String})
    for T in (Int32,Int64,Int128,BigInt)
        @enforce Base.precompile(Tuple{typeof(PeriodicGraphs.extended_gcd),Vector{T}})
    end

    # Symmetries
    @enforce Base.precompile(Tuple{Type{TrivialIdentitySymmetry}})
    @enforce Base.precompile(Tuple{NoSymmetryGroup, Int})
    @enforce Base.precompile(Tuple{typeof(unique),NoSymmetryGroup})
    @enforce Base.precompile(Tuple{typeof(iterate),NoSymmetryGroup})
    @enforce Base.precompile(Tuple{typeof(eltype),Type{NoSymmetryGroup}})
    @enforce Base.precompile(Tuple{typeof(length),NoSymmetryGroup})
    @enforce Base.precompile(Tuple{typeof(one),NoSymmetryGroup})
    @enforce Base.precompile(Tuple{Type{IncludingIdentity},NoSymmetryGroup})

    # Ring statistics
    SmallIntT = PeriodicGraphs.SmallIntType
    CMBitSet = PeriodicGraphs.ConstMiniBitSet
    for T in (UInt32,UInt64)
        @enforce Base.precompile(Tuple{typeof(PeriodicGraphs._constminibitset),T})
        @enforce Base.precompile(Tuple{Type{CMBitSet{T}}})
        @enforce Base.precompile(Tuple{Type{CMBitSet{T}},Int})
        @enforce Base.precompile(Tuple{Type{CMBitSet{T}},Vector{Int}})
        @enforce Base.precompile(Tuple{typeof(union), CMBitSet{T}, CMBitSet{T}})
        @enforce Base.precompile(Tuple{typeof(intersect), CMBitSet{T}, CMBitSet{T}})
        @enforce Base.precompile(Tuple{typeof(setdiff), CMBitSet{T}, CMBitSet{T}})
        @enforce Base.precompile(Tuple{typeof(symdiff), CMBitSet{T}, CMBitSet{T}})
        @enforce Base.precompile(Tuple{typeof(in), Int, CMBitSet{T}})
        @enforce Base.precompile(Tuple{typeof(eltype), Type{CMBitSet{T}}})
        @enforce Base.precompile(Tuple{typeof(length), CMBitSet{T}})
        @enforce Base.precompile(Tuple{typeof(isempty), CMBitSet{T}})
        @enforce Base.precompile(Tuple{typeof(iterate), CMBitSet{T}, T})
        @enforce Base.precompile(Tuple{typeof(iterate), CMBitSet{T}})
        @enforce Base.precompile(Tuple{typeof(Base.isdone), CMBitSet{T}, T})
        @enforce Base.precompile(Tuple{typeof(Base.isdone), CMBitSet{T}})
        @enforce Base.precompile(Tuple{typeof(PeriodicGraphs.hasonly), CMBitSet{T}, Int})
        @enforce Base.precompile(Tuple{typeof(minimum), CMBitSet{T}})
        @enforce Base.precompile(Tuple{typeof(maximum), CMBitSet{T}})
        @enforce Base.precompile(Tuple{typeof(PeriodicGraphs.unsafe_union!),Ptr{CMBitSet{T}},T})
        @enforce Base.precompile(Tuple{Type{PeriodicGraphs.JunctionNode},Int,SmallIntT,CMBitSet{T}})
    end
    @enforce Base.precompile(Tuple{Type{PeriodicGraphs.JunctionNode}})
    @enforce Base.precompile(Tuple{Type{PeriodicGraphs.JunctionNode},Int})
    @enforce Base.precompile(Tuple{typeof(PeriodicGraphs.unsafe_incr!),Ptr{SmallIntT}})
    @enforce Base.precompile(Tuple{typeof(PeriodicGraphs.symdiff_cycles), Vector{Int}, Vector{Int}})
    @enforce Base.precompile(Tuple{typeof(PeriodicGraphs.gaussian_elimination!), PeriodicGraphs.IterativeGaussianElimination, Vector{Int}})


    for i in 0:3
        # Basic vertex and edge construction
        @enforce Base.precompile(Tuple{Type{PeriodicVertex{i}},Int})
        @enforce Base.precompile(Tuple{Type{PeriodicEdge{i}},Int,Int,NTuple{i, Int}})
        @enforce Base.precompile(Tuple{typeof(PeriodicGraphs.isindirectedge), PeriodicEdge{i}})
        @enforce Base.precompile(Tuple{typeof(Graphs.reverse), PeriodicEdge{i}})
        @enforce Base.precompile(Tuple{Type{PeriodicEdge},Int,Int,NTuple{i, Int}})
        @enforce Base.precompile(Tuple{Type{PeriodicEdge},Int,Int,SVector{i, Int}})
        @enforce Base.precompile(Tuple{Type{PeriodicEdge},Int,PeriodicVertex{i}})
        @enforce Base.precompile(Tuple{Type{PeriodicEdge},Tuple{Int,Int,NTuple{i, Int}}})
        @enforce Base.precompile(Tuple{typeof(Base.string),PeriodicEdge{i}})
        @enforce Base.precompile(Tuple{typeof(Base.string),PeriodicGraph{i}})
        @enforce Base.precompile(Tuple{typeof(Base.string),PeriodicVertex{i}})

        # Graph construction
        @enforce Base.precompile(Tuple{Type{PeriodicGraph{i}}})
        @enforce Base.precompile(Tuple{Type{PeriodicGraph{i}},Vector{PeriodicEdge{i}}})
        @enforce Base.precompile(Tuple{Type{PeriodicGraph{i}},String})
        @enforce Base.precompile(Tuple{Type{PeriodicGraph{i}},SubString{String}})
        @enforce Base.precompile(Tuple{Type{PeriodicGraph},PeriodicGraph{i}})
        @enforce Base.precompile(Tuple{typeof(copy),PeriodicGraph{i}})
        for S in (String, SubString{String})
            @enforce Base.precompile(Tuple{typeof(edges_from_string), KeyString{Int,S}, Val{i}})
            @enforce Base.precompile(Tuple{typeof(_parse), Type{PeriodicGraph{i}}, KeyString{Int,S}})
            @enforce Base.precompile(Tuple{typeof(Base.parse), Type{PeriodicGraph{i}}, S})
        end

        # Equality and comparison
        @enforce Base.precompile(Tuple{typeof(==),PeriodicVertex{i},PeriodicVertex{i}})
        @enforce Base.precompile(Tuple{typeof(<),PeriodicVertex{i},PeriodicVertex{i}})
        @enforce Base.precompile(Tuple{typeof(isless),PeriodicVertex{i},PeriodicVertex{i}})
        @enforce Base.precompile(Tuple{typeof(hash),PeriodicVertex{i},UInt})
        @enforce Base.precompile(Tuple{typeof(==),PeriodicEdge{i},PeriodicEdge{i}})
        @enforce Base.precompile(Tuple{typeof(<),PeriodicEdge{i},PeriodicEdge{i}})
        @enforce Base.precompile(Tuple{typeof(hash),PeriodicEdge{i},UInt})
        @enforce Base.precompile(Tuple{typeof(==),PeriodicGraph{i},PeriodicGraph{i}})
        @enforce Base.precompile(Tuple{typeof(hash),PeriodicGraph{i},UInt})

        # Symmetries
        @enforce Base.precompile(Tuple{Type{NoSymmetryGroup}, PeriodicGraph{i}})

        # Graphs.jl API
        @enforce Base.precompile(Tuple{typeof(Graphs.ne),PeriodicGraph{i}})
        @enforce Base.precompile(Tuple{typeof(Graphs.nv),PeriodicGraph{i}})
        @enforce Base.precompile(Tuple{typeof(Graphs.edges),PeriodicGraph{i}})
        @enforce Base.precompile(Tuple{typeof(Graphs.has_edge),PeriodicGraph{i},Int,Int})
        @enforce Base.precompile(Tuple{typeof(Graphs.has_edge),PeriodicGraph{i},Int,PeriodicVertex{i}})
        @enforce Base.precompile(Tuple{typeof(Graphs.has_edge),PeriodicGraph{i},PeriodicEdge{i}})
        @enforce Base.precompile(Tuple{typeof(find_edges),PeriodicGraph{i},Int,Int})
        @enforce Base.precompile(Tuple{typeof(Graphs.outneighbors),PeriodicGraph{i},Int})
        @enforce Base.precompile(Tuple{typeof(Graphs.inneighbors),PeriodicGraph{i},Int})
        @enforce Base.precompile(Tuple{typeof(Graphs.degree),PeriodicGraph{i},Int})
        @enforce Base.precompile(Tuple{typeof(Graphs.degree),PeriodicGraph{i}})
        @enforce Base.precompile(Tuple{typeof(Graphs.add_vertices!),PeriodicGraph{i},Int})
        @enforce Base.precompile(Tuple{typeof(Graphs.add_vertex!),PeriodicGraph{i}})
        @enforce Base.precompile(Tuple{typeof(Graphs.add_edge!),PeriodicGraph{i},PeriodicEdge{i}})
        @enforce Base.precompile(Tuple{typeof(Graphs.rem_edge!),PeriodicGraph{i},PeriodicEdge{i}})
        @enforce Base.precompile(Tuple{typeof(Graphs.rem_vertices!),PeriodicGraph{i},Vector{Int}})
        @enforce Base.precompile(Tuple{typeof(Graphs.rem_vertex!),PeriodicGraph{i},Int})
        @enforce Base.precompile(Tuple{typeof(Graphs.connected_components),PeriodicGraph{i}})
        @enforce Base.precompile(Tuple{typeof(PeriodicGraphs.vertex_permutation),PeriodicGraph{i},Vector{Int}})
        @enforce Base.precompile(Tuple{typeof(Graphs.induced_subgraph),PeriodicGraph{i},Vector{Int}})
        @enforce Base.precompile(Tuple{typeof(iterate),PeriodicNeighborList{i},Int})
        @enforce Base.precompile(Tuple{typeof(iterate),PeriodicNeighborList{i}})
        @enforce Base.precompile(Tuple{typeof(length),PeriodicNeighborList{i}})
        @enforce Base.precompile(Tuple{typeof(eltype),Type{PeriodicNeighborList{i}}})
        @enforce Base.precompile(Tuple{typeof(Graphs.outneighbors),PeriodicGraph{i},PeriodicVertex{i}})
        @enforce Base.precompile(Tuple{typeof(Graphs.inneighbors),PeriodicGraph{i},PeriodicVertex{i}})
        @enforce Base.precompile(Tuple{typeof(Graphs.degree),PeriodicGraph{i},PeriodicVertex{i}})

        # Edge iteration
        @enforce Base.precompile(Tuple{Type{PeriodicEdgeIter{i}},PeriodicGraph{i}})
        @enforce Base.precompile(Tuple{typeof(iterate),PeriodicGraphs.PeriodicEdgeIter{i},Int})
        @enforce Base.precompile(Tuple{typeof(iterate),PeriodicGraphs.PeriodicEdgeIter{i}})
        @enforce Base.precompile(Tuple{typeof(length),PeriodicGraphs.PeriodicEdgeIter{i}})
        @enforce Base.precompile(Tuple{typeof(eltype),Type{PeriodicGraphs.PeriodicEdgeIter{i}}})
        @enforce Base.precompile(Tuple{typeof(collect),PeriodicGraphs.PeriodicEdgeIter{i}})
        @enforce Base.precompile(Tuple{typeof(isempty),PeriodicGraphs.PeriodicEdgeIter{i}})
    
        # Neighborhood
        @enforce Base.precompile(Tuple{typeof(PeriodicGraphs.hash_position),SVector{i,Int}})
        @enforce Base.precompile(Tuple{typeof(PeriodicGraphs.hash_position),PeriodicVertex{i},Int})
        @enforce Base.precompile(Tuple{typeof(PeriodicGraphs.reverse_hash_position),Int,Val{i}})
        @enforce Base.precompile(Tuple{typeof(PeriodicGraphs.reverse_hash_position),Int,Int,Val{i}})
        @enforce Base.precompile(Tuple{typeof(PeriodicGraphs.reverse_hash_position),Int,PeriodicGraph{i}})
        @enforce Base.precompile(Tuple{typeof(graph_width!),PeriodicGraph{i}})
        @enforce Base.precompile(Tuple{typeof(Graphs._neighborhood),PeriodicGraph{i},Int,Int,Graphs.DefaultDistance,typeof(outneighbors)})
        @enforce Base.precompile(Tuple{typeof(Graphs.neighborhood),PeriodicGraph{i},Int,Int})
        @enforce Base.precompile(Tuple{typeof(coordination_sequence),PeriodicGraph{i},Int,Int})

        # Dimensionality
        for T in (Int32,Int64,Int128,BigInt)
            @enforce Base.precompile(Tuple{typeof(PeriodicGraphs.normal_basis),Vector{SVector{i,T}}})
        end
        @enforce Base.precompile(Tuple{typeof(PeriodicGraphs._dimensionality),PeriodicGraph{i}})
        @enforce Base.precompile(Tuple{typeof(dimensionality),PeriodicGraph{i}})
        for j in 0:3
            @enforce Base.precompile(Tuple{typeof(PeriodicGraphs.change_dimension),Type{PeriodicGraph{i}},PeriodicGraph{j}})
        end

        # Other utils
        @enforce Base.precompile(Tuple{typeof(offset_representatives!),PeriodicGraph{i},Vector{Int}})
        @enforce Base.precompile(Tuple{typeof(offset_representatives!),PeriodicGraph{i},Vector{NTuple{i,Int}}})
        @enforce Base.precompile(Tuple{typeof(swap_axes!),PeriodicGraph{i},Vector{Int}})
        @enforce Base.precompile(Tuple{typeof(swap_axes!),PeriodicGraph{i},NTuple{i,Int}})
        @enforce Base.precompile(Tuple{typeof(cellgraph),PeriodicGraph{i}})
        @enforce Base.precompile(Tuple{typeof(periodiccellgraph),PeriodicGraph{i}})

        # Ring statistics
        dag = Vector{PeriodicGraphs.JunctionNode}
        vertexnums = Vector{PeriodicVertex{i}}
        dist = PeriodicGraphs.DistanceRecord{i}
        rea = PeriodicGraphs.RingsEndingAt{Tuple{dist, vertexnums}}
        nosymm = PeriodicGraphs.NoSymmetryGroup
        pair = Tuple{PeriodicVertex{i},PeriodicVertex{i}}
        @enforce Base.precompile(Tuple{typeof(PeriodicGraphs.prepare_phantomdag!), dag, vertexnums, Dict{PeriodicVertex{i},Int}, PeriodicGraph{i}, Int, BitVector, Nothing})
        @enforce Base.precompile(Tuple{typeof(PeriodicGraphs.prepare_phantomdag!), dag, vertexnums, Dict{PeriodicVertex{i},Int}, PeriodicGraph{i}, Int, BitVector, BitVector})
        @enforce Base.precompile(Tuple{typeof(PeriodicGraphs.handle_phantomdag!), PeriodicGraphs.PhantomJunctionNode{i}, PeriodicGraph{i}, dag, Int, Int})
        @enforce Base.precompile(Tuple{typeof(PeriodicGraphs.arcs_list), PeriodicGraph{i}, Int, Int, Nothing, Nothing})
        @enforce Base.precompile(Tuple{typeof(PeriodicGraphs.arcs_list), PeriodicGraph{i}, Int, Int, BitVector, Nothing})
        @enforce Base.precompile(Tuple{typeof(PeriodicGraphs.arcs_list), PeriodicGraph{i}, Int, Int, Nothing, BitVector})
        @enforce Base.precompile(Tuple{typeof(PeriodicGraphs.arcs_list), PeriodicGraph{i}, Int, Int, BitVector, BitVector})
        @enforce Base.precompile(Tuple{Type{PeriodicGraphs.DistanceRecord}, PeriodicGraph{i}, Int})
        @enforce Base.precompile(Tuple{typeof(PeriodicGraphs.known_distance), dist, Int, Int})
        @enforce Base.precompile(Tuple{typeof(PeriodicGraphs.set_distance!), dist, Int, Int, SmallIntT})
        @enforce Base.precompile(Tuple{typeof(PeriodicGraphs.has_been_seen!), dist, Int, Int})
        @enforce Base.precompile(Tuple{typeof(PeriodicGraphs.bfs_smaller!), dist, Int, Int, Int, SmallIntT})
        @enforce Base.precompile(Tuple{typeof(PeriodicGraphs._reorderinit), Vector{Int}, PeriodicVertex{i}, PeriodicVertex{i}})
        @enforce Base.precompile(Tuple{typeof(PeriodicGraphs.init_distance_record!), dist, Int, Int})
        @enforce Base.precompile(Tuple{typeof(PeriodicGraphs.is_distance_smaller!), dist, PeriodicVertex{i}, PeriodicVertex{i}, SmallIntT})
        @enforce Base.precompile(Tuple{typeof(PeriodicGraphs.next_compatible_arc!), Vector{Int}, Vector{UInt64}, Vector{Int}, dag, Bool, dist, vertexnums})
        @enforce Base.precompile(Tuple{typeof(PeriodicGraphs.initial_compatible_arc!), Vector{Int}, Vector{Int}, dag, Bool, dist, vertexnums})
        @enforce Base.precompile(Tuple{Type{PeriodicGraphs.RingsEndingAt}, dag, Int, Tuple{dist, vertexnums}})
        @enforce Base.precompile(Tuple{typeof(iterate), rea, Tuple{Vector{Int}, Vector{Int}, Vector{Int}, UInt64, UInt64, Int, Int}})
        @enforce Base.precompile(Tuple{typeof(iterate), rea, Nothing})
        @enforce Base.precompile(Tuple{typeof(iterate), rea})
        @enforce Base.precompile(Tuple{typeof(PeriodicGraphs.normalize_cycle!), Vector{Int}, PeriodicGraph{i}})
        @enforce Base.precompile(Tuple{typeof(PeriodicGraphs.rings_around), PeriodicGraph{i}, Int, Int, dist, BitVector})
        @enforce Base.precompile(Tuple{typeof(PeriodicGraphs.no_neighboring_nodes), PeriodicGraph{i}, nosymm})
        @enforce Base.precompile(Tuple{typeof(rings), PeriodicGraph{i}, Int, nosymm, dist})
        @enforce Base.precompile(Tuple{typeof(PeriodicGraphs.cages_around), PeriodicGraph{i}, Int})
        @enforce Base.precompile(Tuple{typeof(PeriodicGraphs.sort_cycles), Vector{Vector{Int}}, Vector{pair}, Dict{pair,Int}, PeriodicGraph{i}, Int})
        @enforce Base.precompile(Tuple{typeof(strong_rings), Vector{Vector{Int}}, PeriodicGraph{i}, Int})
        @enforce Base.precompile(Tuple{typeof(strong_rings), PeriodicGraph{i}, Int, nosymm, dist})
        @enforce Base.precompile(Tuple{typeof(getindex), PeriodicGraphs.RingIncluding{i}, Int})
        @enforce Base.precompile(Tuple{typeof(length), PeriodicGraphs.RingIncluding{i}})
        @enforce Base.precompile(Tuple{typeof(eltype), Type{PeriodicGraphs.RingIncluding{i}}})
        @enforce Base.precompile(Tuple{typeof(iterate), PeriodicGraphs.RingIncluding{i}, Int})
        @enforce Base.precompile(Tuple{typeof(iterate), PeriodicGraphs.RingIncluding{i}})
        @enforce Base.precompile(Tuple{Type{RingAttributions}, PeriodicGraph{i}, Bool, Int, nosymm, dist})
        @enforce Base.precompile(Tuple{typeof(getindex), RingAttributions{i}, Int})
        @enforce Base.precompile(Tuple{typeof(length), RingAttributions{i}})
        @enforce Base.precompile(Tuple{typeof(eltype), Type{RingAttributions{i}}})
        @enforce Base.precompile(Tuple{typeof(iterate), RingAttributions{i}, Int})
        @enforce Base.precompile(Tuple{typeof(iterate), RingAttributions{i}})
    end
end

_precompile_()
