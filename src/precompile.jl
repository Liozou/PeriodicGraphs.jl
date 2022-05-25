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
    @enforce Base.precompile(Tuple{Type{IdentitySymmetry}})
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
    for T in (UInt8, UInt16, UInt32, UInt64, UInt128)
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
        @enforce Base.precompile(Tuple{Type{PeriodicGraphs.JunctionNode{T}},Int,SmallIntT,CMBitSet{T}})
        @enforce Base.precompile(Tuple{Type{PeriodicGraphs.JunctionNode{T}}})
        @enforce Base.precompile(Tuple{Type{PeriodicGraphs.JunctionNode{T}},Int})
    end
    @enforce Base.precompile(Tuple{typeof(PeriodicGraphs.unsafe_incr!),Ptr{SmallIntT}})
    @enforce Base.precompile(Tuple{typeof(PeriodicGraphs.symdiff_cycles), Vector{Int}, Vector{Int}})
    @enforce Base.precompile(Tuple{Type{PeriodicGraphs.IterativeGaussianEliminationLength}, Vector{Int}})
    @enforce Base.precompile(Tuple{Type{PeriodicGraphs.IterativeGaussianElimination{Vector{Int32}}}, Vector{Int}})
    @enforce Base.precompile(Tuple{Type{PeriodicGraphs.IterativeGaussianElimination{Vector{Vector{Int32}}}}, Vector{Int}})
    @enforce Base.precompile(Tuple{Type{PeriodicGraphs.IterativeGaussianElimination}, Vector{Int}})
    @enforce Base.precompile(Tuple{typeof(PeriodicGraphs.gaussian_elimination!), PeriodicGraphs.IterativeGaussianEliminationLength, Vector{Int}})
    @enforce Base.precompile(Tuple{typeof(PeriodicGraphs.gaussian_elimination!), PeriodicGraphs.IterativeGaussianEliminationDecomposition, Vector{Int}})
    @enforce Base.precompile(Tuple{typeof(PeriodicGraphs.gaussian_elimination!), PeriodicGraphs.IterativeGaussianEliminationNone, Vector{Int}})
    @enforce Base.precompile(Tuple{typeof(PeriodicGraphs.retrieve_track!), PeriodicGraphs.IterativeGaussianEliminationDecomposition})

    @enforce Base.precompile(Tuple{typeof(maximum), typeof(length), Vector{Vector{Int}}})

    for i in 0:3
        # Basic vertex and edge construction
        @enforce Base.precompile(Tuple{Type{PeriodicVertex{i}},Int})
        @enforce Base.precompile(Tuple{Type{PeriodicEdge{i}},Int,Int,NTuple{i, Int}})
        @enforce Base.precompile(Tuple{typeof(isdirectedge), PeriodicEdge{i}})
        @enforce Base.precompile(Tuple{typeof(directedge), PeriodicEdge{i}})
        @enforce Base.precompile(Tuple{typeof(directedge), Int, Int, SVector{i,Int}})
        @enforce Base.precompile(Tuple{typeof(directedge), Int, Int, NTuple{i,Int}})
        @enforce Base.precompile(Tuple{typeof(directedge), PeriodicVertex{i}, PeriodicVertex{i}})
        @enforce Base.precompile(Tuple{typeof(Graphs.reverse), PeriodicEdge{i}})
        @enforce Base.precompile(Tuple{Type{PeriodicEdge},Int,Int,NTuple{i, Int}})
        @enforce Base.precompile(Tuple{Type{PeriodicEdge},Int,Int,SVector{i, Int}})
        @enforce Base.precompile(Tuple{Type{PeriodicEdge},Int,PeriodicVertex{i}})
        @enforce Base.precompile(Tuple{Type{PeriodicEdge},Tuple{Int,Int,NTuple{i, Int}}})
        @enforce Base.precompile(Tuple{typeof(string),PeriodicEdge{i}})
        @enforce Base.precompile(Tuple{typeof(string),PeriodicGraph{i}})
        @enforce Base.precompile(Tuple{typeof(string),PeriodicVertex{i}})

        # Graph construction
        @enforce Base.precompile(Tuple{Type{PeriodicGraph{i}}})
        @enforce Base.precompile(Tuple{Type{PeriodicGraph{i}},Vector{PeriodicEdge{i}}})
        @enforce Base.precompile(Tuple{Type{PeriodicGraph{i}},String})
        @enforce Base.precompile(Tuple{Type{PeriodicGraph{i}},SubString{String}})
        @enforce Base.precompile(Tuple{Type{PeriodicGraph},PeriodicGraph{i}})
        @enforce Base.precompile(Tuple{typeof(copy),PeriodicGraph{i}})
        for S in (String, SubString{String})
            @enforce Base.precompile(Tuple{Type{KeyString{Int}}, S})
            @enforce Base.precompile(Tuple{typeof(iterate), KeyString{Int,S}, Int})
            @enforce Base.precompile(Tuple{typeof(popfirst!), KeyString{Int,S}})
            @enforce Base.precompile(Tuple{typeof(edges_from_string), KeyString{Int,S}, Val{i}})
            @enforce Base.precompile(Tuple{typeof(_parse), Type{PeriodicGraph{i}}, KeyString{Int,S}})
            @enforce Base.precompile(Tuple{typeof(parse), Type{PeriodicGraph{i}}, S})
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
        @enforce Base.precompile(Tuple{typeof(iterate),OffsetVertexIterator{i},Int})
        @enforce Base.precompile(Tuple{typeof(iterate),OffsetVertexIterator{i}})
        @enforce Base.precompile(Tuple{typeof(length),OffsetVertexIterator{i}})
        @enforce Base.precompile(Tuple{typeof(eltype),Type{OffsetVertexIterator{i}}})
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
            @enforce Base.precompile(Tuple{Type{PeriodicGraph{i}},PeriodicGraph{j}})
        end

        # Other utils
        @enforce Base.precompile(Tuple{typeof(offset_representatives!),PeriodicGraph{i},Vector{Int}})
        @enforce Base.precompile(Tuple{typeof(offset_representatives!),PeriodicGraph{i},Vector{NTuple{i,Int}}})
        @enforce Base.precompile(Tuple{typeof(swap_axes!),PeriodicGraph{i},Vector{Int}})
        @enforce Base.precompile(Tuple{typeof(swap_axes!),PeriodicGraph{i},NTuple{i,Int}})
        @enforce Base.precompile(Tuple{typeof(truncated_graph),PeriodicGraph{i}})
        @enforce Base.precompile(Tuple{typeof(quotient_graph),PeriodicGraph{i}})

        # Ring statistics
        dag{T} = Vector{PeriodicGraphs.JunctionNode{T}}
        vertexnums = Vector{PeriodicVertex{i}}
        dist = PeriodicGraphs.DistanceRecord{i}
        rea{T} = PeriodicGraphs.RingsEndingAt{T,Tuple{dist, vertexnums}}
        nosymm = PeriodicGraphs.NoSymmetryGroup
        pair = Tuple{PeriodicVertex{i},PeriodicVertex{i}}
        for T in (UInt8, UInt16, UInt32, UInt64, UInt128)
            @enforce Base.precompile(Tuple{typeof(PeriodicGraphs.prepare_phantomdag!), dag{T}, vertexnums, Dict{PeriodicVertex{i},Int}, PeriodicGraph{i}, Int, BitVector, Nothing})
            @enforce Base.precompile(Tuple{typeof(PeriodicGraphs.prepare_phantomdag!), dag{T}, vertexnums, Dict{PeriodicVertex{i},Int}, PeriodicGraph{i}, Int, BitVector, BitVector})
            @enforce Base.precompile(Tuple{typeof(PeriodicGraphs.handle_phantomdag!), PeriodicGraphs.PhantomJunctionNode{i}, PeriodicGraph{i}, dag{T}, Int, Int})
        end
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
        for T in (UInt8, UInt16, UInt32, UInt64, UInt128)
            @enforce Base.precompile(Tuple{typeof(PeriodicGraphs.next_compatible_arc!), Vector{Int}, Vector{UInt64}, Vector{Int}, dag{T}, Bool, dist, vertexnums})
            @enforce Base.precompile(Tuple{typeof(PeriodicGraphs.initial_compatible_arc!), Vector{Int}, Vector{Int}, dag{T}, Bool, dist, vertexnums})
            @enforce Base.precompile(Tuple{Type{PeriodicGraphs.RingsEndingAt}, dag{T}, Int, Tuple{dist, vertexnums}})
            @enforce Base.precompile(Tuple{typeof(iterate), rea{T}, Tuple{Vector{Int}, Vector{Int}, Vector{Int}, UInt64, UInt64, Int, Int}})
            @enforce Base.precompile(Tuple{typeof(iterate), rea{T}, Nothing})
            @enforce Base.precompile(Tuple{typeof(iterate), rea{T}})
        end
        @enforce Base.precompile(Tuple{typeof(PeriodicGraphs.normalize_cycle!), Vector{Int}, Int, Val{i}})
        @enforce Base.precompile(Tuple{typeof(PeriodicGraphs.rings_around), PeriodicGraph{i}, Int, Int, dist, BitVector})
        @enforce Base.precompile(Tuple{typeof(PeriodicGraphs.no_neighboring_nodes), PeriodicGraph{i}, nosymm})
        @enforce Base.precompile(Tuple{typeof(rings), PeriodicGraph{i}, Int, nosymm, dist})
        @enforce Base.precompile(Tuple{typeof(PeriodicGraphs.cages_around), PeriodicGraph{i}, Int})
        @enforce Base.precompile(Tuple{Type{PeriodicGraphs.EdgeDict}, PeriodicGraph{i}})
        @enforce Base.precompile(Tuple{typeof(get!), PeriodicGraphs.EdgeDict{i}, PeriodicGraphs.VertexPair{i}})
        @enforce Base.precompile(Tuple{typeof(getindex), PeriodicGraphs.EdgeDict{i}, Int})
        @enforce Base.precompile(Tuple{typeof(PeriodicGraphs.sort_cycles), PeriodicGraph{i}, Vector{Vector{Int}}, Int, Tuple{Vector{Tuple{PeriodicVertex{i},PeriodicVertex{i}}}, Dict{Tuple{PeriodicVertex{i},PeriodicVertex{i}},Int}}})
        @enforce Base.precompile(Tuple{typeof(PeriodicGraphs.sort_cycles), PeriodicGraph{i}, Vector{Vector{Int}}, Int})
        @enforce Base.precompile(Tuple{typeof(PeriodicGraphs.sort_cycles), PeriodicGraph{i}, Vector{Vector{Int}}})
        @enforce Base.precompile(Tuple{typeof(strong_rings), Vector{Vector{Int}}, PeriodicGraph{i}, Int, nosymm})
        @enforce Base.precompile(Tuple{typeof(strong_rings), Vector{Vector{Int}}, PeriodicGraph{i}, Int})
        @enforce Base.precompile(Tuple{typeof(strong_rings), Vector{Vector{Int}}, PeriodicGraph{i}})
        @enforce Base.precompile(Tuple{typeof(strong_rings), PeriodicGraph{i}, Int, nosymm, dist})
        @enforce Base.precompile(Tuple{typeof(strong_rings), PeriodicGraph{i}, Int, nosymm})
        @enforce Base.precompile(Tuple{typeof(strong_rings), PeriodicGraph{i}, Int})
        @enforce Base.precompile(Tuple{typeof(strong_rings), PeriodicGraph{i}})
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
