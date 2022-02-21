macro enforce(expr) # strong @assert
    msg = string(expr)
    return :($(esc(expr)) ? $(nothing) : throw(AssertionError($msg)))
end

function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing

    @enforce Base.precompile(Tuple{Type{PeriodicGraph},String})
    @enforce Base.precompile(Tuple{typeof(PeriodicGraphs.extended_gcd),Vector{Int}})

    for i in 0:3
        @enforce Base.precompile(Tuple{Type{PeriodicEdge{i}},Int,Int,NTuple{i, Int}})
        @enforce Base.precompile(Tuple{Type{PeriodicVertex{i}},Int})

        @enforce Base.precompile(Tuple{Type{PeriodicEdge},Int,Int,NTuple{i, Int}})
        @enforce Base.precompile(Tuple{Type{PeriodicEdge},Int,Int,SVector{i, Int}})
        @enforce Base.precompile(Tuple{Type{PeriodicEdge},Int,PeriodicVertex{i}})
        @enforce Base.precompile(Tuple{Type{PeriodicEdge},Tuple{Int,Int,NTuple{i, Int}}})

        @enforce Base.precompile(Tuple{Type{PeriodicGraph{i}}})
        @enforce Base.precompile(Tuple{Type{PeriodicGraph{i}},Vector{PeriodicEdge{i}}})
        @enforce Base.precompile(Tuple{Type{PeriodicGraph{i}},String})
        @enforce Base.precompile(Tuple{Type{PeriodicGraph},PeriodicGraph{i}})

        @enforce Base.precompile(Tuple{typeof(==),PeriodicVertex{i},PeriodicVertex{i}})
        @enforce Base.precompile(Tuple{typeof(<),PeriodicVertex{i},PeriodicVertex{i}})
        @enforce Base.precompile(Tuple{typeof(isless),PeriodicVertex{i},PeriodicVertex{i}})
        @enforce Base.precompile(Tuple{typeof(hash),PeriodicVertex{i},UInt})
        @enforce Base.precompile(Tuple{typeof(==),PeriodicEdge{i},PeriodicEdge{i}})
        @enforce Base.precompile(Tuple{typeof(<),PeriodicEdge{i},PeriodicEdge{i}})
        @enforce Base.precompile(Tuple{typeof(hash),PeriodicEdge{i},UInt})
        @enforce Base.precompile(Tuple{typeof(==),PeriodicGraph{i},PeriodicGraph{i}})
        @enforce Base.precompile(Tuple{typeof(hash),PeriodicGraph{i},UInt})

        @enforce Base.precompile(Tuple{typeof(Graphs._neighborhood),PeriodicGraph{i},Int,Int,Graphs.DefaultDistance,Function})
        @enforce Base.precompile(Tuple{typeof(Graphs.add_edge!),PeriodicGraph{i},PeriodicEdge{i}})
        @enforce Base.precompile(Tuple{typeof(add_vertex!),PeriodicGraph{i}})
        @enforce Base.precompile(Tuple{typeof(Graphs.add_vertices!),PeriodicGraph{i},Int})
        @enforce Base.precompile(Tuple{typeof(Graphs.rem_edge!),PeriodicGraph{i},PeriodicEdge{i}})
        @enforce Base.precompile(Tuple{typeof(Graphs.rem_vertex!),PeriodicGraph{i},Int})
        @enforce Base.precompile(Tuple{typeof(collect),PeriodicGraphs.PeriodicEdgeIter{i}})
        @enforce Base.precompile(Tuple{typeof(Graphs.connected_components),PeriodicGraph{i}})
        @enforce Base.precompile(Tuple{typeof(coordination_sequence),PeriodicGraph{i},Int,Int})
        @enforce Base.precompile(Tuple{typeof(dimensionality),PeriodicGraph{i}})
        @enforce Base.precompile(Tuple{typeof(find_edges),PeriodicGraph{i},Int,Int})
        @enforce Base.precompile(Tuple{typeof(Graphs.has_edge),PeriodicGraph{i},Int,Int})
        @enforce Base.precompile(Tuple{typeof(Graphs.has_edge),PeriodicGraph{i},PeriodicEdge{i}})
        @enforce Base.precompile(Tuple{typeof(Graphs.induced_subgraph),PeriodicGraph{i},Vector{Int}})
        @enforce Base.precompile(Tuple{typeof(isempty),PeriodicGraphs.PeriodicEdgeIter{i}})
        @enforce Base.precompile(Tuple{typeof(iterate),PeriodicEdgeIter{i}})
        @enforce Base.precompile(Tuple{typeof(Graphs.neighborhood),PeriodicGraph{i},Int,Int})
        @enforce Base.precompile(Tuple{typeof(offset_representatives!),PeriodicGraph{i},Vector{Any}})
        @enforce Base.precompile(Tuple{typeof(offset_representatives!),PeriodicGraph{i},Vector{Int}})
        @enforce Base.precompile(Tuple{typeof(offset_representatives!),PeriodicGraph{i},Vector{NTuple{i,Int}}})
        @enforce Base.precompile(Tuple{typeof(Graphs.reverse),PeriodicEdge{i}})
        @enforce Base.precompile(Tuple{typeof(Base.string),PeriodicEdge{i}})
        @enforce Base.precompile(Tuple{typeof(Base.string),PeriodicGraph{i}})
        @enforce Base.precompile(Tuple{typeof(Base.string),PeriodicVertex{i}})
        @enforce Base.precompile(Tuple{typeof(swap_axes!),PeriodicGraph{i},Vector{Int}})
        @enforce Base.precompile(Tuple{typeof(swap_axes!),PeriodicGraph{i},NTuple{i,Int}})
        @enforce Base.precompile(Tuple{typeof(vertex_permutation),PeriodicGraph{i},Vector{Any}})
        @enforce Base.precompile(Tuple{typeof(periodiccellgraph),PeriodicGraph{i}})
        @enforce Base.precompile(Tuple{typeof(PeriodicGraphs.hash_position),SVector{i,Int}})

        for j in 0:3
            @enforce Base.precompile(Tuple{typeof(PeriodicGraphs.change_dimension),Type{PeriodicGraph{i}},PeriodicGraph{j}})
        end
    end
end

_precompile_()
