function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing

    @assert Base.precompile(Tuple{Type{PeriodicGraph},String})

    for i in 1:3
        @assert Base.precompile(Tuple{Type{PeriodicEdge{i}},Int,Int,NTuple{i, Int}})
        @assert Base.precompile(Tuple{Type{PeriodicVertex{i}},Int})

        @assert Base.precompile(Tuple{Type{PeriodicEdge},Int,Int,NTuple{i, Int}})
        @assert Base.precompile(Tuple{Type{PeriodicEdge},Int,Int,SVector{i, Int}})
        @assert Base.precompile(Tuple{Type{PeriodicEdge},Int,PeriodicVertex{i}})

        @assert Base.precompile(Tuple{Type{PeriodicGraph{i}}})
        @assert Base.precompile(Tuple{Type{PeriodicGraph{i}},Vector{PeriodicEdge{i}}})
        @assert Base.precompile(Tuple{Type{PeriodicGraph{i}},String})

        @assert Base.precompile(Tuple{typeof(==),PeriodicVertex{i},PeriodicVertex{i}})
        @assert Base.precompile(Tuple{typeof(<),PeriodicVertex{i},PeriodicVertex{i}})
        @assert Base.precompile(Tuple{typeof(==),PeriodicEdge{i},PeriodicEdge{i}})
        @assert Base.precompile(Tuple{typeof(<),PeriodicEdge{i},PeriodicEdge{i}})

        @assert Base.precompile(Tuple{typeof(LightGraphs._neighborhood),PeriodicGraph{i},Int,Int,LightGraphs.DefaultDistance,Function})
        @assert Base.precompile(Tuple{typeof(LightGraphs.add_edge!),PeriodicGraph{i},PeriodicEdge{i}})
        @assert Base.precompile(Tuple{typeof(LightGraphs.add_vertices!),PeriodicGraph{i},Int})
        @assert Base.precompile(Tuple{typeof(LightGraphs.rem_edge!),PeriodicGraph{i},PeriodicEdge{i}})
        @assert Base.precompile(Tuple{typeof(LightGraphs.rem_vertex!),PeriodicGraph{i},Int})
        @assert Base.precompile(Tuple{typeof(collect),PeriodicGraphs.PeriodicEdgeIter{i}})
        @assert Base.precompile(Tuple{typeof(LightGraphs.connected_components),PeriodicGraph{i}})
        @assert Base.precompile(Tuple{typeof(coordination_sequence),PeriodicGraph{i},Int,Int})
        @assert Base.precompile(Tuple{typeof(dimensionality),PeriodicGraph{i}})
        @assert Base.precompile(Tuple{typeof(find_edges),PeriodicGraph{i},Int,Int})
        @assert Base.precompile(Tuple{typeof(LightGraphs.has_edge),PeriodicGraph{i},Int,Int})
        @assert Base.precompile(Tuple{typeof(LightGraphs.induced_subgraph),PeriodicGraph{i},Vector{Int}})
        @assert Base.precompile(Tuple{typeof(Base.isempty),PeriodicGraphs.PeriodicEdgeIter{i}})
        @assert Base.precompile(Tuple{typeof(LightGraphs.neighborhood),PeriodicGraph{i},Int,Int})
        @assert Base.precompile(Tuple{typeof(offset_representatives!),PeriodicGraph{i},Vector{Any}})
        @assert Base.precompile(Tuple{typeof(offset_representatives!),PeriodicGraph{i},Vector{Int}})
        @assert Base.precompile(Tuple{typeof(offset_representatives!),PeriodicGraph{i},Vector{NTuple{i,Int}}})
        @assert Base.precompile(Tuple{typeof(LightGraphs.reverse),PeriodicEdge{i}})
        @assert Base.precompile(Tuple{typeof(Base.string),PeriodicEdge{i}})
        @assert Base.precompile(Tuple{typeof(Base.string),PeriodicGraph{i}})
        @assert Base.precompile(Tuple{typeof(Base.string),PeriodicVertex{i}})
        @assert Base.precompile(Tuple{typeof(swap_axes!),PeriodicGraph{i},Vector{Int}})
        @assert Base.precompile(Tuple{typeof(swap_axes!),PeriodicGraph{i},NTuple{i,Int}})
        @assert Base.precompile(Tuple{typeof(vertex_permutation),PeriodicGraph{i},Vector{Any}})

        for j in 1:3
            @assert Base.precompile(Tuple{typeof(PeriodicGraphs.change_dimension),Type{PeriodicGraph{i}},PeriodicGraph{j}})
        end
    end
end

_precompile_()
