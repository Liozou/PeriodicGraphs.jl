function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing

    Base.precompile(Tuple{Type{PeriodicGraph},String})

    for i in 1:3
        Base.precompile(Tuple{Type{PeriodicEdge{i}},Int,Int,NTuple{i, Int}})
        Base.precompile(Tuple{Type{PeriodicVertex{i}},Int})
        
        Base.precompile(Tuple{Type{PeriodicEdge},Int,Int,NTuple{i, Int}})
        Base.precompile(Tuple{Type{PeriodicEdge},Int,Int,SVector{i, Int}})
        Base.precompile(Tuple{Type{PeriodicEdge},Int,PeriodicVertex{i}})

        Base.precompile(Tuple{Type{PeriodicGraph{i}}})
        Base.precompile(Tuple{Type{PeriodicGraph{i}},Vector{PeriodicEdge{i}}})
        Base.precompile(Tuple{Type{PeriodicGraph{i}},String})

        Base.precompile(Tuple{typeof(==),PeriodicVertex{i},PeriodicVertex{i}})
        Base.precompile(Tuple{typeof(<),PeriodicVertex{i},PeriodicVertex{i}})
        Base.precompile(Tuple{typeof(==),PeriodicEdge{i},PeriodicEdge{i}})
        Base.precompile(Tuple{typeof(<),PeriodicEdge{i},PeriodicEdge{i}})

        Base.precompile(Tuple{typeof(LightGraphs._neighborhood),PeriodicGraph{i},Int,Int,LightGraphs.DefaultDistance,Function})
        Base.precompile(Tuple{typeof(add_edge!),PeriodicGraph{i},PeriodicEdge{i}})
        Base.precompile(Tuple{typeof(add_vertices!),PeriodicGraph{i},Int})
        Base.precompile(Tuple{typeof(rem_edge!),PeriodicGraph{i},PeriodicEdge{i}})
        Base.precompile(Tuple{typeof(rem_vertex!),PeriodicGraph{i},Int})
        Base.precompile(Tuple{typeof(connected_components),PeriodicGraph{i}})
        Base.precompile(Tuple{typeof(coordination_sequence),PeriodicGraph{i},Int,Int})
        Base.precompile(Tuple{typeof(dimensionality),PeriodicGraph{i}})
        Base.precompile(Tuple{typeof(has_edge),PeriodicGraph{i},Int,Int})
        Base.precompile(Tuple{typeof(induced_subgraph),PeriodicGraph{i},Vector{Int}})
        Base.precompile(Tuple{typeof(isempty),PeriodicGraphs.PeriodicEdgeIter{i}})
        Base.precompile(Tuple{typeof(neighborhood),PeriodicGraph{i},Int,Int})
        Base.precompile(Tuple{typeof(offset_representatives!),PeriodicGraph{i},Vector{Int}})
        Base.precompile(Tuple{typeof(reverse),PeriodicEdge{i}})
        Base.precompile(Tuple{typeof(string),PeriodicEdge{i}})
        Base.precompile(Tuple{typeof(string),PeriodicGraph{i}})
        Base.precompile(Tuple{typeof(string),PeriodicVertex{i}})
        Base.precompile(Tuple{typeof(swap_axes!),PeriodicGraph{i},Vector{Int}})
        Base.precompile(Tuple{typeof(vertex_permutation),PeriodicGraph{i},Vector{Any}})
    end
end
