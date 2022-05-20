module PeriodicGraphs

using LinearAlgebra
using StaticArrays
using Graphs

import Base: (==), isless, convert, show, showerror, eltype, iterate, zero,
             length, in, ndims, print, cmp, hash


@static if VERSION < v"1.7.0-"
    # copied from Base without the bound checking
    function keepat!(a::Vector, inds)
        i = firstindex(a)
        for k in inds
            if i != k
                @inbounds a[i] = a[k]
            end
            i = nextind(a, i)
        end
        deleteat!(a, i:lastindex(a))
        return a
    end
end

include("definitions/vertices.jl")
include("definitions/edges.jl")
include("definitions/graphs.jl")
include("definitions/symmetries.jl")

include("utils/io.jl")
include("utils/graphsAPI.jl")
include("utils/edgeiter.jl")

include("algorithms/neighborhood.jl")
include("algorithms/dimensionality.jl")
include("algorithms/other.jl")
include("algorithms/rings.jl")

include("precompile.jl")

end
