module PeriodicGraphs

using LinearAlgebra
using StaticArrays
using Graphs

import Base: (==), isless, convert, show, showerror, eltype, iterate, zero,
             length, in, ndims, print, cmp, hash


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
