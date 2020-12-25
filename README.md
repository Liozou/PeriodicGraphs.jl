# PeriodicGraphs

[![Build Status](https://travis-ci.com/Liozou/PeriodicGraphs.jl.svg?branch=master)](https://travis-ci.com/Liozou/PeriodicGraphs.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/Liozou/PeriodicGraphs.jl?svg=true)](https://ci.appveyor.com/project/Liozou/PeriodicGraphs-jl)
[![codecov](https://codecov.io/gh/Liozou/PeriodicGraphs.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/Liozou/PeriodicGraphs.jl)
<!-- [![Aqua QA](https://img.shields.io/badge/Aqua.jl-%F0%9F%8C%A2-aqua.svg)](https://github.com/tkf/Aqua.jl) -->

This module allows to manipulate `N`-dimensional periodic graphs, using the new type `PeriodicGraph{N}`. This is a subtype of `AbstractGraph{Int}` from [LightGraphs.jl](https://github.com/JuliaGraphs/LightGraphs.jl/) and it extends its API.

The main difference with a simple graph is the notion of offset. Each vertex, of type `PeriodicVertex{N}` is uniquely defined using a numeric identifier (a positive integer, like for any simple graph) associated with the offset of the cell in which the designated vertex is, compared to a fixed reference cell. For instance, all vertices in the reference cell have a zero offset, and may be built like so:
```julia
julia> PeriodicVertex{1}(1, (0))
PeriodicVertex{1}(1, (0))

julia> PeriodicVertex{4}(2) # shorthand for the vertices of the reference cell
PeriodicVertex{4}(2, (0,0,0,0))
```

For example, in 2D (`N = 2`), the vertex which is the representative of vertex `3` in the cell with offset `(0, 1)` (that is, on the same plane as the reference cell along the x-axis, and just next to the reference cell positively along the y-axis) is defined as `PeriodicVertex{2}(3, (0,1))`.

An edge, of type `PeriodicEdge{N}`, is defined by its representative starting from the reference cell. Hence, it is uniquely defined by the identifier of the source vertex (in the reference cell) and the destination vertex. It can equivalently be defined by the identifiers of both source and destination vertices, and the offset of the destination, the source being in the reference cell. For example, in 3D, the edge between vertex `1` and vertex `4` in cell `(0, 1, 0)` is `PeriodicEdge{3}(1, PeriodicVertex{3}(4, (0,1,0)))`, or, equivalently, `PeriodicEdge{3}(1, 4, (0,1,0))`. More examples:
```julia
julia> PeriodicEdge{4}(2, PeriodicVertex{4}(3, (0,1,0,0)))
PeriodicEdge{4}(2, 3, (0,1,0,0))

julia> LightGraphs.src(PeriodicEdge{3}(5, 6, (1,0,2)))
5

julia> LightGraphs.dst(PeriodicEdge{3}(5, 6, (1,0,2)))
6

julia> PeriodicGraphs.ofs(PeriodicEdge{3}(5, 6, (1,0,2)))
3-element StaticArrays.SVector{3,Int64} with indices SOneTo(3):
 1
 0
 2

julia> PeriodicEdge{2}(5, 5, (0,0))
ERROR: LoopException: a loop from vertex 5 to itself in the same unit cell is a forbidden edges. Maybe the offset is wrong?
```

Note that loops (that is, edges of the form `(u, u, (0,0,0,...,0))`) are forbidden in the current version of PeriodicGraphs.jl.

Finally, `N`-periodic graphs, represented by the type `PeriodicGraph{N}`, are defined by the number of vertices in the reference cell and the set of edges starting from the reference cell. When the number of vertices is simply the highest number appearing in the list of edges, it can be omitted. Periodic graphs can be built with several methods:
```julia
julia> PeriodicGraph{5}() # create the empty graph
PeriodicGraph{5}(0, PeriodicEdge{5}[])

julia> PeriodicGraph{1}(4) # create a graph with 4 vertices but no edge
PeriodicGraph{1}(4, PeriodicEdge{1}[])

julia> PeriodicGraph(2, PeriodicEdge{4}[(1, 1, (0,0,1,1)), (1, 1, (0,1,0,-1))]) # the dimension can be inferred from the list of edges
PeriodicGraph{4}(2, PeriodicEdge{4}[(1, 1, (0,0,1,1)), (1, 1, (0,1,0,-1))])

julia> PeriodicGraph{3}(PeriodicEdge{3}[(1, 3, (0,1,0)), (2, 2, (0,0,-1)), (1, 2, (1,0,0)), (2, 3, (0,1,1))])
PeriodicGraph3D(3, PeriodicEdge3D[(1, 2, (1,0,0)), (1, 3, (0,1,0)), (2, 2, (0,0,1)), (2, 3, (0,1,1))])

julia> PeriodicGraph("3   1 2  1 0 0   1 3  0 1 0   2 2  0 0 1   2 3  0 1 1") # compact representation of the previous graph
PeriodicGraph3D(3, PeriodicEdge3D[(1, 2, (1,0,0)), (1, 3, (0,1,0)), (2, 2, (0,0,1)), (2, 3, (0,1,1))])

julia> string(ans) # to obtain the compact representation from the graph
"3 1 2 1 0 0 1 3 0 1 0 2 2 0 0 1 2 3 0 1 1"
```

Note that `PeriodicGraph`s are undirected: for this reason, any edge in the graph of the form `(u, v, (x1, x2, ..., xN))` has a reverse edge in the graph of the form `(v, u, (-x1, -x2, ..., -xN))`. For this reason, calling `LightGraphs.edges` on a `PeriodicGraph` yields an iterator to the edges of the graph that will only show the edges under the canonical form `(u, v, ofs)` with either `u < v` or `u == v && ofs > zero(ofs)`.

Neighbors can be obtained using LightGraphs's `neighbors` function, with the same API. It will return the list of neighbors of the given vertex, assuming the vertex is in the reference cell. Each given is given as a `PeriodicVertex`, which contains the offset compared to the reference cell. For instance:
```julia
julia> g = PeriodicGraph("3  1 2 1 0 0  1 3 0 1 0  2 2 0 0 1  2 3 0 1 1")
PeriodicGraph3D(3, PeriodicEdge3D[(1, 2, (1,0,0)), (1, 3, (0,1,0)), (2, 2, (0,0,1)), (2, 3, (0,1,1))])

julia> neighbors(g, 2)
4-element Vector{PeriodicVertex3D}:
 (1, (-1,0,0))
 (2, (0,0,-1))
 (2, (0,0,1))
 (3, (0,1,1))
```

For convenience, aliases are exported for 1D, 2D and 3D (`N = 1`, `N = 2` and `N = 3`) under the names `PeriodicGraph1D`, `PeriodicEdge2D`, `PeriodicVertex3D`, etc.
