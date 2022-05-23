# Basic types and definitions

## Introduction

A [`PeriodicGraph{N}`](@ref) `g` is the representation of an `N`-periodic graph.

Each vertex has a unique representative, indexed from 1 to `n = nv(g)`.
Each vertex `x` of the graph is represented by a [`PeriodicVertex{N}`](@ref) containing
both the representative `v` and the offset `o` between the unit cell containing the vertex
and a reference unit cell, accessible through the syntax `v, o = x`.
The offset is a `N`-uplet of integers, stored as a `SVector{N,Int}`.

Each edge `e` is represented by a [`PeriodicEdge{N}`](@ref) defined by its source vertex
inside the reference unit cell and of representative `src`, and its destination vertex
`dst`, a `PeriodicVertex{N}`, which can be accessed through the syntax `src, dst = e`.

For convenience, aliases are exported for 1D, 2D and 3D (`N = 1`, `N = 2` and `N = 3`)
under the names `PeriodicGraph1D`, `PeriodicEdge2D`, `PeriodicVertex3D`, etc.

## `PeriodicVertex`

A `PeriodicVertex` can be built like so:

```jldoctest
julia> PeriodicVertex1D(5, (-1,))
PeriodicVertex1D(5, (-1))

julia> PeriodicVertex{4}(2) # shorthand for the vertices of the reference cell
PeriodicVertex{4}(2, (0,0,0,0))
```

## `PeriodicEdge`

A `PeriodicEdge`, can be defined from `src` and `dst` or, equivalently, be the identifiers
of both source and destination vertices, and the cell offset between source and destination.
For example, in 3D, the edge between vertex `1` and vertex `4` in cell `(0, 1, 0)` is
`PeriodicEdge3D(1, PeriodicVertex3D(4, (0,1,0)))`, or, equivalently,
`PeriodicEdge3D(1, 4, (0,1,0))`. Since `PeriodicEdge(u, v, ofs)` and
`PeriodicEdge(v, u, .-ofs)` represent the same edge, one of them is called the *direct*
edge, when it has either `u < v` or `u == v && ofs > zero(ofs)`, and the other is the
*indirect* edge. Functions [`isdirectedge`](@ref) and [`directedge`](@ref) are used for
this.

More examples:

```jldoctest
julia> e = PeriodicEdge(3, PeriodicVertex(2, (0,1,0,0))) # dimensions are inferred from the input
PeriodicEdge{4}(3, 2, (0,1,0,0))

julia> isdirectedge(e)
false

julia> directedge(e)
PeriodicEdge{4}(2, 3, (0,-1,0,0))

julia> src, (dstv, ofs) = PeriodicEdge3D(5, 6, (1,0,2));

julia> src
5

julia> dst
6

julia> ofs
3-element StaticArrays.SVector{3,$Int} with indices SOneTo(3):
 1
 0
 2

julia> PeriodicEdge2D(5, 5, (0,0))
ERROR: LoopException: a loop from vertex 5 to itself in the same unit cell is a forbidden edges. Maybe the offset is wrong?
```

Note that loops (that is, edges of the form `(u, u, (0,0,0,...,0))`) are forbidden and will
throw a [`LoopException`](@ref) if created. To bypass this check, use the unexported
[`unsafe_edge`](@ref) function.

## `PeriodicGraph`

Finally, `N`-periodic graphs, represented by the type `PeriodicGraph{N}`, are defined by
the number of vertices in the reference cell and the set of edges starting from the
reference cell. When the number of vertices is simply the highest number appearing in the
list of edges, it can be omitted.
Periodic graphs can be built through several methods:

```jldoctest
julia> PeriodicGraph{5}() # create the empty graph
PeriodicGraph{5}(0, PeriodicEdge{5}[])

julia> PeriodicGraph1D(4) # create a graph with 4 vertices but no edge
PeriodicGraph1D(4, PeriodicEdge1D[])

julia> PeriodicGraph(2, PeriodicEdge{4}[(1, 1, (0,0,1,1)), (1, 1, (0,1,0,-1))]) # the dimension can be inferred
PeriodicGraph{4}(2, PeriodicEdge{4}[(1, 1, (0,0,1,1)), (1, 1, (0,1,0,-1))])

julia> PeriodicGraph3D(PeriodicEdge3D[(1, 3, (0,1,0)), (2, 2, (0,0,-1)), (1, 2, (1,0,0)), (2, 3, (0,1,1))])
PeriodicGraph3D(3, PeriodicEdge3D[(1, 2, (1,0,0)), (1, 3, (0,1,0)), (2, 2, (0,0,1)), (2, 3, (0,1,1))])

julia> parse(PeriodicGraph3D, "3   1 2  1 0 0   1 3  0 1 0   2 2  0 0 1   2 3  0 1 1") # compact representation of the previous graph
PeriodicGraph3D(3, PeriodicEdge3D[(1, 2, (1,0,0)), (1, 3, (0,1,0)), (2, 2, (0,0,1)), (2, 3, (0,1,1))])

julia> string(ans) # to obtain the compact representation from the graph
"3 1 2 1 0 0 1 3 0 1 0 2 2 0 0 1 2 3 0 1 1"

julia> string(PeriodicGraph("2   1 2 0 0  2 3 0 0  3 1 0 0  3 1 1 0  2 1 0 1  2 3 0 1"))
"2 1 2 0 -1 1 2 0 0 1 3 -1 0 1 3 0 0 2 3 0 0 2 3 0 1"
```

## API

### Type definitions

```@docs
PeriodicVertex
PeriodicEdge
PeriodicGraph
```

### Other essential functions

```@docs
PeriodicGraphs.unsafe_edge
PeriodicGraphs.LoopException
isdirectedge
directedge
```
