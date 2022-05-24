# Dimensionality

A periodic graph `g` of type `PeriodicGraph{N}` has *dimension* `N`, which means that its
unit cell is repeated across `N` independent axes.
However, this number may be larger than the actual number of independent axes necessary
to describe the repeating unit of the graph, which we will call the *dimensionality*.
For example, consider an infinite set of identical 2-periodic graphs stacked across a
third dimension: the resulting graph is of dimension 3, but its dimensionality is 2 (or
below) because all the topological information is stored in the 2-periodic subgraph.

The connected components of `g` can be separated and sorted by dimensionality using the
[`dimensionality`](@ref) function:

```@docs
dimensionality
```

To transpose a graph from one dimension to another use the [`change_dimension`](@ref)
function. This can be useful to manipulate a graph of dimensionality `D` as a graph of
actual dimension `D`, which often reduces computational costs.

```@docs
change_dimension
```

For example, the following function extracts the list of 1-dimensional components from a
given `PeriodicGraph`:

```jldoctest extract1D; setup=:(using PeriodicGraphs, Graphs)
julia> function extract_1D_components(g::PeriodicGraph{D}) where D
           d = dimensionality(g)
           components1D = get(d, 1, Vector{Int}[])
           return [change_dimension(PeriodicGraph1D, g[l]) for l in components1D]
       end
extract_1D_components (generic function with 1 method)
```

Let's test it on the following 2-periodic graph, which has one 0D component
(vertices 4 and 5), two 1D components (vertex 3 alone and vertices 6 and 7 together) and
one 2D component (vertices 1 and 2):

```jldoctest extract1D
julia> g = PeriodicGraph2D("2   1 2 0 0  2 1 1 0  2 1 0 1
                                3 3 1 0
                                4 5 0 0
                                6 7 1 0  7 6 1 0
                           ")
PeriodicGraph2D(7, PeriodicEdge2D[(1, 2, (-1,0)), (1, 2, (0,-1)), (1, 2, (0,0)), (3, 3, (1,0)), (4, 5, (0,0)), (6, 7, (-1,0)), (6, 7, (1,0))])

julia> extract_1D_components(g)
2-element Vector{PeriodicGraph1D}:
 PeriodicGraph1D(1, PeriodicEdge1D[(1, 1, (1,))])
 PeriodicGraph1D(2, PeriodicEdge1D[(1, 2, (-1,)), (1, 2, (1,))])
```

The first subgraph corresponds to `g[[3]]` and the second to `g[[6,7]]`, both converted to
`PeriodicGraph1D`.
