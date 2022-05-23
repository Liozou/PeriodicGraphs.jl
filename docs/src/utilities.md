# Utilities

## Hashing

It is sometimes convenient to be able to associate to each vertex of periodic graph `g` an
integer hash, such that the hashes of vertices in the reference unit cell are between
1 and `nv(g)`, the hashes of the vertices in the unit cells around the reference are next,
then the vertices of the unit cells around those, etc. To do so, `PeriodicGraphs.jl` export
the `hash_position` function as follows:

```@docs
hash_position
```

The reciproque function is also exported:

```@docs
reverse_hash_position
```

Both functions are optimized for dimensions 1, 2 and 3, especially for unit cells not too
far from the reference.

## Isomorphic transformations

Several functions transform a periodic graph into another isomorphic to the input, by
renumbering the vertices ([`vertex_permutation`](@ref PeriodicGraphs.vertex_permutation))
or the axes ([`swap_axes!`](@ref)), or by offsetting the chosen representatives for each
vertex ([`offset_representatives!`](@ref)).

```@docs
PeriodicGraphs.vertex_permutation
swap_axes!
offset_representatives!
```

## Reduction to aperiodic graph

Any `PeriodicGraph` can be naturally reduced to an aperiodic graph by removing all offsets
from the edges and either keeping ([`quotient_graph`](@ref)) or removing
([`truncated_graph](@ref)) edges crossing from one unit cell to another.

```@docs
quotient_graph
truncated_graph
```

## Unclassified other utilities

These other convenience functions may be useful for the manipulation of `PeriodicGraph`s:

```@docs
find_edges
coordination_sequence
```
