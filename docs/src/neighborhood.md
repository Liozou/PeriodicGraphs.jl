# Neighborhoods and graph traversals

## Manual

The `AbstractGraph` API for exploring neighbors and traversing graphs has to be somewhat
adapted to account for the specific aspects of `PeriodicGraphs`.
Here, we present a list of notes relative to these aspects, where `g` is a
`PeriodicGraph{N}`:

- `neighbors(g, i)` where `i::Integer` is the list of neighbors of `PeriodicVertex{N}(i)`.
  Performance-wise, this operation is a simple access on the underlying structure of `g` so
  it is very fast, but the returned `Vector{PeriodicVertex{N}}` should not be modified.
- `neighbors(g, x)` where `x::PeriodicVertex{N}` is an iterator over the neighbors of `x`
  but not an array. The iterator object is a [`PeriodicGraphs.OffsetVertexIterator`](@ref).
- `edges(g)` is an iterator over the *direct* edges of `g`. Every edge representative will
  thus be visited exactly once. It is invalidated by any change to `g`, so `g` should not
  be modified while iterating over `edges(g)`. The iterator object is a
  [`PeriodicGraphs.PeriodicEdgeIter`](@ref).
- `Graphs._neighborhood` has a specialization for `PeriodicGraph` if required.

A typical BFS algorithm on `g` can be implemented like so:

```julia
function bfs(g::PeriodicGraph{N}, starting_vertex, depth) where N
    visited = Set{PeriodicVertex{N}}(starting_vertex)
    Q = Tuple{Int,PeriodicVertex{N}}[(0, starting_vertex)]
    for (distance, u) in Q
        distance > depth && break
        # do stuff with u
        for x in neighbors(g, u)
            x ∈ visited && continue
            push!(x, visited)
            # do stuff with x
            push!(Q, (distance+1, x))
        end
        # do more stuff
    end
    # return something
end
```

In some cases however, the `Set`-interface can end up being the bottleneck. In this kind of
situation, it *may* be better to:

1. replace the initialization of `visited` by

   ```julia
   width = PeriodicGraphs.graph_width!(g)
   seen_size = nv(g)*(2*(1 + fld(depth-1, width)) + 1)^N
   visited = falses(seen_size)
   ```

2. replace `x ∈ visited` by `visited[hash_position(x, g)]` and
3. replace `push!(x, visited)` by `visited[hash_position(x, g)] = true`

although some care should be taken since `PeriodicGraphs.graph_width!` is only a heuristic.

Such algorithms can be used to compute the topological invariants like
[`coordination_sequence`](@ref) for example.

## API

```@docs
PeriodicGraphs.PeriodicEdgeIter
PeriodicGraphs.OffsetVertexIterator
PeriodicGraphs.graph_width!
```
