# Ring statistics

## Definitions

In this section, we use the terminology recommended by
[Blatov, O’Keeffe and Proserpio](https://doi.org/10.1016/j.jssc.2005.06.011).
In particular:
- a *cycle* is a sequence of vertices `v₁, v₂, ..., vₙ` such that for all `i` between 2 and
  `n`, `vᵢ₋₁` and `vᵢ` are neighbors, as well as `v₁` and `vₙ`, and no `vᵢ` occurs more
  than once in the sequence. It can be equivalently represented by the sequence of edges
  between two consecutive vertices and between `v₁` and `vₙ`.
- the *sum* of two or more cycles is the set of edges occurring only an odd number of times
  in the set of input cycles. The sum of two cycles is thus the symmetric difference of
  their edges. Note that the sum of cycles may be empty, or may be the juxtaposition of
  several edge-disjoint cycles.
- the *length* of a cycle is its number of edges. It is also equal to its number of
  vertices since no vertex is repeated in the sequence. A cycle `a` is strictly smaller
  than another `b` when the length of `a` is strictly lower than that of `b`.
- a *ring* is a cycle which is not the sum of two strictly smaller cycles. An equivalent
  definition is that a ring is a cycle which does not admit a short-circuit: between two
  vertices of the ring, there is no path of the graph strictly smaller than both branches
  of the ring linking the two vertices.
- a *strong ring* is a cycle which is not the sum of any number of strictly smaller cycles.

## Manual

`PeriodicGraphs.jl` provides an algorithm for the determination of all rings and strong
rings up to a given size in a `PeriodicGraph`.
It can also be used on finite graphs by converting them to `PeriodicGraph{0}`, although the
code is not optimized for this case.

The [`rings`](@ref) (respectively [`strong_rings`](@ref)) function returns the list of all
rings (respectively strong rings) in the graph up to a given size:

```@docs
rings
strong_rings
```

A few notes on the output:
- The output is a pair `(rs, symm)` where `rs` is the list of rings (respectively strong
  rings) and `symm` is an [`AbstractSymmetryGroup`](@ref) which contains symmetries of `rs`.
  `symm` is always a [`NoSymmetryGroup`](@ref) if the optional argument `symmetries` is not
  provided.
- The returned list of rings `rs` is generally unsorted.
- All translations of the same ring are represented by a unique ring, which means that a
  ring crossing through different unit cells will only appear once in the list, even though
  it may appear several times in a single unit cell.
- The `symmetries` optional argument reduces the computational cost of the algorithm.
  The output lists `rs` with and without the optional argument are identical except for the
  order of their elements.

The optional argument `depth` defaults to 15, which means that rings containing up to 33
edges will be considered. This default value is chosen to accomodate the vast majority of
periodic nets encountered as crystal nets, for which the ring size rarely exceeds 20.

Let's take as example an aperiodic graph representing a small house, made of a cube with
a pyramid on top:

```jldoctest house; setup=:(using PeriodicGraphs, Graphs)
julia> house = PeriodicGraph{0}("0 "*
                    "1 2  2 3  3 4  4 1 "* # square base of the house
                    "1 5  2 6  3 7  4 8 "* # 4 vertical pillars
                    "5 6  6 7  7 8  8 5 "* # square ceiling
                    "5 9  6 9  7 9  8 9 "  # pyramidal roof
       );

julia> sort!(first(rings(house)))
14-element Vector{Vector{Int64}}:
 [1, 2, 3, 4]
 [1, 2, 3, 7, 8, 5]
 [1, 2, 6, 5]
 [1, 2, 6, 7, 8, 4]
 [1, 4, 3, 7, 6, 5]
 [1, 4, 8, 5]
 [2, 3, 4, 8, 5, 6]
 [2, 3, 7, 6]
 [3, 4, 8, 7]
 [5, 6, 7, 8]
 [5, 6, 9]
 [5, 8, 9]
 [6, 7, 9]
 [7, 8, 9]

julia> sort!(first(strong_rings(house)))
9-element Vector{Vector{Int64}}:
 [1, 2, 3, 4]
 [1, 2, 6, 5]
 [1, 4, 8, 5]
 [2, 3, 7, 6]
 [3, 4, 8, 7]
 [5, 6, 9]
 [5, 8, 9]
 [6, 7, 9]
 [7, 8, 9]
```

We can see that the house has four weak rings of size 6, six rings of size 4 among
which five are strong, and four strong rings of size 3.

The strong rings are the faces of the house: there are four triangles that make the roof,
four squares that make the walls and one last square for the base of the house. The square
corresponding to the ceiling is actually the sum of the four triangles of the roof, which
is why it is not a strong ring. The four weak rings of size 6 are those that go through
each vertex of the cube except for one of the four opposite pairs.

To explore the ring distributions around individual vertices, the [`RingAttributions`](@ref)
struct factors the ring distribution by vertex. The list of rings including a particular
vertex is factored into a [`RingIncluding`](@ref) struct:

```@docs
RingAttributions
RingIncluding
```

To avoid useless computations, the lists of rings and the rings themselves are returned as
iterables instead of `Vector{PeriodicVertex{N}}` and such, so they should be `collect`ed
if required.

Let's look all the rings on our little house:

```jldoctest house
julia> ras = RingAttributions(house)
RingAttributions{0}(rings per node: [6, 6, 6, 6, 8, 8, 8, 8, 4])

julia> roofpeak = ras[9] # the list of rings including the top of the roof
4-element RingIncluding{0}:
 [5, 8, 9]
 [5, 6, 9]
 [7, 8, 9]
 [6, 7, 9]

julia> roofpeak[2] # the second ring including the top of the roof
3-element PeriodicGraphs.OffsetVertexIterator{0}:
 5
 6
 9

julia> rasstrong = RingAttributions(house, true)
RingAttributions{0}(rings per node: [3, 3, 3, 3, 4, 4, 4, 4, 4])

julia> rasstrong[1] # the base and two walls make the strong rings around vertex 1
3-element RingIncluding{0}:
 [1, 2, 6, 5]
 [1, 4, 8, 5]
 [1, 2, 3, 4]

julia> collect(rasstrong[5]) # two rooftiles and two walls make the strong rings around vertex 5
4-element Vector{PeriodicGraphs.OffsetVertexIterator{0}}:
 [5, 8, 9]
 [5, 6, 9]
 [1, 2, 6, 5]
 [1, 4, 8, 5]
```

## Internal API

Here is a collection of internal utilities used for the algorithms of [`rings`](@ref) and
[`strong_rings`](@ref):

```@docs
PeriodicGraphs.ConstMiniBitSet
PeriodicGraphs.DistanceRecord
PeriodicGraphs.JunctionNode
PeriodicGraphs.PhantomJunctionNode
PeriodicGraphs.arcs_list
PeriodicGraphs.RingsEndingAt
PeriodicGraphs.normalize_cycle!
PeriodicGraphs.symdiff_cycles!
PeriodicGraphs.symdiff_cycles
PeriodicGraphs.IterativeGaussianElimination
PeriodicGraphs.gaussian_elimination!
PeriodicGraphs.intersect_cycles!
PeriodicGraphs.intersect_cycles
PeriodicGraphs.retrieve_track!
PeriodicGraphs.rings_around
PeriodicGraphs.EdgeDict
strong_erings
PeriodicGraphs.convert_to_ering!
```
