# Symmetries

`PeriodicGraphs.jl` comes with a general API for the manipulation of symmetries.
Concrete implementations of the API can be found in `PeriodicGraphEmbeddings.jl`.

## Single symmetry operations

```@docs
AbstractSymmetry
IdentitySymmetry
SimpleSymmetry
SimpleSymmetry(map::T, dict::D) where {T,D}
SimpleSymmetry(map)
```

## Symmetry group

```@docs
AbstractSymmetryGroup
NoSymmetryGroup
IncludingIdentity
```
