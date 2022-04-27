# Symmetry type definitions and interface

# See also PeriodicGraphEmbeddings.jl for a concrete type implementation.

export AbstractGraphSymmetry,
       AbstractGraphSymmetryGroup,
       IdentityGraphSymmetry,
       NoSymmetryGroup,
       IncludingIdentity

"""
    abstract type AbstractGraphSymmetry end

An abstract type representing a symmetry of a graph

## Interface

Subtypes `T` of `AbstractGraphSymmetry` must implement methods for `Base.getindex` so that,
for `symm` of type `T`:
- if `x` is a vertex, `symm[x]` is the image of `x` by the symmetry
- if `i` is the identifier of vertex `x`, `symm[i]` is the identifier of `symm[x]`

`symm` should additionally be callable on any object upon which the symmetry can act.
For example, if `T` represents a symmetry of a 3D embedding of a graph:
- if `x` represents a 3D point, then `symm(x)` should be the image of that point.
- if `x` represents a 3D basis of space, then `symm(x)` should be the image of that basis.
"""
abstract type AbstractGraphSymmetry end

# # Prototype: a typical AbstractGraphSymmetry could be the following PeriodicGraphSymmetry
#
# struct PeriodicGraphSymmetry{D} <: AbstractGraphSymmetry
#     vmap::Vector{PeriodicVertex{D}}
#     rotation::Matrix{Int}
# end
# Base.getindex(symm::PeriodicGraphSymmetry, i::Integer) = symm.vmap[i].v
# function Base.getindex(symm::PeriodicGraphSymmetry{D}, x::PeriodicVertex{D}) where D
#     dst = symm.vmap[x.v]
#     _ofs = muladd(symm.rotation, x.ofs, dst.ofs)
#     PeriodicVertex{D}(dst.v, _ofs)
# end

"""
    IdentityGraphSymmetry <: AbstractGraphSymmetry

Identity graph symmetry, i.e. such that for `s::IdentityGraphSymmetry`, ∀`x``, `s[x] == x`.
"""
struct IdentityGraphSymmetry <: AbstractGraphSymmetry end
Base.getindex(::IdentityGraphSymmetry, x) = x
(::IdentityGraphSymmetry)(x) = x

"""
    AbstractGraphSymmetryGroup{T<:AbstractGraphSymmetry}

An abstract type representing the set of symmetries of a graph.

## Interface

Any `AbstractGraphSymmetryGroup` type must define methods for `Base` functions
`unique`, `iterate`, `length` and `one`
such that, for any `s` of type
`<: AbstractGraphSymmetryGroup`:
- `s(x)` is a representative on the symmetry orbit of vertex `x` such that all vertices on
  the orbit share the same representative. The representative should be an integer `i` such
  that `first(Iterators.rest(vertices(g), i))` is a vertex on the symmetry orbit of `x` in
  graph `g`.
- `unique(s)` is a `<:AbstractVector{Int}` listing such representatives.
- iterating over `s` yields the list of symmetry operations `symm`, each represented as an
  object of type `T` (where `T <: AbstractGraphSymmetry` is the parameter to `typeof(s)`).
  The identity symmetry should not be part of these yielded `symm`, except for the
  specific `IncludingIdentity` subtype of `AbstractGraphSymmetryGroup`.
- `one(s)` is the identity symmetry of type `T`.
"""
abstract type AbstractGraphSymmetryGroup{T<:AbstractGraphSymmetry} end
Base.eltype(::Type{<:AbstractGraphSymmetryGroup{T}}) where {T} = T

"""
    NoSymmetryGroup <: AbstractGraphSymmetryGroup{IdentityGraphSymmetry}

The trivial `AbstractGraphSymmetryGroup` devoid of any symmetry operation.
"""
struct NoSymmetryGroup <: AbstractGraphSymmetryGroup{IdentityGraphSymmetry}
    num::Int
    NoSymmetryGroup(g::PeriodicGraph) = new(nv(g))
end
(::NoSymmetryGroup)(i) = i
Base.unique(x::NoSymmetryGroup) = Base.OneTo(x.num)
Base.iterate(::NoSymmetryGroup) = nothing
Base.length(::NoSymmetryGroup) = 0
Base.one(::NoSymmetryGroup) = IdentityGraphSymmetry()

"""
    IncludingIdentity{S<:AbstractGraphSymmetry,T<:AbstractGraphSymmetryGroup{S}} <: AbstractGraphSymmetryGroup{S}

Wrapper around an `AbstractGraphSymmetry` that explicitly includes the identity operation.
"""
struct IncludingIdentity{S<:AbstractGraphSymmetry,T<:AbstractGraphSymmetryGroup{S}} <: AbstractGraphSymmetryGroup{S}
    symm::T
    IncludingIdentity(s::T) where {T<:AbstractGraphSymmetryGroup} = new{eltype(T),T}(s)
end
IncludingIdentity(s::IncludingIdentity) = s

(s::IncludingIdentity)(i) = (s.symm)(i)
Base.unique(s::IncludingIdentity) = unique(s.symm)
Base.iterate(s::IncludingIdentity) = (one(s.symm), nothing)
function Base.iterate(s::IncludingIdentity, state)
    x = state isa Nothing ? iterate(s.symm) : iterate(s.symm, something(state))
    x isa Nothing && return nothing
    a, b = x
    return a, Some(b)
end
Base.length(s::IncludingIdentity) = 1 + length(s.symm)
function Base.getindex(s::IncludingIdentity, i::Integer)
    i == 1 && return one(s.symm)
    s.symm[i-1] # typeof(s.symm) should implement getindex to make this work
end
Base.one(::IncludingIdentity) = error("Querying the identity symmetry operation of an IncludingIdentity is forbidden.")
