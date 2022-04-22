# Symmetry type definitions and interface

# The main `SymmetryGroup3D` type is defined in PeriodicGraphEmbeddings.jl since it requires a
# 3D embedding

export AbstractSymmetry, PeriodicSymmetry, AbstractSymmetryGroup, NoSymmetryGroup, IncludingIdentity

"""
    abstract type AbstractSymmetry end

An abstract type representing a symmetry of a graph

## Interface

Subtypes of `AbstractSymmetry` must implement a method for `Base.getindex` so that, for
`symm` of type `T <: AbstractSymmetry` and `x` a vertex, `symm[x]` is the image of `x` by
the symmetry.
Additionally, `symm[i]` should be the identifier of vertex `symm[x]` when `i` is the
identifier of vertex `x`.
"""
abstract type AbstractSymmetry end

# # Prototype: a typical AbstractSymmetry could be the following PeriodicSymmetry
#
# struct PeriodicSymmetry{D} <: AbstractSymmetry
#     vmap::Vector{PeriodicVertex{D}}
#     rotation::Matrix{Int}
# end
# Base.getindex(symm::PeriodicSymmetry, i::Integer) = symm.vmap[i].v
# function Base.getindex(symm::PeriodicSymmetry{D}, x::PeriodicVertex{D}) where D
#     dst = symm.vmap[x.v]
#     _ofs = muladd(symm.rotation, x.ofs, dst.ofs)
#     PeriodicVertex{D}(dst.v, _ofs)
# end

struct TrivialIdentitySymmetry <: AbstractSymmetry end
Base.getindex(::TrivialIdentitySymmetry, x) = x

"""
    AbstractSymmetryGroup{T<:AbstractSymmetry}

An abstract type representing the set of symmetries of a graph.

## Interface

Any `AbstractSymmetryGroup` type must define methods for `Base` functions
`unique`, `iterate`, `eltype`, `length` and `one`
such that, for any `s` of type
`<: AbstractSymmetryGroup`:
- `s(i)` is a representative on the symmetry orbit of vertex `i` such that all vertices on
  the orbit share the same representative.
- `unique(s)` is the list of such representatives.
- iterating over `s` yields the list of symmetry operations `symm`, each represented as an
  object of type `T` (where `T <: AbstractSymmetry` is the parameter to `typeof(s)`).
  The identity symmetry should not be part of these yielded `symm`, except for the
  specific `IncludingIdentity` subtype of `AbstractSymmetryGroup`.
- `one(s)` is the identity symmetry of type `T`.
"""
abstract type AbstractSymmetryGroup{T<:AbstractSymmetry} end

"""
    NoSymmetryGroup <: AbstractSymmetryGroup{PeriodicSymmetry}

The trivial `AbstractSymmetryGroup` devoid of any symmetry operation.
"""
struct NoSymmetryGroup <: AbstractSymmetryGroup{TrivialIdentitySymmetry}
    num::Int
    NoSymmetryGroup(g::PeriodicGraph) = new(nv(g))
end
(::NoSymmetryGroup)(i::Integer) = i
Base.unique(x::NoSymmetryGroup) = Base.OneTo(x.num)
Base.iterate(::NoSymmetryGroup) = nothing
Base.eltype(::Type{NoSymmetryGroup}) = PeriodicSymmetry{0}
Base.length(::NoSymmetryGroup) = 0
Base.one(::NoSymmetryGroup) = TrivialIdentitySymmetry()

"""
    IncludingIdentity{T<:AbstractSymmetryGroup} <: AbstractSymmetryGroup

Wrapper around an `AbstractSymmetry` that explicitly includes the identity operation.
"""
struct IncludingIdentity{S<:AbstractSymmetry,T<:AbstractSymmetryGroup{S}} <: AbstractSymmetryGroup{S}
    symm::T
    IncludingIdentity{T}(symm::AbstractSymmetryGroup{S}) where {S,T} = IncludingIdentity{S,T}(symm)
end
IncludingIdentity(s::T) where {T<:AbstractSymmetryGroup} = IncludingIdentity{T}(s)
IncludingIdentity(s::IncludingIdentity) = s

(s::IncludingIdentity)(i::Integer) = (s.symm)(i)
Base.unique(s::IncludingIdentity) = unique(s.symm)
Base.iterate(s::IncludingIdentity) = (one(s.symm), nothing)
function Base.iterate(s::IncludingIdentity{T}, state) where T
    x = state isa Nothing ? iterate(s.symm) : iterate(s.symm, something(state))
    x isa Nothing && return nothing
    a, b = x
    return a, Some(b)
end
Base.eltype(::Type{IncludingIdentity{T}}) where {T} = eltype(T)
Base.length(s::IncludingIdentity) = 1 + length(s.symm)
function Base.getindex(s::IncludingIdentity, i::Integer)
    i == 1 && return one(s.symm)
    s.symm[i-1]
end
Base.one(::IncludingIdentity) = error("Querying the identity symmetry operation of an IncludingIdentity is forbidden.")
