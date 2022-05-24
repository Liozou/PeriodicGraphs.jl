# Symmetry type definitions and interface

# The main `SymmetryGroup3D` type is defined in PeriodicGraphEmbeddings.jl since it requires a
# 3D embedding

export AbstractSymmetry,
       AbstractSymmetryGroup,
       NoSymmetryGroup,
       IncludingIdentity,
       IdentitySymmetry,
       SimpleSymmetry

"""
    abstract type AbstractSymmetry end

An abstract type representing a symmetry.

## Interface

Subtypes `T` of `AbstractSymmetry` should implement a method to make their objects `symm`
callable on any object upon which the symmetry can act.
For example, if `T` represents a symmetry of a 3D embedding of a graph:
- if `x` represents a 3D point, then `symm(x)` should be the image of that point.
- if `x` represents a 3D basis of space, then `symm(x)` should be the image of that basis.

When the symmetry is naturally defined on a discrete set of objects, their image should be
accessible by indexing on the symmetry.
With the same example as before:
- if `x` is a vertex, `symm[x]` is the image of `x` by the symmetry.
- if `i` is the integer identifier of vertex `x`, `symm[i]` is the identifier of `symm[x]`.
"""
abstract type AbstractSymmetry end
(symm::AbstractSymmetry)(x) = symm[x]

"""
    IdentitySymmetry <: AbstractSymmetry

Identity symmetry, i.e. such that for `s::IdentitySymmetry`, ∀`x`, `s[x] == x`.
"""
struct IdentitySymmetry <: AbstractSymmetry end
Base.getindex(::IdentitySymmetry, x) = x
Base.keytype(::Type{IdentitySymmetry}) = Any


"""
    SimpleSymmetry{K,T,D} <: AbstractSymmetry

Symmetry operation defined by a map of type `T` and an optional dict of type `D` accepting
keys of type `K`.
"""
struct SimpleSymmetry{K,T,D} <: AbstractSymmetry
    map::T
    dict::D
    SimpleSymmetry{K,T,D}(map::T, dict::D) where {K,T,D} = new{K,T,D}(map, dict)
end

"""
    SimpleSymmetry(map::T, dict::D) where {T,D}

Create a `SimpleSymmetry{keytype(D),T,D}` object `symm` such that `symm[x]` is
- `map[dict[x]]` if `x isa keytype(D)`.
- `map[x]` otherwise, if `x isa Integer`.
"""
SimpleSymmetry(map::T, dict::D) where {T,D} = SimpleSymmetry{keytype(D),T,D}(map, dict)

"""
    SimpleSymmetry(map)

Create a `SimpleSymmetry` object `symm` such that `symm[x] == map[x]` for all `x`.
"""
SimpleSymmetry(map) = SimpleSymmetry(map, IdentitySymmetry())

function Base.getindex(s::SimpleSymmetry{K}, x::Union{Integer,K}) where K
    if x isa K
        s.map[s.dict[x]]
    else
        s.map[x]
    end
end


"""
    AbstractSymmetryGroup{T<:AbstractSymmetry}

An abstract type representing the set of symmetries of a graph.

## Interface

Any `AbstractSymmetryGroup` type must define methods for `Base` functions
`unique`, `iterate`, `length` and `one`
such that, for any `s` of type `<: AbstractSymmetryGroup{T}`:
- `s(i)` is a representative on the symmetry orbit of `i` such that all elements on
  the orbit share the same representative. The representative should be an integer.
- `unique(s)` is an iterator over such representatives.
- iterating over `s` yields the list of symmetry operations `symm`, each represented as an
  object of type `T` (where `T <: `[`AbstractSymmetry`](@ref) is the parameter to `typeof(s)`).
  The identity symmetry should not be part of these yielded `symm`, except for the
  specific [`IncludingIdentity`](@ref) subtype of `AbstractSymmetryGroup`.
- `one(s)` is the identity symmetry of type `T`.
"""
abstract type AbstractSymmetryGroup{T<:AbstractSymmetry} end
Base.eltype(::Type{<:AbstractSymmetryGroup{T}}) where {T} = T

"""
    NoSymmetryGroup <: AbstractSymmetryGroup{IdentitySymmetry}

The trivial [`AbstractSymmetryGroup`](@ref) devoid of any symmetry operation.
`NoSymmetryGroup(num)` creates a `NoSymmetryGroup` over `num` unique representatives.
"""
struct NoSymmetryGroup <: AbstractSymmetryGroup{IdentitySymmetry}
    num::Int
end
(::NoSymmetryGroup)(i) = i
Base.unique(x::NoSymmetryGroup) = Base.OneTo(x.num)
Base.iterate(::NoSymmetryGroup) = nothing
Base.length(::NoSymmetryGroup) = 0
Base.one(::NoSymmetryGroup) = IdentitySymmetry()

NoSymmetryGroup(g::PeriodicGraph) = NoSymmetryGroup(nv(g))

"""
    IncludingIdentity{S<:AbstractSymmetry,T<:AbstractSymmetryGroup{S}} <: AbstractSymmetryGroup{S}

Wrapper around an [`AbstractSymmetry`](@ref) that explicitly includes the identity
operation.
"""
struct IncludingIdentity{S<:AbstractSymmetry,T<:AbstractSymmetryGroup{S}} <: AbstractSymmetryGroup{S}
    symm::T
    IncludingIdentity(s::T) where {T<:AbstractSymmetryGroup} = new{eltype(T),T}(s)
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
