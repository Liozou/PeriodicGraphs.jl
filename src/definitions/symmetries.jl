# Symmetry type definitions and interface

# The main `Symmetries` type is defined in PeriodicGraphEmbeddings.jl since it requires a
# 3D embedding

export AbstractSymmetries, NoSymmetry, IncludingIdentity

"""
    AbstractSymmetries

An abstract type representing the set of symmetries of a graph.

## Interface

Any `AbstractSymmetries` type must define methods for `Base` functions `getindex`, `unique`,
`iterate`, `eltype` and `length` such that, for `s` of type `T <: AbstractSymmetries`,
- `s[i]` is a representative on the symmetry orbit of vertex `i` such that all vertices on
  the orbit share the same representative.
- `unique(s)` is the list of such representatives.
- `for vmap in s` yields the list of symmetry operations, represented as `vmaps`. The
  identity should not be part of these `vmaps`, except for the `IncludingIdentity` type.
"""
abstract type AbstractSymmetries end

"""
    NoSymmetry <: AbstractSymmetries

The trivial `AbstractSymmetries` devoid of any symmetry operation.
"""
struct NoSymmetry <: AbstractSymmetries
    num::Int
    NoSymmetry(g::PeriodicGraph) = new(nv(g))
end
Base.getindex(::NoSymmetry, i::Integer) = i
Base.unique(x::NoSymmetry) = Base.OneTo(x.num)
Base.iterate(::NoSymmetry, _=nothing) = nothing
Base.eltype(::Type{NoSymmetry}) = Base.OneTo
Base.length(::NoSymmetry) = 0


"""
    IncludingIdentity{T<:AbstractSymmetries} <: AbstractSymmetries

Wrapper around an `AbstractSymmetry` that explicitly includes the identity operation.
"""
struct IncludingIdentity{T<:AbstractSymmetries} <: AbstractSymmetries
    symm::T
end

Base.getindex(s::IncludingIdentity, i::Integer) = s.symm[i]
Base.unique(s::IncludingIdentity) = unique(s.symm)
function Base.iterate(s::IncludingIdentity{T}, state=nothing) where T
    if state isa Nothing
        n = T <: NoSymmetry ? s.symm.num : length(first(s.symm))
        elt = eltype(s.symm)
        if elt isa SubArray{Int,1,Matrix{Int},Tuple{Base.Slice{Base.OneTo{Int}},Int},true}
            # Special case for the `Symmetries` type because SubArray{...}(1:n) errors
            return first(eachcol(reshape(collect(Base.OneTo(n)), n, 1))), 1
        else
            return elt(Base.OneTo(n)), 1
        end
    end
    return iterate(s.symm, state)
end
Base.eltype(::Type{IncludingIdentity{T}}) where {T} = eltype(T)
Base.length(s::IncludingIdentity) = 1 + length(s.symm)
