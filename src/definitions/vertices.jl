# PeriodicVertex definition and basic functions

export PeriodicVertex, PeriodicVertex1D, PeriodicVertex2D, PeriodicVertex3D

"""
    PeriodicVertex{N}

Vertex type for an `N`-periodic graph.

A vertex is uniquely determined by the identifier of its representative in the
a fixed initial cell, and the offset of the cell containing the vertex compared
to the the initial cell.
Vertex identifiers start at 1.
"""
struct PeriodicVertex{N}
    v::Int
    ofs::SVector{N,Int}
end
PeriodicVertex{N}(n::Integer) where {N} = PeriodicVertex{N}(n, zero(SVector{N,Int}))
PeriodicVertex(v::Integer, ofs::Union{SVector{N},NTuple{N}}) where {N} = PeriodicVertex{N}(v, ofs)

function show(io::IO, x::PeriodicVertex{N}) where N
    if get(io, :typeinfo, Any) != PeriodicVertex{N}
        print(io, PeriodicVertex{N})
    end
    print(io, '(', x.v, ", (", join(x.ofs, ','))
    N == 1 && print(io, ',')
    print(io, ')', ')')
end
function convert(::Type{PeriodicVertex{N}}, (dst, offset)::Tuple{Integer,Any}) where N
    PeriodicVertex{N}(dst, offset)
end
function convert(::Type{PeriodicVertex}, (dst, offset)::Tuple{Integer,Union{SVector{N},NTuple{N}}}) where N
    PeriodicVertex{N}(dst, offset)
end
function cmp(x::PeriodicVertex{N}, y::PeriodicVertex{N}) where N
    c = cmp(x.v, y.v)
    iszero(c) || return c
    return cmp(x.ofs, y.ofs)
end
isless(x::PeriodicVertex{N}, y::PeriodicVertex{N}) where {N} = cmp(x, y) < 0
==(x::PeriodicVertex{N}, y::PeriodicVertex{M}) where {N,M} = N == M && iszero(cmp(x, y))
hash(x::PeriodicVertex, h::UInt) = hash(x.ofs, hash(x.v, h))

const PeriodicVertex1D = PeriodicVertex{1}
const PeriodicVertex2D = PeriodicVertex{2}
const PeriodicVertex3D = PeriodicVertex{3}


ZtoN(x::Signed) = -(x<0) + 2*abs(x)
function hash_position((x1,x2,x3)::SVector{3,<:Integer})
    x1 = ZtoN(x1); x2 = ZtoN(x2); x3 = ZtoN(x3)
    _b = x1 >= x2
    b1 = _b & (x1 >= x3)
    b2 = (!_b) & (x2 >= x3)
    b3 = (!b1) & (!b2)
    return b1*(x2*(x1 + 1) + x3 + x1^3) +
           b2*((x2 + 1)*(x1 + x2 + 1) + x3 + x2^3) +
           b3*((x3 + 1)*(2x3 + 1) + x3*(x1 + x3^2) + x2)
end

function hash_position((x1,x2)::SVector{2,<:Integer})
    x1 = ZtoN(x1); x2 = ZtoN(x2)
    b = x1 >= x2
    return b*(x2 + x1^2) + (!b)*(x1 + x2*(x2 + 1) + 1)
end

function hash_position((x,)::SVector{1,<:Integer})
    return ZtoN(x)
end

function hash_position(::SVector{0,<:Integer})
    return 0
end

"""
    hash_position(x::PeriodicVertex{N}, n::Integer) where {N}

Given a `PeriodicVertex{N}` and the number n of vertex identifiers in a graph, compute a
unique hash for the given vertex.

This hash function is a bijection between the set of all the vertices of the graph and
the set of positive integers. Its value is an integer between `1+n\\*(2d-1)^N`
(or `1` if `d == 0`) and `n*(2d+1)^N`, where `d = max.(abs.(x.ofs))`.
This means that when one unit cell A is further than another B (for the Manhattan
distance), all vertices in A will have a larger hash than all vertices in B.
"""
function hash_position(x::PeriodicVertex, n::Integer)
    return x.v + n*hash_position(x.ofs)
end

ndims(::PeriodicVertex{N}) where {N} = N
