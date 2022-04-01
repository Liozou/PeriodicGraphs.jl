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

ndims(::PeriodicVertex{N}) where {N} = N


@inline ZtoN(x::Signed) = abs(2x + signbit(x)) # abs(bitrotate(x, 1)) # abs((x << 1) | (x >>> 63)) # (abs(x) << 1) - (x<0)
@inline function hash_position(x::SVector{3,<:Integer})
    @inbounds begin x1 = x[1]; x2 = x[2]; x3 = x[3] end
    x1 = ZtoN(x1); x2 = ZtoN(x2); x3 = ZtoN(x3)
    _b = x1 >= x2
    b1 = _b & (x1 >= x3)    # x1 ≥ x2 and x1 ≥ x3
    b2 = (!_b) & (x2 >= x3) # x2 > x1 and x2 ≥ x3
    b3 = (!b1) & (!b2)      # x3 > x1 and x3 > x2
    return b1*(x2*(x1 + 1) + x3 + x1^3) +
           b2*((x2 + 1)*(x1 + x2 + 1) + x3 + x2^3) +
           b3*((x3 + 1)*(2x3 + 1) + x3*(x1 + x3^2) + x2)
end

@inline function hash_position(x::SVector{2,<:Integer})
    @inbounds begin x1 = x[1]; x2 = x[2] end
    x1 = ZtoN(x1); x2 = ZtoN(x2)
    return ifelse(x1 >= x2, x2 + x1^2, x1 + x2*(x2 + 1) + 1)
end

@inline function hash_position(x::SVector{1,<:Integer})
    return ZtoN(@inbounds x[1])
end

@inline function hash_position(::SVector{0,<:Integer})
    return 0
end

@static if Int === Int64
    for N in 1:3 # accelerated hash computation
        @eval hash(x::PeriodicVertex{$N}, h::UInt) = hash(x.v << 32 ⊻ hash_position(x.ofs), h)
    end
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
@inline function hash_position(x::PeriodicVertex, n::Integer)
    return x.v + n*hash_position(x.ofs)
end



@inline NtoZ(n) = (((-1)^isodd(n))*n) >> 1

function reverse_hash_position(x::Integer, ::Val{3})
    x == 0 && return zero(SVector{3,Int})
    n = floor(Int, cbrt(x))
    zn = NtoZ(n)
    n2 = n^2
    n3 = n*n2
    y = x - n3
    if y ≤ 2n+n2
        x2, x3 = divrem(y, n+1)
        return SVector{3,Int}(zn, NtoZ(x2), NtoZ(x3))
    end
    y -= n2 + 2n + 1
    if y ≤ n2 - 1 + n
        x1, x3 = divrem(y, n+1)
        return SVector{3,Int}(NtoZ(x1), zn, NtoZ(x3))
    end
    y -= n2 + n
    x1, x2 = divrem(y, n)
    return SVector{3,Int}(NtoZ(x1), NtoZ(x2), zn)
end

function reverse_hash_position(x::Integer, ::Val{2})
    x == 0 && return zero(SVector{2,Int})
    n = floor(Int, sqrt(x))
    zn = NtoZ(n)
    n2 = n^2
    if x ≤ n^2 + n
        return SVector{2,Int}(zn, NtoZ(x - n2))
    end
    return SVector{2,Int}(NtoZ(x - 1 - n2 - n), zn)
end

reverse_hash_position(x::Integer, ::Val{1}) = SVector{1,Int}(NtoZ(x))

reverse_hash_position(::Integer, ::Val{0}) = SVector{0,Int}()

function reverse_hash_position(hash::Integer, n::Integer, ::Val{N}) where N
    x, node = fldmod1(hash, n)
    return PeriodicVertex{N}(node, reverse_hash_position(x-1, Val(N)))
end
