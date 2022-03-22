# PeriodicEdge definition and basic functions

export PeriodicEdge, PeriodicEdge1D, PeriodicEdge2D, PeriodicEdge3D

"""
    LoopException <: Exception

Error type for constructing an invalid `PeriodicEdge{N}` of the form (u, u, zeros(Int,N)).
Loops are not expected in the algorithms implemented in PeriodicGraphs.jl. If you
still want to construct them, use the `unsafe_edge{N}` constructor instead of
`PeriodicEdge{N}``.
"""
struct LoopException <: Exception
    src::Int
end
showerror(io::IO, e::LoopException) = print(io, "LoopException: a loop from vertex $(e.src) to itself in the same unit cell is a forbidden edges. Maybe the offset is wrong?")
@noinline __throw_loopexception(src) = throw(LoopException(src))

"""
    unsafe_edge{N}

Internal constructor for `PeriodicEdge{N}` that bypasses the loop check.
This function is not part of the official API, use with caution.
"""
struct unsafe_edge{N} end

"""
    PeriodicEdge{N} <: Graphs.SimpleGraphs.AbstractSimpleEdge{Int}

Edge type for an `N`-periodic graph.

An edge is uniquely determined by the vertex identifiers of its source and
destination, and the cell offset between the source vertex and the destination vertex.
"""
struct PeriodicEdge{N} <: Graphs.SimpleGraphs.AbstractSimpleEdge{Int}
    src::Int
    dst::PeriodicVertex{N}
    function (::Type{unsafe_edge{N}})(src, dst) where N
        return new{N}(src, dst)
    end
end
unsafe_edge{N}(src, dst, ofs) where {N} = unsafe_edge{N}(src, PeriodicVertex{N}(dst, ofs))
function PeriodicEdge{N}(src, dst::PeriodicVertex{N}) where N
    src == dst.v && all(iszero.(dst.ofs)) && __throw_loopexception(src)
    return unsafe_edge{N}(src, dst)
end
PeriodicEdge(src, dst::PeriodicVertex{N}) where {N} = PeriodicEdge{N}(src, dst)
function PeriodicEdge{N}(src, dst, offset) where N
    return PeriodicEdge{N}(src, PeriodicVertex{N}(dst, offset))
end
function PeriodicEdge(src, dst, offset::Union{SVector{N,T},NTuple{N,T}}) where {N,T}
    PeriodicEdge{N}(src, dst, offset)
end
function PeriodicEdge{N}((src, dst, offset)::Tuple{Any,Any,Union{SVector{N,T},NTuple{N,T}}}) where {N,T<:Integer}
    PeriodicEdge{N}(src, dst, offset)
end
function PeriodicEdge((src, dst, offset)::Tuple{Any,Any,Union{SVector{N,T},NTuple{N,T}}}) where {N,T<:Integer}
    PeriodicEdge{N}(src, dst, offset)
end

function convert(::Type{PeriodicEdge{N}}, (src, dst, offset)::Tuple{Any,Any,Any}) where {N}
    PeriodicEdge{N}(src, dst, offset)
end
function convert(::Type{PeriodicEdge}, (src, dst, offset)::Tuple{Any,Any,Union{SVector{N,T},NTuple{N,T}}}) where {N,T<:Integer}
    PeriodicEdge{N}(src, dst, offset)
end
function convert(::Union{Type{PeriodicEdge},Type{PeriodicEdge{N}}}, (src, v)::Tuple{Any,PeriodicVertex{N}}) where N
    PeriodicEdge{N}(src, v)
end

function show(io::IO, x::PeriodicEdge{N}) where N
    if get(io, :typeinfo, Any) != PeriodicEdge{N}
        print(io, PeriodicEdge{N})
    end
    print(io, '(', x.src, ", ", x.dst.v, ", (", join(x.dst.ofs, ','))
    N == 1 && print(io, ',')
    print(io, ')', ')')
end

const PeriodicEdge1D = PeriodicEdge{1}
const PeriodicEdge2D = PeriodicEdge{2}
const PeriodicEdge3D = PeriodicEdge{3}

function Graphs.reverse(e::PeriodicEdge{N}) where N
    return unsafe_edge{N}(e.dst.v, e.src, .-e.dst.ofs)
end
Graphs.src(e::PeriodicEdge) = e.src
Graphs.dst(e::PeriodicEdge) = e.dst.v

"""
    ofs(v::PeriodicVertex)
    ofs(e::PeriodicEdge)

The offset of a vertex or an edge.
"""
ofs(v::PeriodicVertex) = v.ofs
ofs(e::PeriodicEdge) = ofs(e.dst)
isindirectedge(e::PeriodicEdge) = e.src > e.dst.v || (e.src == e.dst.v && e.dst.ofs < zero(e.dst.ofs))
function cmp(x::PeriodicEdge{N}, y::PeriodicEdge{N}) where N
    c = cmp(x.src, y.src)
    iszero(c) || return c
    return cmp(x.dst, y.dst)
end
isless(x::PeriodicEdge{N}, y::PeriodicEdge{N}) where {N} = cmp(x, y) < 0
==(x::PeriodicEdge{N}, y::PeriodicEdge{M}) where {N,M} = N == M && iszero(cmp(x, y))
hash(x::PeriodicEdge, h::UInt) = hash(x.dst, hash(x.src, h))

ndims(::PeriodicEdge{N}) where {N} = N
