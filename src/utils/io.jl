# PeriodicGraph parsing and export as a string

struct KeyString{T,S<:AbstractString}
    x::S
    start::Base.RefValue{Int}
end
function KeyString{T}(x) where T
    isempty(x) && return KeyString{T,typeof(x)}(x, Ref(1))
    start = firstindex(x)
    idx = start
    this_char = first(x)
    while isspace(this_char)
        state = iterate(x, start)
        state isa Nothing && return KeyString{T,typeof(x)}(x, Ref(idx))
        start = idx
        this_char, idx = state
    end
    KeyString{T,typeof(x)}(x, Ref(start))
end

function iterate(k::KeyString{T}, start) where T
    state = iterate(k.x, start)
    state isa Nothing && return nothing
    next_char, idx = state
    tmp = start
    stop = start
    while next_char != '\0' && !isspace(next_char)
        state = iterate(k.x, idx)
        stop = tmp
        if state isa Nothing
            next_char = '\0'
        else
            tmp = idx
            next_char, idx = state
        end
    end
    ret = tryparse(T, SubString(k.x, start:stop))
    if ret isa T
        tmp = idx
        while isspace(next_char)
            state = iterate(k.x, tmp)
            idx = tmp
            if state isa Nothing
                next_char = '\0'
            else
                next_char, tmp = state
            end
        end
        return (ret, idx)
    end
    throw(ArgumentError("Input string does not represent a graph"))
end

Base.isempty(k::KeyString) = k.start[] > ncodeunits(k.x)
function iterate(k::KeyString)
    start = k.start[]
    start > ncodeunits(k.x) && return nothing
    iterate(k, start)
end

function Base.popfirst!(k::KeyString)
    next = iterate(k)
    next isa Nothing && throw(ArgumentError("Input string does not represent a graph"))
    char, state = next
    k.start[] = state isa Nothing ? (ncodeunits(k.x) + 1) : last(state)
    return char
end

Base.IteratorSize(::Type{<:KeyString}) = Base.SizeUnknown()

function edges_from_string(key::KeyString{Int}, ::Val{N}) where N
    edgs = PeriodicEdge{N}[]
    while !isempty(key)
        src = popfirst!(key)
        dst = popfirst!(key)
        ofs = SVector{N,Int}([popfirst!(key) for _ in 1:N])
        push!(edgs, PeriodicEdge{N}(src, dst, ofs))
    end
    return edgs
end


"""
    PeriodicGraph(key::AbstractString)
    PeriodicGraph{N}(key::AbstractString)

Construct a `PeriodicGraph{N}` from a `key`, which is a string of whitespace-separated
values of the form
`"N src1 dst1 ofs1_1 ofs1_2 ... ofs1_N src2 dst2 ofs2_1 ofs2_2 ... ofs2_N  ...  srcm dstm ofsm_1 ofsm_2 ... ofsm_N"`
where `N` is the number of repeating dimensions of the graph, `m` is the number of edges
and for all `i` between `1` and `m`, the number of edges, the `i`-th edge is described as
- `srci`, the vertex identifier of the source vertex,
- `dsti`, the vertex identifier of the destination vertex and
- `(ofsi_1, ofsi_2, ..., ofsi_N)` the offset of the edge.

This compact representation of a graph can be obtained simply by `print`ing the graph
or with `string`.

!!! note
    Use `parse(PeriodicGraph, key)` or `parse(PeriodicGraph{N}, key)` for a faster
    implementation if `key` was obtained from `string(g)` with `g` a `PeriodicGraph{N}`.

## Examples
```jldoctest
julia> PeriodicGraph("2  1 2 0 0  2 1 1 0  1 1 0 -1")
PeriodicGraph2D(2, PeriodicEdge2D[(1, 1, (0,1)), (1, 2, (-1,0)), (1, 2, (0,0))])

julia> PeriodicGraph3D("3  1 1 0 0 1  1 1 0 1 0  1 1 1 0 0")
PeriodicGraph3D(1, PeriodicEdge3D[(1, 1, (0,0,1)), (1, 1, (0,1,0)), (1, 1, (1,0,0))])

julia> string(ans)
"3 1 1 0 0 1 1 1 0 1 0 1 1 1 0 0"

julia> string(parse(PeriodicGraph3D, ans)) == ans
true
```
"""
function PeriodicGraph{N}(s::AbstractString) where N
    key = KeyString{Int}(s)
    M = popfirst!(key)
    M != N && throw(DimensionMismatch("Cannot construct a $N-dimensional graph from a $M-dimensional key"))
    return PeriodicGraph{N}(edges_from_string(key, Val(N)))
end
function PeriodicGraph(s::AbstractString)
    key = KeyString{Int}(s)
    N = popfirst!(key)
    return PeriodicGraph{N}(edges_from_string(key, Val(N)))
end

function _parse(::Type{PeriodicGraph{N}}, key::KeyString{Int})::PeriodicGraph{N} where N
    return from_edges(edges_from_string(key, Val(N)))
end

"""
    parse(::Type{PeriodicGraph}, key::AbstractString)
    parse(::Type{PeriodicGraph{N}}, key::AbstractString)

Parse a string representation of a `PeriodicGraph` back to a `PeriodicGraph`.

See [`PeriodicGraph(key::AbstractString)`](@ref PeriodicGraph{N}(key::AbstractString)) or
[`PeriodicGraph{N}(key::AbstractString)`](@ref) for details on the string representations
of a `PeriodicGraph`.

!!! warning
    This function assumes that the string is a valid representation of a `PeriodicGraph`
    with all its edges direct, sorted in lexicographical order and unique, as obtained
    from `string(g)` or `print(g)` where `g` is a `PeriodicGraph`.
    No check is performed to ensure this condition.

    If the input string may not obey these conditions, use
    [`PeriodicGraph(key)`](@ref PeriodicGraph{N}(s::AbstractString)) or
    [`PeriodicGraph{N}(key)`](@ref PeriodicGraph{N}(s::AbstractString)) instead.
"""
function Base.parse(::Type{PeriodicGraph{N}}, s::AbstractString) where N
    key = KeyString{Int}(s)
    popfirst!(key) # No verification is done for this function
    return _parse(PeriodicGraph{N}, key)
end
function Base.parse(::Type{PeriodicGraph}, s::AbstractString)
    key = KeyString{Int}(s)
    N = popfirst!(key)
    return _parse(PeriodicGraph{N}, key)
end

function show(io::IO, g::PeriodicGraph{N}) where N
    print(io, PeriodicGraph{N}, '(', nv(g), ',', ' ', collect(edges(g)), ')')
end
function print(io::IO, g::PeriodicGraph{N}) where N
    print(io, N)
    for e in edges(g)
        print(io, ' ', e.src, ' ', e.dst.v, ' ', join(e.dst.ofs, ' '))
    end
end
