#=
mutable struct List{T}
    num::Int
    root::T
    head::T
    tail::List{T}

    List{T}() where {T} = new{T}(0)
    List(x::T) where {T} = new{T}(1, x, x)
    List(root::T, head::T) where {T} = new{T}(0, root, head)

    List(x, l::List{T}) where {T} = new{T}(l.num+1, l.root, x, l)
end
length(l::List{T}) where {T} = l.num
eltype(::List{T}) where {T} = T

iterate(l::List{T}) where {T} = l.num == 0 ? nothing : (l.head, l.tail)
iterate(::List{T}, l::List{T}) where {T} = iterate(l)

function Base.show(io::IO, l::List)
    print(io, "((")
    for (i,x) in enumerate(l)
        print(io, x)
        i != l.num && print(io, ", ")
    end
    print(io, "))")
end

function collect!(ret, x::List{T}, reverse=false, istart=1) where T
    y = x
    for i in (reverse ? ((x.num+istart-1):-1:istart) : (istart:(x.num+istart-1)))
        ret[i] = y.head
        y = y.tail
    end
    return ret
end
Base.collect(x::List{T}) where {T} = collect!(Vector{T}(undef, x.num), x)
=#

struct ConstMiniBitSet <: AbstractSet{Int}
    x::UInt64
    global _constminibitset(x::UInt64) = new(x)
end

ConstMiniBitSet() = _constminibitset(zero(UInt64))
ConstMiniBitSet(i::Integer) = _constminibitset(one(UInt64) << (i % UInt8))
function ConstMiniBitSet(l::AbstractVector{<:Integer})
    x = zero(UInt64)
    for i in l
        x |= one(UInt64) << (i % UInt8)
    end
    return _constminibitset(x)
end

push(x::ConstMiniBitSet, i::Integer) = _constminibitset(x.x | (one(UInt64) << (i % UInt8)))
Base.union(x::ConstMiniBitSet, y::ConstMiniBitSet) = _constminibitset(x.x | y.x)
Base.intersect(x::ConstMiniBitSet, y::ConstMiniBitSet) = _constminibitset(x.x & y.x)
Base.in(x::ConstMiniBitSet, i::Integer) = ((x.x >> (i % UInt8)) & one(UInt64)) % Bool
Base.setdiff(x::ConstMiniBitSet, y::ConstMiniBitSet) = _constminibitset(x.x & ~y.x)
Base.symdiff(x::ConstMiniBitSet, y::ConstMiniBitSet) = _constminibitset(x.x ⊻ y.x)
Base.length(x::ConstMiniBitSet) = count_ones(x.x)
Base.isempty(x::ConstMiniBitSet) = iszero(x.x)
Base.eltype(::Type{ConstMiniBitSet}) = Int

function Base.iterate(::ConstMiniBitSet, x::UInt64)
    first = trailing_zeros(x)
    first == 64 && return nothing
    return (first, x & (~(one(UInt64) << (first % UInt8))))
end
Base.iterate(x::ConstMiniBitSet) = iterate(x, x.x)
Base.isdone(::ConstMiniBitSet, x::UInt64) = iszero(x)
Base.isdone(x::ConstMiniBitSet) = Base.isdone(x, x.x)

minabove(x::ConstMiniBitSet, i::Integer) = trailing_zeros(x.x & ((~one(UInt64)) << (i % UInt8)))
minaboveeq(x::ConstMiniBitSet, i::Integer) = trailing_zeros(x.x & ((~zero(UInt64)) << (i % UInt8)))
# nonemptyintersect(x::ConstMiniBitSet, y::ConstMiniBitSet) = !iszero(x.x & y.y)

Base.minimum(x::ConstMiniBitSet) = first(x)
Base.maximum(x::ConstMiniBitSet) = 64 - leading_zeros(x.x)

function Base.show(io::IO, x::ConstMiniBitSet)
    print(io, ConstMiniBitSet, "([")
    join(io, x, ", ")
    print(io, "])")
end


mutable struct MiniBitSet <: AbstractSet{Int}
    x::UInt64
    global _minibitset(x::UInt64) = new(x)
end

MiniBitSet(x::ConstMiniBitSet) = _minibitset(x.x)
ConstMiniBitSet(x::MiniBitSet) = _constminibitset(x.x)

MiniBitSet() = _minibitset(zero(UInt64))
MiniBitSet(i::Integer) = _minibitset(one(UInt64) << (i % UInt8))
MiniBitSet(l::AbstractVector{<:Integer}) = MiniBitSet(ConstMiniBitSet(l))

Base.push!(x::MiniBitSet, i::Integer) = (x.x |= (one(UInt64) << (i % UInt8)); x)
Base.union!(x::MiniBitSet, y::MiniBitSet) = (x.x |= y.x; x)
Base.intersect!(x::MiniBitSet, y::MiniBitSet) = (x.x &= y.x; x)
Base.in(x::MiniBitSet, i::Integer) = ((x.x >> (i % UInt8)) & one(UInt64)) % Bool
Base.setdiff!(x::MiniBitSet, y::MiniBitSet) = (x.x &= ~y.x; x)
Base.symdiff!(x::MiniBitSet, y::MiniBitSet) = (x.x ⊻= y.x; x)
Base.copymutable(x::MiniBitSet) = _minibitset(x.x)
Base.copy(x::MiniBitSet) = _minibitset(x.x)

Base.iterate(x::MiniBitSet, u::UInt64=x.x) = iterate(ConstMiniBitSet(), u)
Base.isdone(::MiniBitSet, x::UInt64) = iszero(x)
Base.isdone(x::MiniBitSet) = Base.isdone(x, x.x)

minabove(x::MiniBitSet, i::Integer) = minabove(ConstMiniBitSet(x), i)
minaboveeq(x::MiniBitSet, i::Integer) = minaboveeq(ConstMiniBitSet(x), i)
# isnonemptyintersect(x::MiniBitSet, y::MiniBitSet) = !iszero(x.x & y.x)
# isnonemptysymdiff(x::MiniBitSet, y::MiniBitSet) = !iszero(x.x ⊻ y.x)
hasonly(x::MiniBitSet, i::Integer) = iszero(x.x & ~(one(UInt64) << (i % UInt8)))
# hasotherthan(x::MiniBitSet, i::Integer) = !iszero(x.x & ~(one(UInt64) << (i % UInt8)))

Base.length(x::MiniBitSet) = length(ConstMiniBitSet(x))
Base.eltype(x::MiniBitSet) = eltype(ConstMiniBitSet(x))
Base.maximum(x::MiniBitSet) = maximum(ConstMiniBitSet(x))
Base.minimum(x::MiniBitSet) = minimum(ConstMiniBitSet(x))
Base.isempty(x::MiniBitSet) = isempty(ConstMiniBitSet(x))

function Base.show(io::IO, x::MiniBitSet)
    print(io, MiniBitSet, "([")
    join(io, x, ", ")
    print(io, "])")
end


struct JunctionNode
    heads::Vector{Int}
    num::Int # length of the shortest branch
    shortroots::MiniBitSet
    longroots::MiniBitSet
    lastshort::Base.RefValue{Int}
    function JunctionNode(heads::AbstractVector, num, roots::MiniBitSet)
        new(heads, num, roots, MiniBitSet(), Ref(length(heads)))
    end
end
JunctionNode() = JunctionNode(Int[], -1, MiniBitSet())
JunctionNode(head::Integer, num, roots::MiniBitSet) = JunctionNode(Int[head], num, roots)
JunctionNode(root::Integer) = JunctionNode(1, 0, MiniBitSet(root))

function Base.show(io::IO, x::JunctionNode)
    print(io, "JunctionNode([")
    join(io, x.heads, ", ")
    print(io, "])")
end


function arcs_list(g::PeriodicGraph{D}, i, depth) where D
    Q = Tuple{PeriodicVertex{D},JunctionNode}[(PeriodicVertex{D}(i), JunctionNode())]
    sizehint!(Q, ceil(Int, depth^2.9)) # empirical estimate
    append!(Q, (x, JunctionNode(j+1)) for (j,x) in enumerate(neighbors(g, i)))
    if length(Q) > 62
        error("The vertex has a degree too large (> 62) to be represented in a MiniBitSet. Please open an issue.")
    end
    vertexdict = Dict{PeriodicVertex{D},Int}((x[1] => j) for (j,x) in enumerate(Q))
    counter = 1
    for (u, parents) in Iterators.rest(Q, 2)
        parents.num == depth && break
        counter += 1
        previous = Q[parents.heads[1]][1]
        num = 1 + parents.num
        for x in neighbors(g, u)
            x == previous && continue
            idx = get!(vertexdict, x, length(Q)+1)
            if idx == length(Q)+1
                push!(Q, (x, JunctionNode(counter, num, copy(parents.shortroots))))
            elseif idx > counter
                junction = Q[idx][2]
                push!(junction.heads, counter)
                if num == junction.num
                    junction.lastshort[] += 1
                    union!(junction.shortroots, parents.shortroots)
                else
                    union!(junction.longroots, parents.shortroots)
                end
            end
        end
    end
    return Q, vertexdict
end


function next_compatible_arc!(buffer, last_positions, idx_stack, Q, check)
    num = length(idx_stack)
    incompatibleflag = true
    updated_last_positions = last_positions
    haslongarc = length(buffer) - 2*num == 3
    root = buffer[num - haslongarc + 2]
    while incompatibleflag
        trail = trailing_ones(updated_last_positions) % UInt8
        next_toupdate = num - trail
        next_toupdate < 1 && return ~zero(UInt64)
        updated_last_positions = updated_last_positions & (~zero(UInt64) << trail)
        incompatibleflag = false
        idx = idx_stack[next_toupdate] + 1
        head = buffer[check ? next_toupdate : end-next_toupdate]
        parents = Q[head][2]
        for i in next_toupdate:num
            heads = parents.heads
            head = heads[idx]
            _i = i + 1
            last_head = parents.lastshort[]
            if check
                next_parents = Q[head][2]
                while hasonly(next_parents.shortroots, root) || head == buffer[end-_i-haslongarc]
                    if idx == last_head
                        last_positions |= (~zero(UInt64) >> ((63 - num + i) % UInt8))
                        incompatibleflag = true
                        break
                    end
                    idx += 1
                    head = heads[idx]
                    next_parents = Q[head][2]
                end
                incompatibleflag && break
                parents = next_parents
            else
                parents = Q[head][2]
            end
            buffer[check ? _i : end-_i] = head
            idx_stack[i] = idx
            if idx == last_head
                updated_last_positions |= (one(UInt64) << ((num - i) % UInt8))
            end
            idx = 1
        end
    end
    return updated_last_positions
end

function initial_compatible_arc!(buffer, idx_stack, Q, check)
    last_positions = zero(UInt64)
    num = length(idx_stack)
    head = buffer[check ? 1 : end-1]
    haslongarc = length(buffer) - 2*num == 3
    root = buffer[num - haslongarc + 2]
    parents = Q[head][2]
    for i in 1:num
        heads = parents.heads
        head = heads[1]
        _i = i + 1
        idx = 1
        last_head = parents.lastshort[]
        if check
            next_parents = Q[head][2]
            while hasonly(next_parents.shortroots, root) || head == buffer[end-_i-haslongarc]
                if idx == last_head
                    last_positions |= (~zero(UInt64) >> ((63 - num + i) % UInt8))
                    return next_compatible_arc!(buffer, last_positions, idx_stack, Q, true)
                end
                idx += 1
                head = heads[idx]
                next_parents = Q[head][2]
            end
            parents = next_parents
        else
            parents = Q[head][2]
        end
        idx_stack[i] = idx
        buffer[check ? _i : end-_i] = head
        if idx == last_head
            last_positions |= (one(UInt64) << ((num-i) % UInt8))
        end
    end
    return last_positions
end

function _case_0(heads, lastshort, i_stop)
    return Vector{Int}[Int[1, heads[i], i_stop] for i in length(heads):-1:(lastshort+1)]
end

function cycles_ending_at(Q::Vector{Tuple{PeriodicVertex{D},JunctionNode}}, i_stop) where D
    parents = Q[i_stop][2]
    length(parents.heads) ≥ 2 || return Vector{Int}[]
    length(union(ConstMiniBitSet(parents.shortroots), ConstMiniBitSet(parents.longroots))) ≥ 2 || return Vector{Int}[]
    shortest_n = parents.num
    parents_lastshort = parents.lastshort[]
    shortest_n == 0 && return _case_0(parents.heads, parents_lastshort, i_stop)
    ret = Vector{Int}[]
    idx_stack1 = Vector{Int}(undef, shortest_n)
    idx_stack2 = Vector{Int}(undef, shortest_n-1)
    buffer = Vector{Int}(undef, 2*shortest_n + 3)
    buffer[shortest_n+1] = 1
    buffer[end] = i_stop
    for i1 in length(parents.heads):-1:2
        if i1 == parents_lastshort
            Base._deleteend!(buffer, 1)
            Base._deleteend!(idx_stack1, 1)
            buffer[end] = i_stop
        end
        buffer[end-1] = parents.heads[i1]
        last_positions1 = initial_compatible_arc!(buffer, idx_stack1, Q, false)
        while last_positions1 != ~zero(UInt64)
            for i2 in 1:min(i1-1, parents_lastshort)
                head2 = parents.heads[i2]
                i1 > parents_lastshort && buffer[end-2] == head2 && continue # 3-cycle near i_stop
                buffer[1] = head2
                last_positions2 = initial_compatible_arc!(buffer, idx_stack2, Q, true)
                while last_positions2 != ~zero(UInt64)
                    push!(ret, copy(buffer))
                    last_positions2 = next_compatible_arc!(buffer, last_positions2, idx_stack2, Q, true)
                end
            end
            last_positions1 = next_compatible_arc!(buffer, last_positions1, idx_stack1, Q, false)
        end
    end
    return ret
end


#=
function cycles_around(g, node, depth)
    Q, vertexdict = arcs_list(g, node, depth)
    cycles = Vector{PeriodicVertex3D}[]
    for (k, (x, arcs)) in enumerate(Q)
        length(arcs) ≤ 1 && continue
        len = arcs[1].num
        for (i1, arc1) in enumerate(arcs)
            (arc1.num > len || i1 == length(arcs)) && break
            arc1list = Vector{PeriodicVertex{D}}(undef, length(arc1) + 1)
            arc1list[1] = PeriodicVertex{D}(node)
            collect!(arc1list, arc1, true, 2)
            root = arc1.root
            for i2 in (i1+1):length(arcs)
                arc2 = arcs[i2]
                arc2.root == root && continue
            end
        end
    end
end
=#


function simple_bfs(g::PeriodicGraph{D}, i, depth, vertexdict, previous=nothing) where D
    seen = falses(length(vertexdict))
    seen[i isa Int ? i : vertexdict[i]] = true
    for x in neighbors(g, i)
        seen[vertexdict[x]] = true
    end
    Q::Vector{Tuple{PeriodicVertex{D}, Int}} = if previous === nothing
        [(x,1) for x in neighbors(g, i)]
    else
        _Q = Vector{Tuple{PeriodicVertex{D}, Int}}(undef, degree(g, i) - 1)
        _count = 0
        for x in neighbors(g, i)
            x == previous && continue
            _count += 1
            _Q[_count] = (x, 1)
        end
        _Q
    end
    for (u, d) in Q
        d == depth && break
        for x in neighbors(g, u)
            idx = vertexdict[x]
            seen[idx] && continue
            seen[idx] = true
            push!(Q, (x, d+1))
        end
    end
    return Q
end
