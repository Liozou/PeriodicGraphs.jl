struct ConstMiniBitSet{T} <: AbstractSet{Int}
    x::T
    global _constminibitset(x::T) where {T} = new{T}(x)
end

ConstMiniBitSet{T}() where {T} = _constminibitset(zero(T))
ConstMiniBitSet{T}(i::Integer) where {T} = _constminibitset(one(T) << (i % UInt8))
function ConstMiniBitSet{T}(l::AbstractVector{<:Integer}) where T
    x = zero(T)
    for i in l
        x |= one(T) << (i % UInt8)
    end
    return _constminibitset(x)
end

ConstMiniBitSet() = ConstMiniBitSet{UInt64}()
ConstMiniBitSet(i::Integer) = ConstMiniBitSet{UInt64}(i)
ConstMiniBitSet(l::AbstractVector{<:Integer}) = ConstMiniBitSet{UInt64}(l)

push(x::ConstMiniBitSet{T}, i::Integer) where {T} = _constminibitset(x.x | (one(T) << (i % UInt8)))
Base.union(x::ConstMiniBitSet, y::ConstMiniBitSet) = _constminibitset(x.x | y.x)
Base.intersect(x::ConstMiniBitSet, y::ConstMiniBitSet) = _constminibitset(x.x & y.x)
Base.in(x::ConstMiniBitSet{T}, i::Integer) where {T} = ((x.x >> (i % UInt8)) & one(T)) % Bool
Base.setdiff(x::ConstMiniBitSet, y::ConstMiniBitSet) = _constminibitset(x.x & ~y.x)
Base.symdiff(x::ConstMiniBitSet, y::ConstMiniBitSet) = _constminibitset(x.x ⊻ y.x)
Base.length(x::ConstMiniBitSet) = count_ones(x.x)
Base.isempty(x::ConstMiniBitSet) = iszero(x.x)
Base.eltype(::Type{ConstMiniBitSet}) = Int

function Base.iterate(::ConstMiniBitSet{T}, x::T) where T
    first = trailing_zeros(x)
    first == 8*sizeof(T) && return nothing
    return (first, x & (~(one(T) << (first % UInt8))))
end
Base.iterate(x::ConstMiniBitSet) = iterate(x, x.x)
Base.isdone(::ConstMiniBitSet{T}, x::T) where {T} = iszero(x)
Base.isdone(x::ConstMiniBitSet) = Base.isdone(x, x.x)

# minabove(x::ConstMiniBitSet, i::Integer) = trailing_zeros(x.x & ((~one(UInt64)) << (i % UInt8)))
# minaboveeq(x::ConstMiniBitSet, i::Integer) = trailing_zeros(x.x & ((~zero(UInt64)) << (i % UInt8)))
# nonemptyintersect(x::ConstMiniBitSet, y::ConstMiniBitSet) = !iszero(x.x & y.y)
hasonly(x::ConstMiniBitSet{T}, i::Integer) where {T} = iszero(x.x & ~(one(T) << (i % UInt8)))

Base.minimum(x::ConstMiniBitSet) = first(x)
Base.maximum(x::ConstMiniBitSet{T}) where {T} = 8*sizeof(T) - leading_zeros(x.x)

function Base.show(io::IO, x::ConstMiniBitSet)
    print(io, ConstMiniBitSet, "([")
    join(io, x, ", ")
    print(io, "])")
end


# mutable struct MiniBitSet <: AbstractSet{Int}
#     x::UInt64
#     global _minibitset(x::UInt64) = new(x)
# end

# MiniBitSet(x::ConstMiniBitSet) = _minibitset(x.x)
# ConstMiniBitSet(x::MiniBitSet) = _constminibitset(x.x)

# MiniBitSet() = _minibitset(zero(UInt64))
# MiniBitSet(i::Integer) = _minibitset(one(UInt64) << (i % UInt8))
# MiniBitSet(l::AbstractVector{<:Integer}) = MiniBitSet(ConstMiniBitSet(l))

# Base.push!(x::MiniBitSet, i::Integer) = (x.x |= (one(UInt64) << (i % UInt8)); x)
# Base.union!(x::MiniBitSet, y::MiniBitSet) = (x.x |= y.x; x)
# Base.intersect!(x::MiniBitSet, y::MiniBitSet) = (x.x &= y.x; x)
# Base.in(x::MiniBitSet, i::Integer) = ((x.x >> (i % UInt8)) & one(UInt64)) % Bool
# Base.setdiff!(x::MiniBitSet, y::MiniBitSet) = (x.x &= ~y.x; x)
# Base.symdiff!(x::MiniBitSet, y::MiniBitSet) = (x.x ⊻= y.x; x)
# Base.copymutable(x::MiniBitSet) = _minibitset(x.x)
# Base.copy(x::MiniBitSet) = _minibitset(x.x)

# Base.iterate(x::MiniBitSet, u::UInt64=x.x) = iterate(ConstMiniBitSet(), u)
# Base.isdone(::MiniBitSet, x::UInt64) = iszero(x)
# Base.isdone(x::MiniBitSet) = Base.isdone(x, x.x)

# # minabove(x::MiniBitSet, i::Integer) = minabove(ConstMiniBitSet(x), i)
# # minaboveeq(x::MiniBitSet, i::Integer) = minaboveeq(ConstMiniBitSet(x), i)
# # isnonemptyintersect(x::MiniBitSet, y::MiniBitSet) = !iszero(x.x & y.x)
# # isnonemptysymdiff(x::MiniBitSet, y::MiniBitSet) = !iszero(x.x ⊻ y.x)
# # hasonly(x::MiniBitSet, i::Integer) = iszero(x.x & ~(one(UInt64) << (i % UInt8)))
# # hasotherthan(x::MiniBitSet, i::Integer) = !iszero(x.x & ~(one(UInt64) << (i % UInt8)))

# Base.length(x::MiniBitSet) = length(ConstMiniBitSet(x))
# Base.eltype(x::MiniBitSet) = eltype(ConstMiniBitSet(x))
# Base.maximum(x::MiniBitSet) = maximum(ConstMiniBitSet(x))
# Base.minimum(x::MiniBitSet) = minimum(ConstMiniBitSet(x))
# Base.isempty(x::MiniBitSet) = isempty(ConstMiniBitSet(x))

# function Base.show(io::IO, x::MiniBitSet)
#     print(io, MiniBitSet, "([")
#     join(io, x, ", ")
#     print(io, "])")
# end


const SmallIntType = Int8 # used for distances mostly

struct JunctionNode
    shortroots::ConstMiniBitSet{UInt32}
    longroots::ConstMiniBitSet{UInt32}
    lastshort::SmallIntType # (use SmallIntType for better padding)
    num::SmallIntType # length of the shortest branch
    heads::Vector{Int}
    function JunctionNode(heads::AbstractVector, num, roots::ConstMiniBitSet)
        new(roots, ConstMiniBitSet{UInt32}(), length(heads) % SmallIntType, num % SmallIntType, heads)
    end
end
JunctionNode() = JunctionNode(Int[], -one(SmallIntType), ConstMiniBitSet{UInt32}())
JunctionNode(head::Integer, num, roots::ConstMiniBitSet) = JunctionNode(Int[head], num, roots)
function JunctionNode(root::Integer, skip::Bool)
    if skip
        JunctionNode(Int[], root, ConstMiniBitSet{UInt32}())
    else
        JunctionNode(1, zero(SmallIntType), ConstMiniBitSet{UInt32}(root))
    end
end

function Base.show(io::IO, x::JunctionNode)
    print(io, "JunctionNode([")
    join(io, x.heads, ", ")
    print(io, "])")
end


function prepare_avoid!(dag, vertexnums, g, i, hintsize)
    neighs = neighbors(g, i)
    n = length(neighs) + 1
    todelete = falses(n)
    finaldagoffset = zeros(Int, n)
    sizehint!(finaldagoffset, hintsize)
    finaldagoffset[1] = 1
    offsetcounter = 1
    next_toavoid = n + 1
    resize!(dag, n)
    resize!(vertexnums, n)
    #=@inbounds=# for j in 1:n-1
        x = neighs[j]
        if iszero(x.ofs) && avoid[x.v]
            next_toavoid -= 1
            dag[next_toavoid] = JunctionNode(zero(SmallIntType), true)
            todelete[next_toavoid] = true
            vertexnums[next_toavoid] = x
        else
            offsetcounter += 1
            dag[offsetcounter] = JunctionNode(offsetcounter, false)
            finaldagoffset[offsetcounter] = offsetcounter
            vertexnums[offsetcounter] = x
        end
    end
    return todelete, finaldagoffset, offsetcounter, next_toavoid
end

function update_offsets!(dag, finaldagoffset, start, val)
    val == 0 && return nothing
    #=@inbounds=# for j in start:length(dag)
        finaldagoffset[j] -= val
        heads = dag[j].heads
        for (i, h) in enumerate(heads)
            if h ≥ start
                heads[i] = h - val
            end
        end
    end
    nothing
end

const shortrootsoffset = fieldoffset(JunctionNode, findfirst(==(:shortroots), fieldnames(JunctionNode)))
const longrootsoffset = fieldoffset(JunctionNode, findfirst(==(:longroots), fieldnames(JunctionNode)))
const lastshortoffset = fieldoffset(JunctionNode, findfirst(==(:lastshort), fieldnames(JunctionNode)))
@inline unsafe_incr!(ptr::Ptr{T}) where {T} = unsafe_store!(ptr, unsafe_load(ptr) + one(T))
@inline unsafe_union!(ptr, val) = unsafe_store!(ptr, _constminibitset(val | unsafe_load(ptr).x))

"""
    arcs_list(g::PeriodicGraph{D}, i, depth, avoid)

Compute the list of shortest arcs starting from vertex `i` up to length `depth+1` that do
not pass through a vertex in `avoid`.

`avoid` can be set to `nothing` to disregard that condition.

Return `(dag, vertexnums)` where `dag` is a `Vector{JunctionNode}` representing, for each
visited node, the dag of all arcs from that node back to `i`, in a compact representation.
`vertexnums` is a `Vector{PeriodicVertex{D}}` whose `k`-th value is the vertex represented
by number `k` in the `dag`.
"""
function arcs_list(g::PeriodicGraph{D}, i, depth, avoid) where D
    vertexnums = PeriodicVertex{D}[PeriodicVertex{D}(i)]
    dag = JunctionNode[JunctionNode()]
    hintsize = ceil(Int, depth^2.9) # empirical estimate
    sizehint!(vertexnums, hintsize)
    sizehint!(dag, hintsize)
    if length(vertexnums) > 62
        error("The vertex has a degree too large (> 62) to be represented in a MiniBitSet. Please open an issue.")
    end
    hasavoid = !(avoid isa Nothing)
    #= TODO: redo logic for hasavoid:
    - when handling nodes at distance `num`, `dag` contains valid nodes, then dead nodes,
      then nodes at distance `num+1`
    - as long as we are handling valid nodes, the resulting created nodes are either valid
      or explicitly in `avoid` -> store those in a separate list.
    - when reaching the last valid node at distance `num`, `dag` only contains dead node at
      distance `num`, then valid nodes at distance `num+1`. Append to `dag` the separate
      list of `dead` nodes at distance `num+1`
    - when handling dead nodes, only modify `dag` by appending new dead nodes at distance
      `num+1`. No valid node at distance `num+1` can be created like so, they would already
      be in `dag` since they require a valid neighbour at distance `num`, already handled.
    =#
    if hasavoid
        todelete, finaldagoffset, offsetcounter, next_toavoid = prepare_avoid!(dag, vertexnums, g, i, hintsize)
        next_num = length(dag) + 1
    else
        append!(vertexnums, neighbors(g, i))
        append!(dag, JunctionNode(j, false) for j in 2:length(vertexnums))
    end
    vertexdict = Dict{PeriodicVertex{D},Int}((x => j) for (j,x) in enumerate(vertexnums))
    counter = 1
    _depth = depth % SmallIntType
    next_next_toavoid = 0
    removed_avoid = 0
    for parents in Iterators.rest(dag, 2)
        parents.num == _depth && break
        counter += 1
        if hasavoid
            if counter == next_toavoid
                next_next_toavoid = length(dag)+1
                removed_avoid = 0
            end
            if counter == next_num
                next_num = length(dag) + 1
                next_toavoid = next_next_toavoid
                next_next_toavoid = 0
                update_offsets!(dag, finaldagoffset, counter, removed_avoid)
            end
            if isempty(parents.heads)
                @assert next_next_toavoid != 0
                removed_avoid += 1
                todelete[counter] = true
                push!(parents.heads, 1)
                finaldagoffset[counter] = 0
            end
        end
        previous = vertexnums[parents.heads[1]]
        newhead = hasavoid ? finaldagoffset[counter] : counter
        num = parents.num + one(SmallIntType)
        for x in neighbors(g, vertexnums[counter])
            x == previous && continue
            skip = hasavoid && avoid[x.v] && iszero(x.ofs)
            idx = get!(vertexdict, x, length(dag)+1)
            if idx == length(dag)+1
                # push!(dag, skip ? JunctionNode(Int[], num, ConstMiniBitSet{UInt32}()) :
                #                   JunctionNode(counter, num, parents.shortroots))
                push!(vertexnums, x)
                if hasavoid
                    push!(finaldagoffset, skip ? 0 : (offsetcounter += 1))
                    push!(todelete, skip)
                    if newhead == 0
                        push!(dag, JunctionNode(num, true))
                    else
                        push!(dag, JunctionNode(newhead, num, parents.shortroots))
                    end
                else
                    push!(dag, JunctionNode(counter, num, parents.shortroots))
                end
            elseif idx > counter && !skip
                junction = dag[idx]
                hasavoid && (todelete[idx] || newhead == 0) && continue
                push!(junction.heads, newhead)
                # We use pointer and unsafe_store! to modify the immutable type JunctionNode
                # This is only possible because it is stored in a Vector, imposing fixed addresses.
                ptr = pointer(dag, idx)
                if num == junction.num
                    unsafe_incr!(Ptr{SmallIntType}(ptr + lastshortoffset))
                    unsafe_union!(Ptr{ConstMiniBitSet{UInt32}}(ptr + shortrootsoffset), parents.shortroots.x)
                else
                    unsafe_union!(Ptr{ConstMiniBitSet{UInt32}}(ptr + longrootsoffset), parents.shortroots.x)
                end
            end
        end
    end
    if hasavoid
        deleteat!(dag, todelete)
        deleteat!(vertexnums, todelete)
    end
    return dag, vertexnums
end


function next_compatible_arc!(buffer, last_positions, idx_stack, dag, check, dist, vertexnums)
    num_choice = length(idx_stack)
    updated_last_positions = last_positions
    haslongarc = isodd(length(buffer))
    num = (num_choice + 2 - haslongarc) % SmallIntType
    root = buffer[num+1] # unused if !check

    @label loop
    trail = trailing_ones(updated_last_positions) % UInt8
    next_toupdate = num_choice - trail
    next_toupdate < 1 && return ~zero(UInt64)
    updated_last_positions = updated_last_positions & (~zero(UInt64) << trail)
    idx = idx_stack[next_toupdate] + 1
    head = buffer[check ? next_toupdate : end-next_toupdate]
    parents = dag[head]
    for i in next_toupdate:num_choice
        heads = parents.heads
        head = heads[idx]
        _i = i + 1
        last_head = parents.lastshort % Int
        if check
            next_parents = dag[head]
            while hasonly(next_parents.shortroots, root) ||
                    head == buffer[end-_i-haslongarc] ||
                    (dist !== nothing && # checking for rings instead of cycles
                    (is_distance_smaller!(dist, vertexnums[buffer[num+_i]], vertexnums[head], num) ||
                    (haslongarc && is_distance_smaller!(dist, vertexnums[buffer[num+_i+1]], vertexnums[head], num))))
                if idx == last_head
                    updated_last_positions |= (~zero(UInt64) >> ((63 - num_choice + i) % UInt8))
                    @goto loop
                end
                idx += 1
                head = heads[idx]
                next_parents = dag[head]
            end
            parents = next_parents
        else
            parents = dag[head]
        end
        buffer[check ? _i : end-_i] = head
        idx_stack[i] = idx
        if idx == last_head
            updated_last_positions |= (one(UInt64) << ((num_choice - i) % UInt8))
        end
        idx = 1
    end
    return updated_last_positions
end

function initial_compatible_arc!(buffer, idx_stack, dag, check, dist, vertexnums)
    last_positions = zero(UInt64)
    num_choice = length(idx_stack)
    head = buffer[check ? 1 : end-1]
    haslongarc = isodd(length(buffer))
    num = (num_choice + 2 - haslongarc) % SmallIntType
    root = buffer[num+1] # unused if !check
    parents = dag[head]
    for i in 1:num_choice
        heads = parents.heads
        head = heads[1]
        _i = i + 1
        idx = 1
        last_head = parents.lastshort % Int
        if check
            next_parents = dag[head]
            while hasonly(next_parents.shortroots, root) ||
                  head == buffer[end-_i-haslongarc] ||
                  (dist !== nothing && # checking for rings instead of cycles
                  (is_distance_smaller!(dist, vertexnums[buffer[num+_i]], vertexnums[head], num) ||
                  (haslongarc && is_distance_smaller!(dist, vertexnums[buffer[num+_i+1]], vertexnums[head], num))))
                if idx == last_head
                    last_positions |= (~zero(UInt64) >> ((63 - num_choice + i) % UInt8))
                    return next_compatible_arc!(buffer, last_positions, idx_stack, dag, true, dist, vertexnums)
                end
                idx += 1
                head = heads[idx]
                next_parents = dag[head]
            end
            parents = next_parents
        else
            parents = dag[head]
        end
        idx_stack[i] = idx
        buffer[check ? _i : end-_i] = head
        if idx == last_head
            last_positions |= (one(UInt64) << ((num_choice-i) % UInt8))
        end
    end
    return last_positions
end

function _case_0(heads, lastshort, midnode)
    return Vector{Int}[Int[1, heads[i], midnode] for i in length(heads):-1:(lastshort+1)]
end

function cycles_ending_at(dag::Vector{JunctionNode}, midnode, dist=nothing, vertexnums=nothing)
    parents = dag[midnode]
    length(parents.heads) ≥ 2 || return Vector{Int}[]
    length(union(parents.shortroots, parents.longroots)) ≥ 2 || return Vector{Int}[]
    shortest_n = parents.num % Int
    num = parents.num + one(SmallIntType)
    parents_lastshort = parents.lastshort % Int
    iszero(shortest_n) && return _case_0(parents.heads, parents_lastshort, midnode)
    ret = Vector{Int}[]
    idx_stack1 = Vector{Int}(undef, shortest_n)
    idx_stack2 = Vector{Int}(undef, shortest_n-1)
    buffer = Vector{Int}(undef, 2*shortest_n + 3)
    buffer[shortest_n+1] = 1
    buffer[end] = midnode
    for i1 in length(parents.heads):-1:2
        if i1 == parents_lastshort
            Base._deleteend!(buffer, 1)
            Base._deleteend!(idx_stack1, 1)
            buffer[end] = midnode
        end
        buffer[end-1] = parents.heads[i1]
        last_positions1 = initial_compatible_arc!(buffer, idx_stack1, dag, false, dist, vertexnums)
        while last_positions1 != ~zero(UInt64)
            for i2 in 1:min(i1-1, parents_lastshort)
                head2 = parents.heads[i2]
                i1 > parents_lastshort && buffer[end-2] == head2 && continue # 3-cycle near midnode
                buffer[1] = head2
                if dist !== nothing && # checking for rings instead of cycles
                  (is_distance_smaller!(dist, vertexnums[buffer[shortest_n+2]], vertexnums[head2], num) ||
                  (i1 > parents_lastshort && is_distance_smaller!(dist, vertexnums[buffer[shortest_n+3]], vertexnums[head2], num)))
                    continue
                end
                last_positions2 = initial_compatible_arc!(buffer, idx_stack2, dag, true, dist, vertexnums)
                while last_positions2 != ~zero(UInt64)
                    push!(ret, copy(buffer))
                    last_positions2 = next_compatible_arc!(buffer, last_positions2, idx_stack2, dag, true, dist, vertexnums)
                end
            end
            last_positions1 = next_compatible_arc!(buffer, last_positions1, idx_stack1, dag, false, dist, vertexnums)
        end
    end
    return ret
end


"""
    DistanceRecord{D}

Record of the computed distances between vertices of a graph.
"""
struct DistanceRecord{D}
    g::PeriodicGraph{D}
    shortdist::Matrix{SmallIntType}
    distances::Dict{Tuple{Int,Int},SmallIntType}
    Q_lists::Vector{Vector{Tuple{PeriodicVertex{D},SmallIntType}}}
    first_indices::Vector{Int}
    seenshort::BitMatrix
    seens::Vector{Set{Int}}
end
function DistanceRecord(g::PeriodicGraph{D}, depth) where D
    n = nv(g)
    dim = n*(1 + fld(depth, 2))^D
    shortdist = zeros(SmallIntType, n, dim)
    distances = Dict{Tuple{Int,Int},Int}()
    Q_lists = Vector{Vector{Tuple{PeriodicVertex{D},SmallIntType}}}(undef, n)
    first_indices = zeros(Int, n)
    seenshort = falses(n, dim)
    seens = Vector{Set{Int}}(undef, n)
    return DistanceRecord{D}(g, shortdist, distances, Q_lists, first_indices, seenshort, seens)
end

function known_distance(dist::DistanceRecord, i, j)
    a, b = minmax(i,j)
    if checkbounds(Bool, dist.shortdist, a, b)
        @inbounds dist.shortdist[a,b]
    else
        get(dist.distances, (a, b), zero(SmallIntType))
    end
end
function set_distance!(dist::DistanceRecord, i, j, d)
    a, b = minmax(i,j)
    if checkbounds(Bool, dist.shortdist, a, b)
        @inbounds dist.shortdist[a,b] = d % SmallIntType
    else
        dist.distances[(a,b)] = d % SmallIntType
    end
end

function has_been_seen!(dist, i, h)
    if h ≤ size(dist.seenshort)[2]
        @inbounds dist.seenshort[i,h] && return true
        @inbounds dist.seenshort[i,h] = true
    else
        seen = @inbounds dist.seens[i]
        h ∈ seen && return true
        push!(seen, h)
    end
    false
end

struct CyclesEndingAt{T}
    dag::Vector{JunctionNode}
    midnode::Int
    heads::Vector{Int}
    lastshort::Int
    parentsnum::SmallIntType
    record::T
end

"""
    CyclesEndingAt(dag, midnode, record=(nothing,nothing))

Iterable over the cycles of graph `g` around node `i` with `midnode` as vertex furthest
from `i`. `dag` is obtained from `first(arcs_list(g, i, ...))`.

If `record` is set to `(dist, vertexnums)` where `dist == DistanceRecord(g, depth)` and
`dag, vertexnums == first(arcs_list(g, i, depth, ...))`, this will iterate over the rings
instead of all cycles.
"""
@inline function CyclesEndingAt(dag, midnode, record=(nothing,nothing))
    parents = dag[midnode]
    CyclesEndingAt{typeof(record)}(dag, midnode, parents.heads, parents.lastshort % Int, parents.num, record)
end
Base.IteratorSize(::CyclesEndingAt) = Base.SizeUnknown()
Base.eltype(::CyclesEndingAt) = Vector{Int}

function Base.iterate(x::CyclesEndingAt, state=nothing)
    heads = x.heads
    length(heads) ≥ 2 || return nothing
    midnode = x.midnode
    lastshort = x.lastshort
    shortest_n = x.parentsnum % Int
    if iszero(shortest_n)
        length(heads) > lastshort || return nothing
        if state === nothing
            state = (Int[1, length(heads)+1, midnode], Int[], Int[], zero(UInt64), zero(UInt64), 0, 0)
        end
        i0 = @inbounds state[1][2] -= 1
        i0 > lastshort || return nothing
        return state[1], state
    end
    dag = x.dag
    dist, vertexnums = x.record
    num = x.parentsnum + one(SmallIntType)
    state === nothing || @goto next

    length(union(dag[midnode].shortroots, dag[midnode].longroots)) ≥ 2 || return nothing
    idx_stack1 = Vector{Int}(undef, shortest_n)
    idx_stack2 = Vector{Int}(undef, shortest_n-1)
    buffer = Vector{Int}(undef, 2*shortest_n + 3)
    buffer[shortest_n+1] = 1
    buffer[end] = midnode
    i1 = length(heads)
    while i1 > 1
        if i1 == lastshort
            Base._deleteend!(buffer, 1)
            Base._deleteend!(idx_stack1, 1)
            buffer[end] = midnode
        end
        buffer[end-1] = heads[i1]
        last_positions1 = initial_compatible_arc!(buffer, idx_stack1, dag, false, dist, vertexnums)
        while last_positions1 != ~zero(UInt64)
            i2 = min(i1, lastshort+1)
            while i2 > 1
                i2 -= 1
                head2 = heads[i2]
                i1 > lastshort && buffer[end-2] == head2 && continue # 3-cycle near midnode
                buffer[1] = head2
                if dist !== nothing && # checking for rings instead of cycles
                  (is_distance_smaller!(dist, vertexnums[buffer[shortest_n+2]], vertexnums[head2], num) ||
                  (i1 > lastshort && is_distance_smaller!(dist, vertexnums[buffer[shortest_n+3]], vertexnums[head2], num)))
                    continue
                end
                last_positions2 = initial_compatible_arc!(buffer, idx_stack2, dag, true, dist, vertexnums)
                while last_positions2 != ~zero(UInt64)
                    return (buffer, (buffer, idx_stack1, idx_stack2, last_positions1, last_positions2, i1, i2))

                    @label next
                    @inbounds buffer, idx_stack1, idx_stack2, last_positions1, last_positions2, i1, i2 = state
                    last_positions2 = next_compatible_arc!(buffer, last_positions2, idx_stack2, dag, true, dist, vertexnums)
                end
            end
            last_positions1 = next_compatible_arc!(buffer, last_positions1, idx_stack1, dag, false, dist, vertexnums)
        end
        i1 -= 1
    end
    nothing
end

function bfs_smaller!(dist, i, j, start, stop)
    Q = dist.Q_lists[i]
    graph = dist.g
    counter = start
    encountered = false
    n = nv(graph)
    for (u, d) in Iterators.rest(Q, start)
        (encountered || d+one(SmallIntType) ≥ stop) && break
        counter += 1
        for x in neighbors(graph, u)
            h = hash_position(x, n)
            has_been_seen!(dist, i, h) && continue
            push!(Q, (x, d+1))
            set_distance!(dist, i, h, d+one(SmallIntType))
            encountered |= h == j
        end
    end
    dist.first_indices[i] = counter
    if !encountered
        set_distance!(dist, i, j, stop)
        return false
    end
    return true
end

function _reorderinit(first_indices, verti, vertj::PeriodicVertex{D}) where D
    n = length(first_indices)
    counterij = (first_indices[verti.v], first_indices[vertj.v])
    if counterij[1] < counterij[2]
        return vertj.v, hash_position(PeriodicVertex{D}(verti.v, verti.ofs .- vertj.ofs), n), counterij[2]
    end
    return verti.v, hash_position(PeriodicVertex{D}(vertj.v, vertj.ofs .- verti.ofs), n), counterij[1]
end

function init_distance_record!(dist, i, j)
    dist.first_indices[i] = 1
    g = dist.g
    dist.Q_lists[i] = [(x, one(SmallIntType)) for x in neighbors(g, i)]
    _seen = Set{Int}(i)
    encountered = false
    n = nv(g)
    for x in neighbors(g, i)
        _h = hash_position(x, n)
        set_distance!(dist, i, _h, one(SmallIntType))
        if _h ≤ (@inbounds size(dist.seenshort)[2])
            dist.seenshort[i,_h] = true
        else
            push!(_seen, _h)
        end
        encountered |= _h == j
    end
    dist.seens[i] = _seen
    encountered
end

function is_distance_smaller!(dist, verti, vertj, stop)
    i, j, start = _reorderinit(dist.first_indices, verti, vertj)
    known = known_distance(dist, i, j)
    iszero(known) || return known < stop
    if start == 0
        init_distance_record!(dist, i, j) && return one(stop) < stop
        start = 1
    end
    return bfs_smaller!(dist, i, j, start, stop)
end

"""
    normalize_cycle!(cycle)

In-place rotate and possibly reverse `cycle` so that `cycle[1] == minimum(cycle)` and
`cycle[2] < cycle[end]`.
"""
function normalize_cycle!(cycle)
    fst = argmin(cycle)
    lenc = length(cycle)
    endreverse = cycle[mod1(fst+1, lenc)] < cycle[mod1(fst-1, lenc)]
    if fst != 1 || !endreverse
        reverse!(cycle, 1, fst - endreverse)
        reverse!(cycle, fst - endreverse + 1, lenc)
        endreverse && reverse!(cycle)
    end
    cycle
end

function rings_around_old(g::PeriodicGraph{D}, i, depth=15, visited=nothing) where D
    n = nv(g)
    dist = DistanceRecord(g, depth)
    dag, vertexnums = arcs_list(g, i, depth, visited)
    hashes = [hash_position(x, n) for x in vertexnums]
    ret = Vector{Int}[]
    for midnode in 2:length(dag)
        for cycle in CyclesEndingAt(dag, midnode)
        # for cycle in cycles_ending_at(dag, midnode)
            isempty(cycle) && break
            invalidcycle = false
            lenc = length(cycle)
            num = (lenc ÷ 2) % SmallIntType
            if isodd(lenc)
                for i in 1:(num-1)
                    verti = vertexnums[cycle[i]]
                    vertj1 = vertexnums[cycle[num+i]]
                    vertj2 = vertexnums[cycle[num+i+1]]
                    if is_distance_smaller!(dist, verti, vertj1, num) ||
                        is_distance_smaller!(dist, verti, vertj2, num)
                        invalidcycle = true
                        break
                    end
                end
            else
                for i in 1:(num-1)
                    verti = vertexnums[cycle[i]]
                    vertj = vertexnums[cycle[num+i]]
                    if is_distance_smaller!(dist, verti, vertj, num)
                        invalidcycle = true
                        break
                    end
                end
            end
            invalidcycle && continue
            push!(ret, normalize_cycle!(hashes[cycle]))
        end
    end
    return ret
end

"""
    rings_around(g::PeriodicGraph{D}, i, depth=15) where D

Return the list of all rings around node `i` in graph `g` up to length `2*depth+3`.

The returned rings are the list of `hash_position` of the corresponding vertices. To get
back the list of actual `PeriodicVertex` of a returned `ring` in the list, do
```julia
[reverse_hash_position(x, n, Val(D)) for x in ring]
```
where `n == nv(g)`. If the offsets of the corresponding vertices are not needed, simply do
```julia
[mod1(x, n) for x in ring]
```
"""
function rings_around(g::PeriodicGraph{D}, i, depth=15, dist=DistanceRecord(g,depth), visited=nothing) where D
    n = nv(g)
    dag, vertexnums = arcs_list(g, i, depth, visited)
    hashes = [hash_position(x, n) for x in vertexnums]
    ret = Vector{Int}[]
    for midnode in length(dag):-1:2
        for cycle in CyclesEndingAt(dag, midnode, (dist, vertexnums))
        # for cycle in cycles_ending_at(dag, midnode, dist, vertexnums)
            newcycle = normalize_cycle!(hashes[cycle])
            if newcycle == [5, 12, 16, 6, 15, 143, 134, 144, 20, 187, 21, 13]
                @show i, midnode, cycle
                @show hashes[cycle]
            end
            push!(ret, newcycle)
        end
    end
    return ret
end

cycles_around(g::PeriodicGraph, i, depth=15) = rings_around(g, i, depth, nothing)


struct NoSymmetry
    num::Int
    NoSymmetry(g::PeriodicGraph) = new(nv(g))
end
Base.getindex(::NoSymmetry, i::Integer) = i
Base.unique(x::NoSymmetry) = Base.OneTo(x.num)
Base.iterate(::NoSymmetry) = nothing


"""
    no_neighboring_nodes(g, symmetries)

Return a list of nodes representatives (modulo symmetries) such that all the vertices in
`g` either have their representative in the list or are surrounded by nodes whose
representatives are in the list.

This means that all cycles and rings of `g` have at least one node whose representative is
in the list.
"""
function no_neighboring_nodes(g, symmetries)
    n = nv(g)
    # category: 0 => unvisited, 1 => to explore, 2 => unvisited unknown, 3 => not to explore
    categories = zeros(Int8, n)
    categories[1] = 2
    toexplore = Int[]
    u_init = 1
    while u_init !== nothing
        Q = [u_init]
        for u in Q
            newcat = categories[u] == 2 ? UInt8(1) : UInt(2)
            if newcat == 1
                if any(x -> symmetries[x.v] == u, neighbors(g, u))
                    categories[u] = 1
                    push!(toexplore, u)
                    newcat = 2
                else
                    categories[u] = 3
                end
            end
            for x in neighbors(g, u)
                v = symmetries[x.v]
                cat = categories[v]
                if iseven(cat)
                    if cat == 0 && iszero(x.ofs)
                        push!(Q, v)
                    end
                    categories[v] = newcat
                    newcat == 1 && push!(toexplore, v)
                end
            end
        end
        u_init = findfirst(x -> categories[x] == 0, unique(symmetries))
    end
    return toexplore
end

function rings(g::PeriodicGraph{D}, depth=15, symmetries=NoSymmetry(g), dist=DistanceRecord(g,depth)) where D
    toexplore = no_neighboring_nodes(g, symmetries)
    allrings = Set{Vector{Int}}()
    visited = falses(nv(g))
    for x in toexplore
        newrings = rings_around(g, x, depth, dist, visited)
        union!(allrings, newrings)
        visited[x] = true
        for symm in symmetries
            visited[symm[x]] = true
        end
    end
    ret = collect(allrings)
    symmetries isa NoSymmetry && return ret
    sort!(ret; by=length, rev=true) # to avoid resizing the buffer below too many times
    buffer = similar(first(ret))
    for symm in symmetries
        for ring in allrings
            nr = length(ring)
            resize!(buffer, nr)
            #=@inbounds=# for i in 1:nr
                buffer[i] = symm[ring[x]]
            end
            normalize_cycle!(buffer)
            if buffer ∉ allrings
                push!(ret, copy(buffer))
            end
        end
    end
    return ret
end

cycles(g::PeriodicGraph, depth=15, symmetries=NoSymmetry(g)) = rings(g, depth, symmetries, nothing)
