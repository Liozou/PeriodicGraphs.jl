# Ring statistics

export rings, strong_rings, RingAttributions


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
JunctionNode(num::Integer, ::Nothing) = JunctionNode(Int[], num, ConstMiniBitSet{UInt32}())
JunctionNode(root::Integer) = JunctionNode(1, zero(SmallIntType), ConstMiniBitSet{UInt32}(root))
# function JunctionNode(root::Integer, skip::Bool)
#     if skip
#         JunctionNode(Int[], root, ConstMiniBitSet{UInt32}())
#     else
#         JunctionNode(1, zero(SmallIntType), ConstMiniBitSet{UInt32}(root))
#     end
# end

function Base.show(io::IO, x::JunctionNode)
    print(io, "JunctionNode([")
    join(io, x.heads, ", ")
    print(io, "])")
end

struct PhantomJunctionNode{D}
    self::PeriodicVertex{D}
    parent::PeriodicVertex{D}
    num::SmallIntType
end

function prepare_phantomdag!(dag, vertexnums, vertexdict, g::PeriodicGraph{D}, i, avoid, cycleavoid) where D
    phantomdag = PhantomJunctionNode{D}[]
    for x in neighbors(g, i)
        !(cycleavoid isa Nothing) && cycleavoid[x.v] && iszero(x.ofs) && continue
        if avoid[x.v]
            push!(phantomdag, PhantomJunctionNode(x, PeriodicVertex{D}(i), zero(SmallIntType)))
            vertexdict[x] = -length(phantomdag)
        else
            push!(vertexnums, x)
            push!(dag, JunctionNode(length(vertexnums)))
        end
    end
    return phantomdag
end

function handle_phantomdag!(phantomdag, vertexdict, g, last_stop, next_stop)
    counter = last_stop-1
    for node in Iterators.rest(phantomdag, last_stop)
        counter += 1
        counter == next_stop && return
        num = node.num + one(SmallIntType)
        for x in neighbors(g, node.self)
            x == node.parent && continue
            idx = get!(vertexdict, x, -length(phantomdag)-1)
            if idx == -length(phantomdag)-1
                push!(phantomdag, PhantomJunctionNode(x, node.self, num))
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
    arcs_list(g::PeriodicGraph{D}, i, depth, ringavoid=nothing, cycleavoid=nothing)

Compute the list of shortest arcs starting from vertex `i` up to length `depth+1`.
Vertices in `ringavoid` are not included in the returned arcs, but considered still part of
the graph for distance computations.
Vertices in `cycleavoid` are considered removed from the graph completely.

Return `(dag, vertexnums)` where `dag` is a `Vector{JunctionNode}` representing, for each
visited node, the dag of all arcs from that node back to `i`, in a compact representation.
`vertexnums` is a `Vector{PeriodicVertex{D}}` whose `k`-th value is the vertex represented
by number `k` in the `dag`.

If `ringavoid !== nothing`, this list will not include arcs that pass through nodes of the
form `PeriodicVertex{D}(j, ofs)` with `j == i` and `ofs < zero(SVector{Int,D})`: this
allows eagerly pruning cycles that are translations of others. Note that this can result in
missing cycles if those pass through at least three nodes with `j == i`, but that situation
should be rare.
"""
function arcs_list(g::PeriodicGraph{D}, i, depth, ringavoid=nothing, cycleavoid=nothing) where D
    dag = JunctionNode[JunctionNode()]
    vertexnums = PeriodicVertex{D}[PeriodicVertex{D}(i)]
    vertexdict = Dict{PeriodicVertex{D},Int}()
    hintsize = ceil(Int, depth^2.9) # empirical estimate
    sizehint!(vertexnums, hintsize)
    sizehint!(dag, hintsize)
    sizehint!(vertexdict, hintsize)
    if length(vertexnums) > 62
        error("The vertex has a degree too large (> 62) to be represented in a MiniBitSet. Please open an issue.")
    end
    hasavoid = !(ringavoid isa Nothing)
    #= Logic with hasavoid:
    A node `x` is considered dead when they `avoid[x] == true`. Moreover, any node at
    distance `num` from `i` that cannot be reached from `i` through `num` valid vertices
    is also considered dead. Those dead nodes are stored in `phantomdag`, while the rest
    ("valid" nodes) are stored in `dag`.
    - when handling valid nodes at distance `num`, `dag` contains valid nodes at distance
      `num`, then valid nodes at distance `num+1`. `phantomdag` only contains dead nodes at
      distance `num+1`.
    - as long as we are handling valid nodes, the resulting created nodes are either valid
      (-> store them in `dag`) or explicitly in `avoid` (-> store those in `phantomdag`).
    - when reaching the last valid node at distance `num`, `dag` only contains valid nodes
      at distance `num+1` while `phantomdag` contains all dead nodes at distance `num` and
      some at distance `num+1`. Handle all the dead nodes at distance `num` by disregarding
      all alive neighbors and only appending new dead nodes at distance `num+1` to
      `phantomdag`.
    - once all dead nodes at distance `num` have been handled, deal with valid nodes at
      distance `num+1`, etc.
    =#
    if hasavoid
        phantomdag = prepare_phantomdag!(dag, vertexnums, vertexdict, g, i, ringavoid, cycleavoid)
        next_stop = length(dag) + 1
        last_stop_phantomdag = 1
        next_stop_phantomdag = length(phantomdag) + 1
    else
        if cycleavoid isa Nothing
            append!(vertexnums, neighbors(g, i))
        else
            append!(vertexnums, x for x in neighbors(g, i) if !cycleavoid[x.v])
        end
        append!(dag, JunctionNode(j) for j in 2:length(vertexnums))
    end
    for (j, x) in enumerate(vertexnums)
        vertexdict[x] = j
    end
    counter = 1
    _depth = depth % SmallIntType
    for parents in Iterators.rest(dag, 2)
        counter += 1
        if hasavoid 
            if counter == next_stop
                parents.num == _depth && break
                next_stop = length(dag) + 1
                handle_phantomdag!(phantomdag, vertexdict, g, last_stop_phantomdag, next_stop_phantomdag)
                last_stop_phantomdag = next_stop_phantomdag
                next_stop_phantomdag = length(phantomdag) + 1
            end
        else
            parents.num == _depth && break
        end
        previous = vertexnums[parents.heads[1]]
        num = parents.num + one(SmallIntType)
        current_node = vertexnums[counter]
        for x in neighbors(g, current_node)
            x == previous && continue
            !(cycleavoid isa Nothing) && cycleavoid[x.v] && continue
            dead = hasavoid && (ringavoid[x.v] || (x.v == i && x.ofs < zero(SVector{D})))
            idx = get!(vertexdict, x, dead ? -length(phantomdag)-1 : length(dag)+1)
            if idx == length(dag)+1
                push!(vertexnums, x)
                push!(dag, JunctionNode(counter, num, parents.shortroots))
            elseif idx > counter
                junction = dag[idx]
                push!(junction.heads, counter)
                # We use pointer and unsafe_store! to modify the immutable type JunctionNode
                # This is only possible because it is stored in a Vector, imposing fixed addresses.
                ptr = pointer(dag, idx)
                if num == junction.num
                    unsafe_incr!(Ptr{SmallIntType}(ptr + lastshortoffset))
                    unsafe_union!(Ptr{ConstMiniBitSet{UInt32}}(ptr + shortrootsoffset), parents.shortroots.x)
                else
                    unsafe_union!(Ptr{ConstMiniBitSet{UInt32}}(ptr + longrootsoffset), parents.shortroots.x)
                end
            elseif dead && idx == -length(phantomdag)-1
                push!(phantomdag, PhantomJunctionNode(x, current_node, num))
            end
        end
    end
    return dag, vertexnums
end


function arcs_list_minimal(g::PeriodicGraph{D}, i, depth) where D
    dag = JunctionNode[JunctionNode()]
    vertexnums = PeriodicVertex{D}[PeriodicVertex{D}(i)]
    hintsize = ceil(Int, depth^2.9) # empirical estimate
    sizehint!(vertexnums, hintsize)
    sizehint!(dag, hintsize)
    append!(vertexnums, neighbors(g, i))
    append!(dag, JunctionNode(j) for j in 2:length(vertexnums))    
    vertexdict = Dict{PeriodicVertex{D},Int}((x => j) for (j,x) in enumerate(vertexnums))
    # sizehint!(vertexdict, hintsize)
    counter = 1
    _depth = depth % SmallIntType
    for parents in Iterators.rest(dag, 2)
        counter += 1
        parents.num == _depth && break
        previous = vertexnums[parents.heads[1]]
        num = parents.num + one(SmallIntType)
        for x in neighbors(g, vertexnums[counter])
            x == previous && continue
            idx = get!(vertexdict, x, length(dag)+1)
            if idx == length(dag)+1
                push!(vertexnums, x)
                push!(dag, JunctionNode(counter, num, parents.shortroots))
            elseif idx > counter
                junction = dag[idx]
                push!(junction.heads, counter)
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

# function _case_0(heads, lastshort, midnode)
#     return Vector{Int}[Int[1, heads[i], midnode] for i in length(heads):-1:(lastshort+1)]
# end
# function cycles_ending_at(dag::Vector{JunctionNode}, midnode, dist=nothing, vertexnums=nothing)
#     parents = dag[midnode]
#     length(parents.heads) ≥ 2 || return Vector{Int}[]
#     length(union(parents.shortroots, parents.longroots)) ≥ 2 || return Vector{Int}[]
#     shortest_n = parents.num % Int
#     num = parents.num + one(SmallIntType)
#     parents_lastshort = parents.lastshort % Int
#     iszero(shortest_n) && return _case_0(parents.heads, parents_lastshort, midnode)
#     ret = Vector{Int}[]
#     idx_stack1 = Vector{Int}(undef, shortest_n)
#     idx_stack2 = Vector{Int}(undef, shortest_n-1)
#     buffer = Vector{Int}(undef, 2*shortest_n + 3)
#     buffer[shortest_n+1] = 1
#     buffer[end] = midnode
#     for i1 in length(parents.heads):-1:2
#         if i1 == parents_lastshort
#             Base._deleteend!(buffer, 1)
#             Base._deleteend!(idx_stack1, 1)
#             buffer[end] = midnode
#         end
#         buffer[end-1] = parents.heads[i1]
#         last_positions1 = initial_compatible_arc!(buffer, idx_stack1, dag, false, dist, vertexnums)
#         while last_positions1 != ~zero(UInt64)
#             for i2 in 1:min(i1-1, parents_lastshort)
#                 head2 = parents.heads[i2]
#                 i1 > parents_lastshort && buffer[end-2] == head2 && continue # 3-cycle near midnode
#                 buffer[1] = head2
#                 if dist !== nothing && # checking for rings instead of cycles
#                   (is_distance_smaller!(dist, vertexnums[buffer[shortest_n+2]], vertexnums[head2], num) ||
#                   (i1 > parents_lastshort && is_distance_smaller!(dist, vertexnums[buffer[shortest_n+3]], vertexnums[head2], num)))
#                     continue
#                 end
#                 last_positions2 = initial_compatible_arc!(buffer, idx_stack2, dag, true, dist, vertexnums)
#                 while last_positions2 != ~zero(UInt64)
#                     push!(ret, copy(buffer))
#                     last_positions2 = next_compatible_arc!(buffer, last_positions2, idx_stack2, dag, true, dist, vertexnums)
#                 end
#             end
#             last_positions1 = next_compatible_arc!(buffer, last_positions1, idx_stack1, dag, false, dist, vertexnums)
#         end
#     end
#     return ret
# end


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

struct RingsEndingAt{T}
    dag::Vector{JunctionNode}
    midnode::Int
    heads::Vector{Int}
    lastshort::Int
    parentsnum::SmallIntType
    record::T
end

"""
    RingsEndingAt(dag, midnode, record)

Iterable over the rings of graph `g` around node `i` with `midnode` as vertex furthest
from `i`.

`record` should be set to `(dist, vertexnums)` where `dist == DistanceRecord(g, depth)` and
`dag, vertexnums == first(arcs_list(g, i, depth, ...))`, otherwise the iterator will return
many more cycles that may not be rings.

!!! warning
    In order to efficiently cycle through the rings, the iterator reuses a buffer on which
    the rings are written. This means that performing an iteration will change the value
    of the previously returned result: for example, `collect(RingsEndingAt(...))` will
    yield a list containing the same list (unlikely to be an actual ring) repeated over.
    To actually obtain the list of rings, copy the result as they arrive by doing
    `map(copy, RingsEndingAt(...))` or `[copy(x) for x in RingsEndingAt(...)]` for example.

    This also means that the list returned at each iteration should never be modified
    directly: copy it before.
"""
@inline function RingsEndingAt(dag, midnode, record=(nothing,nothing))
    parents = dag[midnode]
    RingsEndingAt{typeof(record)}(dag, midnode, parents.heads, parents.lastshort % Int, parents.num, record)
end
Base.IteratorSize(::RingsEndingAt) = Base.SizeUnknown()
Base.eltype(::RingsEndingAt) = Vector{Int}

function Base.iterate(x::RingsEndingAt, state=nothing)
    heads = x.heads
    length(heads) ≥ 2 || return nothing
    midnode = x.midnode
    lastshort = x.lastshort
    shortest_n = x.parentsnum % Int
    if iszero(shortest_n) # rings ending on a neighbour of the root
        length(heads) > lastshort || return nothing
        if state === nothing
            state = (Int[1, 0, midnode], Int[length(heads)+1], Int[], zero(UInt64), zero(UInt64), 0, 0)
        end
        i0 = @inbounds state[2][1] -= 1
        i0 > lastshort || return nothing
        state[1][2] = heads[i0]
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
            pop!(buffer)
            pop!(idx_stack1)
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

# function rings_around_old(g::PeriodicGraph{D}, i, depth=15, dist=DistanceRecord(g, depth)) where D
#     n = nv(g)
#     dag, vertexnums = arcs_list_minimal(g, i, depth)
#     hashes = [hash_position(x, n) for x in vertexnums]
#     ret = Vector{Int}[]
#     for midnode in 2:length(dag)
#         for cycle in RingsEndingAt(dag, midnode)
#         # for cycle in cycles_ending_at(dag, midnode)
#             isempty(cycle) && break
#             invalidcycle = false
#             lenc = length(cycle)
#             num = (lenc ÷ 2) % SmallIntType
#             if isodd(lenc)
#                 for i in 1:(num-1)
#                     verti = vertexnums[cycle[i]]
#                     vertj1 = vertexnums[cycle[num+i]]
#                     vertj2 = vertexnums[cycle[num+i+1]]
#                     if is_distance_smaller!(dist, verti, vertj1, num) ||
#                         is_distance_smaller!(dist, verti, vertj2, num)
#                         invalidcycle = true
#                         break
#                     end
#                 end
#             else
#                 for i in 1:(num-1)
#                     verti = vertexnums[cycle[i]]
#                     vertj = vertexnums[cycle[num+i]]
#                     if is_distance_smaller!(dist, verti, vertj, num)
#                         invalidcycle = true
#                         break
#                     end
#                 end
#             end
#             invalidcycle && continue
#             push!(ret, normalize_cycle!(hashes[cycle]))
#         end
#     end
#     return ret
# end

"""
    rings_around(g::PeriodicGraph{D}, i, depth=15, dist=DistanceRecord(g,depth), visited=nothing) where D

Return the list of all rings around node `i` in graph `g` up to length `2*depth+3`.

The returned rings are the list of `hash_position` of the corresponding vertices. To get
back the list of actual `PeriodicVertex` of a returned `ring` in the list, do
```julia
[reverse_hash_position(x, n, Val(D)) for x in ring]
```
where `n == nv(g)`. If the offsets of the corresponding vertices are not needed, simply do
```julia
[mod1(x, n) for x in ring]

`visited` is interpreted as the `ringavoid` argument of `arcs_list` unless
`dist === nothing`, in which case it is interpreted as the `cycleavoid` argument.
In particular, unless `dist === nothing`, only one ring will appear in the list even if
some of its translated images also pass through `PeriodicVertex{D}(i)`.
```
"""
function rings_around(g::PeriodicGraph{D}, i, depth=15, dist=DistanceRecord(g,depth), visited=nothing) where D
    n = nv(g)
    ringavoid, cycleavoid = dist isa Nothing ? (nothing, visited) : (visited, nothing)
    dag, vertexnums = arcs_list(g, i, depth, ringavoid, cycleavoid)
    # dag, vertexnums = arcs_list_minimal(g, i, depth)
    hashes = [hash_position(x, n) for x in vertexnums]
    ret = Vector{Int}[]
    for midnode in length(dag):-1:2
        for cycle in RingsEndingAt(dag, midnode, (dist, vertexnums))
        # for cycle in cycles_ending_at(dag, midnode, dist, vertexnums)
            newcycle = normalize_cycle!(hashes[cycle])
            push!(ret, newcycle)
        end
    end
    return ret
end

# cycles_around(g::PeriodicGraph, i, depth=15) = rings_around(g, i, depth, nothing)


"""
    no_neighboring_nodes(g, symmetries::AbstractSymmetryGroup)

Return a list of nodes representatives (modulo symmetries) such that all the vertices in
`g` either have their representative in the list or are surrounded by nodes whose
representatives are in the list.

This means that all cycles and rings of `g` have at least one node whose representative is
in the list.
"""
function no_neighboring_nodes(g, symmetries::AbstractSymmetryGroup)
    n = nv(g)
    # category: 0 => unvisited, 1 => to explore, 2 => unvisited unknown, 3 => not to explore
    categories = zeros(Int8, n)
    categories[1] = 2
    toexplore = Int[]
    next_i_init = 1
    uniques = unique(symmetries)
    Q = Int[1]
    @label loop
    for u in Q
        newcat = categories[u] == 2 ? UInt8(1) : UInt(2)
        if newcat == 1
            if any(x -> symmetries(x.v) == u, neighbors(g, u))
                categories[u] = 1
                push!(toexplore, u)
                newcat = 2
            else
                categories[u] = 3
            end
        end
        for x in neighbors(g, u)
            v = symmetries(x.v)
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
    for i_candidate in next_i_init:length(uniques)
        u_candidate = uniques[i_candidate]
        if categories[u_candidate] == 0
            next_i_init = i_candidate + 1
            empty!(Q)
            push!(Q, u_candidate)
            @goto loop
        end
    end
    return toexplore
end


# function rings_old(g::PeriodicGraph{D}, depth=15, symmetries=NoSymmetryGroup(g), dist=DistanceRecord(g,depth)) where D
#     toexplore = no_neighboring_nodes(g, symmetries)
#     allrings = Set{Vector{Int}}()
#     for x in toexplore
#         # newrings = rings_around_old(g, x, depth)
#         # newrings = rings_around_old(g, x, depth, dist)
#         newrings = rings_around(g, x, depth, dist)
#         union!(allrings, newrings)
#     end
#     return collect(allrings)
# end

function rings(g::PeriodicGraph{D}, depth=15, symmetries::AbstractSymmetryGroup=NoSymmetryGroup(g), dist=DistanceRecord(g,depth)) where D
    toexplore = sort!(no_neighboring_nodes(g, symmetries))
    ret = Vector{Int}[]
    visited = falses(nv(g))
    for x in toexplore
        newrings = rings_around(g, x, depth, dist, visited)
        append!(ret, newrings)
        visited[x] = true
        for symm in symmetries
            visited[symm[x]] = true
        end
    end
    symmetries isa NoSymmetryGroup && return ret
    sort!(ret; by=length, rev=true) # to avoid resizing the buffer below too many times
    buffer = similar(first(ret))
    n = nv(g)
    uniquerings = [[reverse_hash_position(x, n, Val(D)) for x in r] for r in ret]
    for symm in symmetries
        for ring in uniquerings
            nr = length(ring)
            resize!(buffer, nr)
            #=@inbounds=# for i in 1:nr
                buffer[i] = hash_position(symm[ring[i]], n)
            end
            normalize_cycle!(buffer)
            if buffer ∉ uniquerings
                push!(ret, copy(buffer))
            end
        end
    end
    return ret
end

# cycles(g::PeriodicGraph, depth=15, symmetries=NoSymmetryGroup(g)) = rings(g, depth, symmetries, nothing)

function cages_around(::PeriodicGraph{D}, depth) where D
    buffer = MVector{D,Int}(undef)
    buffer[end] = -depth - 1
    top = ~zero(UInt64) >> ((65 - D) % UInt8)
    num = (2*depth+1)^D
    ret = Vector{SVector{D,Int}}(undef, num)
    D > 64 && return ret # that's actually an error
    for i in 1:num
        fst = trailing_ones(top)
        top &= ~zero(UInt64) << (fst % UInt8)
        for i in 1:fst
            buffer[i] = -depth
        end
        if (buffer[fst+1] += 1) == depth
            top |= one(UInt64) << (fst % UInt8)
        end
        ret[i] = buffer
    end
    return ret
end

function find_known_pair!(known_pairs_dict, known_pairs, x)
    pairid = get!(known_pairs_dict, x, length(known_pairs)+1)
    pairid == length(known_pairs)+1 && push!(known_pairs, x)
    pairid
end

function unique_order(cycles)
    I = sortperm(cycles)
    last_c = cycles[I[1]]
    n = length(cycles)
    toremove = falses(n)
    i = 2
    while i ≤ n
        cycle = cycles[I[i]]
        if cycle == last_c
            toremove[i] = true
        else
            last_c = cycle
        end
        i += 1
    end
    deleteat!(I, toremove)
    J = sort(I; by= x->length(cycles[x]))
    return J
end

function sort_cycles(rs, known_pairs, known_pairs_dict, g::PeriodicGraph{D}, depth) where D
    ofss = cages_around(g, 2)
    tot_ofss = length(ofss)
    cycles = Vector{Vector{Int}}(undef, tot_ofss * length(rs))
    origin = zeros(Int, length(cycles))
    buffer = Vector{Int}(undef, 2*depth+3)
    ringbuffer = Vector{PeriodicVertex{D}}(undef, length(buffer))
    n = nv(g)
    zero_ofs = (tot_ofss+1) ÷ 2
    for (i, ring) in enumerate(rs)
        base = (i-1)*tot_ofss
        @simd for k in 1:length(ring)
            ringbuffer[k] = reverse_hash_position(ring[k], n, Val(D))
        end
        origin[base+zero_ofs] = i
        for (i_ofs, ofs) in enumerate(ofss)
            _last_p = ringbuffer[length(ring)]
            last_p = PeriodicVertex{D}(_last_p.v, _last_p.ofs .+ ofs)
            for j in 1:length(ring)
                _new_p = ringbuffer[j]
                new_p = PeriodicVertex{D}(_new_p.v, _new_p.ofs .+ ofs)
                buffer[j] = find_known_pair!(known_pairs_dict, known_pairs, minmax(last_p, new_p))
                last_p = new_p
            end
            cycles[base+i_ofs] = sort!(buffer[1:length(ring)])
        end
    end
    I = unique_order(cycles)
    return cycles[I], origin[I]
end

# Complexity O(len(a) + len(b))
function symdiff_cycles(a::Vector{Int}, b::Vector{Int})
    lenb = length(b)
    lena = length(a)
    c = Vector{Int}(undef, lenb + lena)
    counter_b = 1
    y = lenb == 0 ? typemax(Int) : (@inbounds b[1])
    i = 0
    j = 1
    n = length(a)
    @inbounds while i < n
        i += 1
        x = a[i]
        while y < x
            c[j] = y
            j += 1
            counter_b += 1
            counter_b > lenb && @goto fillwitha
            y = b[counter_b]
        end
        if y == x
            counter_b += 1
            if counter_b > lenb
                i += 1
                @goto fillwitha
            end
            y = b[counter_b]
        else
            c[j] = x
            j += 1
        end
    end
    remaining_towriteb = lenb - counter_b + 1
    remaining_towriteb > 0 && unsafe_copyto!(c, j, b, counter_b, remaining_towriteb)
    @inbounds resize!(c, (j + remaining_towriteb - 1) % UInt)
    return c

    @label fillwitha
    remaining_towritea = lena - i + 1
    remaining_towritea > 0 && unsafe_copyto!(c, j, a, i, remaining_towritea)
    @inbounds resize!(c, (j + remaining_towritea - 1) % UInt)
    return c
end

struct IterativeGaussianElimination
    rings::Vector{Vector{Int}} # The rows of the matrix, in sparse format
    lengths::Vector{Int} # For each row, the length of the corresponding ring
    next::Vector{Int} # Order of the rings (linked list compact format)
    previous::Vector{Int} # verifies previous[next[i+1]] == i
    shortcuts::Vector{Int} # verifies shortcuts[i] = 0 || rings[shortcuts[i]][1] == i
end
function IterativeGaussianElimination(ring::Vector{Int}, sizehint=ring[1])
    r1 = ring[1]
    shortcuts = zeros(Int, sizehint)
    shortcuts[r1] = 1
    IterativeGaussianElimination([ring], [length(ring)], [1,0], [0], shortcuts)
end

# For an IterativeGaussianElimination of i rings of size at most l, symdiff_cycles cost at
# most (l + lenr) + (2l + lenr) + ... + (il + lenr) = O(i^2*l +i*lenr)
# If the size of the k-th ring is at most k*l, then the worst-case complexity is O(i^3*l)
# Calling gaussian_elimination ν times thus leads to a worst-case complexity in O(ν^4*l)
# The reasonning can probably be refined...
function gaussian_elimination!(gauss::IterativeGaussianElimination, r::Vector{Int})
    rings = gauss.rings
    nexts = gauss.next
    lengths = gauss.lengths
    shortcuts = gauss.shortcuts
    lenshort = length(shortcuts)
    len = length(r)
    maxlen = 0
    last_idx = 0
    r1 = r[1]

    short = r1 > lenshort ? 0 : shortcuts[r1]
    idx = short == 0 ? nexts[1] : short
    ridx = rings[idx]
    ridx1 = ridx[1]
    short == 0 || @goto inloop

    # while true
    while r1 > ridx1
        @label whilesuperior
        last_idx = idx
        idx = nexts[idx+1]
        idx == 0 && @goto outloop
        ridx = rings[idx]
        ridx1 = ridx[1]
    end
    r1 == ridx1 || @goto outloop
    # if length(r) < length(ridx) # optional optimization
    #     r, ridx = ridx, r
    #     rings[idx] = ridx
    # end
    @label inloop
    maxlen = max(lengths[idx], maxlen)
    r = symdiff_cycles(r, ridx)
    isempty(r) && return false, maxlen < len
    r1 = r[1]
    short = r1 > lenshort ? 0 : shortcuts[r1]
    short == 0 && @goto whilesuperior
    idx = short
    ridx = rings[idx]
    ridx1 = ridx[1]
    @goto inloop

    # end

    @label outloop
    push!(rings, r)
    nrings = length(rings)
    push!(lengths, len)
    push!(gauss.previous, last_idx)
    nextlast = nexts[last_idx+1]
    if nextlast != 0
        gauss.previous[nextlast] = nrings
    end
    if last_idx == 0
        push!(nexts, nexts[1])
        nexts[1] = nrings
    else
        push!(nexts, idx)
        nexts[last_idx+1] = nrings
    end
    r1 > length(shortcuts) && append!(shortcuts, 0 for _ in 1:(r1-length(shortcuts)))
    shortcuts[r1] = nrings
    true, false
end

function retrieve_vcycle(ecycle, known_pairs)
    fst_pair = known_pairs[ecycle[1]]
    last_pair = known_pairs[ecycle[end]]
    last_v, other_v = fst_pair[1], fst_pair[2]
    if last_v == last_pair[1] || last_v == last_pair[2]
        last_v = fst_pair[2]
        other_v = fst_pair[1]
    end
    n = length(ecycle)
    ret = Vector{typeof(last_v)}(undef, n)
    ret[1] = last_v
    ret[n] = other_v
    for j in 2:(n-1)
        this_pair = known_pairs[ecycle[j]]
        other_v = this_pair[1]
        if last_v == other_v
            last_v = this_pair[2]
        else
            last_v = other_v
        end
        ret[j] = last_v
    end
    return ret
end

function strong_rings(g::PeriodicGraph{D}, depth=15, symmetries::AbstractSymmetryGroup=NoSymmetryGroup(g), dist=DistanceRecord(g,depth)) where D
    rs = rings(g, depth, symmetries, dist)
    known_pairs = Tuple{PeriodicVertex{D},PeriodicVertex{D}}[]
    known_pairs_dict = Dict{Tuple{PeriodicVertex{D},PeriodicVertex{D}},Int}()
    hintsize = 3^D*ne(g)^2
    sizehint!(known_pairs, hintsize)
    sizehint!(known_pairs_dict, hintsize)
    cycles, origin = sort_cycles(rs, known_pairs, known_pairs_dict, g, depth)

    fst_ring = popfirst!(cycles)
    gauss = IterativeGaussianElimination(fst_ring)
    ret = Vector{Int}[]
    # tmpret = Vector{Int}[]
    n = nv(g)
    orig_1 = popfirst!(origin)
    orig_1 != 0 && push!(ret, rs[orig_1])
    # orig_1 != 0 && push!(tmpret, fst_ring)
    for (i, cycle) in enumerate(cycles)
        isfree, onlysmallercycles = gaussian_elimination!(gauss, cycle)
        onlysmallercycles && continue # cycle is a linear combination of smaller cycles
        origin[i] != 0 && push!(ret, rs[origin[i]])
        # origin[i] != 0 && push!(tmpret, cycle)
    end
    return ret#, tmpret, known_pairs, known_pairs_dict, gauss
end


struct RingAttributions{D}
    rings::Vector{Vector{PeriodicVertex{D}}}
    attrs::Vector{Vector{Tuple{Int,Int}}}

    function RingAttributions{D}(n, rs::Vector{Vector{Int}}) where D
        rings = [Vector{PeriodicVertex{D}}(undef, length(r)) for r in rs]
        attrs = [Tuple{Int,Int}[] for _ in 1:n]
        for (i, r) in enumerate(rs)
            ring = rings[i]
            for (j, x) in enumerate(r)
                u = reverse_hash_position(x, n, Val(D))
                ring[j] = u
                push!(attrs[u.v], (i, j))
            end
        end
        return new{D}(rings, attrs)
    end
end

function RingAttributions(g::PeriodicGraph{D}, strong=false, depth=15, symmetries::AbstractSymmetryGroup=NoSymmetryGroup(g), dist=DistanceRecord(g,depth)) where D
    rs = (strong ? strong_rings : rings)(g, depth, symmetries, dist)
    return RingAttributions{D}(nv(g), rs)
end

Base.@propagate_inbounds function Base.getindex(ras::RingAttributions, i::Integer)
    @boundscheck checkbounds(ras.attrs, i)
    RingIncluding(ras, i)
end

struct RingIncluding{D}
    ras::RingAttributions{D}
    i::Int
end
function Base.getindex(ri::RingIncluding{D}, j::Integer) where {D}
    newring_idx, idx = ri.ras.attrs[ri.i][j]
    newring = ri.ras.rings[newring_idx]
    ofs = newring[idx].ofs
    return PeriodicNeighborList{D}(.-ofs, newring)
end
function Base.iterate(ri::RingIncluding{D}, state=1) where {D}
    state > length(ri) && return nothing
    return (@inbounds ri[state]), state+1
end
Base.length(ri::RingIncluding) = length(ri.ras.attrs[ri.i])
Base.eltype(::RingIncluding{D}) where {D} = PeriodicGraphs.PeriodicNeighborList

function Base.show(io::IO, ri::RingIncluding{D}) where D
    print(io, RingIncluding{D}, '(', length(ri), " rings containing vertex ", ri.i, ')')
end
