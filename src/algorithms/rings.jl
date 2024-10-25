# Ring statistics

export rings, strong_rings, strong_erings, RingAttributions, RingIncluding,
       normalize_cycle!

"""
    ConstMiniBitSet{T} <: AbstractSet{Int}

Fixed-size bitset stored on a single word of type `T`, typically a `UInt64` or a `UInt32`.
"""
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
Base.in(i::Integer, x::ConstMiniBitSet{T}) where {T} = ((x.x >> (i % UInt8)) & one(T)) % Bool
Base.setdiff(x::ConstMiniBitSet, y::ConstMiniBitSet) = _constminibitset(x.x & ~y.x)
Base.symdiff(x::ConstMiniBitSet, y::ConstMiniBitSet) = _constminibitset(x.x ⊻ y.x)
Base.length(x::ConstMiniBitSet) = count_ones(x.x)
Base.isempty(x::ConstMiniBitSet) = iszero(x.x)
Base.eltype(::Type{ConstMiniBitSet{T}}) where {T} = Int

function Base.iterate(::ConstMiniBitSet{T}, x::T) where T
    first = trailing_zeros(x)
    first == 8*sizeof(T) && return nothing
    return (first, x & (~(one(T) << (first % UInt8))))
end
Base.iterate(x::ConstMiniBitSet) = iterate(x, x.x)

hasonly(x::ConstMiniBitSet{T}, i::Integer) where {T} = iszero(x.x & ~(one(T) << (i % UInt8)))

Base.minimum(x::ConstMiniBitSet) = first(x)
Base.maximum(x::ConstMiniBitSet{T}) where {T} = 8*sizeof(T) - leading_zeros(x.x) - 1

function Base.show(io::IO, x::ConstMiniBitSet)
    print(io, ConstMiniBitSet, "([")
    join(io, x, ", ")
    print(io, "])")
end


const SmallIntType = Int8 # used for distances mostly

"""
    JunctionNode{T}

Element of the DAG representing the set of arcs linking each vertex `x` to a fixed vertex
`i` of graph `g`. Each `JunctionNode` contains information on the arcs passing through a
particular vertex `x`:
- `num` is the length of the shortest path between `x` and `i`. Since we only collect
  rings, only the shortest paths are of interest, as well as the path of length `num+1`
  which may form an odd-length ring when combined with a shortest path.
- `heads` is a list of neighbors of `x` such that the `lastshort` first are at distance
  `num - 1` from `i`, and the rest are at distance `num` and have a lower value than `x`.
- `shortroots` is the set of roots reachable on a shortest path from `x` to `i`. A root is
  the neighbor of `i` on that path, a.k.a. the second-to-last vertex on that path.
- `longroots` is the set of roots reachable on a path of length `num+1` from `x` to `i`.
"""
struct JunctionNode{T}
    shortroots::ConstMiniBitSet{T}
    longroots::ConstMiniBitSet{T}
    lastshort::SmallIntType # (use SmallIntType for better padding)
    num::SmallIntType # length of the shortest branch
    heads::Vector{Int}
    function JunctionNode{T}(heads::AbstractVector, num, roots::ConstMiniBitSet{T}) where T
        new(roots, ConstMiniBitSet{T}(), length(heads) % SmallIntType, num % SmallIntType, heads)
    end
end
JunctionNode{T}() where {T} = JunctionNode{T}(Int[], -one(SmallIntType), ConstMiniBitSet{T}())
JunctionNode{T}(head::Integer, num, roots::ConstMiniBitSet{T}) where {T} = JunctionNode{T}(Int[head], num, roots)
JunctionNode{T}(root::Integer) where {T} = JunctionNode{T}(1, zero(SmallIntType), ConstMiniBitSet{T}(root))

function Base.show(io::IO, x::JunctionNode{T}) where {T}
    print(io, "JunctionNode{", T, "}([")
    join(io, x.heads, ", ")
    print(io, "])")
end

"""
    PhantomJunctionNode{D}

Element of the phantom DAG.

Similarly to the DAG of `JunctionNode`, the phantom DAG tracks arcs linking vertices `x` to
a fixed vertex `i`, except that the vertices `x` are those that should be ignored in the
returned list of rings.
Thus, only the shortest distance between `x` and `i` needs to be recorded, since the arcs
themselves will be discarded eventually.
"""
struct PhantomJunctionNode{D}
    self::PeriodicVertex{D}
    parent::PeriodicVertex{D} # storing it allows skipping a dict call
    num::SmallIntType
end

function prepare_phantomdag!(dag::Vector{JunctionNode{T}}, vertexnums, vertexdict, g::PeriodicGraph{D}, i, avoid, cycleavoid) where {T,D}
    phantomdag = PhantomJunctionNode{D}[]
    for x in neighbors(g, i)
        !(cycleavoid isa Nothing) && cycleavoid[x.v] && iszero(x.ofs) && continue
        if avoid[x.v] || (x.v == i && x.ofs < zero(SVector{D,Int}))
            push!(phantomdag, PhantomJunctionNode(x, PeriodicVertex{D}(i), zero(SmallIntType)))
            vertexdict[x] = -length(phantomdag)
        else
            push!(vertexnums, x)
            push!(dag, JunctionNode{T}(length(vertexnums)))
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

const lastshortoffsetpos = findfirst(==(:lastshort), fieldnames(JunctionNode))
const shortrootsoffsetpos = findfirst(==(:shortroots), fieldnames(JunctionNode))
const longrootsoffsetpos = findfirst(==(:longroots), fieldnames(JunctionNode))
@inline function lastshortoffsetof(ptr::Ptr{JunctionNode{T}}) where {T}
    Ptr{SmallIntType}(ptr + fieldoffset(JunctionNode{T}, lastshortoffsetpos))
end
@inline function shortrootsoffsetof(ptr::Ptr{JunctionNode{T}}) where {T}
    Ptr{ConstMiniBitSet{T}}(ptr + fieldoffset(JunctionNode{T}, shortrootsoffsetpos))
end
@inline function longrootsoffsetof(ptr::Ptr{JunctionNode{T}}) where {T}
    Ptr{ConstMiniBitSet{T}}(ptr + fieldoffset(JunctionNode{T}, longrootsoffsetpos))
end
@inline unsafe_incr!(ptr::Ptr{T}) where {T} = unsafe_store!(ptr, unsafe_load(ptr) + one(T))
@inline unsafe_union!(ptr, val) = unsafe_store!(ptr, _constminibitset(val | unsafe_load(ptr).x))

"""
    arcs_list(g::PeriodicGraph{D}, i, depth::T, ringavoid=nothing, cycleavoid=nothing) where {D,T}

Compute the list of shortest arcs starting from vertex `i` up to length `depth+1`.
Vertices in `ringavoid` are not included in the returned arcs, but considered still part of
the graph for distance computations.
Vertices in `cycleavoid` are considered removed from the graph completely.

Return `(dag, vertexnums)` where `dag` is a `Vector{JunctionNode{T}}` representing, for each
visited node, the DAG of all arcs from that node back to `i`, in a compact representation.
`vertexnums` is a `Vector{PeriodicVertex{D}}` whose `k`-th value is the vertex represented
by number `k` in `dag`.

If `ringavoid !== nothing`, `dag` will also not include arcs that pass through nodes of the
form `PeriodicVertex{D}(j, ofs)` with `j == i` and `ofs < zero(SVector{Int,D})`: this
allows eagerly pruning cycles that are translations of others. Note that this can result in
missing cycles if those pass through at least three nodes with `j == i`, but that situation
should be exceptionally rare.

!!! note
    The type `T` of `depth` is used as type parameter to the `JunctionNode`, just to avoid
    having a dedicated argument (since `depth` should be at most 62 for the rest of the
    algorithm to work). This size controls the maximal degree a vertex of `g` should have :
    for example, `T == UInt32` indicates that all vertices must have degree at most 31.
"""
function arcs_list(g::PeriodicGraph{D}, i, depth::T, ringavoid=nothing, cycleavoid=nothing) where {D,T}
    dag = JunctionNode{T}[JunctionNode{T}()]
    vertexnums = PeriodicVertex{D}[PeriodicVertex{D}(i)]
    vertexdict = Dict{PeriodicVertex{D},Int}()
    hintsize = ceil(Int, depth^2.9) # empirical estimate
    sizehint!(vertexnums, hintsize)
    sizehint!(dag, hintsize)
    sizehint!(vertexdict, hintsize)
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
        append!(dag, JunctionNode{T}(j) for j in 2:length(vertexnums))
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
                parents.num ≥ _depth && break
                next_stop = length(dag) + 1
                handle_phantomdag!(phantomdag, vertexdict, g, last_stop_phantomdag, next_stop_phantomdag)
                last_stop_phantomdag = next_stop_phantomdag
                next_stop_phantomdag = length(phantomdag) + 1
            end
        else
            parents.num ≥ _depth && break
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
                push!(dag, JunctionNode{T}(counter, num, parents.shortroots))
            elseif idx > counter
                junction = dag[idx]
                push!(junction.heads, counter)
                # We use pointer and unsafe_store! to modify the immutable type JunctionNode
                # This is only possible because it is stored in a Vector, imposing fixed addresses.
                ptr = pointer(dag, idx)
                if num == junction.num
                    unsafe_incr!(lastshortoffsetof(ptr))
                    unsafe_union!(shortrootsoffsetof(ptr), parents.shortroots.x)
                else
                    unsafe_union!(longrootsoffsetof(ptr), parents.shortroots.x)
                end
            elseif dead && idx == -length(phantomdag)-1
                push!(phantomdag, PhantomJunctionNode(x, current_node, num))
            end
        end
    end
    return dag, vertexnums
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

function init_distance_record!(dist::DistanceRecord{D}, i, j) where D
    dist.first_indices[i] = 1
    g::PeriodicGraph{D} = dist.g
    dist.Q_lists[i] = Tuple{PeriodicVertex{D},SmallIntType}[(x, one(SmallIntType)) for x in neighbors(g, i)]
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


function next_compatible_arc!(buffer, last_positions, idx_stack, dag, check, dist, vertexnums)
    num_choice = length(idx_stack)
    updated_last_positions = last_positions
    haslongarc = isodd(length(buffer))
    num = (num_choice + 2) % SmallIntType
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
    num = (num_choice + 2) % SmallIntType
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


struct RingsEndingAt{T,R}
    dag::Vector{JunctionNode{T}}
    midnode::Int
    heads::Vector{Int}
    lastshort::Int
    parentsnum::SmallIntType
    record::R
end

"""
    RingsEndingAt(dag, midnode, record)

Iterable over the rings of graph `g` around node `i` with `midnode` as vertex furthest
from `i`. If there are two such vertices (odd ring), `midnode` is the higher of the two.

`record` should be set to `(dist, vertexnums)` where `dist == DistanceRecord(g, depth)` and
`dag, vertexnums == first(arcs_list(g, i, depth, ...))`, otherwise the iterator will return
many more cycles that may not be rings.

!!! warning
    In order to efficiently cycle through the rings, the iterator reuses a buffer on which
    the rings are written. This means that performing an iteration will change the value
    of the previously returned result: for example, `collect(RingsEndingAt(...))` will
    yield a list containing the same sublist (unlikely to be an actual ring) repeated over.
    To actually obtain the list of rings, copy the result as they arrive by doing
    `map(copy, RingsEndingAt(...))` or `[copy(x) for x in RingsEndingAt(...)]` for example.

    This also means that the list returned at each iteration should never be modified
    directly: copy it before.
"""
@inline function RingsEndingAt(dag, midnode, record=(nothing,nothing))
    parents = dag[midnode]
    RingsEndingAt(dag, midnode, parents.heads, parents.lastshort % Int, parents.num, record)
end
Base.IteratorSize(::RingsEndingAt) = Base.SizeUnknown()
Base.eltype(::Type{RingsEndingAt{T,R}}) where {T,R} = Vector{Int}

#= Iteration strategy for rings passing through `i` and `midnode` (furthest from `i`).
Each ring is made from the composition of two arcs, let's call them arc1 and arc2. We take
arc1 in descending order of the `heads` of `midnode` in the DAG until it reaches the first
head, where we stop. This means arc1 starts with the long rings (odd-size) until reaching
the `lastshort` head, from which there will only be short rings.
For each arc1, arc2 is taken in descending order of the `heads` starting at the head after
that of arc1, or, if that head correspond to a long arc, at `lastshort`.

`shortest_n` is the length of the shortest path from `midnode` to `i`.

Each `state` of the iterator is made of 7 elements:
- the `buffer` stores the ring under exploration. Its length is `2*shortest_n + 3` if arc1
  is still in the long arcs, then it becomes `2*shortest_n + 2`
  Its last element is always `midnode` and its `shortest_n+1` is 1, designating vertex `i`.
  arc2, which is always a short arc, is stored in the `shortest_n` first element, and arc1
  between `shortest_n+2` and the end.
- `idx_stack1` is the stack of indices of the head chosen in each `heads` list to form
  arc1. In other words, `buffer[end-i] == heads[j]` implies `idx_stack1[i] == j`, where
  `heads` is that of the `JunctionNode` corresponding to `buffer[end-i+1]`.
- `idx_stack1` is the same for arc2.
- `last_positions1` is a bitset where element `k` is true if `idx_stack1[i]` is at its
  highest value.
- `last_positions2` is the equivalent for arc2.
- `i1` is the index of the first element of arc1 among the heads of `midnode`. As
  mentionned before, it progresses in decreasing order.
- `i2` is the equivalent for arc2.
=#
function Base.iterate(x::RingsEndingAt{T,R}, state=nothing) where {T,R}
    heads = x.heads
    length(heads) ≥ 2 || return nothing
    midnode = x.midnode
    lastshort = x.lastshort
    shortest_n = x.parentsnum % Int
    if iszero(shortest_n) # `midnode` is a neighbor of the root: handle things differently
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

    # At this point, state === nothing, so this the first iteration.
    # First, check that there can be any ring by checking there are at least 2 compatible roots:
    length(union(dag[midnode].shortroots, dag[midnode].longroots)) ≥ 2 || return nothing

    idx_stack1 = Vector{Int}(undef, shortest_n)
    idx_stack2 = Vector{Int}(undef, shortest_n-1) # idx_stack2 can only contain short arcs
    buffer = Vector{Int}(undef, 2*shortest_n + 3)
    buffer[shortest_n+1] = 1 # Fixed position
    buffer[end] = midnode # Will evolve at most once, when going from long to short arcs
    i1 = length(heads)
    while i1 > 1
        if i1 == lastshort # This corresponds to the transition from long to short arcs.
            pop!(buffer)
            pop!(idx_stack1)
            buffer[end] = midnode
        end
        buffer[end-1] = heads[i1]
        last_positions1 = initial_compatible_arc!(buffer, idx_stack1, dag, false, dist, vertexnums)
        @label start_iter_i1
        if last_positions1 == ~zero(UInt64)
            @goto next_iter_i1
        elseif dist !== nothing && i1 > lastshort
            sroots = dag[midnode].shortroots
            root1 = buffer[shortest_n+2]
            while root1 ∈ sroots || is_distance_smaller!(dist, vertexnums[root1], vertexnums[midnode], num)
                last_positions1 = next_compatible_arc!(buffer, last_positions1, idx_stack1, dag, false, dist, vertexnums)
                last_positions1 == ~zero(UInt64) && @goto next_iter_i1
                root1 = buffer[shortest_n+2]
            end
        end
        # At this position, buffer[shortest_n+2:end] contains the first arc1 corresponding to i1

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
        @goto start_iter_i1

        @label next_iter_i1
        i1 -= 1
    end
    nothing
end

"""
    normalize_cycle!(cycle::Vector{Int}, n, v::Val{D}) where D

In-place rotate and possibly reverse `cycle`, a `Vector{Int}` whose elements are the
`hash_position` of vertices of `g` so that the result is the same for all such vectors that
represent the same cycle, possibly translated to a different unit cell or rotated.

The graph `g::PeriodicGraph{D}` is represented as `n = nv(g)` and `v = Val(D)`
"""
function normalize_cycle!(cycle::Vector{Int}, n, v::Val{D}) where D
    minval, fst = findmin(Base.Fix2(mod1, n), cycle)
    ofsfst = reverse_hash_position(fld1(cycle[fst], n) - 1, v)
    lenc = length(cycle)
    for i in fst+1:lenc
        _f, _m = fldmod1(cycle[i], n)
        _m == minval || continue
        ofs = reverse_hash_position(_f - 1, v)
        if ofs < ofsfst
            ofsfst = ofs
            fst = i
        end
    end
    before_f, before_m = fldmod1(cycle[mod1(fst-1, lenc)], n)
    after_f, after_m = fldmod1(cycle[mod1(fst+1, lenc)], n)
    endreverse = before_m > after_m ||
                 (before_m == after_m && reverse_hash_position(before_f, v) > reverse_hash_position(after_f, v))
    if fst != 1 || !endreverse
        reverse!(cycle, 1, fst - endreverse)
        reverse!(cycle, fst - endreverse + 1, lenc)
        endreverse && reverse!(cycle)
    end
    if !iszero(ofsfst)
        for i in 1:lenc
            rev = reverse_hash_position(cycle[i], n, v)
            cycle[i] = hash_position(PeriodicVertex{D}(rev.v, rev.ofs .- ofsfst), n)
        end
    end
    cycle, ofsfst
end

function normalize_cycle!(cycle::Vector{PeriodicVertex{D}}) where D
    minval, fst = findmin(cycle)
    lenc = length(cycle)
    endreverse = cycle[mod1(fst-1, lenc)] > cycle[mod1(fst+1, lenc)]
    if fst != 1 || !endreverse
        reverse!(cycle, 1, fst - endreverse)
        reverse!(cycle, fst - endreverse + 1, lenc)
        endreverse && reverse!(cycle)
    end
    ofsfst = minval.ofs
    if !iszero(ofsfst)
        for (i, c) in enumerate(cycle)
            cycle[i] = PeriodicVertex{D}(c.v, c.ofs - ofsfst)
        end
    end
    cycle, ofsfst
end


"""
    rings_around(g::PeriodicGraph{D}, i, depth=15, dist::DistanceRecord=DistanceRecord(g,depth), visited=nothing) where D

Return the list of all rings around node `i` in graph `g` up to length `2*depth+3`.

The returned rings are the list of `hash_position` of the corresponding vertices. To get
back the list of actual `PeriodicVertex` of a returned `ring` in the list, do
```julia
[reverse_hash_position(x, g) for x in ring]
```

If the offsets of the corresponding vertices are not needed, simply do
```julia
[mod1(x, n) for x in ring]   # n should be nv(g)
```

`visited` is interpreted as the `ringavoid` argument of `arcs_list` unless
`dist === nothing`, in which case it is interpreted as the `cycleavoid` argument.
In particular, unless `dist === nothing`, only one ring will appear in the list even if
some of its translated images also pass through `PeriodicVertex{D}(i)`.
"""
function rings_around(g::PeriodicGraph{D}, i, depth=15, dist::DistanceRecord=DistanceRecord(g,depth), visited=nothing) where D
    n = nv(g)
    ringavoid, cycleavoid = dist isa Nothing ? (nothing, visited) : (visited, nothing)
    dag, vertexnums = arcs_list(g, i, depth, ringavoid, cycleavoid)
    hashes = [hash_position(x, n) for x in vertexnums]
    ret = Vector{Int}[]
    for midnode in length(dag):-1:2
        for cycle in RingsEndingAt(dag, midnode, (dist, vertexnums))
        # for cycle in cycles_ending_at(dag, midnode, dist, vertexnums)
            newcycle, _ = normalize_cycle!(hashes[cycle], n, Val(D))
            push!(ret, newcycle)
        end
        yield()
    end
    return ret
end

# cycles_around(g::PeriodicGraph, i, depth=15) = rings_around(g, i, depth, nothing)


"""
    no_neighboring_nodes(g, symmetries::AbstractSymmetryGroup=NoSymmetryGroup(g))

Return a list of nodes representatives (modulo symmetries) such that all the vertices in
`g` either have their representative in the list or are surrounded by nodes whose
representatives are in the list.

This means that all cycles and rings of `g` have at least one node whose representative is
in the list.
"""
function no_neighboring_nodes(g, symmetries::AbstractSymmetryGroup=NoSymmetryGroup(g))
    n = nv(g)
    # category: 0 => unvisited, 1 => to explore, 2 => unknown (queued), 3 => not to explore
    categories = zeros(Int8, n)
    toexplore = Int[]
    next_i_init = 2
    uniques = unique(symmetries)
    u_candidate = uniques[1]
    Q = Int[]
    @label loop
    push!(Q, u_candidate)
    categories[u_candidate] = 2
    for u in Q
        newcat = categories[u] == 2 ? UInt8(1) : UInt8(2)
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
                if iszero(cat)
                    push!(Q, v)
                end
                categories[v] = newcat
                newcat == 1 && push!(toexplore, v)
            end
        end
    end
    for i_candidate in next_i_init:length(uniques)
        u_candidate = uniques[i_candidate]
        if iseven(categories[u_candidate])
            next_i_init = i_candidate + 1
            empty!(Q)
            @goto loop
        end
    end
    return toexplore
end


struct RingSymmetry{D,S<:AbstractSymmetry} <: AbstractSymmetry
    symm::S
    nvg::Int
end
function Base.getindex(rs::RingSymmetry{D}, ring::Vector{Int}) where D
    buffer = [hash_position(rs.symm[reverse_hash_position(r, rs.nvg, Val(D))], rs.nvg) for r in ring]
    first(normalize_cycle!(buffer, rs.nvg, Val(D)))
end
function Base.getindex(rs::RingSymmetry{D}, ring::Vector{PeriodicVertex{D}}) where D
    buffer = [rs.symm[x] for x in ring]
    first(normalize_cycle!(buffer))
end

struct RingSymmetryGroup{D,S,T<:AbstractSymmetryGroup{S},U} <: AbstractSymmetryGroup{RingSymmetry{D,S}}
    dict::Dict{Vector{Int},Int}
    uniques::U
    nvg::Int
    symms::T
    function RingSymmetryGroup{D}(dict, uniques::U, nvg, symms::T) where {D,U,T<:AbstractSymmetryGroup}
        new{D,eltype(T),T,U}(dict, uniques, nvg, symms)
    end
end

(rsg::RingSymmetryGroup)(x::AbstractVector{Int}) = rsg.dict[x]
Base.unique(rsg::RingSymmetryGroup) = rsg.uniques
function Base.iterate(rsg::RingSymmetryGroup{D,S}, state=nothing) where {D,S}
    x = state isa Nothing ? iterate(rsg.symms) : iterate(rsg.symms, something(state))
    x isa Nothing && return nothing
    symm, b = x
    return RingSymmetry{D,S}(symm, rsg.nvg), Some(b)
end
Base.length(rsg::RingSymmetryGroup) = length(rsg.symms)
Base.one(rsg::RingSymmetryGroup{D,S}) where {D,S} = RingSymmetry{D,S}(one(rsg.symms), rsg.nvg)

"""
    rings(g::PeriodicGraph{D}, [depth::Integer=15,] symmetries::AbstractSymmetryGroup=NoSymmetryGroup(g), dist::DistanceRecord=DistanceRecord(g,depth)) where D

Compute the list of rings in `g`, up to length `2*depth+3`. Return the list of `Vector{Int}`
where each sublist is a ring whose vertices are the `reverse_hash_position`s of the sublist
elements. Also return an [`AbstractSymmetryGroup`](@ref) acting on the returned rings.

A ring is a cycle of the graph for which there is no shortcut, i.e. no path in the graph
between two vertices of the cycle that is shorter than either path connecting the vertices
in the cycle.

If provided, `symmetries` should represent the symmetries of the graph as a
[`AbstractSymmetryGroup`](@ref) object respecting its documented interface.

A [`PeriodicGraphs.DistanceRecord`](@ref) `dist` can be optionally provided
to track the distances between pairs of vertices in the graph.
"""
function rings(g::PeriodicGraph{D}, depth::Integer=15, symmetries::AbstractSymmetryGroup=NoSymmetryGroup(g), dist::DistanceRecord=DistanceRecord(g,depth)) where D
    # The following errors are irrecoverable without changing the structure of the
    # algorithm to allow larger bitsets. Open an issue if required.
    maxdeg = maximum(degree(g); init=0)
    if maxdeg > 126
        error("A vertex has a degree too large (> 127) to be represented in a ConstMiniBitSet. Please open an issue.")
    elseif depth > 62
        error("Required depth is > 62 which cannot be handled currently. Please open an issue.")
    end
    toexplore = sort!(no_neighboring_nodes(g, symmetries))
    ret = Vector{Int}[]
    visited = falses(nv(g))
    for x in toexplore
        newrings = if maxdeg < 7
            rings_around(g, x, depth % UInt8, dist, visited)
        elseif maxdeg < 15
            rings_around(g, x, depth % UInt16, dist, visited)
        elseif maxdeg < 31
            rings_around(g, x, depth % UInt32, dist, visited)
        elseif maxdeg < 63
            rings_around(g, x, depth % UInt64, dist, visited)
        else
            rings_around(g, x, depth % UInt128, dist, visited)
        end
        append!(ret, newrings)
        visited[x] = true
        for symm in symmetries
            visited[symm[x]] = true
        end
    end
    symmetries isa NoSymmetryGroup && return ret, NoSymmetryGroup(length(ret))
    sort!(ret; by=length, rev=true) # to avoid resizing the buffer below too many times
    buffer = similar(first(ret))
    n = nv(g)
    ringdict = Dict{Vector{Int},Int}()
    symmuniques = Int[]
    toremove = Int[]
    N = length(ret)
    for idx in 1:N
        hring = ret[idx]
        nr = length(hring)
        resize!(buffer, nr)
        corrected_idx = idx - length(toremove)
        newidx = get!(ringdict, hring, corrected_idx)
        if newidx == corrected_idx
            push!(symmuniques, corrected_idx)
        else
            push!(toremove, idx)
        end
        ring = [reverse_hash_position(x, g) for x in hring]
        for symm in symmetries
            #=@inbounds=# for i in 1:nr
                buffer[i] = hash_position(symm[ring[i]], n)
            end
            normalize_cycle!(buffer, n, Val(D))
            if !haskey(ringdict, buffer)
                cp = copy(buffer)
                push!(ret, cp)
                ringdict[cp] = newidx
            end
        end
        yield()
    end
    deleteat!(ret, toremove)
    return ret, RingSymmetryGroup{D}(ringdict, symmuniques, n, symmetries)
end
function rings(g::PeriodicGraph{D}, symmetries::AbstractSymmetryGroup, dist::DistanceRecord=DistanceRecord(g,15)) where D
    rings(g, 15, symmetries, dist)
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
cages_around(::PeriodicGraph{0}, _) = [SVector{0,Int}()]

const VertexPair{D} = Tuple{PeriodicVertex{D},PeriodicVertex{D}}

"""
    EdgeDict{D}

Map from pairs of `PeriodicVertex{D}` to the identifier of the corresponding edge.

`kp::EdgeDict{D}` should be queried by either `get!(kp, minmax(v1, v2))` where `v1` and
`v2` are `PeriodicVertex{D}` to obtain the identifier of the edge and store a new
identifier if absent, or by `kp[i]` where `i` is an `Integer` to obtain the pair of
vertices corresponding to identifier `i`.

`kp` is built by calling `EdgeDict(g)` where `g` is a `PeriodicGraph`.
"""
struct EdgeDict{D}
    direct::Vector{VertexPair{D}}
    reverse::Dict{VertexPair{D},Int}

    function EdgeDict(g::PeriodicGraph{D}) where D
        known_pairs = VertexPair{D}[]
        known_pairs_dict = Dict{VertexPair{D},Int}()
        hintsize = (3^D*ne(g)^2) ÷ 5
        sizehint!(known_pairs, hintsize)
        sizehint!(known_pairs_dict, hintsize)
        return new{D}(known_pairs, known_pairs_dict)
    end
end

Base.length(kp::EdgeDict) = length(kp.direct)
function Base.get!(kp::EdgeDict{D}, x::VertexPair{D}) where D
    n = length(kp) + 1
    pairid = get!(kp.reverse, x, n)
    pairid == n && push!(kp.direct, x)
    pairid
end

Base.getindex(kp::EdgeDict, i::Integer) = kp.direct[i]
Base.getindex(kp::EdgeDict, x::VertexPair) = kp.reverse[x]

"""
    unique_order(cycles)

Return a sublist `I` of indices of `cycles` such that `cycles[I]` is the list of unique
elements of `cycles`, sorted by length first and by value next.
"""
function unique_order(cycles)
    isempty(cycles) && return Int[]
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


"""
    convert_to_ering!(buffer::Vector{Int}, ring::Vector{PeriodicVertex{D}}, kp::EdgeDict{D}, ofs, len) where D

Return the edge ring corresponding to `OffsetVertexIterator{D}(ring[1:len], ofs)` inside
`buffer`, resized to length `len`.

See also [`PeriodicGraphs.convert_to_ering`](@ref).
"""
function convert_to_ering!(buffer::Vector{Int}, ring::Vector{PeriodicVertex{D}}, kp::EdgeDict{D}, ofs, len) where D
    resize!(buffer, len)
    _last_p = ring[len]
    last_p = PeriodicVertex{D}(_last_p.v, _last_p.ofs .+ ofs)
    for j in 1:len
        _new_p = ring[j]
        new_p = PeriodicVertex{D}(_new_p.v, _new_p.ofs .+ ofs)
        buffer[j] = get!(kp, minmax(last_p, new_p))
        last_p = new_p
    end
    sort!(buffer)
end

"""
    convert_to_ering(ring::Vector{PeriodicVertex{D}}, kp::EdgeDict{D}, ofs=zero(SVector{D,Int}), len=length(ring)) where D

Return the edge ring corresponding to `OffsetVertexIterator{D}(ring[1:len], ofs)`.

See also [`PeriodicGraphs.convert_to_ering!`](@ref) to provide a buffer.
"""
function convert_to_ering(ring::Vector{PeriodicVertex{D}}, kp::EdgeDict{D}, ofs=zero(SVector{D,Int}), len=length(ring)) where D
    return convert_to_ering!(Int[], ring, kp, ofs, len)
end


"""
    sort_cycles(g::PeriodicGraph{D}, rs, depth=maximum(length, rs), kp=EdgeDict(g)) where D

Given a list `rs` of rings of `g` with at least one vertex in the origin unit cell, return
a list of edge-rings that cross a unit cell up to distance 2 from the origin, as well as
the indices of the corresponding rings in `rs`.

An edge-ring is a `Vector{Int}` where each element is the index of a pair of
`PeriodicVertex{D}` stored in `kp`, representing an edge between the two vertices.
Two consecutive edges in an edge-ring must share exactly one of their two vertices.
"""
function sort_cycles(g::PeriodicGraph{D}, rs, depth=maximum(length, rs; init=0), kp=EdgeDict(g)) where D
    ofss = cages_around(g, 2)
    tot_ofss = length(ofss)
    cycles = Vector{Vector{Int}}(undef, tot_ofss * length(rs))
    origin = zeros(Int, length(cycles))
    ringbuffer = Vector{PeriodicVertex{D}}(undef, 2*depth+3)
    zero_ofs = (tot_ofss+1) ÷ 2
    for (i, ring) in enumerate(rs)
        base = (i-1)*tot_ofss
        @simd for k in 1:length(ring)
            ringbuffer[k] = reverse_hash_position(ring[k], g)
        end
        origin[base+zero_ofs] = i
        for (i_ofs, ofs) in enumerate(ofss)
            cycles[base+i_ofs] = convert_to_ering(ringbuffer, kp, ofs, length(ring))
        end
    end
    I = unique_order(cycles)
    return cycles[I], origin[I]
end

# Complexity O(len(a) + len(b))
"""
    symdiff_cycles!(c::Vector{T}, a::Vector{T}, b::Vector{T}) where T

Like [`PeriodicGraphs.symdiff_cycles`](@ref) but stores the result in `c`.

`c` will be resized accordingly so its initial length does not matter.
"""
function symdiff_cycles!(c::Vector{T}, a::Vector{T}, b::Vector{T}) where T
    lenb = length(b)
    lena = length(a)
    if lena < lenb
        a, b = b, a
        lena, lenb = lenb, lena
    end
    counter_a = 0
    @inbounds while true
        if lenb == counter_a
            newlenc = (lena - counter_a) % UInt
            resize!(c, newlenc)
            unsafe_copyto!(c, 1, a, counter_a+1, newlenc)
            return c
        end
        counter_a += 1
        a[counter_a] == b[counter_a] || break
    end
    counter_a -= 1
    newlen = (lenb + lena - counter_a) % UInt
    length(c) < newlen && @inbounds resize!(c, newlen)
    counter_b = counter_a + 1
    y = lenb == 0 ? typemax(Int) : (@inbounds b[1+counter_a])
    j = 1
    n = length(a)
    @inbounds while counter_a < n
        counter_a += 1
        x = a[counter_a]
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
                counter_a += 1
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
    remaining_towritea = lena - counter_a + 1
    remaining_towritea > 0 && unsafe_copyto!(c, j, a, counter_a, remaining_towritea)
    @inbounds resize!(c, (j + remaining_towritea - 1) % UInt)
    return c
end

"""
    symdiff_cycles(a, b)

Symmetric difference between two sorted lists `a` and `b`.
Return the sorted list of elements belonging to `a` or to `b` but not to both.

Use [`PeriodicGraphs.symdiff_cycles!`](@ref) to provide a pre-allocated destination.

## Example
```jldoctest
julia> PeriodicGraphs.symdiff_cycles([3,4,5], [4,5,6])
2-element Vector{Int64}:
 3
 6
```
"""
symdiff_cycles(a, b) = symdiff_cycles!(Vector{Int}(undef, length(b) + length(a) - 1), a, b)


"""
    IterativeGaussianElimination{T}

Struct containing the list of sparse columns of the matrix under gaussian elimination on
the 𝔽₂ finite field.

To be used with [`PeriodicGraphs.gaussian_elimination!`](@ref) as one of the three concrete
types:
- `PeriodicGraphs.IterativeGaussianEliminationNone` for simple gaussian elimination,
- `PeriodicGraphs.IterativeGaussianEliminationLength` to detect when a new column can be
  expressed as a sum of strictly smaller columns of the matrix.
- `PeriodicGraphs.IterativeGaussianEliminationDecomposition` to detect when a new column
  can be expressed as a sum of other columns of the matrix and keep track of which.
"""
struct IterativeGaussianElimination{T}
    rings::Vector{Vector{Int}} # The rows of the matrix, in sparse format
    shortcuts::Vector{Int32} # verifies shortcuts[i] = 0 || rings[shortcuts[i]][1] == i
    buffer1::Vector{Int}
    buffer2::Vector{Int}
    track::T # If Vector{UInt8} : lengths of the encountered cycles. Otherwise, tracks the added cycles.
end

const IterativeGaussianEliminationLength = IterativeGaussianElimination{Vector{UInt8}}
const IterativeGaussianEliminationDecomposition = IterativeGaussianElimination{Tuple{Vector{Int32},Vector{Vector{Int32}}}}
const IterativeGaussianEliminationNone = IterativeGaussianElimination{Nothing}

function IterativeGaussianElimination{T}() where T
    track = T == Tuple{Vector{Int32},Vector{Vector{Int32}}} ? (Int32[],Vector{Int32}[]) : T()
    IterativeGaussianElimination{T}(Vector{Int}[], Int32[], Int[], Int[], track)
end
function IterativeGaussianElimination{T}(ring::Vector{Int}, sizehint=ring[1]) where T
    r1 = ring[1]
    shortcuts = zeros(Int, max(r1, sizehint))
    shortcuts[r1] = 1
    buffer = Vector{Int}(undef, length(r1))
    track = T == Vector{UInt8} ? [length(ring) % UInt8] :
            T == Tuple{Vector{Int32},Vector{Vector{Int32}}} ? (Int32[], [Int32[]]) : T()
    IterativeGaussianElimination{T}([ring], shortcuts, buffer, similar(buffer), track)
end

IterativeGaussianElimination(ring, sizehint=ring[1]) = IterativeGaussianEliminationNone(ring, sizehint)
IterativeGaussianElimination() = IterativeGaussianEliminationNone()

"""
    gaussian_elimination(gauss::IterativeGaussianElimination, r::Vector{Int})

Return `notindependent, info` where `notindependent` is `true` if `r` can be expressed as a
sum of vectors stored in `gauss`.

See [`PeriodicGraphs.gaussian_elimination!`](@ref) to store `r` in `gauss` if not, and for
more details dependending on the type of `gauss`.

Call `PeriodicGraphs.gaussian_elimination!(gauss, r, notindependent, info)` to obtain the
result of `PeriodicGraphs.gaussian_elimination!(gauss, r)` without duplicating computation.
"""
function gaussian_elimination(gauss::IterativeGaussianElimination{T}, r::Vector{Int}) where T
    rings = gauss.rings
    shortcuts = gauss.shortcuts
    buffer1::Vector{Int} = gauss.buffer1
    buffer2::Vector{Int} = gauss.buffer2
    lenshort = length(shortcuts)
    if gauss isa IterativeGaussianEliminationLength
        len = length(r) % UInt8
        lengths = gauss.track
    elseif gauss isa IterativeGaussianEliminationDecomposition
        track::Vector{Int32} = first(gauss.track)
        empty!(track)
    end
    r1 = r[1]
    maxlen = 0

    idx::Int32 = r1 > lenshort ? zero(Int32) : shortcuts[r1]
    if !iszero(idx)
        ridx = rings[idx]
        if gauss isa IterativeGaussianEliminationLength
            maxlen = lengths[idx]
        elseif gauss isa IterativeGaussianEliminationDecomposition
            push!(track, idx)
        end
        symdiff_cycles!(buffer1, r, ridx)
        isempty(buffer1) && return true, (r1, maxlen, buffer1)
        r1 = buffer1[1]
        idx = r1 > lenshort ? zero(Int32) : shortcuts[r1]
    else
        buffer1 = r
    end
    while !iszero(idx)
        ridx = rings[idx]
        if gauss isa IterativeGaussianEliminationLength
            maxlen = max(lengths[idx], maxlen)
        elseif gauss isa IterativeGaussianEliminationDecomposition
            push!(track, idx)
        end
        symdiff_cycles!(buffer2, buffer1, ridx)
        isempty(buffer2) && return true, (r1, maxlen, buffer1)
        r1 = buffer2[1]
        idx = r1 > lenshort ? zero(Int32) : shortcuts[r1]
        buffer2, buffer1 = buffer1, buffer2
    end

    return false, (r1, maxlen, buffer1)
end


# For an IterativeGaussianElimination of i rings of size at most l, symdiff_cycles cost at
# most (l + lenr) + (2l + lenr) + ... + (il + lenr) = O(i^2*l +i*lenr)
# If the size of the k-th ring is at most k*l, then the worst-case complexity is O(i^3*l)
# Calling gaussian_elimination ν times thus leads to a worst-case complexity in O(ν^4*l)
# The reasonning can probably be refined...
"""
    gaussian_elimination!(gauss::IterativeGaussianElimination, r::Vector{Int})

Test whether `r` can be expressed as a sum of vectors stored in `gauss`, and store `r` if
not. "sum" refers to the symmetric difference of boolean vectors, represented in sparse
format as the ordered list of non-zero indices.

If `gauss isa IterativeGaussianEliminationLength`, return whether `r` can be expressed as a
sum of strictly smaller vectors.

Otherwise, return `true` when `r` is a sum of any previously encoutered vectors.
If `gauss` isa `IterativeGaussianEliminationDecomposition`, query `retrieve_track(gauss)`
to obtain the sorted list of indices of such previously encountered vectors.

See also [`gaussian_elimination`](@ref) to test `r` without storing it.
"""
function gaussian_elimination!(gauss::IterativeGaussianElimination{T}, r::Vector{Int}) where T
    notindependent, info = gaussian_elimination(gauss, r)
    gaussian_elimination!(gauss, r, notindependent, info)
end

function gaussian_elimination!(gauss::IterativeGaussianElimination{T}, r::Vector{Int}, notindependent, (r1, maxlen, buffer1)) where T
    rings = gauss.rings
    shortcuts = gauss.shortcuts
    if gauss isa IterativeGaussianEliminationLength
        len = length(r) % UInt8
    elseif gauss isa IterativeGaussianEliminationDecomposition
        track::Vector{Int32} = first(gauss.track)
    end

    if notindependent
        gauss isa IterativeGaussianEliminationLength && return maxlen < len
        if gauss isa IterativeGaussianEliminationDecomposition
            push!(rings, buffer1) # dummy, used to keep track of dependent rings
            push!(last(gauss.track), track) # also dummy
        end
        return true
    end

    push!(rings, copy(buffer1))
    r1 > length(shortcuts) && append!(shortcuts, zero(Int32) for _ in 1:(r1-length(shortcuts)))
    shortcuts[r1] = length(rings) % Int32
    if gauss isa IterativeGaussianEliminationLength
        # the ring was not a sum of strictly smaller ones (since it's not a sum of previous ones at all)
        push!(gauss.track, len)
    elseif gauss isa IterativeGaussianEliminationDecomposition
        push!(last(gauss.track), copy(sort!(track)))
    end
    return false # the new ring was independent from the other ones
end

"""
    retrieve_track!([ret::Vector{Int32}, buffer::Vector{Int32},] gauss::IterativeGaussianEliminationDecomposition)

To be called consecutive to a call to `gaussian_elimination!(gauss, x)` that returned `true`.
In that case, `x` was found to be the sum of previously encountered vectors: return the
(reverse-sorted) list of their indices.

!!! warning
    Calling `retrieve_track!` after a call to `gaussian_elimination!` that returned `false`
    will produce an invalid result.
    Calling it twice will also produce an invalid result.
"""
function retrieve_track!(ret::Vector{Int32}, buffer::Vector{Int32}, gauss::IterativeGaussianEliminationDecomposition)
    track = sort!(first(gauss.track))
    tracks = last(gauss.track)
    empty!(ret)
    while !isempty(track)
        x = pop!(track)
        push!(ret, x)
        symdiff_cycles!(buffer, track, tracks[x])
        track, buffer = buffer, track
    end
    return ret
end
retrieve_track!(gauss::IterativeGaussianEliminationDecomposition) = retrieve_track!(Int32[], Int32[], gauss)


"""
    intersect_cycles!(c::Vector{T}, a::Vector{T}, b::Vector{T}) where T

Like [`PeriodicGraphs.intersect_cycles`](@ref) but stores the result in `c`.

`c` will be resized accordingly so its initial length does not matter as long as it is
at least as large as the resulting list.
"""
function intersect_cycles!(c::Vector{T}, a::Vector{T}, b::Vector{T}) where T
    isempty(b) && (empty!(c); return c)
    lenb = length(b)
    counter_b = 1
    xb = b[counter_b]
    j = 0
    for xa in a
        while xb < xa
            counter_b += 1
            counter_b > lenb && @goto ret
            xb = b[counter_b]
        end
        if xa == xb
            j += 1
            c[j] = xa
            counter_b += 1
            counter_b > lenb && @goto ret
            xb = b[counter_b]
        end
    end
    @label ret
    resize!(c, j)
    return c
end

"""
    intersect_cycles(a, b)

Intersection between two sorted lists `a` and `b`.
Return the sorted list of elements belonging to both `a` and `b`.

Use [`PeriodicGraphs.intersect_cycles!`](@ref) to provide a pre-allocated destination.

## Example
```jldoctest
julia> PeriodicGraphs.intersect_cycles([3,4,5], [4,5,6])
2-element Vector{Int64}:
 4
 5
```
"""
intersect_cycles(a, b) = intersect_cycles!(Vector{Int}(undef, min(length(a), length(b))), a, b)


"""
    union_cycles!(c::Vector{T}, a::Vector{T}, b::Vector{T}) where T

Like [`PeriodicGraphs.union_cycles`](@ref) but stores the result in `c`.

`c` will be resized accordingly so its initial length does not matter as long as it is
at least as large as the resulting list.
"""
function union_cycles!(c::Vector{T}, a::Vector{T}, b::Vector{T}) where T
    lenb = length(b)
    lena = length(a)
    counter_a = 0
    counter_b = 1
    y = lenb == 0 ? typemax(Int) : (@inbounds b[1])
    j = 1
    @inbounds while counter_a < lena
        counter_a += 1
        x = a[counter_a]
        while y < x
            c[j] = y
            j += 1
            counter_b += 1
            counter_b > lenb && @goto fillenda
            y = b[counter_b]
        end
        c[j] = x
        j += 1
        if y == x
            counter_b += 1
            if counter_b > lenb
                counter_a += 1
                @goto fillenda
            end
            y = b[counter_b]
        end
    end
    remaining_towriteb = lenb - counter_b + 1
    remaining_towriteb > 0 && unsafe_copyto!(c, j, b, counter_b, remaining_towriteb)
    @inbounds resize!(c, (j + remaining_towriteb - 1) % UInt)
    return c

    @label fillenda
    remaining_towritea = lena - counter_a + 1
    remaining_towritea > 0 && unsafe_copyto!(c, j, a, counter_a, remaining_towritea)
    @inbounds resize!(c, (j + remaining_towritea - 1) % UInt)
    return c
end

"""
    union_cycles(a, b)

Union between two sorted lists `a` and `b`.
Return the sorted list of elements belonging to `a` or `b` or both.

Use [`PeriodicGraphs.union_cycles!`](@ref) to provide a pre-allocated destination.

## Example
```jldoctest
julia> PeriodicGraphs.union_cycles([3,4,5], [4,5,6])
4-element Vector{Int64}:
 3
 4
 5
 6
```
"""
union_cycles(a, b) = union_cycles!(Vector{Int}(undef, length(a) + length(b)), a, b)

# function retrieve_vcycle(ecycle, known_pairs)
#     fst_pair = known_pairs[ecycle[1]]
#     last_pair = known_pairs[ecycle[end]]
#     last_v, other_v = fst_pair[1], fst_pair[2]
#     if last_v == last_pair[1] || last_v == last_pair[2]
#         last_v = fst_pair[2]
#         other_v = fst_pair[1]
#     end
#     n = length(ecycle)
#     ret = Vector{typeof(last_v)}(undef, n)
#     ret[1] = last_v
#     ret[n] = other_v
#     for j in 2:(n-1)
#         this_pair = known_pairs[ecycle[j]]
#         other_v = this_pair[1]
#         if last_v == other_v
#             last_v = this_pair[2]
#         else
#             last_v = other_v
#         end
#         ret[j] = last_v
#     end
#     return ret
# end

"""
    strong_erings([rs::Vector{Vector{Int}},], g::PeriodicGraph{D}, [depth=15,] ringsymms::AbstractSymmetryGroup=NoSymmetryGroup(length(rs))) where D

Compute the list of strong edge rings in `g`, up to length `2*depth+3`.
See [`strong_rings`](@ref) and [`rings`](@ref) for the meaning of the optional arguments.

Return a quadruplet of values:
- the two first values are the list of rings and their symmetry group, identical to the
  result of [`strong_rings`](@ref), unless `rs` is provided (see below).
- the third is the list of edge rings: each edge of the periodic graph is mapped to an
  integer and each ring is represented by the sorted list of its edges.
- the last is the mapping from edges to integers, given as an [`EdgeDict`](@ref).

If `rs` is provided, the first returned value is the list of indices `keep` of `rs` such
that `rs[keep]` is the list of strong rings.
"""
function strong_erings(rs::Vector{Vector{Int}}, g::PeriodicGraph{D}, depth=15, ringsymms::AbstractSymmetryGroup=NoSymmetryGroup(length(rs))) where D
    kp = EdgeDict(g)
    ecycles, origin = sort_cycles(g, rs, depth, kp)
    if isempty(ecycles)
        if ringsymms isa RingSymmetryGroup
            return Int[], RingSymmetryGroup{D}(Dict{Vector{Int},Int}(), Base.OneTo(0), nv(g), ringsymms.symms), Vector{Int}[], kp
        end
        return Int[], NoSymmetryGroup(0), Vector{Int}[], kp
    end
    fst_ring = first(ecycles)
    gauss = IterativeGaussianEliminationLength(fst_ring)
    ret = Int[]
    if ringsymms isa RingSymmetryGroup
        ringdict = Dict{Vector{Int},Int}()
        ringmap = zeros(Int, last(unique(ringsymms)))
    end
    orig_1 = first(origin)
    if orig_1 != 0
        r1 = rs[orig_1]
        push!(ret, 1)
        if ringsymms isa RingSymmetryGroup
            _rsymm1 = ringsymms(r1)
            rsymm1 = rs[_rsymm1]
            ringmap[_rsymm1] = 1
            ringdict[rsymm1] = 1
            if r1 != rsymm1
                ringdict[r1] = 1
            end
        end
    end
    counter = 1
    for (i, ecycle) in Iterators.drop(enumerate(ecycles), 1)
        onlysmallercycles = gaussian_elimination!(gauss, ecycle)
        onlysmallercycles && continue # cycle is a linear combination of smaller cycles
        if origin[i] != 0
            ri = rs[origin[i]]
            push!(ret, i)
            if ringsymms isa RingSymmetryGroup
                _rsymmi = ringsymms(ri)
                rsymmi = rs[_rsymmi]
                if ringmap[_rsymmi] == 0
                    counter += 1
                    ringmap[_rsymmi] = counter
                    ringdict[rsymmi] = counter
                end
                if ri != rsymmi
                    ringdict[ri] = ringmap[_rsymmi]
                end
            end
        end
    end
    keepat!(origin, ret)
    keepat!(ecycles, ret)
    if ringsymms isa RingSymmetryGroup
        return origin, RingSymmetryGroup{D}(ringdict, Base.OneTo(counter), nv(g), ringsymms.symms), ecycles, kp
    end
    return origin, NoSymmetryGroup(length(ret)), ecycles, kp
end

function strong_erings(g::PeriodicGraph, depth::Integer=15, symmetries::AbstractSymmetryGroup=NoSymmetryGroup(g), dist::DistanceRecord=DistanceRecord(g,depth))
    rs, symmg = rings(g, depth, symmetries, dist)
    keep, symms, erings, kp = strong_erings(rs, g, depth, symmg)
    return rs[keep], symms, erings, kp
end
function strong_erings(g::PeriodicGraph, symmetries::AbstractSymmetryGroup, dist::DistanceRecord=DistanceRecord(g,15))
    strong_erings(g, 15, symmetries, dist)
end

function strong_rings(rs::Vector{Vector{Int}}, g::PeriodicGraph{D}, depth=15, ringsymms::AbstractSymmetryGroup=NoSymmetryGroup(length(rs))) where D
    keep, symms = strong_erings(rs, g, depth, ringsymms)
    return rs[keep], symms
end

"""
    strong_rings([rs::Vector{Vector{Int}},] g::PeriodicGraph{D}, [depth::Integer=15,] symmetries::AbstractSymmetryGroup=NoSymmetryGroup(g), dist::DistanceRecord=DistanceRecord(g,depth)) where D

Compute the list of strong rings in `g`, up to length `2*depth+3`. Return them with their
symmetry group. Each ring is represented by the list of [`hash_position`](@ref) of its
vertices.

The optional first argument `rs` is the list of rings which can be provided if previously
computed.

See [`rings`](@ref) for the meaning of the other arguments.

A strong ring is a cycle of the graph which cannot be decomposed into a sum of any number
of smaller cycles. By comparison, a ring is a cycle which cannot be decomposed into a sum
of two smaller cycles. In particular, all strong rings are rings.

See also [`strong_erings`](@ref) to obtain the rings as a list of integers representing the
edges of the ring, instead of a list of integers representing its vertices.
"""
function strong_rings(g::PeriodicGraph, depth::Integer=15, symmetries::AbstractSymmetryGroup=NoSymmetryGroup(g), dist::DistanceRecord=DistanceRecord(g,depth))
    rs, symmg = rings(g, depth, symmetries, dist)
    return strong_rings(rs, g, depth, symmg)
end
function strong_rings(g::PeriodicGraph, symmetries::AbstractSymmetryGroup, dist::DistanceRecord=DistanceRecord(g,15))
    strong_rings(g, 15, symmetries, dist)
end


# trick for mutually-dependent types, see definition after that of RingAttributions
struct _RingIncluding{D,T} <: AbstractVector{OffsetVertexIterator{D}}
    ras::T
    i::Int
end

"""
    RingAttributions{D} <: AbstractVector{RingIncluding{D}}

Represent a set of rings of a `PeriodicGraph{D}`.

For `ra` of type `RingAttributions{D}`, `ra[i]` is a [`RingIncluding{D}`](@ref) object
representing the set of rings including `PeriodicVertex{D}(i)`.
"""
struct RingAttributions{D} <: AbstractVector{_RingIncluding{D,RingAttributions{D}}}
    rings::Vector{Vector{PeriodicVertex{D}}}
    attrs::Vector{Vector{Tuple{Int,Int}}}

    function RingAttributions{D}(n, rs::Vector{Vector{T}}) where {D,T<:Union{Integer,PeriodicVertex{D}}}
        keeprs = T == PeriodicVertex{D}
        rings = keeprs ? rs : [Vector{PeriodicVertex{D}}(undef, length(r)) for r in rs]
        attrs = [Tuple{Int,Int}[] for _ in 1:n]
        for (i, r) in enumerate(rs)
            ring = rings[i]
            for (j, x) in enumerate(r)
                u = T == PeriodicVertex{D} ? x : reverse_hash_position(x, n, Val(D))
                keeprs || (ring[j] = u)
                push!(attrs[u.v], (i, j))
            end
        end
        return new{D}(rings, attrs)
    end
end

"""
    RingAttributions(g::PeriodicGraph{D}, [strong::Bool=false,] [depth::Integer=15,] symmetries::AbstractSymmetryGroup=NoSymmetryGroup(g), dist::DistanceRecord=DistanceRecord(g,depth)) where D

Return the rings of `g` sorted into a `RingAttributions`.

If `strong` is set, only the strong rings are kept.
See [`rings`](@ref) for the meaning of the other arguments.
"""
function RingAttributions(g::PeriodicGraph{D}, strong::Bool=false, depth::Integer=15, symmetries::AbstractSymmetryGroup=NoSymmetryGroup(g), dist::DistanceRecord=DistanceRecord(g,depth)) where D
    rs, _ = (strong ? strong_rings : rings)(g, depth, symmetries, dist)
    return RingAttributions{D}(nv(g), rs)
end
function RingAttributions(g::PeriodicGraph{D}, depth::Integer, symmetries::AbstractSymmetryGroup=NoSymmetryGroup(g), dist::DistanceRecord=DistanceRecord(g,depth)) where D
    RingAttributions(g, false, depth, symmetries, dist)
end
function RingAttributions(g::PeriodicGraph{D}, strong::Bool, symmetries::AbstractSymmetryGroup, dist::DistanceRecord=DistanceRecord(g,15)) where D
    RingAttributions(g, strong, 15, symmetries, dist)
end
function RingAttributions(g::PeriodicGraph{D}, symmetries::AbstractSymmetryGroup, dist::DistanceRecord=DistanceRecord(g,15)) where D
    RingAttributions(g, false, 15, symmetries, dist)
end

Base.@propagate_inbounds function Base.getindex(ras::RingAttributions, i::Integer)
    @boundscheck checkbounds(ras.attrs, i)
    RingIncluding(ras, i)
end
Base.size(ras::RingAttributions) = size(ras.attrs)
Base.IndexStyle(::Type{RingAttributions{D}}) where {D} = Base.IndexLinear()

function Base.show(io::IO, ::MIME"text/plain", ras::RingAttributions)
    println(io, typeof(ras), "(rings per node: ", length.(ras.attrs), ')')
end

"""
    RingIncluding{D} <: AbstractVector{OffsetVertexIterator{D}}

The list of rings of a `PeriodicGraph{D}` including a particular vertex
`PeriodicVertex{D}(i)`.

The object is iterable and indexable by an integer: for `ri` of type `RingIncluding{D}`,
`ri[j]` is an iterable over the vertices of the `j`-th ring including vertex `i`.
"""
const RingIncluding{D} = _RingIncluding{D,RingAttributions{D}}
RingIncluding(ras::RingAttributions{D}, i) where {D} = RingIncluding{D}(ras, i)

function Base.getindex(ri::RingIncluding{D}, j::Integer) where {D}
    newring_idx, idx = ri.ras.attrs[ri.i][j]
    newring = ri.ras.rings[newring_idx]
    ofs = newring[idx].ofs
    return OffsetVertexIterator{D}(.-ofs, newring)
end
Base.size(ri::RingIncluding) = size(ri.ras.attrs[ri.i])
Base.IndexStyle(::Type{RingIncluding{D}}) where {D} = Base.IndexLinear()
