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

# minabove(x::ConstMiniBitSet, i::Integer) = trailing_zeros(x.x & ((~one(UInt64)) << (i % UInt8)))
# minaboveeq(x::ConstMiniBitSet, i::Integer) = trailing_zeros(x.x & ((~zero(UInt64)) << (i % UInt8)))
# nonemptyintersect(x::ConstMiniBitSet, y::ConstMiniBitSet) = !iszero(x.x & y.y)
hasonly(x::ConstMiniBitSet, i::Integer) = iszero(x.x & ~(one(UInt64) << (i % UInt8)))

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

# minabove(x::MiniBitSet, i::Integer) = minabove(ConstMiniBitSet(x), i)
# minaboveeq(x::MiniBitSet, i::Integer) = minaboveeq(ConstMiniBitSet(x), i)
# isnonemptyintersect(x::MiniBitSet, y::MiniBitSet) = !iszero(x.x & y.x)
# isnonemptysymdiff(x::MiniBitSet, y::MiniBitSet) = !iszero(x.x ⊻ y.x)
# hasonly(x::MiniBitSet, i::Integer) = iszero(x.x & ~(one(UInt64) << (i % UInt8)))
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
    shortroots::ConstMiniBitSet
    longroots::ConstMiniBitSet
    lastshort::Int
    num::Int # length of the shortest branch
    heads::Vector{Int}
    function JunctionNode(heads::AbstractVector, num, roots::ConstMiniBitSet)
        new(roots, ConstMiniBitSet(), length(heads), num, heads)
    end
end
JunctionNode() = JunctionNode(Int[], -1, ConstMiniBitSet())
JunctionNode(head::Integer, num, roots::ConstMiniBitSet) = JunctionNode(Int[head], num, roots)
JunctionNode(root::Integer) = JunctionNode(1, 0, ConstMiniBitSet(root))

function Base.show(io::IO, x::JunctionNode)
    print(io, "JunctionNode([")
    join(io, x.heads, ", ")
    print(io, "])")
end

const shortrootsoffset = fieldoffset(JunctionNode, findfirst(==(:shortroots), fieldnames(JunctionNode)))
const longrootsoffset = fieldoffset(JunctionNode, findfirst(==(:longroots), fieldnames(JunctionNode)))
const lastshortoffset = fieldoffset(JunctionNode, findfirst(==(:lastshort), fieldnames(JunctionNode)))

@inline unsafe_incr!(ptr) = unsafe_store!(ptr, unsafe_load(ptr)+1)
@inline unsafe_union!(ptr, val) = unsafe_store!(ptr, _constminibitset(val | unsafe_load(ptr).x))

function arcs_list(g::PeriodicGraph{D}, i, depth) where D
    vertexnums = PeriodicVertex{D}[PeriodicVertex{D}(i)]
    dag = JunctionNode[JunctionNode()]
    hintsize = ceil(Int, depth^2.9) # empirical estimate
    sizehint!(vertexnums, hintsize)
    sizehint!(dag, hintsize)
    append!(vertexnums, neighbors(g, i))
    if length(vertexnums) > 62
        error("The vertex has a degree too large (> 62) to be represented in a MiniBitSet. Please open an issue.")
    end
    append!(dag, JunctionNode(j) for j in 2:length(vertexnums))
    vertexdict = Dict{PeriodicVertex{D},Int}((x => j) for (j,x) in enumerate(vertexnums))
    counter = 1
    for parents in Iterators.rest(dag, 2)
        parents.num == depth && break
        counter += 1
        previous = vertexnums[parents.heads[1]]
        num = 1 + parents.num
        for x in neighbors(g, vertexnums[counter])
            x == previous && continue
            idx = get!(vertexdict, x, length(dag)+1)
            if idx == length(dag)+1
                push!(dag, JunctionNode(counter, num, parents.shortroots))
                push!(vertexnums, x)
            elseif idx > counter
                junction = dag[idx]
                push!(junction.heads, counter)
                # We use pointer and unsafe_store! to modify the immutable type JunctionNode
                # This is only possible because it is stored in a Vector, imposing fixed addresses.
                ptr = pointer(dag, idx)
                if num == junction.num
                    unsafe_incr!(Ptr{Int}(ptr + lastshortoffset))
                    unsafe_union!(Ptr{ConstMiniBitSet}(ptr + shortrootsoffset), parents.shortroots.x)
                else
                    unsafe_union!(Ptr{ConstMiniBitSet}(ptr + longrootsoffset), parents.shortroots.x)
                end
            end
        end
    end
    return dag, vertexnums, vertexdict
end


function next_compatible_arc!(buffer, last_positions, idx_stack, dag, check, dist, vertexnums)
    num_choice = length(idx_stack)
    incompatibleflag = true
    updated_last_positions = last_positions
    haslongarc = isodd(length(buffer))
    num = num_choice + 2 - haslongarc
    @assert buffer[num] == 1
    @assert num == length(buffer) ÷ 2
    root = buffer[num+1] # unused if !check
    while incompatibleflag
        trail = trailing_ones(updated_last_positions) % UInt8
        next_toupdate = num_choice - trail
        next_toupdate < 1 && return ~zero(UInt64)
        updated_last_positions = updated_last_positions & (~zero(UInt64) << trail)
        incompatibleflag = false
        idx = idx_stack[next_toupdate] + 1
        head = buffer[check ? next_toupdate : end-next_toupdate]
        parents = dag[head]
        for i in next_toupdate:num_choice
            heads = parents.heads
            head = heads[idx]
            _i = i + 1
            last_head = parents.lastshort
            if check
                next_parents = dag[head]
                while hasonly(next_parents.shortroots, root) ||
                      head == buffer[end-_i-haslongarc] ||
                      (dist !== nothing && # checking for rings instead of cycles
                      (is_distance_smaller!(dist, vertexnums[buffer[num+_i]], vertexnums[head], num) ||
                      (haslongarc && is_distance_smaller!(dist, vertexnums[buffer[num+_i+1]], vertexnums[head], num))))
                    if idx == last_head
                        updated_last_positions |= (~zero(UInt64) >> ((63 - num_choice + i) % UInt8))
                        incompatibleflag = true
                        break
                    end
                    idx += 1
                    head = heads[idx]
                    next_parents = dag[head]
                end
                incompatibleflag && break
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
    end
    return updated_last_positions
end

function initial_compatible_arc!(buffer, idx_stack, dag, check, dist, vertexnums)
    last_positions = zero(UInt64)
    num_choice = length(idx_stack)
    head = buffer[check ? 1 : end-1]
    haslongarc = isodd(length(buffer))
    num = num_choice + 2 - haslongarc
    @assert buffer[num] == 1
    @assert num == length(buffer) ÷ 2
    root = buffer[num+1] # unused if !check
    parents = dag[head]
    for i in 1:num_choice
        heads = parents.heads
        head = heads[1]
        _i = i + 1
        idx = 1
        last_head = parents.lastshort
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
    shortest_n = parents.num
    parents_lastshort = parents.lastshort
    shortest_n == 0 && return _case_0(parents.heads, parents_lastshort, midnode)
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
                  (is_distance_smaller!(dist, vertexnums[buffer[shortest_n+2]], vertexnums[head2], shortest_n+1) ||
                  (i1 > parents_lastshort && is_distance_smaller!(dist, vertexnums[buffer[shortest_n+3]], vertexnums[head2], shortest_n+1)))
                    continue
                end
                last_positions2 = initial_compatible_arc!(buffer, idx_stack2, dag, true, dist, vertexnums)
                while last_positions2 != ~zero(UInt64)
                    # if buffer[shortest_n+1:end] == [1, 4, 12, 24]
                    #     @show getindex.((args[2],), buffer)
                    #     @show buffer
                    # end
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
    distances::Dict{Tuple{Int,Int},Int}
    Q_lists::Vector{Vector{Tuple{PeriodicVertex{D},Int}}}
    first_indices::Vector{Int}
    seens::Vector{Set{Int}}
end
function DistanceRecord(g::PeriodicGraph{D}) where D
    distances = Dict{Tuple{Int,Int},Int}()
    n = nv(g)
    Q_lists = Vector{Vector{Tuple{PeriodicVertex{D},Int}}}(undef, n)
    first_indices = zeros(Int, n)
    seens = Vector{Set{Int}}(undef, n)
    return DistanceRecord{D}(g, distances, Q_lists, first_indices, seens)
end

function known_distance(dist::DistanceRecord, i, j)
    get(dist.distances, minmax(i,j), 0)
end
function set_distance!(dist::DistanceRecord, i, j, d)
    dist.distances[minmax(i,j)] = d
end

function bfs_smaller!(dist, i, j, start, stop)
    Q = dist.Q_lists[i]
    seen = dist.seens[i]
    graph = dist.g
    counter = start
    encountered = false
    n = nv(graph)
    for (u, d) in Iterators.rest(Q, start)
        (encountered || d+1 ≥ stop) && break
        counter += 1
        for x in neighbors(graph, u)
            h = hash_position(x, n)
            h ∈ seen && continue
            push!(seen, h)
            push!(Q, (x, d+1))
            set_distance!(dist, i, h, d+1)
            encountered |= h == j
        end
    end
    dist.first_indices[i] = counter
    if !encountered
        # set_distance!(dist, i, j, stop)
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
    dist.Q_lists[i] = [(x, 1) for x in neighbors(g, i)]
    _seen = Set{Int}(i)
    encountered = false
    n = nv(g)
    for x in neighbors(g, i)
        _h = hash_position(x, n)
        set_distance!(dist, i, _h, 1)
        push!(_seen, _h)
        encountered |= _h == j
    end
    dist.seens[i] = _seen
    encountered
end

function is_distance_smaller!(dist, verti, vertj, stop)
    i, j, start = _reorderinit(dist.first_indices, verti, vertj)
    known = known_distance(dist, i, j)
    known == 0 || return known < stop
    if start == 0
        init_distance_record!(dist, i, j) && return 1 < stop
        start = 1
    end
    return bfs_smaller!(dist, i, j, start, stop)
end


function rings_around_old(g::PeriodicGraph{D}, i, depth=20) where D
    n = nv(g)
    dist = DistanceRecord(g)
    dag, vertexnums, _ = arcs_list(g, i, depth)
    hashes = [hash_position(x, n) for x in vertexnums]
    ret = Vector{Int}[]
    for midnode in 2:length(dag)
        for cycle in cycles_ending_at(dag, midnode)#, dist, vertexnums)
            invalidcycle = false
            lenc = length(cycle)
            num = lenc ÷ 2
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
            hashcycle = hashes[cycle]
            fst = argmin(hashcycle)
            endreverse = hashcycle[mod1(fst+1, lenc)] < hashcycle[mod1(fst-1, lenc)]
            if fst != 1 || !endreverse
                reverse!(hashcycle, 1, fst - endreverse)
                reverse!(hashcycle, fst - endreverse + 1, lenc)
                endreverse && reverse!(hashcycle)
            end
            push!(ret, hashcycle)
        end
    end
    return ret
end


function rings_around(g::PeriodicGraph{D}, i, depth=20) where D
    n = nv(g)
    dist = DistanceRecord(g)
    dag, vertexnums, _ = arcs_list(g, i, depth)
    hashes = [hash_position(x, n) for x in vertexnums]
    ret = Vector{Int}[]
    for midnode in 2:length(dag)
        for cycle in cycles_ending_at(dag, midnode, dist, vertexnums)
            hashcycle = hashes[cycle]
            fst = argmin(hashcycle)
            lenc = length(cycle)
            endreverse = hashcycle[mod1(fst+1, lenc)] < hashcycle[mod1(fst-1, lenc)]
            if fst != 1 || !endreverse
                reverse!(hashcycle, 1, fst - endreverse)
                reverse!(hashcycle, fst - endreverse + 1, lenc)
                endreverse && reverse!(hashcycle)
            end
            push!(ret, hashcycle)
        end
    end
    return ret
end


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
