import Base.GMP.MPZ
import LinearAlgebra: det
using StaticArrays

export dimensionality, split_catenation

# Algorithm 3 from George Havas, Bohdan S. Majewski, and Keith R. Matthews,
# "Extended GCD and Hermite Normal Form Algorithms via Lattice Basis Reduction"

function gcd_reduce1!(k, i, m, λ, D, a, B, tmpa, tmpD, tmpB, tmpλ)
    q = @inbounds begin
        if !iszero(a[i])
            div(a[k], a[i], RoundNearest)
        elseif 2*λ[i,k] > D[i+1]
            div(λ[i,k], D[i+1], RoundNearest)
        else
            zero(BigInt)
        end
    end

    @inbounds if !iszero(q)
        MPZ.mul!(tmpa, q, a[i])
        MPZ.sub!(a[k], tmpa)
        @simd for j in 1:m
            MPZ.mul!(tmpB[j], q, B[j,i])
            MPZ.sub!(B[j,k], tmpB[j])
        end
        MPZ.mul!(tmpD, q, D[i+1])
        MPZ.sub!(λ[i,k], tmpD)
        @simd for j in 1:i-1
            MPZ.mul!(tmpλ[j], q, λ[j,i])
            MPZ.sub!(λ[j,k], tmpλ[j])
        end
    end
    nothing
end

function gcd_swap!(k, m, λ, D, a, B, tmpa, tmpD, tmpB, tmpλ, tmpt)
    @inbounds begin
        a[k], a[k-1] = a[k-1], a[k]
        B[:,k], B[:,k-1] = B[:,k-1], B[:,k]
        @simd for j in 1:k-2
            λ[j,k], λ[j,k-1] = λ[j,k-1], λ[j,k]
        end
        @simd for i in k+1:m
            MPZ.mul!(tmpB[i], λ[k-1,i], D[k+1])
            MPZ.mul!(tmpλ[i], λ[k,i], λ[k-1,k])
            MPZ.sub!(tmpt[i], tmpB[i], tmpλ[i])
            MPZ.mul!(tmpB[i], λ[k-1,i], λ[k-1,k])
            MPZ.mul!(tmpλ[i], λ[k,i], D[k-1])
            MPZ.add!(tmpλ[i], tmpB[i])
            MPZ.cdiv_q!(λ[k-1,i], tmpλ[i], D[k])
            MPZ.cdiv_q!(λ[k,i], tmpt[i], D[k])
        end
        MPZ.mul!(tmpD, D[k-1], D[k+1])
        MPZ.pow_ui!(tmpa, λ[k-1,k], UInt(2))
        MPZ.add!(tmpD, tmpa)
        MPZ.cdiv_q!(D[k], tmpD, D[k])
    end
    nothing
end

"""
    extended_gcd(s)::Tuple{BigInt, Vector{BigInt}}

Given a list of integers, return a tuple `(d, coefs)` where `d` is the gcd of these
integers and `coefs` is a list of integers such that `dot(s, coefs) == d`.
"""
function extended_gcd(s)
    m = length(s)
    B = BigInt[i==j for i in 1:m, j in 1:m] # LinearAlgebra.I defaults to fill
    λ = BigInt[0 for _ in 1:m, _ in 1:m] # zeros defaults to fill
    # λ = zeros(Int, m, m)
    D = BigInt[1 for _ in 1:m+1]
    a = BigInt.(abs.(s))

    tmpa = BigInt()
    tmpD = BigInt()
    tmpB = [BigInt() for _ in 1:m]
    tmpλ = [BigInt() for _ in 1:m]
    tmpt = [BigInt() for _ in 1:m]

    k = 2
    @inbounds while k <= m
        gcd_reduce1!(k, k-1, m, λ, D, a, B, tmpa, tmpD, tmpB, tmpλ)
        if !iszero(a[k-1]) || ((iszero(a[k-1]) & iszero(a[k])) &&
                               4*(D[k-1]*D[k+1] + λ[k-1,k]^2) < 3*D[k]^2)
            gcd_swap!(k, m, λ, D, a, B, tmpa, tmpD, tmpB, tmpλ, tmpt)
            k -= (k > 2)
        else
            for i in k-2:-1:1
                gcd_reduce1!(k, i, m, λ, D, a, B, tmpa, tmpD, tmpB, tmpλ)
            end
            k += 1
        end
    end
    @inbounds if a[end] < 0
        MPZ.neg!(a[end])
        foreach(MPZ.neg!, @view B[:,end])
    end
    i = 1
    @inbounds for x in s
        if x < 0
            MPZ.neg!(B[i,end])
        end
        i+=1
    end
    return @inbounds (a[end], B[:,end])
end


function subdiagonalenum(::Val{N}) where N
    ret = Tuple{Int,Int}[]
    for j in 1:N-1
        for i in j+1:N
            push!(ret, (i,j))
        end
    end
    return ret
end

"""
    normal_basis(l::AbstractVector{<:StaticVector{N,T}}) where {N,T<:Integer}

Given a list of integer vectors of dimension `N`, return a tuple `(mat, D)` where
`D` is the dimension of the space spanned by the input vectors, and `mat` is an
invertible matrix whose `D` first columns form a basis of this spanned space, which
does not depend on the exact input.

If `D ≠ N`, the remaining columns are set so that `mat` be invertible. These additional
columns will only contain one coefficient equal to 1, and all others to 0. No other
assumption should be made about these columns; in particular, they may depend on the input.

!!! warning
    If modifiable, the input list will be modified in-place during this process.
"""
function normal_basis(l::AbstractVector{<:StaticVector{N,T}}) where {N,T<:Integer}
    n = length(l)
    # coords = [MVector{N,T}(undef) for _ in 1:n]
    n == 0 && return (MMatrix{N,N,T}(LinearAlgebra.I), 0)
    basis = zero(MMatrix{N,N,T}) # The artificial integer basis
    zerocol = Int[]
    @inbounds for j in 1:N
        d, coefs = Tuple{T,Vector{T}}(extended_gcd(x[j] for x in l))
        for i in 1:n
            basis[:,j] .+= coefs[i] .* l[i]
        end
        @assert basis[j,j] ≥ 0
        if iszero(basis[j,j])
            push!(zerocol, j)
            continue
        end
        for i in 1:n
            coord = l[i][j] .÷ d
            # coords[i][j] = coord
            l[i] .-= coord .* basis[:,j]
        end
    end
    @assert all(iszero, l)
    # normalization = MMatrix{N,N,T}(LinearAlgebra.I)
    @inbounds for (i,j) in subdiagonalenum(Val(N))
        # normalization[i,j] = (basis[i,j] ÷ basis[i,i]) - signbit(ret[i,j])
        # intbasis[:,j] .-= normalization[i,j] * basis[:,i]
        iszero(basis[i,i]) && continue
        x = signbit(basis[i,j])
        basis[:,j] .-= (((basis[i,j] + x) ÷ basis[i,i]) - x) * basis[:,i]
    end
    # In this last step, we simply push the empty columns at the end of the basis
    iend = N
    m = N - length(zerocol)
    @inbounds while !isempty(zerocol)
        if zerocol[end] == iend
            basis[iend,iend] = 1
            iend -= 1
            pop!(zerocol)
            continue
        end
        istart = popfirst!(zerocol)
        @simd for k in 1:N
            basis[k,istart] = basis[k,iend]
        end
        basis[:,iend] = 0
        basis[istart,iend] = 1
        iend -= 1
    end
    return (basis, m)
end

function normal_basis(l::AbstractVector{SVector{N,T}}) where {N,T<:Integer}
    return normal_basis([MVector{N,T}(x) for x in l])
end


function _dimensionality(g::PeriodicGraph{0})
    return Dict{Int,Vector{Tuple{Vector{Int},Vector{SVector{0,Int}}}}}(0 => [(collect(vertices(g)), zeros(Int, 0))])
end
function _dimensionality(g::PeriodicGraph{N}) where N
    ret = Dict{Int,Vector{Tuple{Vector{Int},Vector{SVector{N,Int}}}}}()
    n = nv(g)
    visited = falses(n)
    nullofs = zero(SVector{N,Int})
    for i in vertices(g)
        visited[i] && continue
        recordedperiodicities = Set{SVector{N,Int}}()
        component = Dict{Int,SVector{N,Int}}(i => nullofs)
        seen = Set{PeriodicVertex{N}}([PeriodicVertex{N}(i)])
        Q = PeriodicVertex{N}[(PeriodicVertex{N}(i))]
        for src in Q
            @assert !visited[src.v]
            for dst in outneighbors(g, src.v)
                dst = PeriodicVertex{N}(dst.v, dst.ofs .+ src.ofs)
                dst ∈ seen && continue
                push!(seen, dst)
                lastperiodicity = get!(component, dst.v, dst.ofs)
                thisperiodicity = dst.ofs .- lastperiodicity
                c = cmp(thisperiodicity, nullofs)
                if iszero(c) # First time we encounter dst.v
                    push!(Q, dst)
                elseif c > 0
                    push!(recordedperiodicities, thisperiodicity)
                end
            end
        end
        newcomponent = collect(keys(component))
        visited[newcomponent] .= true
        vec::Vector{SVector{N,Int}} = collect(recordedperiodicities)
        dim::Int = rank(reduce(hcat, vec; init=reshape(nullofs, N, 1)))
        x = get!(ret, dim, Tuple{Vector{Int},Vector{SVector{N,Int}}}[])
        push!(x, (newcomponent, vec))
    end
    return ret
end


"""
    dimensionality(g::PeriodicGraph{N}) where N

Determine the actual dimension of each connected component of `g`.
Return a dictionary where each entry `n => [(l1,m1), (l2,m2), ...]` means that
`li` is a list of vertices that form a connected component of dimension `n`, and that
component is present `mi` times per unit cell.

In other words, the connected component `li` has a periodicity that can only be expressed
in a unit cell `mi` times larger than the current one. See also [`split_catenation`](@ref).

## Examples
```jldoctest
julia> dimensionality(PeriodicGraph("2   1 1  1 0   2 2  -1 1   3 3  1 0   3 3  0 1"))
Dict{$Int, Vector{Tuple{Vector{$Int}, $Int}}} with 2 entries:
  2 => [([3], 1)]
  1 => [([1], 1), ([2], 1)]

julia> dimensionality(PeriodicGraph("1   1 1  2"))
Dict{$Int, Vector{Tuple{Vector{$Int}, $Int}}} with 1 entry:
  1 => [([1], 2)]
```
"""
function dimensionality(g::PeriodicGraph)
    dim = _dimensionality(g)
    ret = Dict{Int,Vector{Tuple{Vector{Int},Int}}}()
    @inbounds for (d, x) in dim
        retd = Vector{Tuple{Vector{Int},Int}}(undef, length(x))
        for (k, (l, vec)) in enumerate(x)
            nfoldmat, d2 = normal_basis(vec)
            @assert d2 == d
            nfold = abs(det(nfoldmat))
            retd[k] = (l, nfold)
        end
        ret[d] = retd
    end
    return ret
end


"""
    PeriodicGraph{D}(graph::PeriodicGraph{N}) where {D,N}

Return a graph that has the same structural information as the input `graph` but
embedded in `D` dimensions instead of `N`. It will fail if the dimensionality of
the graph is greater than `D`.

!!! note
    The behaviour is undefined if `D < N` and there are multiple non-identical
    connected components. In particular, the function is expected to fail if these
    components do not share the same orientation.
"""
function PeriodicGraph{D}(graph::PeriodicGraph{N}, dims=_dimensionality(graph)) where {D,N}
    # We handle the trivial cases (where D ≥ N) first
    N == D && return graph
    n = ne(graph)
    if D > N
        edgs = Vector{PeriodicEdge{D}}(undef, n)
        @inbounds for (i, e) in enumerate(edges(graph))
            _ofs = zero(MVector{D,Int})
            _ofs[1:N] .= e.dst.ofs
            edgs[i] = PeriodicEdge{D}(e.src, e.dst.v, SVector{D,Int}(_ofs))
        end
        return PeriodicGraph{D}(nv(graph), edgs)
    end
    # If D ≤ N, we have to find a new D-dimensional basis on which to project the edges
    d = maximum(keys(dims))
    if d > D
        throw(DimensionMismatch("Graph of dimensionality $d cannot be reduced to dimension $D < $d."))
    end
    l::Vector{SVector{N,Int}} = @inbounds (first(dims[d])[2])
    # The vectors of l span the periodicity of the graph
    basis, _d = normal_basis(l)
    @assert _d == d
    _invmat::Matrix{Rational{Int}} = inv(Rational.(basis))[1:d,:]
    maxden = maximum(denominator.(_invmat); dims=2)
    invmat::Matrix{Int} = maximum(maxden; init=0) > 1 ? (maxden .* _invmat) : _invmat
    newedges = Vector{PeriodicEdge{D}}(undef, n)
    if d == D
        for (i, e) in enumerate(edges(graph))
            newedges[i] = PeriodicEdge{D}(e.src, e.dst.v, invmat * e.dst.ofs)
        end
    else # d < D
        for (i, e) in enumerate(edges(graph))
            eofs = zero(MVector{D,Int})
            eofs[1:d] = invmat * e.dst.ofs
            newedges[i] = PeriodicEdge{D}(e.src, e.dst.v, eofs)
        end
    end
    return PeriodicGraph{D}(nv(graph), newedges)
end


function _explore_one_component!(expected, encountered, visited, graph::PeriodicGraph{N}, start) where N
    nullofs = zero(SVector{N,Int})
    recordedperiodicities = Set{SVector{N,Int}}()
    component = Dict{Int,Tuple{Int,SVector{N,Int}}}(start.v => (1, start.ofs))
    seen = Set{PeriodicVertex{N}}([start])
    newedges = PeriodicEdge{N}[]
    Q = PeriodicVertex{N}[start]
    for (j, src) in enumerate(Q)
        for dst in outneighbors(graph, src)
            newvertex = length(component) + 1
            k, lastperiodicity = get!(component, dst.v, (newvertex, dst.ofs))
            thisperiodicity = dst.ofs .- lastperiodicity
            push!(newedges, PeriodicEdge(j, k, thisperiodicity))
            dst ∈ seen && continue
            push!(seen, dst)
            visited[dst.v] = true
            c = cmp(thisperiodicity, nullofs)
            if iszero(c) # First time we encounter dst.v
                @assert k == newvertex
                push!(Q, dst)
            elseif c > 0
                push!(recordedperiodicities, thisperiodicity)
            end
        end
    end
    vec::Vector{SVector{N,Int}} = collect(recordedperiodicities)
    catenationmat, d = normal_basis(vec)
    invmcatenationmat = inv(SMatrix{N,N,Rational{Int},N*N}(catenationmat))
    for (j, (src, (dst, ofs))) in enumerate(newedges)
        newedges[j] = PeriodicEdge{N}(src, dst, Int.(invmcatenationmat*ofs))
    end

    idx = minimum(keys(component))
    @assert expected[idx] == 0
    expected[idx] = Int(det(catenationmat))
    encountered[idx] = (PeriodicGraph{N}(newedges), [OffsetVertexIterator(nullofs, Q)], catenationmat, d)
    idx
end

"""
    split_catenation(graph::PeriodicGraph{N}) where N

Return a list of tuples `(subgraph, vmaps, mat, dim)` where `subgraph` is a
connected component of the input `graph`.
Each sublist `vmap` in `vmaps` is a distinct mapping of the vertices of the graph into
those of the corresponding `subgraph`. The union of all the vertices in all the `vmaps`,
repeated with the according periodicity, exactly covers all the vertices of the input
graph, without repetition.
`dim` is the [`dimensionality`](@ref) of these connected components, and `mat` is the
transformation matrix between the input cell and the component's cell.

Each sublist contains connected components that share the same vertex indices.

## Example
```jldoctest
julia> g = PeriodicGraph("2    1 1  3 0   1 1  0 1   2 3  1 0   3 2  1 0");

julia> dimensionality(g)
Dict{$Int, Vector{Tuple{Vector{$Int}, $Int}}} with 2 entries:
  2 => [([1], 3)]
  1 => [([2, 3], 2)]

julia> splits = split_catenation(g);

julia> last.(splits) # One 2-dimensional catenation (vertex 1), one 1-dimensional (vertices 2 and 3)
2-element Vector{$Int}:
 2
 1

julia> subgraph, vmaps, mat, dim = last(splits);  # the 1-dimensional connected components

julia> vmaps  # first sublist takes vertex 2 from the reference unit cell, second one takes vertex 3.
2-element Vector{PeriodicGraphs.OffsetVertexIterator{2}}:
 [(2, (0,0)), (3, (-1,0))]
 [(2, (1,0)), (3, (0,0))]

julia> subgraph  # the subgraph consists in only the two relevant two vertices
PeriodicGraph2D(2, PeriodicEdge2D[(1, 2, (0,0)), (1, 2, (1,0))])

julia> mat  # the new cell is twice as large as the initial one
2×2 SMatrix{2, 2, $Int, 4} with indices SOneTo(2)×SOneTo(2):
 2  0
 0  1
```
"""
function split_catenation(graph::PeriodicGraph{N}) where N
    n = nv(graph)
    visited = falses(n)
    expected = zeros(Int, n)
    encountered = Vector{Tuple{PeriodicGraph{N}, Vector{OffsetVertexIterator{N}}, SMatrix{N,N,Int,N*N}, Int}}(undef, n)

    for i in 1:n
        visited[i] && continue
        visited[i] = true
        _explore_one_component!(expected, encountered, visited, graph, PeriodicVertex{N}(i))
    end

    keep = Int[]
    for (idx, ex) in enumerate(expected)
        ex == 0 && continue
        push!(keep, idx)
        _, vmaps, mat, _ = encountered[idx]
        if length(vmaps) < ex
            refmap = first(vmaps).nlist
            invmat = inv(SMatrix{N,N,Rational{Int},N*N}(mat))
            dumb_indices = findall(i -> begin onei = zeros(Int, N); onei[i] = 1; mat * onei == onei end, 1:N)
            for h in 1:1000
                x = reverse_hash_position(h, Val{N}())
                all(i -> iszero(x[i]), dumb_indices) || continue
                attempt_x = true
                for ofsvmap in vmaps
                    ofs = ofsvmap.ofs
                    if all(isinteger, invmat * (x .- ofs))
                        attempt_x = false
                        break
                    end
                end
                attempt_x || continue
                newvmap = OffsetVertexIterator(x, refmap)
                # attempt to make it so that one of the vertices is in the reference cell if possible
                if all(x -> !iszero(x.ofs), newvmap)
                    for h2 in 0:100
                        x2 = reverse_hash_position(h2, Val{N}())
                        altx = -x - mat * x2
                        if any(t -> t.ofs == altx, refmap)
                            newvmap = OffsetVertexIterator(-altx, refmap)
                            break
                        end
                    end
                end
                push!(vmaps, newvmap)
                length(vmaps) == ex && break
            end
        end
        length(vmaps) < ex && error("Could not complete the splitting of connected components")
    end
    encountered[keep]
end
