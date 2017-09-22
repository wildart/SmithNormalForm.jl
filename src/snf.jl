# Smith Normal Form

function divisable(y::R, x::R ) where {R}
  x == zero(R) && return y == zero(R)
  return div(y,x)*x == y
end

function divide(y::R, x::R) where {R}
    if x != -one(R)
        return div(y,x)
    else
        return y * x
    end
end

function rcountnz(X, i)
    n = size(X, 1)
    c = 0
    for j in 1:n
        if X[j, i] != 0
           c += 1
        end
    end
    return c
end

function ccountnz(X, j)
    m = size(X, 2)
    c = 0
    for i in 1:m
        if X[j, i] != 0
           c += 1
        end
    end
    return c
end

function rswap!(M::AbstractArray{R,2}, r1::Int, r2::Int) where {R}
    r1 == r2 && return M
    m = size(M, 2)
    @inbounds for i in 1:m
        tmp = M[r1, i]
        M[r1, i] = M[r2, i]
        M[r2, i] = tmp
    end
    return M
end

function cswap!(M::AbstractArray{R,2}, c1::Int, c2::Int) where {R}
    c1 == c2 && return M
    n = size(M, 1)
    @inbounds for j in 1:n
        tmp = M[j, c1]
        M[j, c1] = M[j, c2]
        M[j, c2] = tmp
    end
    return M
end

function rowelimination(D::AbstractArray{R,2}, a::R, b::R, c::R, d::R, i::Int, j::Int) where {R}
    m = size(D, 2)
    @inbounds for k in 1:m
        t = D[i,k]
        s = D[j,k]
        D[i,k] = a*t + b*s
        D[j,k] = c*t + d*s
    end
    return D
end

function colelimination(D::AbstractArray{R,2}, a::R, b::R, c::R, d::R, i::Int, j::Int) where {R}
    n = size(D, 1)
    @inbounds for k in 1:n
        t = D[k,i]
        s = D[k,j]
        D[k,i] = a*t + b*s
        D[k,j] = c*t + d*s
    end
    return D
end

function rowpivot(U::AbstractArray{R,2},
                  Uinv::AbstractArray{R,2},
                  D::AbstractArray{R,2},
                  i, j; inverse=true) where {R}
    for k in findn(D[:,j]) |> reverse
        a = D[i,j]
        b = D[k,j]

        i == k && continue

        s,t,g = bezout(a, b)
        x = divide(a, g)
        y = divide(b, g)

        rowelimination(D, s, t, -y, x, i, k)
        inverse && rowelimination(Uinv, s, t, -y, x, i, k)
        colelimination(U, x, y, -t, s, i, k)
    end
end

function colpivot(V::AbstractArray{R,2},
                  Vinv::AbstractArray{R,2},
                  D::AbstractArray{R,2},
                  i, j; inverse=true) where {R}
    for k in findn(D[i,:])|> reverse
        a = D[i,j]
        b = D[i,k]

        j == k && continue

        s, t, g = bezout(a, b)
        x = divide(a, g)
        y = divide(b, g)

        colelimination(D, s, t, -y, x, j, k)
        inverse && colelimination(Vinv, s, t, -y, x, j, k)
        rowelimination(V, x, y, -t, s, j, k)
    end
end

function smithpivot(U::AbstractArray{R,2},
                    Uinv::AbstractArray{R,2},
                    V::AbstractArray{R,2},
                    Vinv::AbstractArray{R,2},
                    D::AbstractArray{R,2},
                    i, j; inverse=true) where {R}

    pivot = D[i,j]
    @assert pivot != zero(R) "Pivot cannot be zero"
    while ccountnz(D,i) > 1 || rcountnz(D,j) > 1
        colpivot(V, Vinv, D, i, j, inverse=inverse)
        rowpivot(U, Uinv, D, i, j, inverse=inverse)
    end
end

function init(M::AbstractSparseMatrix{R,Ti}; inverse=true) where {R, Ti}
    D = copy(M)
    rows, cols = size(M)

    U = spzeros(R, rows, rows)
    for i in 1:rows
        U[i,i] = one(R)
    end
    Uinv = inverse ? copy(U) : spzeros(R, 0, 0)

    V = spzeros(R, cols, cols)
    for i in 1:cols
        V[i,i] = one(R)
    end
    Vinv = inverse ? copy(V) : spzeros(R, 0, 0)

    return U, V, D, Uinv, Vinv
end

function init(M::AbstractMatrix{R}; inverse=true) where {R}
    D = copy(M)
    rows, cols = size(M)

    U = zeros(R, rows, rows)
    for i in 1:rows
        U[i,i] = one(R)
    end
    Uinv = inverse ? copy(U) : zeros(R, 0, 0)

    V = zeros(R, cols, cols)
    for i in 1:cols
        V[i,i] = one(R)
    end
    Vinv = inverse ? copy(V) : zeros(R, 0, 0)

    return U, V, D, Uinv, Vinv
end

function snf(M::AbstractMatrix{R}; inverse=true, debug=false) where {R}
    rows, cols = size(M)
    U, V, D, Uinv, Vinv = init(M, inverse=inverse)

    t = 1
    for j in 1:cols
        debug && println("Working on column $j out of ", size(D,2))
        debug && display(D)

        rcountnz(D,j) == 0 && continue

        prow = 1
        if D[t,t] != zero(R)
            prow = t
        else
            # Good pivot row for j-th column is the one
            # that have a smallest number of elements
            idxs = findn(D[:,j])
            rsize, i = mapslices(x->count(!iszero, x), D[idxs, :], 2) |> findmin
            prow = idxs[i]
        end

        debug && "Pivot Row selected: t = $t, pivot_row = $prow" |> println
        debug && display(D)
        debug && "Swapping rows. ( $t, $prow )" |> println

        rswap!(D, t, prow)
        inverse && rswap!(Uinv, t, prow)
        cswap!(U, t, prow)

        debug && println("Performing the pivot step at ($t, $j)")
        debug && display(collect(D))

        smithpivot(U, Uinv, V, Vinv, D, t, j, inverse=inverse)

        cswap!(D, t, j)
        inverse && cswap!(Vinv, t, j)
        rswap!(V, t, j)

        t += 1

        debug && println("D: $(size(D))")
        debug && display(collect(D))
        debug && println("U: $(size(U))")
        debug && display(collect(U))
        debug && println("V: $(size(V))")
        debug && display(collect(V))
        debug && println("Vinv:")
        debug && display(collect(Vinv))
    end

    # Make sure that d_i is divisible be d_{i+1}.
    r = minimum(size(D))
    pass = true
    while pass
        pass = false
        for i in 1:r-1
            divisable(D[i+1,i+1], D[i,i]) && continue
            pass = true
            D[i+1,i] = D[i+1,i+1]

            colelimination(Vinv, one(R), one(R), zero(R), one(R), i, i+1)
            rowelimination(V, one(R), zero(R), -one(R), one(R), i, i+1)

            smithpivot(U, Uinv, V, Vinv, D, i, i, inverse=inverse)
        end
    end

    return U, V, D, Uinv, Vinv
end
