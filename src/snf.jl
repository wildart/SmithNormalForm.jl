# Smith Normal Form

function divisable{R}(y::R, x::R )
  x == zero(R) && return y == zero(R)
  return div(y,x)*x == y
end

function divide{R}(y::R, x::R)
    if x != -one(R)
        return div(y,x)
    else
        return y * x
    end
end

function rswap{R}(M::AbstractArray{R,2}, r1, r2)
    r1 == r2 && return
    for i in eachindex(M[r1,:])
        tmp = M[r1, i]
        M[r1, i] = M[r2, i]
        M[r2, i] = tmp
    end
    return
end

function cswap{R}(M::AbstractArray{R,2}, c1, c2)
    c1 == c2 && return
    for j in eachindex(M[:,c1])
        tmp = M[j, c1]
        M[j, c1] = M[j, c2]
        M[j, c2] = tmp
    end
    return
end

function rowelimination{R}(D::AbstractArray{R,2}, a::R, b::R, c::R, d::R, i::Int, j::Int)
    I = view(D, i, :)
    @inbounds for k in eachindex(I)
        t = D[i,k]
        s = D[j,k]
        D[i,k] = a*t + b*s
        D[j,k] = c*t + d*s
    end
    return
end

function colelimination{R}(D::AbstractArray{R,2}, a::R, b::R, c::R, d::R, i::Int, j::Int)
    I = view(D, :, i)
    @inbounds for k in eachindex(I)
        t = D[k,i]
        s = D[k,j]
        D[k,i] = a*t + b*s
        D[k,j] = c*t + d*s
    end
    return
end

function rowpivot{R}(U::AbstractArray{R,2},
                     Uinv::AbstractArray{R,2},
                     D::AbstractArray{R,2},
                     i, j)
    for k in findn(D[:,j]) |> reverse
        a = D[i,j]
        b = D[k,j]

        i == k && continue

        s,t,g = bezout(a, b)
        x = divide(a, g)
        y = divide(b, g)

        rowelimination(D, s, t, -y, x, i, k)
        rowelimination(Uinv, s, t, -y, x, i, k)
        colelimination(U, x, y, -t, s, i, k)
    end
end

function colpivot{R}(V::AbstractArray{R,2},
                     Vinv::AbstractArray{R,2},
                     D::AbstractArray{R,2},
                     i, j)
    for k in findn(D[i,:])|> reverse
        a = D[i,j]
        b = D[i,k]

        j == k && continue

        s, t, g = bezout(a, b)
        x = divide(a, g)
        y = divide(b, g)

        colelimination(D, s, t, -y, x, j, k)
        colelimination(Vinv, s, t, -y, x, j, k)
        rowelimination(V, x, y, -t, s, j, k)
    end
end

function smithpivot{R}(U::AbstractArray{R,2},
                       Uinv::AbstractArray{R,2},
                       V::AbstractArray{R,2},
                       Vinv::AbstractArray{R,2},
                       D::AbstractArray{R,2},
                       i, j)

    pivot = D[i,j]
    @assert pivot != zero(R) "Pivot cannot be zero"

    while countnz(D[i,:]) > 1 || countnz(D[:,j]) > 1
        colpivot(V, Vinv, D, i, j)
        rowpivot(U, Uinv, D, i, j)
    end
end

function snf{R}(M::AbstractArray{R,2})
    D = copy(M)
    rows, cols = size(M)

    U = issparse(M) ? spzeros(R, rows, rows) : zeros(R, rows, rows)
    for i in 1:rows
        U[i,i] = one(R)
    end
    Uinv = copy(U)

    V = issparse(M) ? spzeros(R, cols, cols) : zeros(R, cols, cols)
    for i in 1:cols
        V[i,i] = one(R)
    end
    Vinv = copy(V)

    t = 1
    for j in 1:cols
        # println("Working on column $j out of ", size(D,2))
        # display(D)

        countnz(D[:,j]) == 0 && continue

        prow = 1
        if D[t,t] != zero(R)
            prow = t
        else
            # Good pivot row for j-th column is the one
            # that have a smallest number of elements
            idxs = findn(D[:,j])
            rsize, i = mapslices(countnz, D[idxs, :], 2) |> findmin
            prow = idxs[i]
        end

        # "Pivot Row selected: t = $t, pivot_row = $prow" |> println
        # display(D)
        # "Swapping rows. ( $t, $prow )" |> println

        rswap(D, t, prow)
        rswap(Uinv, t, prow)
        cswap(U, t, prow)

        # println("Performing the pivot step at ($t, $j)")
        # display(collect(D))

        smithpivot(U, Uinv, V, Vinv, D, t, j)

        cswap(D, t, j)
        cswap(Vinv, t, j)
        rswap(V, t, j)

        t += 1

        # println("D: $(size(D))")
        # display(collect(D))
        # println("U: $(size(U))")
        # display(collect(U))
        # println("V: $(size(V))")
        # display(collect(V))
        # println("Vinv:")
        # display(collect(Vinv))
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

            smithpivot(U, Uinv, V, Vinv, D, i, i)
        end
    end

    return U, Uinv, V, Vinv, D
end
