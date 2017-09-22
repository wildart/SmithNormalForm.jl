VERSION >= v"0.4.0-dev+6521" && __precompile__()

module SmithNormalForm

import Base: show, getindex

export snf, snfact

include("bezout.jl")
include("snf.jl")

struct Smith{T,S<:AbstractMatrix} <: Factorization{T}
    P::S
    Pinv::S
    Q::S
    Qinv::S
    A::AbstractVector{T}
    Smith{T,S}(P::AbstractMatrix{T}, Pinv::AbstractMatrix{T},
               Q::AbstractMatrix{T}, Qinv::AbstractMatrix{T},
               A::AbstractVector{T}) where {T,S} = new(P, Pinv, Q, Qinv, A)
end
Smith{T}(P::AbstractMatrix{T}, Q::AbstractMatrix{T}, A::AbstractVector{T}) =
    Smith{T,typeof(A)}(P, similar(P, 0, 0), Q, similar(Q, 0, 0), A)

function snfact{T}(X::AbstractMatrix{T}; inverse=true)
    P, Q, A, Pinv, Qinv = snf(X, inverse=inverse)
    return Smith{T, typeof(X)}(P, Pinv, Q, Qinv, diag(A))
end

function getindex{T,S}(F::Smith{T,S}, d::Symbol)
    if d == :P
        return F.P
    elseif d == :Q
        return F.Q
    elseif d == :SNF
        return F.A
    elseif d == :A
        return if S <: AbstractSparseMatrix
            spdiagm(F.A)*speye(T, size(F.P,1), size(F.Q,1))
        else
            diagm(F.A)*eye(T, size(F.P,1), size(F.Q,1))
        end
    elseif d == :Pinv
        @assert size(F.Pinv,1) != 0 "Inverse are not available, set `inverse` parameter when factorizing."
        return F.Pinv
    elseif d == :Qinv
        @assert size(F.Pinv,1) != 0 "Inverse are not available, set `inverse` parameter when factorizing."
        return F.Qinv
    else
        throw(KeyError(d))
    end
end

function show(io::IO, F::Smith)
    println(io, "Smith normal form:")
    show(io, F[:SNF])
end

end # module
