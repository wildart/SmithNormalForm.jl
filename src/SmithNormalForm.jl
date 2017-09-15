VERSION >= v"0.4.0-dev+6521" && __precompile__()

module SmithNormalForm

import Base: show, getindex

export snf, smithfact

include("bezout.jl")
include("snf.jl")

struct Smith{T,S<:AbstractMatrix} <: Factorization{T}
    P::S
    Pinv::S
    Q::S
    Qinv::S
    A::S
    Smith{T,S}(P::AbstractMatrix{T}, Pinv::AbstractMatrix{T},
               Q::AbstractMatrix{T}, Qinv::AbstractMatrix{T},
               A::AbstractMatrix{T}) where {T,S} = new(P, Pinv, Q, Qinv, A)
end
Smith(P::AbstractMatrix{T}, Q::AbstractMatrix{T}, A::AbstractMatrix{T}) where {T} =
    Smith{T,typeof(A)}(P, inv(P), Q, inv(Q), A)

function smithfact(X::AbstractMatrix{T}) where {T}
    P, Pinv, Q, Qinv, A = snf(X)
    return Smith{T, typeof(X)}(P, Pinv, Q, Qinv, A)
end

function getindex(F::Smith{T,S}, d::Symbol) where {T,S}
    if d == :P
        return F.P
    elseif d == :Q
        return F.Q
    elseif d == :SNF || d == :A
        return F.A
    else
        throw(KeyError(d))
    end
end

function show(io::IO, F::Smith)
    println(io, "Smith normal form:")
    show(io, F[:SNF])
end

end # module
