module SmithNormalForm

# import Base: show, getindex
using LinearAlgebra
using SparseArrays
using Base.CoreLogging

import Base: show
import LinearAlgebra: diagm

export snf, smith, Smith

include("bezout.jl")
include("snf.jl")


struct Smith{P,Q<:AbstractMatrix{P},V<:AbstractVector{P}} <: Factorization{P}
    S::Q
    Sinv::Q
    T::Q
    Tinv::Q
    SNF::V
    Smith{P,Q,V}(S::AbstractMatrix{P}, Sinv::AbstractMatrix{P},
                 T::AbstractMatrix{P}, Tinv::AbstractMatrix{P},
                 D::AbstractVector{P}) where {P,Q,V} = new(S, Sinv, T, Tinv, D)
end
Smith(S::AbstractMatrix{P}, T::AbstractMatrix{P}, SNF::AbstractVector{P}) where {P} =
    Smith{P,typeof(S),typeof(SNF)}(S, similar(S, 0, 0), T, similar(T, 0, 0), SNF)

"""
    smith(X::AbstractMatrix{P}; inverse::Bool=true) --> Smith{P,Q,V}

Return a Smith normal form of an integer matrix `X` as a `Smith` structure (of element type
`P`, matrix type `Q`, and invariant factor type `V`).

The Smith normal form is well-defined for any matrix ``m×n`` matrix `X` with elements in a
principal domain (PID; e.g., integers) and provides a decomposition of `X` into ``m×m``,
``m×n``, and `S`, `Λ`, and `T` as `X = SΛT`, where `Λ` is a diagonal matrix with entries 
("invariant factors") Λᵢ ≥ Λᵢ₊₁ ≥ 0 with nonzero entries divisible in the sense Λᵢ | Λᵢ₊₁.
The invariant factors can be obtained from [`diag(::Smith)`](@ref).

`S` and `T` are invertible matrices; if the keyword argument `inverse` is true (default),
the inverse matrices are computed and returned as part of the `Smith` factorization.
"""
function smith(X::AbstractMatrix{P}; inverse::Bool=true) where {P}
    S, T, D, Sinv, Tinv = snf(X, inverse=inverse)
    SNF = diag(D)
    return Smith{P, typeof(X), typeof(SNF)}(S, Sinv, T, Tinv, SNF)
end

"""Retrive SNF matrix from the factorization"""
function diagm(F::Smith{P,Q}) where {P,Q}
    base = issparse(F.SNF) ? spzeros(P, size(F.S,1), size(F.T,1)) : zeros(P, size(F.S,1), size(F.T,1))
    for i in 1:length(F.SNF)
        base[i,i] = F.SNF[i]
    end
    return base
end

function show(io::IO, F::Smith)
    println(io, "Smith normal form:")
    show(io, diagm(F))
end

end # module
