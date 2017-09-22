using SmithNormalForm
using Base.Test

srand(18743874873);

@testset "Bezout" begin
    a = 12
    b = 42

    s,t,g = SmithNormalForm.bezout(a,b)

    @test s*a + t*b == g
end

@testset "Smith Normal Form" begin
    M=[ 1  0  0  0  0  0 ;
        0  0  0  0  0  0 ;
       -1  0  1  0  0  1 ;
        0 -1  0 -1  0  0 ;
        0  0  0  1  1 -1 ;
        0  1 -1  0 -1  0 ]

    P, Q, A, Pinv, Qinv = snf(M)
    @test !issparse(A)
    @test P*A*Q == M
    @test inv(P) == Pinv
    @test inv(Q) == Qinv
    @testset "dense diag" for i in 1:4
        @test A[i,i] == -1
    end

    P, Q, A, Pinv, Qinv = snf(sparse(M))
    @test issparse(A)
    @test P*A*Q == M
    @test inv(collect(P)) == collect(Pinv)
    @test inv(collect(Q)) == collect(Qinv)
    @testset "sparse diag" for i in 1:4
        @test A[i,i] == -1
    end

    P, Q, A, Pinv, Qinv = snf(sparse(M), inverse=false)
    @test issparse(A)
    @test size(Pinv) == (0,0)
    @test size(Qinv) == (0,0)
    @test P*A*Q == M
    @testset "sparse diag" for i in 1:4
        @test A[i,i] == -1
    end

    n = 3
    m = 5
    M = rand(1:100, n, m)
    F = snfact(M)
    @test !issparse(F[:P])
    @test size(F[:P]) == (n,n)
    @test size(F[:Q]) == (m,m)
    @test size(F[:Pinv]) == (n,n)
    @test size(F[:Qinv]) == (m,m)
    @test all(F[:P]*F[:A]*F[:Q] .== M)
    @test all(F[:P]*F[:Pinv] .== eye(eltype(F[:P]), size(F[:P])...))
    @test all(F[:Q]*F[:Qinv] .== eye(eltype(F[:Q]), size(F[:Q])...))

    F = snfact(sparse(M), inverse=false)
    @test issparse(F[:P])
    @test size(F[:P]) == (n,n)
    @test size(F[:Q]) == (m,m)
    @test all(F[:P]*F[:A]*F[:Q] .== M)
    @test_throws AssertionError F[:Pinv]
    @test_throws AssertionError F[:Qinv]
end
