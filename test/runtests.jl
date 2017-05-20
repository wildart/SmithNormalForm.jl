using SmithNormalForm
using Base.Test

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

    U, Uinv, V, Vinv, D = snf(M)
    @test !issparse(D)
    @test U*D*V == M
    @testset "dense diag" for i in 1:4
        @test D[i,i] == -1
    end

    U, Uinv, V, Vinv, D = snf(sparse(M))
    @test issparse(D)
    @test U*D*V == M
    @testset "sparse diag" for i in 1:4
        @test D[i,i] == -1
    end
end
