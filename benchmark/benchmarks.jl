using SmithNormalForm
using BenchmarkTools

srand(2903872398473)

@benchgroup "SNF" ["matrix", "decomposition"] begin
    rows, cols = 7, 13
    TMP = rand(1:10, rows*cols)
    D = reshape(TMP, rows, cols)
    i, j = 2, 3
    a, b, c, d = 1, 2, 3, 4

    @bench "rowelementation" SmithNormalForm.rowelementation(D, a, b, c, d, i, j)
    @bench "colelementation" SmithNormalForm.colelementation(D, a, b, c, d, i, j)

    Z = sprand(Int, 100, 120, 0.3)
    Z[Z .< 0] = -1
    Z[Z .> 0] = 1
    @bench "snf" snf(Z)
end
