using SmithNormalForm
using BenchmarkTools

srand(2903872398473)
rows, cols = 7, 13
TMP = rand(1:10, rows*cols)
D = reshape(TMP, rows, cols)
i, j = 2, 3
a, b, c, d = 1, 2, 3, 4

suite = BenchmarkGroup()

suite["internal"] = BenchmarkGroup(["swap", "elimination"])
suite["snf"] = BenchmarkGroup(["decomposition"])

suite["internal"]["col_swap"] = @benchmarkable SmithNormalForm.cswap!(D, i, j) samples=100000
suite["internal"]["row_swap"] = @benchmarkable SmithNormalForm.rswap!(D, i, j) samples=100000

suite["internal"]["col_elim"] = @benchmarkable SmithNormalForm.colelimination(D, a, b, c, d, i, j)
suite["internal"]["row_elim"] = @benchmarkable SmithNormalForm.rowelimination(D, a, b, c, d, i, j)

Z = rand(1:100, 100, 120)
suite["snf"]["dense"] = @benchmarkable snf(Z)

Z = sprand(Int, 100, 120, 0.3)
Z[Z .< 0] = -1
Z[Z .> 0] = 1
suite["snf"]["sparse"] = @benchmarkable snf(Z)
