# Smith Normal Form

[![Build Status](https://travis-ci.org/wildart/SmithNormalForm.jl.svg?branch=master)](https://travis-ci.org/wildart/SmithNormalForm.jl)
[![Coverage Status](https://coveralls.io/repos/wildart/SmithNormalForm.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/wildart/SmithNormalForm.jl?branch=master)


The [Smith normal form](https://en.wikipedia.org/wiki/Smith_normal_form) decomposition implementation in Julia.

## Installation
```julia
julia> Pkg.clone("https://github.com/wildart/SmithNormalForm.jl.git")
```

## Example

```julia
julia> using SmithNormalForm

julia> X = rand(1:100, 3, 5)
3×5 Array{Int64,2}:
 14  82  85  39  56
 70  51   4  68  23
 58  27  87   1  52

julia> F = snfact(X)
Smith normal form:
[1, -1, 1]

julia> F[:P]
3×3 Array{Int64,2}:
   1      0  0
 640      1  0
 807  29427  1

julia> F[:Q]
5×5 Array{Int64,2}:
        14          82          85         39          56
      8890       52429       54396      24892       35817
 261594790  1542762036  1600642584  732465412  1053941719
  87198265   514254012   533547528  244155142   351313913
         0  -109493940  -113601877          0     -536735

julia> F[:A]
3×5 Array{Int64,2}:
 1   0  0  0  0
 0  -1  0  0  0
 0   0  1  0  0

julia> F[:SNF]
3-element Array{Int64,1}:
  1
 -1
  1

julia> F[:P]*F[:A]*F[:Q]
3×5 Array{Int64,2}:
 14  82  85  39  56
 70  51   4  68  23
 58  27  87   1  52
```

*Note:* Currently works only for integer domain.
