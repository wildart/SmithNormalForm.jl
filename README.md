# Smith Normal Form

[![Build Status](https://travis-ci.org/wildart/SmithNormalForm.jl.svg?branch=master)](https://travis-ci.org/wildart/SmithNormalForm.jl)
[![Coverage Status](https://coveralls.io/repos/wildart/SmithNormalForm.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/wildart/SmithNormalForm.jl?branch=master)


The [Smith normal form](https://en.wikipedia.org/wiki/Smith_normal_form) decomposition implementation in Julia. *Note:* Works only over integer domain.

## Installation
```
pkg> add https://github.com/wildart/SmithNormalForm.jl.git#0.2.0
```

## Example

```julia
julia> using SmithNormalForm, LinearAlgebra

julia> X = rand(1:100, 3, 5)
3×5 Array{Int64,2}:
 14  82  85  39  56
 70  51   4  68  23
 58  27  87   1  52

julia> F = smith(X)
Smith normal form:
[1, -1, 1]

julia> F.S
3×3 Array{Int64,2}:
   1      0  0
 640      1  0
 807  29427  1

julia> F.T
5×5 Array{Int64,2}:
        14          82          85         39          56
      8890       52429       54396      24892       35817
 261594790  1542762036  1600642584  732465412  1053941719
  87198265   514254012   533547528  244155142   351313913
         0  -109493940  -113601877          0     -536735

julia> diagm(F)
3×5 Array{Int64,2}:
 1   0  0  0  0
 0  -1  0  0  0
 0   0  1  0  0

julia> F.SNF
3-element Array{Int64,1}:
  1
 -1
  1

julia> F.S*diagm(F)*F.T
3×5 Array{Int64,2}:
 14  82  85  39  56
 70  51   4  68  23
 58  27  87   1  52
```
