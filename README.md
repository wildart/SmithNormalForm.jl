# Smith Normal Form

[![Build Status](https://travis-ci.org/wildart/SmithNormalForm.jl.svg?branch=master)](https://travis-ci.org/wildart/SmithNormalForm.jl)
[![Coverage Status](https://coveralls.io/repos/wildart/SmithNormalForm.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/wildart/SmithNormalForm.jl?branch=master)
[![codecov.io](http://codecov.io/github/wildart/SmithNormalForm.jl/coverage.svg?branch=master)](http://codecov.io/github/wildart/SmithNormalForm.jl?branch=master)

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
 28  17  72  30   4
 95  39  78  58  86
 65  64  59  97   3

julia> F = smithfact(X)
Smith normal form:
[1 0 … 0 0; 0 1 … 0 0; 0 0 … 0 0]

julia> F[:P]
3×3 Array{Int64,2}:
    1      0  0
 4391      1  0
 -544  38656  1

julia> F[:Q]
5×5 Array{Int64,2}:
          28           17            72           30           4
     -122853       -74608       -316074      -131672      -17478
 -4749020865  -2884056160  -12218195771  -5089929249  -675631747
          -7            0             0           -7          -1
   758381427    431684520    1828815281    809408475   107921130

julia> F[:P]*F[:SNF]*F[:Q]
3×5 Array{Int64,2}:
 28  17  72  30   4
 95  39  78  58  86
 65  64  59  97   3
```

*Note:* Currently works only for integer domain.
