VERSION >= v"0.4.0-dev+6521" && __precompile__()

module SmithNormalForm

export snf

include("bezout.jl")
include("snf.jl")

end # module
