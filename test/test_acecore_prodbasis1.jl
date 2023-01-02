
using ACEcore, BenchmarkTools
using ACEcore: SimpleProdBasis
using CxxWrap
using Test

@info("Loading C++ module")
include("../source/cxxwrap/libace.jl")

##

M = 20 
spec = ACEcore.Testing.generate_SO2_spec(5, M)
A = randn(ComplexF64, 2*M+1)

## 

@info("Test consistency of SimpleProdBasis with C++ version")
basis1 = SimpleProdBasis(spec)
AA1 = basis1(A)
AA1_cpp = LibACE.SimpleProdBasis_evaluate(
    StdVector([StdVector(Int32.(s)) for s in spec]),
    StdVector(A))
println(@test AA1 ≈ AA1_cpp)
