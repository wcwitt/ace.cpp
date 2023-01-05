
using ACEcore
using ACEcore: SimpleProdBasis, SparseSymmProd
using CxxWrap
using Test

@info("Loading C++ module")
include("../source/cxxwrap/libace.jl")

##

M = 20 
spec = ACEcore.Testing.generate_SO2_spec(5, M)  # 5, M
A = randn(ComplexF64, 2*M+1)

## 

@info("Test consistency of SparseSymmetricProduct with SimpleProdBasis")
basis1 = SimpleProdBasis(spec)
AA1 = basis1(A)

basis2 = SparseSymmProd(spec; T = ComplexF64)
AA2 = basis2(A)

println(@test AA1 ≈ AA2)

## 

@info("Test consistency of SimpleProdBasis with C++ version")
AA1_cpp = LibACE.SimpleProdBasis_evaluate(
    StdVector([StdVector(Int32.(s)) for s in spec]),
    StdVector(A))
println(@test AA1 ≈ AA1_cpp)

@info("Test consistency of SparseSymmProd with C++ version")
AA2_cpp = LibACE.SparseSymmProd_evaluate(
    StdVector([StdVector(Int32.(s.-1)) for s in spec]),
    StdVector(A))
println(@test AA2 ≈ AA2_cpp)
