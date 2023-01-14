# To run this we need ACEcore.jl needs to be the version in *this* repo, which 
# can be added via `add ../..` or `dev ../..` or something like it. 
#

##

using ACE1, ACEcore, BenchmarkTools
using ACE1: PolyTransform, transformed_jacobi

using CxxWrap
using Test
@info("Loading C++ module")
include("../source/cxxwrap/libace.jl")

##
# create a real potential 
# tune the degree to get different size graphs, I like the following: 
#
# maxdeg  | basis size 
#   12    |  ≈ 7k
#   15    |  ≈ 33k
#   18    |  ≈ 135k 
#   21    |  ≈ 466k 

@info("Basic test of PIPotential construction and evaluation")
maxdeg = 12
order = 3
r0 = 1.0
rcut = 3.0
trans = PolyTransform(1, r0)
Pr = transformed_jacobi(maxdeg, trans, rcut; pcut = 2)
D = ACE1.SparsePSHDegree()
P1 = ACE1.BasicPSH1pBasis(Pr; species = :X, D = D)
basis = ACE1.PIBasis(P1, order, D, maxdeg)

@show length(basis)

## 
# convert to the new format 

orders = basis.inner[1].orders
iAA2iA = basis.inner[1].iAA2iA
new_spec = [ iAA2iA[i, 1:orders[i]][:] for i = 1:length(orders) ]
bAA = ACEcore.SparseSymmProdDAG(new_spec, T = ComplexF64)
len_A = maximum(iAA2iA)
len_AA = length(bAA)

## 
# non-batched benchmark 
A1 = randn(ComplexF64, len_A)
AA1 = zeros(ComplexF64, len_AA)
rAA1 = real.(AA1)
rA1 = real.(A1)

println("non-batched benchmark:")
print(" complex: "); @btime ACEcore.evaluate!($AA1, $bAA, $A1) samples=1 evals=10
print("    real: "); @btime ACEcore.evaluate!($rAA1, $bAA, $rA1) #samples=1 evals=10
println() 

# repeat with cpp
AA1_cpp = zeros(ComplexF64, len_AA)
rAA1_cpp = real.(AA1)
nodes_cpp = zeros(Int32, 2*length(bAA.nodes))
for i in 1:length(bAA.nodes)
    nodes_cpp[2*i-1] = bAA.nodes[i][1]-1
    nodes_cpp[2*i] = bAA.nodes[i][2]-1
end
print(" complex (cpp): "); @btime LibACE.evaluate_complex(AA1_cpp, nodes_cpp, A1, len_AA, len_A)
print("    real (cpp): "); @btime LibACE.evaluate_real(rAA1_cpp, nodes_cpp, rA1, len_AA, len_A)
println(@test AA1 ≈ AA1_cpp)
println(@test rAA1 ≈ rAA1_cpp)

###

Nbatch = 32 
A = randn(ComplexF64, Nbatch, len_A)
AA = zeros(ComplexF64, Nbatch, len_AA)
rAA = real.(AA)
rA = real.(A)

println("batched benchmark with batch-size $Nbatch:")
print(" complex: "); @btime ACEcore.evaluate!($AA, $bAA, $A)
print("    real: "); @btime ACEcore.evaluate!($rAA, $bAA, $rA)
println() 

# repeat with cpp
AA_cpp = zeros(ComplexF64, Nbatch*len_AA)
rAA_cpp = zeros(Float64, Nbatch*len_AA)
A_cpp = Array(vec(A))
rA_cpp = Array(vec(rA))
print(" complex (cpp): "); @btime LibACE.evaluate_batch_complex(AA_cpp, nodes_cpp, A_cpp, len_AA, len_A, Nbatch)
print("    real (cpp): "); @btime LibACE.evaluate_batch_real(rAA_cpp, nodes_cpp, rA_cpp, len_AA, len_A, Nbatch)
AA_cpp = reshape(AA_cpp,(Nbatch,len_AA))
rAA_cpp = reshape(rAA_cpp,(Nbatch,len_AA))
println(@test AA ≈ AA_cpp)
println(@test rAA ≈ rAA_cpp)

#println("repeated single-input benchmark")
#println("Of course this is not entirely fair since the same memory will reused.")
#print(" complex: "); @btime (for _=1:$Nbatch; ACEcore.evaluate!($AA1, $bAA, $A1); end)
#print("    real: "); @btime (for _=1:$Nbatch; ACEcore.evaluate!($rAA1, $bAA, $rA1); end)

###
#
#using LinearAlgebra: mul!
#println() 
#println("For comparison: real matrix-vector multiplication AA * c")
#vals = zeros(size(AA, 1))
#c = zeros(len_AA)
#@btime mul!($vals, $rAA, $c)
#println("This shows that we MUST incorporate the contraction into the AA evaluation.")
