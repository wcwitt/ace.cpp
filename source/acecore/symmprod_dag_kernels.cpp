#include "symmprod_dag_kernels.hpp"

void evaluate(
    std::vector<std::complex<double>> &AA,
    const SparseSymmProdDAG &dag,
    std::vector<std::complex<double>> A)
{
    auto nodes = dag.nodes;
    // @assert length(AA) >= dag.numstore
    // @assert length(A) >= dag.num1

    // Stage-1: copy the 1-particle basis into AA
    for (int i=0; i<dag.num1; ++i) {
        AA[i] = A[i];
    }
 
    // Stage-2: go through the dag and store the intermediate results we need
    for (int i=dag.num1; i<nodes.size(); ++i) {
        auto [n1, n2] = nodes[i];
        AA[i] = AA[n1] * AA[n2];
    }
 
    return;
}
