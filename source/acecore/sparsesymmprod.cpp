#include "symmprod_dag.hpp"
#include "symmprod_dag_kernels.hpp"
#include "sparsesymmprod.hpp"

std::vector<std::complex<double>> SparseSymmProd_evaluate(
    std::vector<std::vector<int>> spec,
    std::vector<std::complex<double>> A)
{
    auto dag = BuildSparseSymmProdDAG(spec);

    std::vector<std::complex<double>> AAdag(dag.nodes.size());
    evaluate(AAdag, dag, A);

    std::vector<std::complex<double>> AA(spec.size());
    // serial projection
    for (int i=0; i<AA.size(); ++i)
        AA[i] = AAdag[dag.projection[i]];

    return AA;
}
