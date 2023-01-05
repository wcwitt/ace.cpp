#pragma once

#include <complex>
#include <vector>

#include "symmprod_dag.hpp"

struct SparseSymmProd
{
    SparseSymmProdDAG dag;
    std::vector<int> proj;
};

std::vector<std::complex<double>> SparseSymmProd_evaluate(
    std::vector<std::vector<int>> specv,
    std::vector<std::complex<double>> A);
