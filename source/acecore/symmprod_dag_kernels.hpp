#pragma once

#include <complex>
#include <vector>

#include "symmprod_dag.hpp"

// --------------- interface functions


// --------------- old evaluation

void evaluate(
    std::vector<std::complex<double>> &AA,
    const SparseSymmProdDAG &dag,
    std::vector<std::complex<double>> A);
