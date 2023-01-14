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

// ----- new for testing

void evaluate_real_stl(
    std::vector<double> AA,
    const std::vector<int> nodes,
    const std::vector<double> A);

void evaluate_real(
    double AA[],
    const int nodes[],
    const double A[],
    const int AA_size,
    const int A_size);

void evaluate_complex(
    std::complex<double> AA[],
    const int nodes[],
    std::complex<double> A[],
    const int AA_size,
    const int A_size);

void evaluate_batch_real(
    double AA[],
    const int nodes[],
    const double A[],
    const int AA_size,
    const int A_size,
    const int batch_size);

void evaluate_batch_complex(
    std::complex<double> AA[],
    const int nodes[],
    const std::complex<double> A[],
    const int AA_size,
    const int A_size,
    const int batch_size);
