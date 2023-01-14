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

void evaluate_real_stl(
    std::vector<double> AA,
    const std::vector<int> nodes,
    const std::vector<double> A)
{
    // Stage-1: copy the 1-particle basis into AA
    for (int i=0; i<A.size(); ++i) {
        AA[i] = A[i];
    }

    // Stage-2: go through the dag and store the intermediate results we need
    for (int i=A.size(); i<AA.size(); ++i) {
        AA[i] = AA[nodes[2*i]] * AA[nodes[2*i+1]];
    }

    return;
}

void evaluate_real(
    double AA[],
    const int nodes[],
    const double A[],
    const int AA_size,
    const int A_size)
{
    // Stage-1: copy the 1-particle basis into AA
    for (int i=0; i<A_size; ++i) {
        AA[i] = A[i];
    }

    // Stage-2: go through the dag and store the intermediate results we need
    for (int i=A_size; i<AA_size; ++i) {
        AA[i] = AA[nodes[2*i]] * AA[nodes[2*i+1]];
    }

    return;
}

void evaluate_complex(
    std::complex<double> AA[],
    const int nodes[],
    std::complex<double> A[],
    const int AA_size,
    const int A_size)
{
    // Stage-1: copy the 1-particle basis into AA
    for (int i=0; i<A_size; ++i) {
        AA[i] = A[i];
    }

    // Stage-2: go through the dag and store the intermediate results we need
    for (int i=A_size; i<AA_size; ++i) {
        AA[i] = AA[nodes[2*i]] * AA[nodes[2*i+1]];
    }

    return;
}

void evaluate_batch_real(
    double AA[],
    const int nodes[],
    const double A[],
    const int AA_size,
    const int A_size,
    const int batch_size)
{
    // Stage-1: copy the 1-particle basis into AA
    for (int i=0; i<A_size; ++i) {
        double* AA_i = &AA[i*batch_size];
        const double* A_i = &A[i*batch_size];
        #pragma omp simd
        for (int j=0; j<batch_size; ++j) {
            AA_i[j] = A_i[j];
        }
    }

    // Stage-2: go through the dag and store the intermediate results we need
    for (int i=A_size; i<AA_size; ++i) {
        const int n0 = nodes[2*i];
        const int n1 = nodes[2*i+1];
        double* AA_i = &AA[i*batch_size];
        const double* AA_0 = &AA[n0*batch_size];
        const double* AA_1 = &AA[n1*batch_size];
        #pragma omp simd
        for (int j=0; j<batch_size; ++j) {
            AA_i[j] = AA_0[j] * AA_1[j];
        }
    }
    return;
}

void evaluate_batch_complex(
    std::complex<double> AA[],
    const int nodes[],
    const std::complex<double> A[],
    const int AA_size,
    const int A_size,
    const int batch_size)
{
    // Stage-1: copy the 1-particle basis into AA
    for (int i=0; i<A_size; ++i) {
        std::complex<double>* AA_i = &AA[i*batch_size];
        const std::complex<double>* A_i = &A[i*batch_size];
        #pragma omp simd
        for (int j=0; j<batch_size; ++j) {
            AA_i[j] = A_i[j];
        }
    }

    // Stage-2: go through the dag and store the intermediate results we need
    for (int i=A_size; i<AA_size; ++i) {
        const int n0 = nodes[2*i];
        const int n1 = nodes[2*i+1];
        std::complex<double>* AA_i = &AA[i*batch_size];
        const std::complex<double>* AA_0 = &AA[n0*batch_size];
        const std::complex<double>* AA_1 = &AA[n1*batch_size];
        #pragma omp simd
        for (int j=0; j<batch_size; ++j) {
            AA_i[j] = AA_0[j] * AA_1[j];
        }
    }
    return;
}
