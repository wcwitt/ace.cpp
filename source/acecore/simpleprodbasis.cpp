#include "simpleprodbasis.hpp"

std::vector<std::complex<double>> SimpleProdBasis_evaluate(
    std::vector<std::vector<int>> specv,
    std::vector<std::complex<double>> A)
{
    std::vector<std::complex<double>> AA(specv.size());
    for (int i=0; i<specv.size(); ++i) {
        AA[i] = A[specv[i][0]-1];
        for (int j=1; j<specv[i].size(); ++j) {
            AA[i] *= A[specv[i][j]-1];
        }
    }
    return AA;
}

//std::vector<std::complex<double>> SimpleProdBasis_evaluate(
//    std::vector<std::vector<int>> specm,
//    std::vector<int> orders,
//    std::vector<double> A)
//{
//    std::vector<std::complex<double>> AA(specm.size());
//    for (int i=0; i<specm.size(); ++i) {
//        AA[i] = A[specm[i][0]];
//        for (int j=1; j<=orders[i]; ++j) {
//            AA[i] *= A[specm[i][j]];
//        }
//    }
//    return AA;
//}
