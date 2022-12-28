#include <cmath>
#include <vector>

#include "spherical_bessel.hpp"
#include "spherical_bessel_zeros.hpp"

double spherical_bessel_radial(int n, int l, double r, double r_c) {
    double z_nl = spherical_bessel_zeros[l][n-1];
    double c_nlm;
    if (l==0) {
        c_nlm = std::pow(z_nl/r_c,1.5)*std::pow(0.5*z_nl-0.25*std::sin(2*z_nl),-0.5);
    } else {
        c_nlm = std::pow(-0.5*r_c*r_c*r_c*std::sph_bessel(l-1,z_nl)*std::sph_bessel(l+1,z_nl),-0.5);
    }
    return c_nlm * std::sph_bessel(l, z_nl*r/r_c);
}

std::tuple<std::vector<int>,
           std::vector<int>,
           std::vector<int>,
           std::vector<double>>
    determine_basis(double r_max, double e_max) {
    // todo: manage these hardcoded values in better way
    int num_l = 10;
    int num_n = 50;
    // end todo
    double k_max = std::sqrt(2*e_max);
    std::vector<int> A_n;
    std::vector<int> A_l;
    std::vector<int> A_m;
    std::vector<double> A_k;
    for (int l=0; l<num_l; ++l) {
        for (int n=1; n<=num_n; ++n) {
            double z_nl = spherical_bessel_zeros[l][n-1];
            if (z_nl/r_max < k_max) {
                for (int m=-l; m<=l; ++m) {
                    A_n.push_back(n);
                    A_l.push_back(l);
                    A_m.push_back(m);
                    A_k.push_back(z_nl/r_max);
                }
            }
        }
    }
    return {A_n, A_l, A_m, A_k};
}
