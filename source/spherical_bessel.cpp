#include <cmath>

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
