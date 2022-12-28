#include <cmath>

#include "spherical_harmonic.hpp"

std::complex<double> naive_sph_harm(int l, int m, double polar, double azimuth) {
    if (m >=0) {
        return std::pow(-1,m)
            * std::sqrt((2*l+1)/(4*M_PI)*std::tgamma(1+l-m)/std::tgamma(1+l+m))
            * std::assoc_legendre(l,m,std::cos(polar))
            * std::exp(std::complex<double>(0,1)*double(m)*azimuth);
    } else {
        return std::pow(-1,m)
            * std::sqrt((2*l+1)/(4*M_PI)*std::tgamma(1+l-m)/std::tgamma(1+l+m))
            * std::pow(-1,m)*std::tgamma(1+l+m)/std::tgamma(1+l-m)*std::assoc_legendre(l,-m,std::cos(polar))
            * std::exp(std::complex<double>(0,1)*double(m)*azimuth);
    }
}
