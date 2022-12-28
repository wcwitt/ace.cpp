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

std::complex<double> naive_sph_harm_xyz(int l, int m, double x, double y, double z) {
    double r = std::sqrt(x*x+y*y+z*z);
    double polar = std::acos(z/r);
    double azimuth = ((y>0)-(y<0)) * std::acos(x/std::sqrt(x*x+y*y));
    return naive_sph_harm(l, m, polar, azimuth);
}
