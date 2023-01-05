#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/stl.h>

#include "partitions.hpp"
#include "spherical_bessel.hpp"
#include "spherical_harmonic.hpp"

PYBIND11_MODULE(ace, m) {
    m.doc() = "pybind11 example plugin"; // optional module docstring

    m.def("determine_basis", &determine_basis, "Docstring.");
    m.def("naive_sph_harm", &naive_sph_harm, "A naive spherical harmonic.");
    m.def("naive_sph_harm_xyz", &naive_sph_harm_xyz, "A naive spherical harmonic.");
    m.def("partitions", &partitions);
    m.def("spherical_bessel_radial", &spherical_bessel_radial, "Docstring.");
}
