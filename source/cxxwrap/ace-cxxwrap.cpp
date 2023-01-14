#include "jlcxx/jlcxx.hpp"
#include "jlcxx/stl.hpp"

#include "acecore/simpleprodbasis.hpp"
#include "acecore/sparsesymmprod.hpp"
#include "acecore/symmprod_dag_kernels.hpp"
#include "partitions.hpp"
#include "spherical_harmonic.hpp"
#include "spherical_bessel.hpp"

JLCXX_MODULE define_julia_module(jlcxx::Module& m)
{
    m.method("naive_sph_harm", &naive_sph_harm);
    m.method("naive_sph_harm_xyz", &naive_sph_harm_xyz);
    m.method("spherical_bessel_radial", &spherical_bessel_radial);
    m.method("determine_basis", &determine_basis);
    m.method("partitions", &partitions);
    
    m.method("SimpleProdBasis_evaluate", &SimpleProdBasis_evaluate);
    m.method("SparseSymmProd_evaluate", &SparseSymmProd_evaluate);
    m.method("evaluate_real_stl", &evaluate_real_stl);
    m.method("evaluate_real", &evaluate_real);
    m.method("evaluate_complex", &evaluate_complex);
    m.method("evaluate_batch_real", &evaluate_batch_real);
    m.method("evaluate_batch_complex", &evaluate_batch_complex);
}
