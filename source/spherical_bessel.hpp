#include <tuple>
#include <vector>

double spherical_bessel_radial(int n, int l, double r, double r_c);

std::tuple<std::vector<int>,
           std::vector<int>,
           std::vector<int>,
           std::vector<double>>
    determine_basis(double r_max, double e_max);
