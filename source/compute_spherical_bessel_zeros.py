import numpy as np
from scipy.optimize import bisect
from scipy.special import jn_zeros, spherical_jn

num_l = 10
num_n = 50

def spherical_jn_zero(l, n):
    a, b = jn_zeros(l, n)[-1], jn_zeros(l+1, n)[-1]
    return bisect(lambda x: spherical_jn(l,x), a, b, xtol=1e-15)

with open('spherical_bessel_zeros.hpp', 'w') as f:
    f.write('#include <array>\n\n')
    f.write('const std::array<std::array<double,{}>,{}> spherical_bessel_zeros =\n'.format(num_n,num_l))
    f.write('    {{\n')
    for l in range(num_l):
        f.write('        {{   // l = {}\n'.format(l))
        for n in range(1,num_n+1):
            f.write('            {:20.15e}'.format(spherical_jn_zero(l,n)))
            if n < num_n:
                f.write(',')
            f.write('\n')
        f.write('        }')
        if l < num_l-1:
            f.write(',')
        f.write('\n')
    f.write('    }};')
