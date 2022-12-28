import numpy as np
import os
import pytest
import sys
from scipy.optimize import bisect
from scipy.special import jn_zeros
from scipy.special import spherical_jn

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../build/'))
from ace import spherical_bessel_radial

def test_spherical_bessel_radial():

    def scipy_spherical_bessel_radial(n, l, r, r_c):
    
        a, b = jn_zeros(l, n)[-1], jn_zeros(l+1, n)[-1]
        z_nl = bisect(lambda x: spherical_jn(l, x), a, b, xtol=1e-15)
        if l==0:
            c_nlm = (z_nl/r_c)**(1.5)*(0.5*z_nl-0.25*np.sin(2*z_nl))**(-0.5)
        else:
            c_nlm = (-0.5*r_c**3*spherical_jn(l-1,z_nl)*spherical_jn(l+1,z_nl))**(-0.5)
        return c_nlm * spherical_jn(l, z_nl*r/r_c)

    for l in range(0, 7):
        for n in range(1, 15):
            for _ in range(50):
                r = np.random.uniform(0.1, 15)
                r_c = r + np.random.uniform(0, 15)
                x1 = scipy_spherical_bessel_radial(n, l, r, r_c)
                x2 = spherical_bessel_radial(n, l, r, r_c)
                assert x1 == pytest.approx(x2)
