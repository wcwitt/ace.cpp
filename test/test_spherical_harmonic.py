import numpy as np
import os
import pytest
import sys
from scipy.special import sph_harm as scipy_sph_harm

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../build/'))
from ace import naive_sph_harm, naive_sph_harm_xyz

def test_naive_sph_harm():

    for l in range(0,11):
        for m in range(-l,l+1):
            for _ in range(100):
                polar = np.random.uniform(0, np.pi)
                azimuth = np.random.uniform(0, 2*np.pi)
                sh1 = naive_sph_harm(l, m, polar, azimuth)
                sh2 = scipy_sph_harm(m, l, azimuth, polar)
                assert sh1 == pytest.approx(sh2)

def test_naive_sph_harm_xyz():

    def scipy_sph_harm_xyz(l, m, x, y, z):
        r = np.sqrt(x*x + y*y + z*z)
        polar = np.arccos(z/r)
        azimuth = np.arctan2(y,x) + (y<0)*2*np.pi
        return scipy_sph_harm(m,l,azimuth,polar)

    for l in range(0,11):
        for m in range(-l,l+1):
            for _ in range(100):
                xyz = np.random.uniform(-1,1,3)
                sh1 = naive_sph_harm_xyz(l, m, xyz[0], xyz[1], xyz[2])
                sh2 = scipy_sph_harm_xyz(l, m, xyz[0], xyz[1], xyz[2])
                assert sh1 == pytest.approx(sh2)

