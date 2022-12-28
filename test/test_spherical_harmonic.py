import numpy as np
import os
import pytest
import sys
from scipy.special import sph_harm as scipy_sph_harm

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../build/'))
from ace import naive_sph_harm

def test_naive_sph_harm():
    for l in range(0,11):
        for m in range(-l,l+1):
            for _ in range(100):
                polar = np.random.uniform(0, np.pi)
                azimuth = np.random.uniform(0, 2*np.pi)
                sh1 = naive_sph_harm(l, m, polar, azimuth)
                sh2 = scipy_sph_harm(m, l, azimuth, polar)
                assert sh1 == pytest.approx(sh2)
