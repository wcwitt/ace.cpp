import os
import pytest
import sys

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../build/'))
from ace import partitions

def test_partitions():

    num_partitions = []
    for i in range(1,8):
        num_partitions.append(len(partitions(range(1,i+1))))
    assert num_partitions == [1, 2, 5, 15, 52, 203, 877]
