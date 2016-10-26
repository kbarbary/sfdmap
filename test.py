#!/usr/bin/env py.test

import numpy as np
from numpy.testing import assert_allclose
from astropy.coordinates import SkyCoord

import sfdmap


def test_single():
    # value from http://irsa.ipac.caltech.edu/applications/DUST/
    # with input "204.0 -30.0 J2000"
    true_ebv = 0.0477

    # Use interpolate=False to match IRSA value
    ebv = sfdmap.ebv(204.0, -30.0, interpolate=False)
    assert_allclose(ebv, true_ebv, rtol=0.01)

def test_skycoord():
    coords = SkyCoord(204.0, -30.0, frame='icrs', unit='deg')

    # value from http://irsa.ipac.caltech.edu/applications/DUST/
    # with input "204.0 -30.0 J2000"
    true_ebv = 0.0477

    ebv = sfdmap.ebv(coords, interpolate=False)
    assert_allclose(ebv, true_ebv, rtol=0.01)
