#!/usr/bin/env python

import time
import numpy as np
import sfdmap
from astropy.coordinates import SkyCoord

interpolate = True

# get one coordinate to trigger loading the file for the first time
# The first time, the file might or might not be cached in RAM.
for temp in ("possibly cold", "hot"):
    print("access one coordinate (FITS file read; {})".format(temp))
    t = time.time()
    sfdmap.ebv(0., 0.)
    t = time.time() - t
    print("{:5d} : {:7.2f} ms".format(1, t*1000.))

m = sfdmap.SFDMap()
m.ebv(0., 0.)  # trigger file load

# time single coordinate access
print("access one coordinate at a time (after initial read)")
for n in [1, 10, 100]:
    t = time.time()
    for _ in range(n):
        m.ebv(0., 0., interpolate=interpolate)
    t = time.time() - t
    print("{:5d} : {:7.2f} ms".format(n, t*1000.))

# time single coordinate access
print("access array of coordinates (after initial read)")
for n in [1, 10, 100]:
    ra, dec = np.zeros(n), np.zeros(n)
    t = time.time()
    m.ebv(ra, dec, interpolate=interpolate)
    t = time.time() - t
    print("{:5d} : {:7.2f} ms".format(n, t*1000.))
