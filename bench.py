#!/usr/bin/env python

import time
import numpy as np
import sfdmap
from astropy.coordinates import SkyCoord

# time single coordinate access
print("access one coordinate at a time")
for n in [1, 10, 100]:
    t = time.time()
    m = sfdmap.SFDMap()
    for _ in range(n):
        m.ebv(0., 0.)
    t = time.time() - t
    print("{:5d} : {:7.2f} ms".format(n, t*1000.))

# time single coordinate access
print("access array of coordinates")
for n in [1, 10, 100]:
    ra, dec = np.zeros(n), np.zeros(n)
    t = time.time()
    m = sfdmap.SFDMap()
    m.ebv(ra, dec)
    t = time.time() - t
    print("{:5d} : {:7.2f} ms".format(n, t*1000.))
