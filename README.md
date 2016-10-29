# sfdmap

A minimal, fast, MIT-licensed Python module for getting E(B-V) values from
[Schlegel, Finkbeiner & Davis (1998)](http://adsabs.harvard.edu/abs/1998ApJ...500..525S) dust map FITS files.

## Install

This package requires only numpy and a FITS reader of some sort
(either fitsio or astropy).

```
pip install sfdmap  # doesn't yet work!
```

The FITS files comprising the map must be downloaded separately. Among other
locations, they are available from http://github.com/kbarbary/sfddata. On UNIX
systems, run the following to download the maps (~92 MB download size):

```
wget https://github.com/kbarbary/sfddata/archive/master.tar.gz
tar xzf master.tar.gz
```

A directory "sfddata-master" will be created. Move or rename as you like.

## Usage

Initialize map:

```python
>>> import sfdmap

>>> m = sfdmap.SFDMap('/path/to/dustmap/files')
>>> m = sfdmap.SFDMap()  # get directory from SFD_DIR environment variable
```

Get E(B-V) value at RA, Dec = 0., 0. (ICRS frame)

```python
>>> m.ebv(0., 0.)
0.031814847141504288
```

Get E(B-V) at three locations (first argument is RA, second is Dec):

```python
>>> m.ebv([0., 5., 10.], [0., 1.5, 2.1])
array([ 0.03181879,  0.02374864,  0.01746732])
```

By default the coordinates are assumed to be in degrees in the ICRS coordinate
system (e.g., "J2000"). This can be changed with the `frame` and `unit`
keywords:

```python
>>> m.ebv(1.68140, -1.0504884, frame='galactic', unit='radian')
0.031820329230751863
```

The dust map values are linearly interpolated by default. Change this with the
`interpolate` keyword:

```python
>>> m.ebv(1.68140, -1.0504884, frame='galactic', unit='radian', interpolate=False)
0.031526423990726471
```

You can pass an astropy `SkyCoord` instance:

```python
>>> from astropy.coordinates import SkyCoord

>>> coords = SkyCoord([0., 5., 10.], [0., 1.5, 2.1], frame='icrs', unit='degree')

>>> m.ebv(coords)
array([ 0.03181879,  0.02374864,  0.01746732])
```

There is a convenience function in the module so that you don't have to
construct a `SFDMap` instance if you just want to query the map once:

```python
>>> sfdmap.ebv(0., 0.)  # get map directory from SFD_DIR environment variable
0.031818788521008

>>> sfdmap.ebv(0., 0., mapdir='/path/to/dust/files')
0.031818788521008
```

### Performance Note

Note that while passing an astropy `SkyCoord` object works, if you
have coordinates in the ICRS or FK5 (epoch 2000) (e.g., J2000) systems, it
is far faster to pass latitute and longitude directly. This is particularly
true for small numbers of coordinates or scalar coordinates:

```python
In [1]: from astropy.coordinates import SkyCoord

In [2]: import sfdmap

In [3]: m = sfdmap.SFDMap()

In [4]: m.ebv(0., 0.)  # evaluate once to trigger reading the FITS file
Out[4]: 0.03181878852100873

In [5]: coord = SkyCoord(0., 0., unit='degree')

In [6]: %timeit m.ebv(coord)  # time with SkyCoord object
100 loops, best of 3: 18.1 ms per loop

In [7]: %timeit m.ebv(0., 0., unit='degree')  # pass ra, dec directly
10000 loops, best of 3: 80 µs per loop
```

## Alternatives

There are a couple other packages that support getting dust
values from this map. Both these packages have a bigger scope than
this one and include several other 2-d or 3-d galactic dust maps.
Check them out if you want to compare between different maps or need 3-d
maps. Below, I note a few relevant differences from this package.

#### [mwdust](http://github.com/jobovy/mwdust) [[docs](https://pypi.python.org/pypi/mwdust)]

`mwdust.SFD` gives the extinction in a given band rather than E(B-V).
The API is geared towards 3-d maps, so a distance must be given. Python 3
is not currently supported. The license is BSD.

#### [dustmaps](http://github.com/gregreen) [[docs](http://dustmaps.readthedocs.io/en/latest/)]

`dustmaps.sfd.SFDQuery` uses astropy's SkyCoord for all coordinate conversions
and therefore has suboptimal performance for small numbers of coordinates (see "Performance Note" above). The license is GPLv2.
