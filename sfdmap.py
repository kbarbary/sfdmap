# Licensed under an MIT "Expat" license - See LICENSE
"""Get E(B-V) values from the Schlegel, Finkbeiner & Davis (1998) dust map."""

import os

import numpy as np

# require a FITS reader of some sort.
try:
    import fitsio
    getdata = fitsio.read
except ImportError:
    try:
        from astropy.io.fits import getdata
    except ImportError:
        raise ImportError("could not import fitsio or astropy.io.fits")

__all__ = ['SFDMap', 'ebv']
__version__ = "0.1.0"


def _isiterable(obj):
    """Returns `True` if the given object is iterable."""

    try:
        iter(obj)
        return True
    except TypeError:
        return False

# -----------------------------------------------------------------------------
# Coordinate conversion
#
# astropy's coordinate conversions have a gigantic overhead of about
# 20ms.  This kills performance in situations where you need to get a
# single position at a time. We can do better (well under 0.1ms!) for
# common static systems (IRCS, FK5J2000) by including a bit of code
# from astropy.coordinates here.

# Create rotation matrix about a given axis (x, y, z)
def zrotmat(angle):
    s = np.sin(angle)
    c = np.cos(angle)
    return np.array(((c, s, 0),
                     (-s, c, 0),
                     (0, 0, 1)))

def yrotmat(angle):
    s = np.sin(angle)
    c = np.cos(angle)
    return np.array(((c, 0, -s),
                     (0, 1, 0),
                     (s, 0, c)))

def xrotmat(angle):
    s = np.sin(angle)
    c = np.cos(angle)
    return np.array(((1, 0, 0),
                     (0, c, s),
                     (0, -s, c)))

# constant ICRS --> FK5J2000 (See USNO Circular 179, section 3.5)
eta0 = np.deg2rad(-19.9 / 3600000.)
xi0 = np.deg2rad(9.1 / 3600000.)
da0 = np.deg2rad(-22.9 / 3600000.)
ICRS_TO_FK5J2000 = np.dot(np.dot(xrotmat(-eta0), yrotmat(xi0)), zrotmat(da0))

# FK5J2000 --> Gal
# Note that galactic pole and zeropoint of l are somewhat arbitrary
# and not officially defined (as far as I know). The values below are
# from astropy.coordinates, which includes the following comment:
# | "This gives better consistency with other codes than using the values
# |  from Reid & Brunthaler 2004 and the best self-consistency between FK5
# |  -> Galactic and FK5 -> FK4 -> Galactic. The lon0 angle was found by
# |  optimizing the self-consistency."
ngp_fk5j2000_ra = np.deg2rad(192.8594812065348)
ngp_fk5j2000_dec = np.deg2rad(27.12825118085622)
lon0_fk5j2000 = np.deg2rad(122.9319185680026)
FK5J2000_TO_GAL = np.dot(np.dot(zrotmat(np.pi - lon0_fk5j2000),
                                yrotmat(np.pi/2. - ngp_fk5j2000_dec)),
                         zrotmat(ngp_fk5j2000_ra))

# ICRS --> Gal: simply chain through FK5J2000 
ICRS_TO_GAL = np.dot(FK5J2000_TO_GAL, ICRS_TO_FK5J2000)


# (lon, lat) -> [x, y, z] unit vector
def coords2cart(lon, lat):
    coslat = np.cos(lat)
    return np.array((coslat * np.cos(lon), coslat * np.sin(lon), np.sin(lat)))


# [x, y, z] unit vector -> (lon, lat)
def cart2coords(xyz):
    x, y, z = xyz
    return np.arctan2(y, x), np.arctan2(z, np.sqrt(x*x + y*y))


def _icrs_to_gal(lon, lat):
    return cart2coords(np.dot(ICRS_TO_GAL, coords2cart(lon, lat)))


def _fk5j2000_to_gal(lon, lat):
    return cart2coords(np.dot(FK5J2000_TO_GAL, coords2cart(lon, lat)))


# -----------------------------------------------------------------------------
# bilinear interpolation (because scipy.ndimage.map_coordinates is really slow
# for this)

def _bilinear_interpolate(data, y, x):
    yfloor = np.floor(y)
    xfloor = np.floor(x)
    yw = y - yfloor
    xw = x - xfloor

    # pixel locations
    y0 = yfloor.astype(np.int)
    y1 = y0 + 1
    x0 = xfloor.astype(np.int)
    x1 = x0 + 1

    # clip locations out of range
    ny, nx = data.shape
    y0 = np.maximum(y0, 0)
    y1 = np.minimum(y1, ny-1)
    x0 = np.maximum(x0, 0)
    x1 = np.minimum(x1, nx-1)

    return ((1.0-xw) * (1.0-yw) * data[y0, x0] +
            xw       * (1.0-yw) * data[y0, x1] + 
            (1.0-xw) * yw       * data[y1, x0] +
            xw       * yw       * data[y1, x1])   


# -----------------------------------------------------------------------------

class _Hemisphere(object):
    """Represents one of the hemispheres (in a single fle)"""

    def __init__(self, fname, scaling):
        self.data, header = getdata(fname, header=True)
        self.data *= scaling
        self.crpix1 = header['CRPIX1']
        self.crpix2 = header['CRPIX2']
        self.lam_scal = header['LAM_SCAL']
        self.sign = header['LAM_NSGP']  # north = 1, south = -1
        
    def ebv(self, l, b, interpolate):
        # Project from galactic longitude/latitude to lambert pixels.
        # (See SFD98 or SFD data FITS header).
        x = (self.crpix1 - 1.0 +
             self.lam_scal * np.cos(l) *
             np.sqrt(1.0 - self.sign * np.sin(b)))
        y = (self.crpix2 - 1.0 -
             self.sign * self.lam_scal * np.sin(l) *
             np.sqrt(1.0 - self.sign * np.sin(b)))
            
        # Get map values at these pixel coordinates.
        if interpolate:
            return _bilinear_interpolate(self.data, y, x)
        else:
            x = np.round(x).astype(np.int)
            y = np.round(y).astype(np.int)

            # some valid coordinates are right on the border (e.g., x/y = 4096)
            x = np.clip(x, 0, self.data.shape[1]-1)
            y = np.clip(y, 0, self.data.shape[0]-1) 
            return self.data[y, x]


class SFDMap(object):
    """Map of E(B-V) from Schlegel, Finkbeiner and Davis (1998).

    Use this class for repeated retrieval of E(B-V) values when
    there is no way to retrieve all the values at the same time: It keeps
    a reference to the FITS data from the maps so that each FITS image
    is read only once.

    Parameters
    ----------

    mapdir : str, optional

        Directory in which to find dust map FITS images, named
        ``SFD_dust_4096_ngp.fits`` and ``SFD_dust_4096_sgp.fits`` by
        default. If not specified, the value of the ``SFD_DIR``
        environment variable is used, otherwise an empty string is
        used.

    north, south : str, optional

        Names of north and south galactic pole FITS files. Defaults are
        ``SFD_dust_4096_ngp.fits`` and ``SFD_dust_4096_sgp.fits``
        respectively.

    scaling : float, optional
        Scale all E(B-V) map values by this factor. Default is 0.86,
        corresponding to recalibration from Schlafly & Finkbeiner (2011).
    """

    def __init__(self, mapdir=None, north="SFD_dust_4096_ngp.fits",
                 south="SFD_dust_4096_sgp.fits", scaling=0.86):

        if mapdir is None:
            mapdir = os.environ.get('SFD_DIR', '')
        mapdir = os.path.expanduser(mapdir)
        mapdir = os.path.expandvars(mapdir)
        self.mapdir = mapdir

        # don't load maps initially
        self.fnames = {'north': north, 'south': south}
        self.hemispheres = {'north': None, 'south': None}

        self.scaling = scaling

    def ebv(self, *args, **kwargs):
        """Get E(B-V) value(s) at given coordinate(s).

        Parameters
        ----------

        coordinates or ra, dec: SkyCoord or numpy.ndarray

            If one argument is passed, assumed to be a
            `astropy.coordinates.SkyCoords` instance. In this case,
            the ``frame`` and ``unit`` keyword arguments are
            ignored. If two arguments are passed, they are treated as
            ``latitute, longitude`` (can be scalars or arrays).  In
            the two argument case, the frame and unit is taken from
            the keywords.

        frame : {'icrs', 'fk5j2000', 'galactic'}, optional
            Coordinate frame, if two arguments are passed. Default is
            ``'icrs'``.

        unit : {'degree', 'radian'}, optional

            Unit of coordinates, if two arguments are passed. Default
            is ``'degree'``.

        interpolate : bool, optional

            Interpolate between the map values using bilinear interpolation.
            Default is True.

        Returns
        -------

        `~numpy.ndarray`

            Specific extinction E(B-V) at the given locations.

        """

        # collect kwargs
        frame = kwargs.get('frame', 'icrs')
        unit = kwargs.get('unit', 'degree')
        interpolate = kwargs.get('interpolate', True)

        # compatibility: treat single argument 2-tuple as (RA, Dec)
        if ((len(args) == 1) and (type(args[0]) is tuple) and
            (len(args[0]) == 2)):
            args = args[0]

        if len(args) == 1:
            # treat object as an astropy.coordinates.SkyCoords
            try:
                coordinates = args[0].galactic
                l = coordinates.l.radian
                b = coordinates.b.radian
            except AttributeError:
                raise ValueError("single argument must be "
                                 "astropy.coordinates.SkyCoord")

        elif len(args) == 2:
            lat, lon = args

            # convert to radians
            if unit in ('deg', 'degree'):
                lat = np.deg2rad(lat)
                lon = np.deg2rad(lon)
            elif unit in ('rad', 'radian'):
                pass
            else:
                raise ValueError("unit not understood")

            # convert to galactic
            if frame == 'icrs':
                l, b = _icrs_to_gal(lat, lon)
            elif frame == 'galactic':
                l, b = lat, lon
            elif frame in ('fk5j2000', 'j2000'):
                l, b = _fk5j2000_to_gal(lat, lon)
            else:
                raise ValueError("frame not understood")

        else:
            raise ValueError("too many arguments")


        # Check if l, b are scalar. If so, convert to 1-d arrays.
        return_scalar = False
        if not _isiterable(l):
            return_scalar = True
            l, b = np.array([l]), np.array([b])

        # Initialize return array
        values = np.empty_like(l)

        # Treat north (b>0) separately from south (b<0).
        for pole, mask in (('north', b >= 0), ('south', b < 0)):
            if not np.any(mask):
                continue

            # Initialize hemisphere if it hasn't already been done.
            if self.hemispheres[pole] is None:
                fname = os.path.join(self.mapdir, self.fnames[pole])
                self.hemispheres[pole] = _Hemisphere(fname, self.scaling)

            values[mask] = self.hemispheres[pole].ebv(l[mask], b[mask],
                                                      interpolate)

        if return_scalar:
            return values[0]
        else:
            return values

    def __repr__(self):
        return ("SFDMap(mapdir={!r}, north={!r}, south={!r}, scaling={!r})"
                .format(self.mapdir, self.fnames['north'],
                        self.fnames['south'], self.scaling))


def ebv(*args, **kwargs):
    """Convenience function, equivalent to SFDMap().ebv(*args)"""
    
    m = SFDMap(mapdir=kwargs.get('mapdir', None),
               north=kwargs.get('north', "SFD_dust_4096_ngp.fits"),
               south=kwargs.get('south', "SFD_dust_4096_sgp.fits"),
               scaling=kwargs.get('scaling', 0.86))
    return m.ebv(*args, **kwargs)
