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
# astropy's coordinate conversions have a gigantic overhead of 30-40ms. This
# kills performance in situations where you need to get a single position
# at a time. We can do better by including the core code for the most
# common coordinate conversions (IRCS/FK5J2000) from astropy here.

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

    def __init__(self, fname):
        self.data, header = getdata(fname, header=True)
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
            return self.data[y, x]


class SFDMap(object):
    """Map of E(B-V) from Schlegel, Finkbeiner and Davis (1998).

    This class is useful for repeated retrieval of E(B-V) values when
    there is no way to retrieve all the values at the same time: It keeps
    a reference to the FITS data from the maps so that each FITS image
    is read only once.  Note that there is still a large overhead due to
    coordinate conversion: When possible, pass arrays of coordinates to
    `SFD98Map.get_ebv` or `get_ebv_from_map`.

    Parameters
    ----------
    mapdir : str, optional
        Directory in which to find dust map FITS images, which must be
        named ``SFD_dust_4096_[ngp,sgp].fits``. If not specified,
        the value of the SFD_MAP_DIR configuration item is used. By
        default, this is ``'.'``.  The value of ``SFD_MAP_DIR`` can be set
        in the configuration file, typically located in
        ``$HOME/.astropy/config/sncosmo.cfg``.

    See Also
    --------
    get_ebv_from_map

    Examples
    --------

    >>> m = SFDMap(mapdir='/path/to/SFD98/images')

    Get E(B-V) value at RA, Dec = 0., 0.:

    >>> m.ebv((0., 0.))
    0.031814847141504288

    Get E(B-V) at RA, Dec = (0., 0.) and (1., 0.):

    >>> m.ebv([0., 1.], [0., 0.])
    array([ 0.03181485,  0.03275469])
    """

    def __init__(self, mapdir=None, north="SFD_dust_4096_ngp.fits",
                 south="SFD_dust_4096_sgp.fits"):

        # Get mapdir
        if mapdir is None:
            mapdir = os.environ.get('SFD_DIR', '')
        mapdir = os.path.expanduser(mapdir)
        mapdir = os.path.expandvars(mapdir)

        self.fnames = {'north': os.path.join(mapdir, north),
                       'south': os.path.join(mapdir, south)}

        # don't load maps initially
        self.hemispheres = {'north': None, 'south': None}

        
    def ebv(self, *args, frame='icrs', unit='degree', interpolate=True):
        """Get E(B-V) value(s) at given coordinate(s).

        Parameters
        ----------
        coordinates or ra, dec:

            If one argument is passed, assumed to be a
            `astropy.coordinates.SkyCoords` instance.  If two
            arguments, treated as ``RA, Dec`` (can be scalars or
            arrays).  In the two argument case, the frame and unit is
            take from the keywords.

        frame : {'icrs', 'fk5j2000', 'galactic'}
            Coordinate frame, if two arguments are passed.

        unit : {'degree', 'radian'}

            Unit of coordinates, if two arguments are passed. Default is ``'degree'``.

        interpolate : bool, optional

            Interpolate between the map values using bilinear interpolation.
            Default is True.

        Returns
        -------

        float or `~numpy.ndarray`

            Specific extinction E(B-V) at the given locations.

        """

        if len(args) == 1:
            # treat object as an astropy.coordinates.SkyCoords
            try:
                coordinates = args[0].galactic
                l = coordinates.l.radian
                b = coordinates.b.radian
            except AttributeError:
                raise ValueError("single argument must be astropy.coordinates.SkyCoord")

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
        for label, mask in (('north', b >= 0), ('south', b < 0)):
            if not np.any(mask):
                continue

            # Initialize hemisphere if it hasn't already been done.
            if self.hemispheres[label] is None:
                self.hemispheres[label] = _Hemisphere(self.fnames[label])
                       
            values[mask] = self.hemispheres[label].ebv(l[mask], b[mask],
                                                       interpolate)

        if return_scalar:
            return values[0]
        else:
            return values


def ebv(*args, mapdir=None, frame='icrs', unit='degree', interpolate=True):
    return SFDMap(mapdir=mapdir).ebv(*args, frame=frame, unit=unit,
                                     interpolate=interpolate)
