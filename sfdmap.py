# Licensed under an MIT "Expat" license - See LICENSE
"""Get E(B-V) values from the Schlegel, Finkbeiner & Davis (1998) dust map."""

import os

import numpy as np
from scipy.ndimage import map_coordinates
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.io import fits
from astropy.utils import isiterable

__all__ = ['SFDMap', 'ebv']


class _Hemisphere(object):
    def __init__(self, fname):
        hdulist = fits.open(fname)
        header = hdulist[0].header
        self.crpix1 = header['CRPIX1']
        self.crpix2 = header['CRPIX2']
        self.lam_scal = header['LAM_SCAL']
        self.sign = header['LAM_NSGP']  # north = 1, south = -1
        self.data = hdulist[0].data
        hdulist.close()

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
            return map_coordinates(self.data, [y, x], order=1)
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

        
    def ebv(self, *args, frame='ircs', interpolate=True):
        """Get E(B-V) value(s) at given coordinate(s).

        Parameters
        ----------
        coordinates or ra, dec:

            If one argument is passed, assumed to be a
            `astropy.coordinates.SkyCoords` instance.  If two
            arguments, treated as ``RA, Dec`` in degrees in the ICRS
            (e.g., "J2000") system. RA and Dec can each be float or
            list or numpy array.

        interpolate : bool, optional

            Interpolate between the map values using
            `scipy.ndimage.map_coordinates`. Default is ``True``.

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
                raise ValueError("single argument must be SkyCoord")

        elif len(args) == 2:
            lat, lon = args
            lat = np.deg2rad(lat)
            lon = np.deg2rad(lon)
            if frame == 'galactic':
                l, b = lat, lon
            else:
                c = SkyCoord(ra=lat, dec=lon, frame='icrs',
                             unit=u.rad).galactic
                l = c.l.radian
                b = c.b.radian

        else:
            raise ValueError("too many arguments")


        # Check if l, b are scalar. If so, convert to 1-d arrays.
        return_scalar = False
        if not isiterable(l):
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


def ebv(*args, mapdir=None, interpolate=True):
    """Get E(B-V) value(s) from Schlegel, Finkbeiner, and Davis (1998)
    extinction maps.

    Parameters
    ----------

    coordinates : astropy Coordinates object or tuple

        If tuple, treated as (RA, Dec) in degrees in the ICRS (e.g.,
        "J2000") system. RA and Dec can each be float or list or numpy
        array.

    mapdir : str, optional

        Directory in which to find dust map FITS images, which must be
        named ``SFD_dust_4096_[ngp,sgp].fits``. If `None` (default),
        the value of the ``sfd98_dir`` configuration item is used. By
        default, this is ``'.'``.  The value of ``sfd98_dir`` can be set
        in the configuration file, typically located in
        ``$HOME/.astropy/config/sncosmo.cfg``.

    interpolate : bool
        Interpolate between the map values using
        `scipy.ndimage.map_coordinates`. Default is ``True``.
    order : int
        Interpolation order used for interpolate=True. Default is 1.

    Returns
    -------
    ebv : float or `~numpy.ndarray`
        Specific extinction E(B-V) at the given locations.

    Examples
    --------

    Get E(B-V) value at RA, Dec = (0., 0.):

    >>> get_ebv_from_map((0., 0.), mapdir='/path/to/dir')
    0.031814847141504288

    Get E(B-V) at RA, Dec = (0., 0.) and (1., 0.):

    >>> get_ebv_from_map(([0., 1.], [0., 0.]), mapdir='/path/to/dir')
    array([ 0.03181485,  0.03275469])
    """

    return SFDMap(mapdir=mapdir).ebv(*args, interpolate=interpolate)
