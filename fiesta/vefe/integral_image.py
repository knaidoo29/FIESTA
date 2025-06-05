import numpy as np


def get_integral_image2D(fgrid):
    """
    Computes the 2D integral image.

    Parameters
    ----------
    fgrid : array_like
        The 2D field grid.
    
    Returns
    -------
    igrid : array_like
        The 2D integral image.
    """
    shape = np.shape(fgrid)
    _fgrid = np.zeros((shape[0]+1, shape[1]+1))
    _fgrid[1:,1:] = fgrid
    igrid = np.cumsum(_fgrid, axis=0, dtype=np.float64)
    igrid = np.cumsum(igrid, axis=1, dtype=np.float64)
    return igrid


def get_integral_image3D(fgrid):
    """
    Computes the 3D integral image.

    Parameters
    ----------
    fgrid : array_like
        The 3D field grid.
    
    Returns
    -------
    igrid : array_like
        The 3D integral image.
    """
    shape = np.shape(fgrid)
    _fgrid = np.zeros((shape[0]+1, shape[1]+1, shape[2]+1))
    _fgrid[1:,1:,1:] = fgrid
    igrid = np.cumsum(_fgrid, axis=0, dtype=np.float64)
    igrid = np.cumsum(igrid, axis=1, dtype=np.float64)
    igrid = np.cumsum(igrid, axis=2, dtype=np.float64)
    return igrid