import numpy as np
from numba import njit


@njit
def weight_cic(xp, xg, dl):
    """
    Cloud-in-Cell (CIC) weights in 1D.
    
    Parameters:
    xp : float
        Coordinate of the particle.
    xg : float
        Coordinate of the grid point.
    dl : float
        Size of the grid cell.

    Returns:
    w : float
        Weight computed based on the distance between the particle and the grid point.
    """
    # Compute the absolute distance between the particle and the grid point
    dx = abs(xp - xg)
    
    # Calculate the weight using the Cloud-in-Cell method
    if dx <= dl:
        w = (1. - (dx / dl)) / dl
    else:
        w = 0.0
    
    return w


@njit
def weight_tsc(xp, xg, dl):
    """
    Triangular-shaped-cloud (TSC) weights in 1D.
    
    Parameters:
    xp : float
        Coordinate of the particle.
    xg : float
        Coordinate of the grid point.
    dl : float
        Size of the grid cell.
    
    Returns:
    w : float
        Weight computed based on the triangular-shaped cloud function.
    """
    # Compute the absolute distance between the particle and the grid point
    dx = abs(xp - xg)
    
    # Calculate the weight using the Triangular-Shaped-Cloud method
    if dx <= dl / 2.:
        w = 3. / 4. - (dx / dl) ** 2.
    elif dx <= 3. * dl / 2.:
        w = (1./2.) * ((3. / 2. - (dx / dl)) ** 2.)
    else:
        w = 0.0
    
    # Normalize the weight by the cell size
    w = w / dl
    
    return w


@njit
def weight_pcs(xp, xg, dl):
    """
    Piecewise-Cubic-Spline (PCS) weights in 1D.
    
    Parameters:
    xp : float
        Coordinate of the particle.
    xg : float
        Coordinate of the grid point.
    dl : float
        Size of the grid cell.
    
    Returns:
    w : float
        Weight computed based on the Piecewise-Cubic-Spline function.
    """
    # Compute the absolute distance between the particle and the grid point
    dx = abs(xp - xg)
    
    # Calculate the weight using the Piecewise-Cubic-Spline method
    if dx <= dl:
        w = (1. / 6.) * (4. - 6. * (dx / dl) ** 2. + 3. * (dx / dl) ** 3.)
    elif dx <= 2. * dl:
        w = (1. / 6.) * ((2. - (dx / dl)) ** 3.)
    else:
        w = 0.0
    
    # Normalize the weight by the cell size
    w = w / dl
    
    return w
