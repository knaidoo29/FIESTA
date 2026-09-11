import numpy as np
import shift

from typing import List, Union
from .. import coords
from . import part2grid


def mpi_part2grid2D(
    x: np.ndarray,
    y: np.ndarray,
    boxsize: Union[float, List[float]],
    ngrid: Union[int, List[int]],
    MPI: object,
    f: np.ndarray = None,
    method: str = "TSC",
    periodic: Union[bool, List[bool]] = True,
    origin: Union[float, List[float]] = 0.0,
    check_distributed: bool = True,
    dtype: np.dtype = np.float64
) -> np.ndarray:
    """
    Assign particle values to a distributed two-dimensional grid.

    Particles are distributed between MPI ranks according to their
    x coordinate and assigned to the local grid using the selected
    particle-to-grid assignment scheme.

    Ghost cells are added along the MPI-decomposed x direction to
    capture assignment-kernel contributions that cross rank boundaries.
    These contributions are subsequently exchanged with neighbouring
    MPI ranks and accumulated onto the physical grid.

    Parameters
    ----------
    x, y : ndarray
        One-dimensional arrays containing particle coordinates.
    boxsize : float or sequence of float
        Physical size of the global grid. A scalar specifies the same
        length along both axes.
    ngrid : int or sequence of int
        Number of global grid cells. A scalar specifies the same number
        along both axes.
    MPI : object
        MPI utility object controlling the slab decomposition and
        neighbouring-rank communication.
    f : ndarray, optional
        Values associated with each particle. If ``None``, unit particle
        weights are assumed.
    method : {'NGP', 'CIC', 'TSC', 'PCS'}, optional
        Particle-to-grid assignment scheme.
    periodic : bool or sequence of bool, optional
        Periodic boundary conditions along the x and y axes.
    origin : float or sequence of float, optional
        Minimum coordinate of the global grid along each axis.
    check_distributed : bool, optional
        If ``True``, particles are redistributed between MPI ranks
        according to their x coordinate before grid assignment.
    dtype : numpy dtype, optional
        Data type used for the returned grid. Default is ``np.float64``.

    Returns
    -------
    fgrid : ndarray
        Local two-dimensional MPI slab containing the assigned field,
        with shape ``(nx_local, nygrid)``.

    Notes
    -----
    The MPI decomposition is performed along the x direction.

    CIC and TSC require one ghost cell on either side of each local
    slab. PCS requires two ghost cells because its one-dimensional
    assignment stencil spans four neighbouring grid cells. NGP requires
    no ghost cells.

    Periodicity along the x direction is handled through the MPI halo
    exchange rather than by the local particle-to-grid routine.
    """
    method = method.upper()

    # Does a particle field exist?
    fexist = f is not None

    # ------------------------------------------------------------
    # Check/distribute particles
    # ------------------------------------------------------------

    if check_distributed:

        _fexist = int(fexist)

        ftotal = MPI.sum(_fexist)
        ftotal = MPI.broadcast(ftotal)

        # Only distribute f if it exists on every rank.
        fexist = ftotal == MPI.size

        if x is None:
            data = None
        else:
            if fexist:
                data = coords.coord2points([x, y, f])
            else:
                data = coords.coord2points([x, y])

        data = coords.distribute_points_by_x(
            data,
            boxsize,
            ngrid,
            origin,
            MPI,
        )

        x = data[:, 0]
        y = data[:, 1]

        if fexist:
            f = data[:, 2]

    # ------------------------------------------------------------
    # Grid geometry
    # ------------------------------------------------------------

    if np.isscalar(boxsize):
        xlength = boxsize
        ylength = boxsize
    else:
        xlength = boxsize[0]
        ylength = boxsize[1]

    if np.isscalar(origin):
        xmin = origin
        ymin = origin
    else:
        xmin = origin[0]
        ymin = origin[1]

    if np.isscalar(ngrid):
        nxgrid = ngrid
        nygrid = ngrid
    else:
        nxgrid = ngrid[0]
        nygrid = ngrid[1]

    if np.isscalar(periodic):
        periodx = periodic
        periody = periodic
    else:
        periodx = periodic[0]
        periody = periodic[1]

    # ------------------------------------------------------------
    # Assignment stencil / halo width
    # ------------------------------------------------------------

    if method == "NGP":
        nhalo = 0
    elif method in ("CIC", "TSC"):
        nhalo = 1
    elif method == "PCS":
        nhalo = 2
    else:
        raise ValueError(
            "method must be 'NGP', 'CIC', 'TSC' or 'PCS'"
        )

    # ------------------------------------------------------------
    # Local MPI slab
    # ------------------------------------------------------------

    xedges, xgrid = shift.cart.mpi_grid1D(
        xlength,
        nxgrid,
        MPI,
        origin=xmin,
    )

    xmin = xedges[0]
    xmax = xedges[-1]

    dx = xedges[1] - xedges[0]

    nxgrid = len(xgrid)

    # Extend local mesh by required ghost cells.
    if nhalo > 0:
        xmin -= nhalo * dx
        xmax += nhalo * dx
        nxgrid += 2 * nhalo

    xlength = xmax - xmin

    # ------------------------------------------------------------
    # Grid assignment
    # ------------------------------------------------------------

    if method == "NGP":

        if fexist:
            fgrid = part2grid.part2grid_ngp_2d(
                x, y, f,
                xlength, ylength,
                xmin, ymin,
                nxgrid, nygrid,
                dtype=dtype,
            )
        else:
            fgrid = part2grid.part2grid_ngp_2d_unit(
                x, y,
                xlength, ylength,
                xmin, ymin,
                nxgrid, nygrid,
                dtype=dtype,
            )

    elif method == "CIC":

        if fexist:
            fgrid = part2grid.part2grid_cic_2d(
                x, y, f,
                xlength, ylength,
                xmin, ymin,
                nxgrid, nygrid,
                False, periody,
                dtype=dtype,
            )
        else:
            fgrid = part2grid.part2grid_cic_2d_unit(
                x, y,
                xlength, ylength,
                xmin, ymin,
                nxgrid, nygrid,
                False, periody,
                dtype=dtype,
            )

    elif method == "TSC":

        if fexist:
            fgrid = part2grid.part2grid_tsc_2d(
                x, y, f,
                xlength, ylength,
                xmin, ymin,
                nxgrid, nygrid,
                False, periody,
                dtype=dtype,
            )
        else:
            fgrid = part2grid.part2grid_tsc_2d_unit(
                x, y,
                xlength, ylength,
                xmin, ymin,
                nxgrid, nygrid,
                False, periody,
                dtype=dtype,
            )

    else:  # PCS

        if fexist:
            fgrid = part2grid.part2grid_pcs_2d(
                x, y, f,
                xlength, ylength,
                xmin, ymin,
                nxgrid, nygrid,
                False, periody,
                dtype=dtype,
            )
        else:
            fgrid = part2grid.part2grid_pcs_2d_unit(
                x, y,
                xlength, ylength,
                xmin, ymin,
                nxgrid, nygrid,
                False, periody,
                dtype=dtype,
            )

    # ------------------------------------------------------------
    # Exchange ghost-cell contributions
    # ------------------------------------------------------------

    if nhalo > 0:

        fgrid_send_up = MPI.send_up(
            np.ascontiguousarray(fgrid[-nhalo:])
        )

        fgrid_send_down = MPI.send_down(
            np.ascontiguousarray(fgrid[:nhalo])
        )

        # Remove ghost cells.
        fgrid = fgrid[nhalo:-nhalo]

        if periodx or MPI.rank > 0:
            fgrid[:nhalo] += fgrid_send_up

        if periodx or MPI.rank < MPI.size - 1:
            fgrid[-nhalo:] += fgrid_send_down

    return fgrid

def mpi_part2grid3D(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    boxsize: Union[float, List[float]],
    ngrid: Union[int, List[int]],
    MPI: object,
    f: np.ndarray = None,
    method: str = "TSC",
    periodic: Union[bool, List[bool]] = True,
    origin: Union[float, List[float]] = 0.0,
    check_distributed: bool = True,
    dtype: np.dtype = np.float64
) -> np.ndarray:
    """
    Assign particle values to a distributed three-dimensional grid.

    Particles are distributed between MPI ranks according to their
    x coordinate and assigned to the local grid using the selected
    particle-to-grid assignment scheme.

    Ghost cells are added along the MPI-decomposed x direction to
    capture assignment-kernel contributions that cross rank boundaries.
    These contributions are subsequently exchanged with neighbouring
    MPI ranks and accumulated onto the physical grid.

    Parameters
    ----------
    x, y, z : ndarray
        One-dimensional arrays containing particle coordinates.
    boxsize : float or sequence of float
        Physical size of the global grid. A scalar specifies the same
        length along all three axes.
    ngrid : int or sequence of int
        Number of global grid cells. A scalar specifies the same number
        along all three axes.
    MPI : object
        MPI utility object controlling the slab decomposition and
        neighbouring-rank communication.
    f : ndarray, optional
        Values associated with each particle. If ``None``, unit particle
        weights are assumed.
    method : {'NGP', 'CIC', 'TSC', 'PCS'}, optional
        Particle-to-grid assignment scheme.
    periodic : bool or sequence of bool, optional
        Periodic boundary conditions along the x, y and z axes.
    origin : float or sequence of float, optional
        Minimum coordinate of the global grid along each axis.
    check_distributed : bool, optional
        If ``True``, particles are redistributed between MPI ranks
        according to their x coordinate before grid assignment.
    dtype : numpy dtype, optional
        Data type used for the returned grid. Default is ``np.float64``.

    Returns
    -------
    fgrid : ndarray
        Local three-dimensional MPI slab containing the assigned field,
        with shape ``(nx_local, nygrid, nzgrid)``.

    Notes
    -----
    The MPI decomposition is performed along the x direction.

    CIC and TSC require one ghost cell on either side of each local
    slab. PCS requires two ghost cells because its one-dimensional
    assignment stencil spans four neighbouring grid cells. NGP requires
    no ghost cells.

    Periodicity along the x direction is handled through the MPI halo
    exchange rather than by the local particle-to-grid routine.
    """
    method = method.upper()

    # Does a particle field exist?
    fexist = f is not None

    # ------------------------------------------------------------
    # Check/distribute particles
    # ------------------------------------------------------------

    if check_distributed:

        _fexist = int(fexist)

        ftotal = MPI.sum(_fexist)
        ftotal = MPI.broadcast(ftotal)

        # Only distribute f if it exists on every rank.
        fexist = ftotal == MPI.size

        if x is None:
            data = None
        else:
            if fexist:
                data = coords.coord2points([x, y, z, f])
            else:
                data = coords.coord2points([x, y, z])

        data = coords.distribute_points_by_x(
            data,
            boxsize,
            ngrid,
            origin,
            MPI,
        )

        x = data[:, 0]
        y = data[:, 1]
        z = data[:, 2]

        if fexist:
            f = data[:, 3]

    # ------------------------------------------------------------
    # Grid geometry
    # ------------------------------------------------------------

    if np.isscalar(boxsize):
        xlength = boxsize
        ylength = boxsize
        zlength = boxsize
    else:
        xlength = boxsize[0]
        ylength = boxsize[1]
        zlength = boxsize[2]

    if np.isscalar(origin):
        xmin = origin
        ymin = origin
        zmin = origin
    else:
        xmin = origin[0]
        ymin = origin[1]
        zmin = origin[2]

    if np.isscalar(ngrid):
        nxgrid = ngrid
        nygrid = ngrid
        nzgrid = ngrid
    else:
        nxgrid = ngrid[0]
        nygrid = ngrid[1]
        nzgrid = ngrid[2]

    if np.isscalar(periodic):
        periodx = periodic
        periody = periodic
        periodz = periodic
    else:
        periodx = periodic[0]
        periody = periodic[1]
        periodz = periodic[2]

    # ------------------------------------------------------------
    # Assignment stencil / halo width
    # ------------------------------------------------------------

    if method == "NGP":
        nhalo = 0
    elif method in ("CIC", "TSC"):
        nhalo = 1
    elif method == "PCS":
        nhalo = 2
    else:
        raise ValueError(
            "method must be 'NGP', 'CIC', 'TSC' or 'PCS'"
        )

    # ------------------------------------------------------------
    # Local MPI slab
    # ------------------------------------------------------------

    xedges, xgrid = shift.cart.mpi_grid1D(
        xlength,
        nxgrid,
        MPI,
        origin=xmin,
    )

    xmin = xedges[0]
    xmax = xedges[-1]

    dx = xedges[1] - xedges[0]

    nxgrid = len(xgrid)

    if nhalo > 0:
        xmin -= nhalo * dx
        xmax += nhalo * dx
        nxgrid += 2 * nhalo

    xlength = xmax - xmin

    # ------------------------------------------------------------
    # Grid assignment
    # ------------------------------------------------------------

    if method == "NGP":

        if fexist:
            fgrid = part2grid.part2grid_ngp_3d(
                x,
                y,
                z,
                f,
                xlength,
                ylength,
                zlength,
                xmin,
                ymin,
                zmin,
                nxgrid,
                nygrid,
                nzgrid,
                dtype=dtype,
            )
        else:
            fgrid = part2grid.part2grid_ngp_3d_unit(
                x,
                y,
                z,
                xlength,
                ylength,
                zlength,
                xmin,
                ymin,
                zmin,
                nxgrid,
                nygrid,
                nzgrid,
                dtype=dtype,
            )

    elif method == "CIC":

        if fexist:
            fgrid = part2grid.part2grid_cic_3d(
                x,
                y,
                z,
                f,
                xlength,
                ylength,
                zlength,
                xmin,
                ymin,
                zmin,
                nxgrid,
                nygrid,
                nzgrid,
                False,
                periody,
                periodz,
                dtype=dtype,
            )
        else:
            fgrid = part2grid.part2grid_cic_3d_unit(
                x,
                y,
                z,
                xlength,
                ylength,
                zlength,
                xmin,
                ymin,
                zmin,
                nxgrid,
                nygrid,
                nzgrid,
                False,
                periody,
                periodz,
                dtype=dtype,
            )

    elif method == "TSC":

        if fexist:
            fgrid = part2grid.part2grid_tsc_3d(
                x,
                y,
                z,
                f,
                xlength,
                ylength,
                zlength,
                xmin,
                ymin,
                zmin,
                nxgrid,
                nygrid,
                nzgrid,
                False,
                periody,
                periodz,
                dtype=dtype,
            )
        else:
            fgrid = part2grid.part2grid_tsc_3d_unit(
                x,
                y,
                z,
                xlength,
                ylength,
                zlength,
                xmin,
                ymin,
                zmin,
                nxgrid,
                nygrid,
                nzgrid,
                False,
                periody,
                periodz,
                dtype=dtype,
            )

    else:  # PCS

        if fexist:
            fgrid = part2grid.part2grid_pcs_3d(
                x,
                y,
                z,
                f,
                xlength,
                ylength,
                zlength,
                xmin,
                ymin,
                zmin,
                nxgrid,
                nygrid,
                nzgrid,
                False,
                periody,
                periodz,
                dtype=dtype,
            )
        else:
            fgrid = part2grid.part2grid_pcs_3d_unit(
                x,
                y,
                z,
                xlength,
                ylength,
                zlength,
                xmin,
                ymin,
                zmin,
                nxgrid,
                nygrid,
                nzgrid,
                False,
                periody,
                periodz,
                dtype=dtype,
            )

    # ------------------------------------------------------------
    # Exchange ghost-cell contributions
    # ------------------------------------------------------------

    if nhalo > 0:

        fgrid_send_up = MPI.send_up(fgrid[-nhalo:])
        fgrid_send_down = MPI.send_down(fgrid[:nhalo])

        fgrid = fgrid[nhalo:-nhalo]

        if periodx or MPI.rank > 0:
            fgrid[:nhalo] += fgrid_send_up

        if periodx or MPI.rank < MPI.size - 1:
            fgrid[-nhalo:] += fgrid_send_down

    return fgrid