import numpy as np
import pytest

pytest.importorskip("mpi4py")
import shift
from shift.mpiutils import MPI

from fiesta.maths import dfdx, dfdy, dfdz, mpi_dfdx, mpi_dfdy, mpi_dfdz


def collect_arrays(mpi, local_array):
    collected = mpi.collect(local_array, outlist=True)
    if mpi.rank == 0:
        return np.concatenate(collected)
    return None


def test_mpi_dfdx_matches_serial_derivative_nonperiodic():
    mpi = MPI()
    npoints = max(16, mpi.size * 4)
    xedges, xgrid = shift.cart.mpi_grid1D(1.0, npoints, mpi, origin=0.0)
    f_local = np.sin(xgrid)

    result_local = mpi_dfdx(xgrid, f_local, mpi, periodic=False)
    assert result_local.shape == xgrid.shape
    assert np.all(np.isfinite(result_local))

    x_global = collect_arrays(mpi, xgrid)
    result_global = collect_arrays(mpi, result_local)

    if mpi.rank == 0:
        expected = dfdx(x_global, np.sin(x_global), periodic=False)
        assert np.allclose(result_global, expected, atol=1e-8, rtol=1e-6)


def test_mpi_dfdx_matches_serial_derivative_periodic():
    mpi = MPI()
    npoints = max(16, mpi.size * 4)
    xedges, xgrid = shift.cart.mpi_grid1D(1.0, npoints, mpi, origin=0.0)
    f_local = np.sin(2.0 * np.pi * xgrid)

    result_local = mpi_dfdx(xgrid, f_local, mpi, periodic=True)
    assert result_local.shape == xgrid.shape
    assert np.all(np.isfinite(result_local))

    x_global = collect_arrays(mpi, xgrid)
    result_global = collect_arrays(mpi, result_local)

    if mpi.rank == 0:
        expected = dfdx(x_global, np.sin(2.0 * np.pi * x_global), periodic=True)
        assert np.allclose(result_global, expected, atol=1e-8, rtol=1e-6)


def test_mpi_dfdy_matches_serial_derivative_nonperiodic():
    mpi = MPI()
    ygrid = np.linspace(0.0, 1.0, 10)
    xgrid = np.linspace(0.0, 1.0, 5)
    _, Y = np.meshgrid(xgrid, ygrid, indexing="xy")
    f = np.sin(Y)

    result = mpi_dfdy(ygrid, f, mpi, periodic=False)
    expected = dfdy(ygrid, f, periodic=False)

    assert result.shape == expected.shape
    assert np.all(np.isfinite(result))
    assert np.allclose(result, expected, atol=1e-8, rtol=1e-6)


def test_mpi_dfdy_matches_serial_derivative_periodic():
    mpi = MPI()
    ygrid = np.linspace(0.0, 2.0 * np.pi, 10, endpoint=False)
    xgrid = np.linspace(0.0, 1.0, 5)
    _, Y = np.meshgrid(xgrid, ygrid, indexing="xy")
    f = np.sin(Y)

    result = mpi_dfdy(ygrid, f, mpi, periodic=True)
    expected = dfdy(ygrid, f, periodic=True)

    assert result.shape == expected.shape
    assert np.all(np.isfinite(result))
    assert np.allclose(result, expected, atol=1e-8, rtol=1e-6)


def test_mpi_dfdz_matches_serial_derivative_nonperiodic():
    mpi = MPI()
    zgrid = np.linspace(0.0, 1.0, 10)
    xgrid = np.linspace(0.0, 1.0, 3)
    ygrid = np.linspace(0.0, 1.0, 4)
    _, _, Z = np.meshgrid(xgrid, ygrid, zgrid, indexing="ij")
    f = np.sin(Z)

    result = mpi_dfdz(zgrid, f, mpi, periodic=False)
    expected = dfdz(zgrid, f, periodic=False)

    assert result.shape == expected.shape
    assert np.all(np.isfinite(result))
    assert np.allclose(result, expected, atol=1e-8, rtol=1e-6)


def test_mpi_dfdz_matches_serial_derivative_periodic():
    mpi = MPI()
    zgrid = np.linspace(0.0, 2.0 * np.pi, 10, endpoint=False)
    xgrid = np.linspace(0.0, 1.0, 3)
    ygrid = np.linspace(0.0, 1.0, 4)
    _, _, Z = np.meshgrid(xgrid, ygrid, zgrid, indexing="ij")
    f = np.sin(Z)

    result = mpi_dfdz(zgrid, f, mpi, periodic=True)
    expected = dfdz(zgrid, f, periodic=True)

    assert result.shape == expected.shape
    assert np.all(np.isfinite(result))
    assert np.allclose(result, expected, atol=1e-8, rtol=1e-6)
