include "grid.f90"
include "part2gridw.f90"


subroutine p2g_ngp_2d(x, y, f, boxsize, ngrid, npart, fgrid)
  ! Assigns points using the nearest grid point mass assignment.
  !
  ! Parameters
  ! ----------
  ! x : array
  !   X coordinates.
  ! y : array
  !   Y coordinates.
  ! boxsize : float
  !   Box size.
  ! ngrid : int
  !   Size of the grid to apply the mass assignment.
  ! npart : int
  !   Number of particles.
  ! f : array
  !   Value for the particles.
  !
  ! Returns
  ! -------
  ! fgrid : array
  !   Nearest grid point value assignment.

  implicit none
  integer, parameter :: dp = kind(1.d0)

  ! define variables

  integer, intent(in) :: ngrid, npart
  real(kind=dp), intent(in) :: x(npart), y(npart), f(npart), boxsize
  real(kind=dp), intent(out) :: fgrid(ngrid*ngrid)

  real(kind=dp) :: dl, xp, yp, fp
  integer :: i, ix, iy, q11

  ! assigning grid value 0. before assignment

  do i = 1, ngrid*ngrid
    fgrid(i) = 0.
  end do

  dl = boxsize / real(ngrid)

  ! nearest grid point assignment.

  do i = 1, npart

    xp = x(i)
    yp = y(i)
    fp = f(i)
    ix = int(xp / dl)
    iy = int(yp / dl)
    q11 = iy + ngrid*ix + 1
    fgrid(q11) = fgrid(q11) + fp/(dl*dl)

  end do

end subroutine p2g_ngp_2d


subroutine p2g_cic_2d_periodic(x, y, f, boxsize, ngrid, npart, fgrid)

  ! Assigns points using the cloud in cell mass assignment.
  !
  ! Parameters
  ! ----------
  ! x : array
  !   X coordinates.
  ! y : array
  !   Y coordinates.
  ! boxsize : float
  !   Box size.
  ! ngrid : int
  !   Size of the grid to apply the mass assignment.
  ! npart : int
  !   Number of particles.
  ! f : array
  !   Value for the particles.
  !
  ! Returns
  ! -------
  ! fgrid : array
  !   Cloud in cell value assignment.

  implicit none
  integer, parameter :: dp = kind(1.d0)

  ! define variables

  integer, intent(in) :: ngrid, npart
  real(kind=dp), intent(in) :: x(npart), y(npart), f(npart), boxsize
  real(kind=dp), intent(out) :: fgrid(ngrid*ngrid)

  real(kind=dp) :: xp, yp, fp, dl, xg, yg, wx1, wx2, wy1, wy2
  integer :: ix1, ix2, iy1, iy2
  integer :: q11, q12, q21, q22
  integer :: i

  ! assigning grid value 0. before assignment

  do i = 1, ngrid*ngrid
    fgrid(i) = 0.
  end do

  dl = boxsize / real(ngrid)

  ! cloud in cell assignment.

  do i = 1, npart

    xp = x(i)
    yp = y(i)
    fp = f(i)

    if (xp - dl/2. < 0.) then
      xp = xp + boxsize
    end if
    if (yp - dl/2. < 0.) then
      yp = yp + boxsize
    end if

    ix1 = int((xp - dl/2.) / dl)
    call xgrid(ix1, dl, xg)
    call weight_cic(xp, xg, dl, wx1)

    ix2 = ix1 + 1
    call xgrid(ix2, dl, xg)
    call weight_cic(xp, xg, dl, wx2)

    if (ix2 .eq. ngrid) then
      ix2 = ix2 - ngrid
    end if

    iy1 = int((yp - dl/2.) / dl)
    call xgrid(iy1, dl, yg)
    call weight_cic(yp, yg, dl, wy1)

    iy2 = iy1 + 1
    call xgrid(iy2, dl, yg)
    call weight_cic(yp, yg, dl, wy2)

    if (iy2 .eq. ngrid) then
      iy2 = iy2 - ngrid
    end if

    q11 = iy1 + ngrid*ix1 + 1
    q21 = iy2 + ngrid*ix1 + 1
    q12 = iy1 + ngrid*ix2 + 1
    q22 = iy2 + ngrid*ix2 + 1

    fgrid(q11) = fgrid(q11) + fp*wx1*wy1
    fgrid(q12) = fgrid(q12) + fp*wx1*wy2
    fgrid(q21) = fgrid(q21) + fp*wx2*wy1
    fgrid(q22) = fgrid(q22) + fp*wx2*wy2

  end do

end subroutine p2g_cic_2d_periodic


subroutine p2g_cic_2d_nonperiodic(x, y, f, boxsize, ngrid, npart, fgrid)

  ! Assigns points using the cloud in cell mass assignment.
  !
  ! Parameters
  ! ----------
  ! x : array
  !   X coordinates.
  ! y : array
  !   Y coordinates.
  ! boxsize : float
  !   Box size.
  ! ngrid : int
  !   Size of the grid to apply the mass assignment.
  ! npart : int
  !   Number of particles.
  ! f : array
  !   Value for the particles.
  !
  ! Returns
  ! -------
  ! fgrid : array
  !   Cloud in cell value assignment.

  implicit none
  integer, parameter :: dp = kind(1.d0)

  ! define variables

  integer, intent(in) :: ngrid, npart
  real(kind=dp), intent(in) :: x(npart), y(npart), f(npart), boxsize
  real(kind=dp), intent(out) :: fgrid(ngrid*ngrid)

  real(kind=dp) :: xp, yp, fp, dl, xg, yg, wx1, wx2, wy1, wy2
  integer :: ix1, ix2, iy1, iy2
  integer :: q11, q12, q21, q22
  integer :: i

  ! assigning grid value 0. before assignment

  do i = 1, ngrid*ngrid
    fgrid(i) = 0.
  end do

  dl = boxsize / real(ngrid)

  ! cloud in cell assignment.

  do i = 1, npart

    xp = x(i)
    yp = y(i)
    fp = f(i)

    if (xp - dl/2. < 0.) then
      ix1 = -1
    else if (xp > boxsize - dl/2.) then
      ix1 = ngrid - 1
    else
      ix1 = int((xp - dl/2.) / dl)
    end if

    ix2 = ix1 + 1

    call xgrid(ix1, dl, xg)
    call weight_cic(xp, xg, dl, wx1)
    call xgrid(ix2, dl, xg)
    call weight_cic(xp, xg, dl, wx2)

    if (yp - dl/2. < 0.) then
      iy1 = -1
    else if (yp > boxsize - dl/2.) then
      iy1 = ngrid - 1
    else
      iy1 = int((yp - dl/2.) / dl)
    end if

    iy2 = iy1 + 1

    call xgrid(iy1, dl, yg)
    call weight_cic(yp, yg, dl, wy1)
    call xgrid(iy2, dl, yg)
    call weight_cic(yp, yg, dl, wy2)

    if ((ix1 .GE. 0) .AND. (ix1 .LE. ngrid-1) .AND. (iy1 .GE. 0) .AND. (iy1 .LE. ngrid-1)) then
      q11 = iy1 + ngrid*ix1 + 1
      fgrid(q11) = fgrid(q11) + fp*wx1*wy1
    end if
    if ((ix1 .GE. 0) .AND. (ix1 .LE. ngrid-1) .AND. (iy2 .GE. 0) .AND. (iy2 .LE. ngrid-1)) then
      q21 = iy2 + ngrid*ix1 + 1
      fgrid(q12) = fgrid(q12) + fp*wx1*wy2
    end if
    if ((ix2 .GE. 0) .AND. (ix2 .LE. ngrid-1) .AND. (iy1 .GE. 0) .AND. (iy1 .LE. ngrid-1)) then
      q12 = iy1 + ngrid*ix2 + 1
      fgrid(q21) = fgrid(q21) + fp*wx2*wy1
    end if
    if ((ix2 .GE. 0) .AND. (ix2 .LE. ngrid-1) .AND. (iy2 .GE. 0) .AND. (iy2 .LE. ngrid-1)) then
      q22 = iy2 + ngrid*ix2 + 1
      fgrid(q22) = fgrid(q22) + fp*wx2*wy2
    end if

  end do

end subroutine p2g_cic_2d_nonperiodic


subroutine p2g_tsc_2d_periodic(x, y, f, boxsize, ngrid, npart, fgrid)

  ! Assigns points using the triangular shaped cloud mass assignment in 2D.
  !
  ! Parameters
  ! ----------
  ! x : array
  !   X coordinates.
  ! y : array
  !   Y coordinates.
  ! boxsize : float
  !   Box size.
  ! ngrid : int
  !   Size of the grid to apply the mass assignment.
  ! npart : int
  !   Number of particles.
  ! f : array
  !   Value for the particles.
  !
  ! Returns
  ! -------
  ! fgrid : array
  !   Triangular shaped cloud value assignment.

  implicit none
  integer, parameter :: dp = kind(1.d0)

  ! define variables

  integer, intent(in) :: ngrid, npart
  real(kind=dp), intent(in) :: x(npart), y(npart), f(npart), boxsize
  real(kind=dp), intent(out) :: fgrid(ngrid*ngrid)

  real(kind=dp) :: xp, yp, fp, xg, yg, dl, wx1, wx2, wx3, wy1, wy2, wy3
  integer :: i, ix1, ix2, ix3, iy1, iy2, iy3
  integer :: q11, q12, q13, q21, q22, q23, q31, q32, q33

  ! assigning grid value 0. before assignment

  do i = 1, ngrid*ngrid
    fgrid(i) = 0.
  end do

  dl = boxsize / real(ngrid)

  ! triangular shaped cloud assignment.

  do i = 1, npart
    xp = x(i)
    yp = y(i)
    fp = f(i)

    if (xp - dl/2. < 0.) then
      xp = xp + boxsize
    end if

    if (yp - dl/2. < 0.) then
      yp = yp + boxsize
    end if

    ix2 = int((xp - dl/2.) / dl)
    call xgrid(ix2, dl, xg)
    call weight_tsc(xp, xg, dl, wx2)

    iy2 = int((yp - dl/2.) / dl)
    call xgrid(iy2, dl, yg)
    call weight_tsc(yp, yg, dl, wy2)

    ix1 = ix2 - 1
    call xgrid(ix1, dl, xg)
    call weight_tsc(xp, xg, dl, wx1)

    if (ix1 .eq. -1) then
      ix1 = ix1 + ngrid
    end if

    iy1 = iy2 - 1
    call xgrid(iy1, dl, yg)
    call weight_tsc(yp, yg, dl, wy1)

    if (iy1 .eq. -1) then
      iy1 = iy1 + ngrid
    end if

    ix3 = ix2 + 1
    call xgrid(ix3, dl, xg)
    call weight_tsc(xp, xg, dl, wx3)

    if (ix3 .eq. ngrid) then
      ix3 = ix3 - ngrid
    end if

    iy3 = iy2 + 1
    call xgrid(iy3, dl, yg)
    call weight_tsc(yp, yg, dl, wy3)

    if (iy3 .eq. ngrid) then
      iy3 = iy3 - ngrid
    end if

    q11 = iy1 + ngrid*ix1 + 1
    q21 = iy2 + ngrid*ix1 + 1
    q31 = iy3 + ngrid*ix1 + 1
    q12 = iy1 + ngrid*ix2 + 1
    q22 = iy2 + ngrid*ix2 + 1
    q32 = iy3 + ngrid*ix2 + 1
    q13 = iy1 + ngrid*ix3 + 1
    q23 = iy2 + ngrid*ix3 + 1
    q33 = iy3 + ngrid*ix3 + 1

    fgrid(q11) = fgrid(q11) + wx1*wy1*fp
    fgrid(q21) = fgrid(q21) + wx1*wy2*fp
    fgrid(q31) = fgrid(q31) + wx1*wy3*fp
    fgrid(q12) = fgrid(q12) + wx2*wy1*fp
    fgrid(q22) = fgrid(q22) + wx2*wy2*fp
    fgrid(q32) = fgrid(q32) + wx2*wy3*fp
    fgrid(q13) = fgrid(q13) + wx3*wy1*fp
    fgrid(q23) = fgrid(q23) + wx3*wy2*fp
    fgrid(q33) = fgrid(q33) + wx3*wy3*fp

  end do

end subroutine p2g_tsc_2d_periodic


subroutine p2g_tsc_2d_nonperiodic(x, y, f, boxsize, ngrid, npart, fgrid)

  ! Assigns points using the triangular shaped cloud mass assignment in 2D.
  !
  ! Parameters
  ! ----------
  ! x : array
  !   X coordinates.
  ! y : array
  !   Y coordinates.
  ! boxsize : float
  !   Box size.
  ! ngrid : int
  !   Size of the grid to apply the mass assignment.
  ! npart : int
  !   Number of particles.
  ! f : array
  !   Value for the particles.
  !
  ! Returns
  ! -------
  ! fgrid : array
  !   Triangular shaped cloud value assignment.

  implicit none
  integer, parameter :: dp = kind(1.d0)

  ! define variables

  integer, intent(in) :: ngrid, npart
  real(kind=dp), intent(in) :: x(npart), y(npart), f(npart), boxsize
  real(kind=dp), intent(out) :: fgrid(ngrid*ngrid)

  real(kind=dp) :: xp, yp, fp, xg, yg, dl, wx1, wx2, wx3, wy1, wy2, wy3
  integer :: i, ix1, ix2, ix3, iy1, iy2, iy3
  integer :: q11, q12, q13, q21, q22, q23, q31, q32, q33

  ! assigning grid value 0. before assignment

  do i = 1, ngrid*ngrid
    fgrid(i) = 0.
  end do

  dl = boxsize / real(ngrid)

  ! triangular shaped cloud assignment.

  do i = 1, npart
    xp = x(i)
    yp = y(i)
    fp = f(i)

    if (xp - dl/2. < 0.) then
      ix2 = -1
    else if (xp > boxsize - dl/2.) then
      ix2 = ngrid - 1
    else
      ix2 = int((xp - dl/2.) / dl)
    end if

    ix1 = ix2 - 1
    ix3 = ix2 + 1

    call xgrid(ix1, dl, xg)
    call weight_tsc(xp, xg, dl, wx1)
    call xgrid(ix2, dl, xg)
    call weight_tsc(xp, xg, dl, wx2)
    call xgrid(ix3, dl, xg)
    call weight_tsc(xp, xg, dl, wx3)

    if (yp - dl/2. < 0.) then
      iy2 = -1
    else if (yp > boxsize - dl/2.) then
      iy2 = ngrid - 1
    else
      iy2 = int((yp - dl/2.) / dl)
    end if

    iy1 = iy2 - 1
    iy3 = iy2 + 1

    call xgrid(iy1, dl, yg)
    call weight_tsc(yp, yg, dl, wy1)
    call xgrid(iy2, dl, yg)
    call weight_tsc(yp, yg, dl, wy2)
    call xgrid(iy3, dl, yg)
    call weight_tsc(yp, yg, dl, wy3)

    if ((ix1 .GE. 0) .AND. (ix1 .LE. ngrid-1) .AND. (iy1 .GE. 0) .AND. (iy1 .LE. ngrid-1)) then
      q11 = iy1 + ngrid*ix1 + 1
      fgrid(q11) = fgrid(q11) + wx1*wy1*fp
    end if
    if ((ix1 .GE. 0) .AND. (ix1 .LE. ngrid-1) .AND. (iy2 .GE. 0) .AND. (iy2 .LE. ngrid-1)) then
      q21 = iy2 + ngrid*ix1 + 1
      fgrid(q21) = fgrid(q21) + wx1*wy2*fp
    end if
    if ((ix1 .GE. 0) .AND. (ix1 .LE. ngrid-1) .AND. (iy3 .GE. 0) .AND. (iy3 .LE. ngrid-1)) then
      q31 = iy3 + ngrid*ix1 + 1
      fgrid(q31) = fgrid(q31) + wx1*wy3*fp
    end if
    if ((ix2 .GE. 0) .AND. (ix2 .LE. ngrid-1) .AND. (iy1 .GE. 0) .AND. (iy1 .LE. ngrid-1)) then
      q12 = iy1 + ngrid*ix2 + 1
      fgrid(q12) = fgrid(q12) + wx2*wy1*fp
    end if
    if ((ix2 .GE. 0) .AND. (ix2 .LE. ngrid-1) .AND. (iy2 .GE. 0) .AND. (iy2 .LE. ngrid-1)) then
      q22 = iy2 + ngrid*ix2 + 1
      fgrid(q22) = fgrid(q22) + wx2*wy2*fp
    end if
    if ((ix2 .GE. 0) .AND. (ix2 .LE. ngrid-1) .AND. (iy3 .GE. 0) .AND. (iy3 .LE. ngrid-1)) then
      q32 = iy3 + ngrid*ix2 + 1
      fgrid(q32) = fgrid(q32) + wx2*wy3*fp
    end if
    if ((ix3 .GE. 0) .AND. (ix3 .LE. ngrid-1) .AND. (iy1 .GE. 0) .AND. (iy1 .LE. ngrid-1)) then
      q13 = iy1 + ngrid*ix3 + 1
      fgrid(q13) = fgrid(q13) + wx3*wy1*fp
    end if
    if ((ix3 .GE. 0) .AND. (ix3 .LE. ngrid-1) .AND. (iy2 .GE. 0) .AND. (iy2 .LE. ngrid-1)) then
      q23 = iy2 + ngrid*ix3 + 1
      fgrid(q23) = fgrid(q23) + wx3*wy2*fp
    end if
    if ((ix3 .GE. 0) .AND. (ix3 .LE. ngrid-1) .AND. (iy3 .GE. 0) .AND. (iy3 .LE. ngrid-1)) then
      q33 = iy3 + ngrid*ix3 + 1
      fgrid(q33) = fgrid(q33) + wx3*wy3*fp
    end if
  end do

end subroutine p2g_tsc_2d_nonperiodic
