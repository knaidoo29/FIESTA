include "grid.f90"


subroutine bilinear_periodic(fgrid, x, y, boxsize, ngrid, npart, f)

  ! Bilinear interpolation of field defined on a grid.
  !
  ! Parameters
  ! ----------
  ! fgrid : array
  !   Field values on the grid.
  ! x : array
  !   X coordinates where we need interpolated values.
  ! y : array
  !   Y coordinates where we need interpolated values.
  ! boxsize : float
  !   Size of the box.
  ! ngrid : int
  !   Size of the grid alone one axis.
  ! npart : int
  !   Number of particles.
  !
  ! Returns
  ! -------
  ! f : array
  !   Interpolated field values.

  implicit none

  ! define variables

  integer, intent(in) :: ngrid, npart
  real, intent(in) :: fgrid(ngrid*ngrid), x(npart), y(npart), boxsize
  real, intent(out) :: f(npart)

  real :: dl, xp, yp, xg1, xg2, yg1, yg2, f11, f12, f21, f22, f1, f2
  integer :: ix1, ix2, iy1, iy2, q11, q12, q21, q22
  integer :: i

  dl = boxsize / real(ngrid)

  ! bilinear interpolation.

  do i = 1, npart

    xp = x(i)
    yp = y(i)

    if (xp - dl/2. < 0.) then
      xp = xp + boxsize
    end if
    if (yp - dl/2. < 0.) then
      yp = yp + boxsize
    end if

    ix1 = int((xp - dl/2.) / dl)
    call xgrid(ix1, dl, xg1)

    ix2 = ix1 + 1
    call xgrid(ix2, dl, xg2)

    if (ix2 .eq. ngrid) then
      ix2 = ix2 - ngrid
    end if

    iy1 = int((yp - dl/2.) / dl)
    call xgrid(iy1, dl, yg1)

    iy2 = iy1 + 1
    call xgrid(iy2, dl, yg2)

    if (iy2 .eq. ngrid) then
      iy2 = iy2 - ngrid
    end if

    ! surround points in the grid of a single point for interpolation.

    q11 = iy1 + ngrid*ix1 + 1
    q12 = iy1 + ngrid*ix2 + 1
    q21 = iy2 + ngrid*ix1 + 1
    q22 = iy2 + ngrid*ix2 + 1

    f11 = fgrid(q11)
    f12 = fgrid(q12)
    f21 = fgrid(q21)
    f22 = fgrid(q22)

    f1 = ((xg2 - xp)/(xg2 - xg1))*f11 + ((xp - xg1)/(xg2 - xg1))*f12
    f2 = ((xg2 - xp)/(xg2 - xg1))*f21 + ((xp - xg1)/(xg2 - xg1))*f22
    f(i) = ((yg2 - yp)/(yg2 - yg1))*f1 + ((yp - yg1)/(yg2 - yg1))*f2

  end do

end subroutine bilinear_periodic


subroutine bilinear_nonperiodic(fgrid, x, y, boxsize, ngrid, npart, f)

  ! Bilinear interpolation of field defined on a grid.
  !
  ! Parameters
  ! ----------
  ! fgrid : array
  !   Field values on the grid.
  ! x : array
  !   X coordinates where we need interpolated values.
  ! y : array
  !   Y coordinates where we need interpolated values.
  ! boxsize : float
  !   Size of the box.
  ! ngrid : int
  !   Size of the grid alone one axis.
  ! npart : int
  !   Number of particles.
  !
  ! Returns
  ! -------
  ! f : array
  !   Interpolated field values.

  implicit none

  ! define variables

  integer, intent(in) :: ngrid, npart
  real, intent(in) :: fgrid(ngrid*ngrid), x(npart), y(npart), boxsize
  real, intent(out) :: f(npart)

  real :: dl, xp, yp, xg1, xg2, yg1, yg2, f11, f12, f21, f22, f1, f2
  integer :: ix1, ix2, iy1, iy2, q11, q12, q21, q22
  integer :: i

  dl = boxsize / real(ngrid)

  ! bilinear interpolation.

  do i = 1, npart

    xp = x(i)
    yp = y(i)

    if (xp - dl/2. < 0.) then
      ix1 = -1
      ix2 = 0
      call xgrid(ix1, dl, xg1)
      call xgrid(ix2, dl, xg2)
      ix1 = 0
    else if (xp > boxsize - dl/2.) then
      ix1 = ngrid - 1
      ix2 = ngrid
      call xgrid(ix1, dl, xg1)
      call xgrid(ix2, dl, xg2)
      ix2 = ngrid - 1
    else
      ix1 = int((xp - dl/2.) / dl)
      ix2 = ix1 + 1
      call xgrid(ix1, dl, xg1)
      call xgrid(ix2, dl, xg2)
    end if

    if (yp - dl/2. < 0.) then
      iy1 = -1
      iy2 = 0
      call xgrid(iy1, dl, yg1)
      call xgrid(iy2, dl, yg2)
      iy1 = 0
    else if (yp > boxsize - dl/2.) then
      iy1 = ngrid - 1
      iy2 = ngrid
      call xgrid(iy1, dl, yg1)
      call xgrid(iy2, dl, yg2)
      iy2 = ngrid - 1
    else
      iy1 = int((yp - dl/2.) / dl)
      iy2 = iy1 + 1
      call xgrid(iy1, dl, yg1)
      call xgrid(iy2, dl, yg2)
    end if

    ! surround points in the grid of a single point for interpolation.

    q11 = iy1 + ngrid*ix1 + 1
    q12 = iy1 + ngrid*ix2 + 1
    q21 = iy2 + ngrid*ix1 + 1
    q22 = iy2 + ngrid*ix2 + 1

    f11 = fgrid(q11)
    f12 = fgrid(q12)
    f21 = fgrid(q21)
    f22 = fgrid(q22)

    f1 = ((xg2 - xp)/(xg2 - xg1))*f11 + ((xp - xg1)/(xg2 - xg1))*f12
    f2 = ((xg2 - xp)/(xg2 - xg1))*f21 + ((xp - xg1)/(xg2 - xg1))*f22
    f(i) = ((yg2 - yp)/(yg2 - yg1))*f1 + ((yp - yg1)/(yg2 - yg1))*f2

  end do

end subroutine bilinear_nonperiodic
