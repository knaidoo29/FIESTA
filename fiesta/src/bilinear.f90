include "grid.f90"


subroutine bilinear_periodic(fgrid, x, y, xbox, ybox, ngridx, ngridy, npart, f)

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
  ! xbox, ybox : float
  !   Size of the box.
  ! ngridx, ngridy : int
  !   Size of the grid alone each axes.
  ! npart : int
  !   Number of particles.
  !
  ! Returns
  ! -------
  ! f : array
  !   Interpolated field values.

  implicit none

  integer, parameter :: dp = kind(1.d0)

  ! define variables

  integer, intent(in) :: ngridx, ngridy, npart
  real(kind=dp), intent(in) :: fgrid(ngridx*ngridy), x(npart), y(npart), xbox, ybox
  real(kind=dp), intent(out) :: f(npart)

  real(kind=dp) :: dx, dy, xp, yp, xg1, xg2, yg1, yg2, f11, f12, f21, f22, f1, f2, minx
  integer :: ix1, ix2, iy1, iy2, q11, q12, q21, q22
  integer :: i

  dx = xbox / real(ngridx)
  dy = ybox / real(ngridy)

  minx = 0.

  ! bilinear interpolation.

  do i = 1, npart

    xp = x(i)
    yp = y(i)

    if (xp - dx/2. < 0.) then
      xp = xp + xbox
    end if
    if (yp - dy/2. < 0.) then
      yp = yp + ybox
    end if

    ix1 = int((xp - dx/2.) / dx)
    call xgrid(ix1, dx, minx, xg1)

    ix2 = ix1 + 1
    call xgrid(ix2, dx, minx, xg2)

    if (ix2 .eq. ngridx) then
      ix2 = ix2 - ngridx
    end if

    iy1 = int((yp - dy/2.) / dy)
    call xgrid(iy1, dy, minx, yg1)

    iy2 = iy1 + 1
    call xgrid(iy2, dy, minx, yg2)

    if (iy2 .eq. ngridy) then
      iy2 = iy2 - ngridy
    end if

    ! surround points in the grid of a single point for interpolation.

    q11 = iy1 + ngridy*ix1 + 1
    q12 = iy1 + ngridy*ix2 + 1
    q21 = iy2 + ngridy*ix1 + 1
    q22 = iy2 + ngridy*ix2 + 1

    f11 = fgrid(q11)
    f12 = fgrid(q12)
    f21 = fgrid(q21)
    f22 = fgrid(q22)

    f1 = ((xg2 - xp)/(xg2 - xg1))*f11 + ((xp - xg1)/(xg2 - xg1))*f12
    f2 = ((xg2 - xp)/(xg2 - xg1))*f21 + ((xp - xg1)/(xg2 - xg1))*f22
    f(i) = ((yg2 - yp)/(yg2 - yg1))*f1 + ((yp - yg1)/(yg2 - yg1))*f2

  end do

end subroutine bilinear_periodic


subroutine bilinear_nonperiodic(fgrid, x, y, xbox, ybox, ngridx, ngridy, npart, f)

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
  ! xbox, ybox : float
  !   Size of the box.
  ! ngridx, ngridy : int
  !   Size of the grid alone each axis.
  ! npart : int
  !   Number of particles.
  !
  ! Returns
  ! -------
  ! f : array
  !   Interpolated field values.

  implicit none
  integer, parameter :: dp = kind(1.d0)

  ! define variables

  integer, intent(in) :: ngridx, ngridy, npart
  real(kind=dp), intent(in) :: fgrid(ngridx*ngridy), x(npart), y(npart), xbox, ybox
  real(kind=dp), intent(out) :: f(npart)

  real(kind=dp) :: dx, dy, xp, yp, xg1, xg2, yg1, yg2, f11, f12, f21, f22, f1, f2, minx
  integer :: ix1, ix2, iy1, iy2, q11, q12, q21, q22
  integer :: i

  minx = 0.

  dx = xbox / real(ngridx)
  dy = ybox / real(ngridy)

  ! bilinear interpolation.

  do i = 1, npart

    xp = x(i)
    yp = y(i)

    if (xp - dx/2. < 0.) then
      ix1 = -1
      ix2 = 0
      call xgrid(ix1, dx, minx, xg1)
      call xgrid(ix2, dx, minx, xg2)
      ix1 = 0
    else if (xp > xbox - dx/2.) then
      ix1 = ngridx - 1
      ix2 = ngridx
      call xgrid(ix1, dx, minx, xg1)
      call xgrid(ix2, dx, minx, xg2)
      ix2 = ngridx - 1
    else
      ix1 = int((xp - dx/2.) / dx)
      ix2 = ix1 + 1
      call xgrid(ix1, dx, minx, xg1)
      call xgrid(ix2, dx, minx, xg2)
    end if

    if (yp - dy/2. < 0.) then
      iy1 = -1
      iy2 = 0
      call xgrid(iy1, dy, minx, yg1)
      call xgrid(iy2, dy, minx, yg2)
      iy1 = 0
    else if (yp > ybox - dy/2.) then
      iy1 = ngridy - 1
      iy2 = ngridy
      call xgrid(iy1, dy, minx, yg1)
      call xgrid(iy2, dy, minx, yg2)
      iy2 = ngridy - 1
    else
      iy1 = int((yp - dy/2.) / dy)
      iy2 = iy1 + 1
      call xgrid(iy1, dy, minx, yg1)
      call xgrid(iy2, dy, minx, yg2)
    end if

    ! surround points in the grid of a single point for interpolation.

    q11 = iy1 + ngridy*ix1 + 1
    q12 = iy1 + ngridy*ix2 + 1
    q21 = iy2 + ngridy*ix1 + 1
    q22 = iy2 + ngridy*ix2 + 1

    f11 = fgrid(q11)
    f12 = fgrid(q12)
    f21 = fgrid(q21)
    f22 = fgrid(q22)

    f1 = ((xg2 - xp)/(xg2 - xg1))*f11 + ((xp - xg1)/(xg2 - xg1))*f12
    f2 = ((xg2 - xp)/(xg2 - xg1))*f21 + ((xp - xg1)/(xg2 - xg1))*f22
    f(i) = ((yg2 - yp)/(yg2 - yg1))*f1 + ((yp - yg1)/(yg2 - yg1))*f2

  end do

end subroutine bilinear_nonperiodic


subroutine bilinear_axisperiodic(fgrid, x, y, xbox, ybox, perix, periy, ngridx, ngridy, npart, f)

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
  ! xbox, ybox : float
  !   Size of the box.
  ! ngridx, ngridy : int
  !   Size of the grid alone each axis.
  ! npart : int
  !   Number of particles.
  !
  ! Returns
  ! -------
  ! f : array
  !   Interpolated field values.

  implicit none
  integer, parameter :: dp = kind(1.d0)

  ! define variables

  integer, intent(in) :: ngridx, ngridy, npart, perix, periy
  real(kind=dp), intent(in) :: fgrid(ngridx*ngridy), x(npart), y(npart), xbox, ybox
  real(kind=dp), intent(out) :: f(npart)

  real(kind=dp) :: dx, dy, xp, yp, xg1, xg2, yg1, yg2, f11, f12, f21, f22, f1, f2, minx
  integer :: ix1, ix2, iy1, iy2, q11, q12, q21, q22
  integer :: i

  minx = 0.

  dx = xbox / real(ngridx)
  dy = ybox / real(ngridy)

  ! bilinear interpolation.

  do i = 1, npart

    xp = x(i)
    yp = y(i)

    if (perix == 1) then
      if (xp - dx/2. < 0.) then
        xp = xp + xbox
      end if
      ix1 = int((xp - dx/2.) / dx)
      call xgrid(ix1, dx, minx, xg1)

      ix2 = ix1 + 1
      call xgrid(ix2, dx, minx, xg2)

      if (ix2 .eq. ngridx) then
        ix2 = ix2 - ngridx
      end if
    else
      if (xp - dx/2. < 0.) then
        ix1 = -1
        ix2 = 0
        call xgrid(ix1, dx, minx, xg1)
        call xgrid(ix2, dx, minx, xg2)
        ix1 = 0
      else if (xp > xbox - dx/2.) then
        ix1 = ngridx - 1
        ix2 = ngridx
        call xgrid(ix1, dx, minx, xg1)
        call xgrid(ix2, dx, minx, xg2)
        ix2 = ngridx - 1
      else
        ix1 = int((xp - dx/2.) / dx)
        ix2 = ix1 + 1
        call xgrid(ix1, dx, minx, xg1)
        call xgrid(ix2, dx, minx, xg2)
      end if
    end if
    if (perix == 1) then
      if (yp - dy/2. < 0.) then
        yp = yp + ybox
      end if
      iy1 = int((yp - dy/2.) / dy)
      call xgrid(iy1, dy, minx, yg1)

      iy2 = iy1 + 1
      call xgrid(iy2, dy, minx, yg2)

      if (iy2 .eq. ngridy) then
        iy2 = iy2 - ngridy
      end if
    else
      if (yp - dy/2. < 0.) then
        iy1 = -1
        iy2 = 0
        call xgrid(iy1, dy, minx, yg1)
        call xgrid(iy2, dy, minx, yg2)
        iy1 = 0
      else if (yp > ybox - dy/2.) then
        iy1 = ngridy - 1
        iy2 = ngridy
        call xgrid(iy1, dy, minx, yg1)
        call xgrid(iy2, dy, minx, yg2)
        iy2 = ngridy - 1
      else
        iy1 = int((yp - dy/2.) / dy)
        iy2 = iy1 + 1
        call xgrid(iy1, dy, minx, yg1)
        call xgrid(iy2, dy, minx, yg2)
      end if
    end if

    ! surround points in the grid of a single point for interpolation.

    q11 = iy1 + ngridy*ix1 + 1
    q12 = iy1 + ngridy*ix2 + 1
    q21 = iy2 + ngridy*ix1 + 1
    q22 = iy2 + ngridy*ix2 + 1

    f11 = fgrid(q11)
    f12 = fgrid(q12)
    f21 = fgrid(q21)
    f22 = fgrid(q22)

    f1 = ((xg2 - xp)/(xg2 - xg1))*f11 + ((xp - xg1)/(xg2 - xg1))*f12
    f2 = ((xg2 - xp)/(xg2 - xg1))*f21 + ((xp - xg1)/(xg2 - xg1))*f22
    f(i) = ((yg2 - yp)/(yg2 - yg1))*f1 + ((yp - yg1)/(yg2 - yg1))*f2

  end do

end subroutine bilinear_axisperiodic
