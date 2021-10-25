include "grid.f90"

subroutine trilinear_periodic(fgrid, x, y, z, boxsize, ngrid, npart, f)

  ! Trilinear interpolation of field defined on a grid.
  !
  ! Parameters
  ! ----------
  ! fgrid : array
  !   Field values on the grid.
  ! x : array
  !   X coordinates where we need interpolated values.
  ! y : array
  !   Y coordinates where we need interpolated values.
  ! z : array
  !   Z coordinates where we need interpolated values.
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
  integer, parameter :: dp = kind(1.d0)

  ! define variables

  integer, intent(in) :: ngrid, npart
  real(kind=dp), intent(in) :: fgrid(ngrid*ngrid*ngrid)
  real(kind=dp), intent(in) :: x(npart), y(npart), z(npart), boxsize
  real(kind=dp), intent(out) :: f(npart)

  real(kind=dp) :: dl, xp, yp, zp, xg1, xg2, yg1, yg2, zg1, zg2, xd, yd, zd
  real(kind=dp) :: f111, f112, f121, f122, f211, f212, f221, f222
  real(kind=dp) :: f11, f12, f21, f22, f1, f2
  integer :: q111, q112, q121, q122, q211, q212, q221, q222
  integer :: i, ix1, ix2, iy1, iy2, iz1, iz2

  dl = boxsize / real(ngrid)

  do i = 1, npart

    xp = x(i)
    yp = y(i)
    zp = z(i)

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

    iz1 = int((zp - dl/2.) / dl)
    call xgrid(iz1, dl, zg1)

    iz2 = iz1 + 1
    call xgrid(iz2, dl, zg2)

    if (iz2 .eq. ngrid) then
      iz2 = iz2 - ngrid
    end if

    ! surround points in the grid of a single point for interpolation.

    q111 = iz1 + ngrid*(iy1 + ngrid*ix1) + 1
    q112 = iz1 + ngrid*(iy1 + ngrid*ix2) + 1
    q121 = iz1 + ngrid*(iy2 + ngrid*ix1) + 1
    q122 = iz1 + ngrid*(iy2 + ngrid*ix2) + 1
    q211 = iz2 + ngrid*(iy1 + ngrid*ix1) + 1
    q212 = iz2 + ngrid*(iy1 + ngrid*ix2) + 1
    q221 = iz2 + ngrid*(iy2 + ngrid*ix1) + 1
    q222 = iz2 + ngrid*(iy2 + ngrid*ix2) + 1

    f111 = fgrid(q111)
    f112 = fgrid(q112)
    f121 = fgrid(q121)
    f122 = fgrid(q122)
    f211 = fgrid(q211)
    f212 = fgrid(q212)
    f221 = fgrid(q221)
    f222 = fgrid(q222)

    xd = (xp - xg1) / (xg2 - xg1)
    yd = (yp - yg1) / (yg2 - yg1)
    zd = (zp - zg1) / (zg2 - zg1)

    f11 = f111*(1-xd) + f112*xd
    f21 = f211*(1-xd) + f212*xd
    f12 = f121*(1-xd) + f122*xd
    f22 = f221*(1-xd) + f222*xd

    f1 = f11*(1-yd) + f12*yd
    f2 = f21*(1-yd) + f22*yd

    f(i) = f1*(1-zd) + f2*zd

  end do

end subroutine trilinear_periodic


subroutine trilinear_nonperiodic(fgrid, x, y, z, boxsize, ngrid, npart, f)

  ! Trilinear interpolation of field defined on a grid.
  !
  ! Parameters
  ! ----------
  ! fgrid : array
  !   Field values on the grid.
  ! x : array
  !   X coordinates where we need interpolated values.
  ! y : array
  !   Y coordinates where we need interpolated values.
  ! z : array
  !   Z coordinates where we need interpolated values.
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
  integer, parameter :: dp = kind(1.d0)

  ! define variables

  integer, intent(in) :: ngrid, npart
  real(kind=dp), intent(in) :: fgrid(ngrid*ngrid*ngrid)
  real(kind=dp), intent(in) :: x(npart), y(npart), z(npart), boxsize
  real(kind=dp), intent(out) :: f(npart)

  real(kind=dp) :: dl, xp, yp, zp, xg1, xg2, yg1, yg2, zg1, zg2, xd, yd, zd
  real(kind=dp) :: f111, f112, f121, f122, f211, f212, f221, f222
  real(kind=dp) :: f11, f12, f21, f22, f1, f2
  integer :: q111, q112, q121, q122, q211, q212, q221, q222
  integer :: i, ix1, ix2, iy1, iy2, iz1, iz2

  dl = boxsize / real(ngrid)

  do i = 1, npart

    xp = x(i)
    yp = y(i)
    zp = z(i)

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

    if (zp - dl/2. < 0.) then
      iz1 = -1
      iz2 = 0
      call xgrid(iz1, dl, zg1)
      call xgrid(iz2, dl, zg2)
      iz1 = 0
    else if (zp > boxsize - dl/2.) then
      iz1 = ngrid - 1
      iz2 = ngrid
      call xgrid(iz1, dl, zg1)
      call xgrid(iz2, dl, zg2)
      iz2 = ngrid - 1
    else
      iz1 = int((zp - dl/2.) / dl)
      iz2 = iz1 + 1
      call xgrid(iz1, dl, zg1)
      call xgrid(iz2, dl, zg2)
    end if

    ! surround points in the grid of a single point for interpolation.

    q111 = iz1 + ngrid*(iy1 + ngrid*ix1) + 1
    q112 = iz1 + ngrid*(iy1 + ngrid*ix2) + 1
    q121 = iz1 + ngrid*(iy2 + ngrid*ix1) + 1
    q122 = iz1 + ngrid*(iy2 + ngrid*ix2) + 1
    q211 = iz2 + ngrid*(iy1 + ngrid*ix1) + 1
    q212 = iz2 + ngrid*(iy1 + ngrid*ix2) + 1
    q221 = iz2 + ngrid*(iy2 + ngrid*ix1) + 1
    q222 = iz2 + ngrid*(iy2 + ngrid*ix2) + 1

    f111 = fgrid(q111)
    f112 = fgrid(q112)
    f121 = fgrid(q121)
    f122 = fgrid(q122)
    f211 = fgrid(q211)
    f212 = fgrid(q212)
    f221 = fgrid(q221)
    f222 = fgrid(q222)

    xd = (xp - xg1) / (xg2 - xg1)
    yd = (yp - yg1) / (yg2 - yg1)
    zd = (zp - zg1) / (zg2 - zg1)

    f11 = f111*(1-xd) + f112*xd
    f21 = f211*(1-xd) + f212*xd
    f12 = f121*(1-xd) + f122*xd
    f22 = f221*(1-xd) + f222*xd

    f1 = f11*(1-yd) + f12*yd
    f2 = f21*(1-yd) + f22*yd

    f(i) = f1*(1-zd) + f2*zd

  end do

end subroutine trilinear_nonperiodic
