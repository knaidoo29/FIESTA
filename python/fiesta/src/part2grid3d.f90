include "grid.f90"
include "part2gridw.f90"


subroutine p2g_ngp_3d(x, y, z, f, boxsize, ngrid, npart, fgrid)

  ! Assigns points using the nearest grid point mass assignment.
  !
  ! Parameters
  ! ----------
  ! x : array
  !   X coordinates.
  ! y : array
  !   Y coordinates.
  ! z : array
  !   Z coordinates.
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
  real(kind=dp), intent(in) :: x(npart), y(npart), z(npart), f(npart), boxsize
  real(kind=dp), intent(out) :: fgrid(ngrid*ngrid*ngrid)

  real(kind=dp) :: dl, xp, yp, zp, fp
  integer :: i, ix, iy, iz, q111

  ! assigning grid value 0. before assignment

  do i = 1, ngrid*ngrid*ngrid
    fgrid(i) = 0.
  end do

  dl = boxsize / real(ngrid)

  ! nearest grid point assignment.

  do i = 1, npart

    xp = x(i)
    yp = y(i)
    zp = z(i)
    fp = f(i)
    ix = int(xp / dl)
    iy = int(yp / dl)
    iz = int(zp / dl)
    q111 = iz + ngrid*(iy + ngrid*ix) + 1
    fgrid(q111) = fgrid(q111) + fp/(dl*dl*dl)
  end do

end subroutine p2g_ngp_3d


subroutine p2g_cic_3d_periodic(x, y, z, f, boxsize, ngrid, npart, fgrid)

  ! Assigns points using the cloud in cell mass assignment.
  !
  ! Parameters
  ! ----------
  ! x : array
  !   X coordinates.
  ! y : array
  !   Y coordinates.
  ! z : array
  !   Z coordinates
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
  real(kind=dp), intent(in) :: x(npart), y(npart), z(npart), f(npart), boxsize
  real(kind=dp), intent(out) :: fgrid(ngrid*ngrid*ngrid)

  real(kind=dp) :: xp, yp, zp, fp, dl, xg, yg, zg, wx1, wx2, wy1, wy2, wz1, wz2
  integer :: ix1, ix2, iy1, iy2, iz1, iz2
  integer :: q111, q112, q121, q122, q211, q212, q221, q222
  integer :: i

  ! assigning grid value 0. before assignment

  do i = 1, ngrid*ngrid*ngrid
    fgrid(i) = 0.
  end do

  dl = boxsize / real(ngrid)

  ! cloud in cell assignment.

  do i = 1, npart

    xp = x(i)
    yp = y(i)
    zp = z(i)
    fp = f(i)

    if (xp - dl/2. < 0.) then
      xp = xp + boxsize
    end if
    if (yp - dl/2. < 0.) then
      yp = yp + boxsize
    end if
    if (zp - dl/2. < 0.) then
      zp = zp + boxsize
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

    iz1 = int((zp - dl/2.) / dl)
    call xgrid(iz1, dl, zg)
    call weight_cic(zp, zg, dl, wz1)

    iz2 = iz1 + 1
    call xgrid(iz2, dl, zg)
    call weight_cic(zp, zg, dl, wz2)

    if (iz2 .eq. ngrid) then
      iz2 = iz2 - ngrid
    end if

    q111 = iz1 + ngrid*(iy1 + ngrid*ix1) + 1
    q112 = iz1 + ngrid*(iy1 + ngrid*ix2) + 1
    q121 = iz1 + ngrid*(iy2 + ngrid*ix1) + 1
    q122 = iz1 + ngrid*(iy2 + ngrid*ix2) + 1
    q211 = iz2 + ngrid*(iy1 + ngrid*ix1) + 1
    q212 = iz2 + ngrid*(iy1 + ngrid*ix2) + 1
    q221 = iz2 + ngrid*(iy2 + ngrid*ix1) + 1
    q222 = iz2 + ngrid*(iy2 + ngrid*ix2) + 1

    fgrid(q111) = fgrid(q111) + fp*wx1*wy1*wz1
    fgrid(q112) = fgrid(q112) + fp*wx2*wy1*wz1
    fgrid(q121) = fgrid(q121) + fp*wx1*wy2*wz1
    fgrid(q122) = fgrid(q122) + fp*wx2*wy2*wz1
    fgrid(q211) = fgrid(q211) + fp*wx1*wy1*wz2
    fgrid(q212) = fgrid(q212) + fp*wx2*wy1*wz2
    fgrid(q221) = fgrid(q221) + fp*wx1*wy2*wz2
    fgrid(q222) = fgrid(q222) + fp*wx2*wy2*wz2

  end do

end subroutine p2g_cic_3d_periodic


subroutine p2g_cic_3d_nonperiodic(x, y, z, f, boxsize, ngrid, npart, fgrid)

  ! Assigns points using the cloud in cell mass assignment.
  !
  ! Parameters
  ! ----------
  ! x : array
  !   X coordinates.
  ! y : array
  !   Y coordinates.
  ! z : array
  !   Z coordinates
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
  real(kind=dp), intent(in) :: x(npart), y(npart), z(npart), f(npart), boxsize
  real(kind=dp), intent(out) :: fgrid(ngrid*ngrid*ngrid)

  real(kind=dp) :: xp, yp, zp, fp, dl, xg, yg, zg, wx1, wx2, wy1, wy2, wz1, wz2
  integer :: ix1, ix2, iy1, iy2, iz1, iz2
  integer :: q111, q112, q121, q122, q211, q212, q221, q222
  integer :: i

  ! assigning grid value 0. before assignment

  do i = 1, ngrid*ngrid*ngrid
    fgrid(i) = 0.
  end do

  dl = boxsize / real(ngrid)

  ! cloud in cell assignment.

  do i = 1, npart

    xp = x(i)
    yp = y(i)
    zp = z(i)
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

    if (zp - dl/2. < 0.) then
      iz1 = -1
    else if (zp > boxsize - dl/2.) then
      iz1 = ngrid - 1
    else
      iz1 = int((zp - dl/2.) / dl)
    end if

    iz2 = iz1 + 1

    call xgrid(iz1, dl, zg)
    call weight_cic(zp, zg, dl, wz1)
    call xgrid(iz2, dl, zg)
    call weight_cic(zp, zg, dl, wz2)

    if ((iz1 .GE. 0) .AND. (iz1 .LE. ngrid-1)) then
      if ((ix1 .GE. 0) .AND. (ix1 .LE. ngrid-1) .AND. (iy1 .GE. 0) .AND. (iy1 .LE. ngrid-1)) then
        q111 = iz1 + ngrid*(iy1 + ngrid*ix1) + 1
        fgrid(q111) = fgrid(q111) + fp*wx1*wy1*wz1
      end if
      if ((ix2 .GE. 0) .AND. (ix2 .LE. ngrid-1) .AND. (iy1 .GE. 0) .AND. (iy1 .LE. ngrid-1)) then
        q112 = iz1 + ngrid*(iy1 + ngrid*ix2) + 1
        fgrid(q112) = fgrid(q112) + fp*wx2*wy1*wz1
      end if
      if ((ix1 .GE. 0) .AND. (ix1 .LE. ngrid-1) .AND. (iy2 .GE. 0) .AND. (iy2 .LE. ngrid-1)) then
        q121 = iz1 + ngrid*(iy2 + ngrid*ix1) + 1
        fgrid(q121) = fgrid(q121) + fp*wx1*wy2*wz1
      end if
      if ((ix2 .GE. 0) .AND. (ix2 .LE. ngrid-1) .AND. (iy2 .GE. 0) .AND. (iy2 .LE. ngrid-1)) then
        q122 = iz1 + ngrid*(iy2 + ngrid*ix2) + 1
        fgrid(q122) = fgrid(q122) + fp*wx2*wy2*wz1
      end if
    end if
    if ((iz2 .GE. 0) .AND. (iz2 .LE. ngrid-1)) then
      if ((ix1 .GE. 0) .AND. (ix1 .LE. ngrid-1) .AND. (iy1 .GE. 0) .AND. (iy1 .LE. ngrid-1)) then
        q211 = iz2 + ngrid*(iy1 + ngrid*ix1) + 1
        fgrid(q211) = fgrid(q211) + fp*wx1*wy1*wz2
      end if
      if ((ix2 .GE. 0) .AND. (ix2 .LE. ngrid-1) .AND. (iy1 .GE. 0) .AND. (iy1 .LE. ngrid-1)) then
        q212 = iz2 + ngrid*(iy1 + ngrid*ix2) + 1
        fgrid(q212) = fgrid(q212) + fp*wx2*wy1*wz2
      end if
      if ((ix1 .GE. 0) .AND. (ix1 .LE. ngrid-1) .AND. (iy2 .GE. 0) .AND. (iy2 .LE. ngrid-1)) then
        q221 = iz2 + ngrid*(iy2 + ngrid*ix1) + 1
        fgrid(q221) = fgrid(q221) + fp*wx1*wy2*wz2
      end if
      if ((ix2 .GE. 0) .AND. (ix2 .LE. ngrid-1) .AND. (iy2 .GE. 0) .AND. (iy2 .LE. ngrid-1)) then
        q222 = iz2 + ngrid*(iy2 + ngrid*ix2) + 1
        fgrid(q222) = fgrid(q222) + fp*wx2*wy2*wz2
      end if
    end if

  end do

end subroutine p2g_cic_3d_nonperiodic


subroutine p2g_tsc_3d_periodic(x, y, z, f, boxsize, ngrid, npart, fgrid)

  ! Assigns points using the triangular shaped cloud mass assignment in 2D.
  !
  ! Parameters
  ! ----------
  ! x : array
  !   X coordinates.
  ! y : array
  !   Y coordinates.
  ! z : array
  !   Z coordinates.
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
  real(kind=dp), intent(in) :: x(npart), y(npart), z(npart), f(npart), boxsize
  real(kind=dp), intent(out) :: fgrid(ngrid*ngrid*ngrid)

  real(kind=dp) :: xp, yp, zp, fp, xg, yg, zg, dl
  real(kind=dp) :: wx1, wx2, wx3, wy1, wy2, wy3, wz1, wz2, wz3
  integer :: i, ix1, ix2, ix3, iy1, iy2, iy3, iz1, iz2, iz3
  integer :: q111, q112, q113, q121, q122, q123, q131, q132, q133
  integer :: q211, q212, q213, q221, q222, q223, q231, q232, q233
  integer :: q311, q312, q313, q321, q322, q323, q331, q332, q333

  ! assigning grid value 0. before assignment

  do i = 1, ngrid*ngrid*ngrid
    fgrid(i) = 0.
  end do

  dl = boxsize / real(ngrid)

  ! triangular shaped cloud assignment.

  do i = 1, npart
    xp = x(i)
    yp = y(i)
    zp = z(i)
    fp = f(i)

    if (xp - dl/2. < 0.) then
      xp = xp + boxsize
    end if

    if (yp - dl/2. < 0.) then
      yp = yp + boxsize
    end if

    if (zp - dl/2. < 0.) then
      zp = zp + boxsize
    end if

    ix2 = int((xp - dl/2.) / dl)
    call xgrid(ix2, dl, xg)
    call weight_tsc(xp, xg, dl, wx2)

    iy2 = int((yp - dl/2.) / dl)
    call xgrid(iy2, dl, yg)
    call weight_tsc(yp, yg, dl, wy2)

    iz2 = int((zp - dl/2.) / dl)
    call xgrid(iz2, dl, zg)
    call weight_tsc(zp, zg, dl, wz2)

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

    iz1 = iz2 - 1
    call xgrid(iz1, dl, zg)
    call weight_tsc(zp, zg, dl, wz1)

    if (iz1 .eq. -1) then
      iz1 = iz1 + ngrid
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

    iz3 = iz2 + 1
    call xgrid(iz3, dl, zg)
    call weight_tsc(zp, zg, dl, wz3)

    if (iz3 .eq. ngrid) then
      iz3 = iz3 - ngrid
    end if

    q111 = iz1 + ngrid*(iy1 + ngrid*ix1) + 1
    q121 = iz1 + ngrid*(iy2 + ngrid*ix1) + 1
    q131 = iz1 + ngrid*(iy3 + ngrid*ix1) + 1
    q112 = iz1 + ngrid*(iy1 + ngrid*ix2) + 1
    q122 = iz1 + ngrid*(iy2 + ngrid*ix2) + 1
    q132 = iz1 + ngrid*(iy3 + ngrid*ix2) + 1
    q113 = iz1 + ngrid*(iy1 + ngrid*ix3) + 1
    q123 = iz1 + ngrid*(iy2 + ngrid*ix3) + 1
    q133 = iz1 + ngrid*(iy3 + ngrid*ix3) + 1

    q211 = iz2 + ngrid*(iy1 + ngrid*ix1) + 1
    q221 = iz2 + ngrid*(iy2 + ngrid*ix1) + 1
    q231 = iz2 + ngrid*(iy3 + ngrid*ix1) + 1
    q212 = iz2 + ngrid*(iy1 + ngrid*ix2) + 1
    q222 = iz2 + ngrid*(iy2 + ngrid*ix2) + 1
    q232 = iz2 + ngrid*(iy3 + ngrid*ix2) + 1
    q213 = iz2 + ngrid*(iy1 + ngrid*ix3) + 1
    q223 = iz2 + ngrid*(iy2 + ngrid*ix3) + 1
    q233 = iz2 + ngrid*(iy3 + ngrid*ix3) + 1

    q311 = iz3 + ngrid*(iy1 + ngrid*ix1) + 1
    q321 = iz3 + ngrid*(iy2 + ngrid*ix1) + 1
    q331 = iz3 + ngrid*(iy3 + ngrid*ix1) + 1
    q312 = iz3 + ngrid*(iy1 + ngrid*ix2) + 1
    q322 = iz3 + ngrid*(iy2 + ngrid*ix2) + 1
    q332 = iz3 + ngrid*(iy3 + ngrid*ix2) + 1
    q313 = iz3 + ngrid*(iy1 + ngrid*ix3) + 1
    q323 = iz3 + ngrid*(iy2 + ngrid*ix3) + 1
    q333 = iz3 + ngrid*(iy3 + ngrid*ix3) + 1

    fgrid(q111) = fgrid(q111) + wx1*wy1*wz1*fp
    fgrid(q121) = fgrid(q121) + wx1*wy2*wz1*fp
    fgrid(q131) = fgrid(q131) + wx1*wy3*wz1*fp
    fgrid(q112) = fgrid(q112) + wx2*wy1*wz1*fp
    fgrid(q122) = fgrid(q122) + wx2*wy2*wz1*fp
    fgrid(q132) = fgrid(q132) + wx2*wy3*wz1*fp
    fgrid(q113) = fgrid(q113) + wx3*wy1*wz1*fp
    fgrid(q123) = fgrid(q123) + wx3*wy2*wz1*fp
    fgrid(q133) = fgrid(q133) + wx3*wy3*wz1*fp

    fgrid(q211) = fgrid(q211) + wx1*wy1*wz2*fp
    fgrid(q221) = fgrid(q221) + wx1*wy2*wz2*fp
    fgrid(q231) = fgrid(q231) + wx1*wy3*wz2*fp
    fgrid(q212) = fgrid(q212) + wx2*wy1*wz2*fp
    fgrid(q222) = fgrid(q222) + wx2*wy2*wz2*fp
    fgrid(q232) = fgrid(q232) + wx2*wy3*wz2*fp
    fgrid(q213) = fgrid(q213) + wx3*wy1*wz2*fp
    fgrid(q223) = fgrid(q223) + wx3*wy2*wz2*fp
    fgrid(q233) = fgrid(q233) + wx3*wy3*wz2*fp

    fgrid(q311) = fgrid(q311) + wx1*wy1*wz3*fp
    fgrid(q321) = fgrid(q321) + wx1*wy2*wz3*fp
    fgrid(q331) = fgrid(q331) + wx1*wy3*wz3*fp
    fgrid(q312) = fgrid(q312) + wx2*wy1*wz3*fp
    fgrid(q322) = fgrid(q322) + wx2*wy2*wz3*fp
    fgrid(q332) = fgrid(q332) + wx2*wy3*wz3*fp
    fgrid(q313) = fgrid(q313) + wx3*wy1*wz3*fp
    fgrid(q323) = fgrid(q323) + wx3*wy2*wz3*fp
    fgrid(q333) = fgrid(q333) + wx3*wy3*wz3*fp

  end do

end subroutine p2g_tsc_3d_periodic


subroutine p2g_tsc_3d_nonperiodic(x, y, z, f, boxsize, ngrid, npart, fgrid)

  ! Assigns points using the triangular shaped cloud mass assignment in 2D.
  !
  ! Parameters
  ! ----------
  ! x : array
  !   X coordinates.
  ! y : array
  !   Y coordinates.
  ! z : array
  !   Z coordinates.
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
  real(kind=dp), intent(in) :: x(npart), y(npart), z(npart), f(npart), boxsize
  real(kind=dp), intent(out) :: fgrid(ngrid*ngrid*ngrid)

  real(kind=dp) :: xp, yp, zp, fp, xg, yg, zg, dl
  real(kind=dp) :: wx1, wx2, wx3, wy1, wy2, wy3, wz1, wz2, wz3
  integer :: i, ix1, ix2, ix3, iy1, iy2, iy3, iz1, iz2, iz3
  integer :: q111, q112, q113, q121, q122, q123, q131, q132, q133
  integer :: q211, q212, q213, q221, q222, q223, q231, q232, q233
  integer :: q311, q312, q313, q321, q322, q323, q331, q332, q333

  ! assigning grid value 0. before assignment

  do i = 1, ngrid*ngrid*ngrid
    fgrid(i) = 0.
  end do

  dl = boxsize / real(ngrid)

  ! triangular shaped cloud assignment.

  do i = 1, npart
    xp = x(i)
    yp = y(i)
    zp = z(i)
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

    if (zp - dl/2. < 0.) then
      iz2 = -1
    else if (zp > boxsize - dl/2.) then
      iz2 = ngrid - 1
    else
      iz2 = int((zp - dl/2.) / dl)
    end if

    iz1 = iz2 - 1
    iz3 = iz2 + 1

    call xgrid(iz1, dl, zg)
    call weight_tsc(zp, zg, dl, wz1)
    call xgrid(iz2, dl, zg)
    call weight_tsc(zp, zg, dl, wz2)
    call xgrid(iz3, dl, zg)
    call weight_tsc(zp, zg, dl, wz3)

    if ((iz1 .GE. 0) .AND. (iz1 .LE. ngrid-1)) then
      if ((ix1 .GE. 0) .AND. (ix1 .LE. ngrid-1) .AND. (iy1 .GE. 0) .AND. (iy1 .LE. ngrid-1)) then
        q111 = iz1 + ngrid*(iy1 + ngrid*ix1) + 1
        fgrid(q111) = fgrid(q111) + wx1*wy1*wz1*fp
      end if
      if ((ix1 .GE. 0) .AND. (ix1 .LE. ngrid-1) .AND. (iy2 .GE. 0) .AND. (iy2 .LE. ngrid-1)) then
        q121 = iz1 + ngrid*(iy2 + ngrid*ix1) + 1
        fgrid(q121) = fgrid(q121) + wx1*wy2*wz1*fp
      end if
      if ((ix1 .GE. 0) .AND. (ix1 .LE. ngrid-1) .AND. (iy3 .GE. 0) .AND. (iy3 .LE. ngrid-1)) then
        q131 = iz1 + ngrid*(iy3 + ngrid*ix1) + 1
        fgrid(q131) = fgrid(q131) + wx1*wy3*wz1*fp
      end if
      if ((ix2 .GE. 0) .AND. (ix2 .LE. ngrid-1) .AND. (iy1 .GE. 0) .AND. (iy1 .LE. ngrid-1)) then
        q112 = iz1 + ngrid*(iy1 + ngrid*ix2) + 1
        fgrid(q112) = fgrid(q112) + wx2*wy1*wz1*fp
      end if
      if ((ix2 .GE. 0) .AND. (ix2 .LE. ngrid-1) .AND. (iy2 .GE. 0) .AND. (iy2 .LE. ngrid-1)) then
        q122 = iz1 + ngrid*(iy2 + ngrid*ix2) + 1
        fgrid(q122) = fgrid(q122) + wx2*wy2*wz1*fp
      end if
      if ((ix2 .GE. 0) .AND. (ix2 .LE. ngrid-1) .AND. (iy3 .GE. 0) .AND. (iy3 .LE. ngrid-1)) then
        q132 = iz1 + ngrid*(iy3 + ngrid*ix2) + 1
        fgrid(q132) = fgrid(q132) + wx2*wy3*wz1*fp
      end if
      if ((ix3 .GE. 0) .AND. (ix3 .LE. ngrid-1) .AND. (iy1 .GE. 0) .AND. (iy1 .LE. ngrid-1)) then
        q113 = iz1 + ngrid*(iy1 + ngrid*ix3) + 1
        fgrid(q113) = fgrid(q113) + wx3*wy1*wz1*fp
      end if
      if ((ix3 .GE. 0) .AND. (ix3 .LE. ngrid-1) .AND. (iy2 .GE. 0) .AND. (iy2 .LE. ngrid-1)) then
        q123 = iz1 + ngrid*(iy2 + ngrid*ix3) + 1
        fgrid(q123) = fgrid(q123) + wx3*wy2*wz1*fp
      end if
      if ((ix3 .GE. 0) .AND. (ix3 .LE. ngrid-1) .AND. (iy3 .GE. 0) .AND. (iy3 .LE. ngrid-1)) then
        q133 = iz1 + ngrid*(iy3 + ngrid*ix3) + 1
        fgrid(q133) = fgrid(q133) + wx3*wy3*wz1*fp
      end if
    end if

    if ((iz2 .GE. 0) .AND. (iz2 .LE. ngrid-1)) then
      if ((ix1 .GE. 0) .AND. (ix1 .LE. ngrid-1) .AND. (iy1 .GE. 0) .AND. (iy1 .LE. ngrid-1)) then
        q211 = iz2 + ngrid*(iy1 + ngrid*ix1) + 1
        fgrid(q211) = fgrid(q211) + wx1*wy1*wz2*fp
      end if
      if ((ix1 .GE. 0) .AND. (ix1 .LE. ngrid-1) .AND. (iy2 .GE. 0) .AND. (iy2 .LE. ngrid-1)) then
        q221 = iz2 + ngrid*(iy2 + ngrid*ix1) + 1
        fgrid(q221) = fgrid(q221) + wx1*wy2*wz2*fp
      end if
      if ((ix1 .GE. 0) .AND. (ix1 .LE. ngrid-1) .AND. (iy3 .GE. 0) .AND. (iy3 .LE. ngrid-1)) then
        q231 = iz2 + ngrid*(iy3 + ngrid*ix1) + 1
        fgrid(q231) = fgrid(q231) + wx1*wy3*wz2*fp
      end if
      if ((ix2 .GE. 0) .AND. (ix2 .LE. ngrid-1) .AND. (iy1 .GE. 0) .AND. (iy1 .LE. ngrid-1)) then
        q212 = iz2 + ngrid*(iy1 + ngrid*ix2) + 1
        fgrid(q212) = fgrid(q212) + wx2*wy1*wz2*fp
      end if
      if ((ix2 .GE. 0) .AND. (ix2 .LE. ngrid-1) .AND. (iy2 .GE. 0) .AND. (iy2 .LE. ngrid-1)) then
        q222 = iz2 + ngrid*(iy2 + ngrid*ix2) + 1
        fgrid(q222) = fgrid(q222) + wx2*wy2*wz2*fp
      end if
      if ((ix2 .GE. 0) .AND. (ix2 .LE. ngrid-1) .AND. (iy3 .GE. 0) .AND. (iy3 .LE. ngrid-1)) then
        q232 = iz2 + ngrid*(iy3 + ngrid*ix2) + 1
        fgrid(q232) = fgrid(q232) + wx2*wy3*wz2*fp
      end if
      if ((ix3 .GE. 0) .AND. (ix3 .LE. ngrid-1) .AND. (iy1 .GE. 0) .AND. (iy1 .LE. ngrid-1)) then
        q213 = iz2 + ngrid*(iy1 + ngrid*ix3) + 1
        fgrid(q213) = fgrid(q213) + wx3*wy1*wz2*fp
      end if
      if ((ix3 .GE. 0) .AND. (ix3 .LE. ngrid-1) .AND. (iy2 .GE. 0) .AND. (iy2 .LE. ngrid-1)) then
        q223 = iz2 + ngrid*(iy2 + ngrid*ix3) + 1
        fgrid(q223) = fgrid(q223) + wx3*wy2*wz2*fp
      end if
      if ((ix3 .GE. 0) .AND. (ix3 .LE. ngrid-1) .AND. (iy3 .GE. 0) .AND. (iy3 .LE. ngrid-1)) then
        q233 = iz2 + ngrid*(iy3 + ngrid*ix3) + 1
        fgrid(q233) = fgrid(q233) + wx3*wy3*wz2*fp
      end if
    end if

    if ((iz3 .GE. 0) .AND. (iz3 .LE. ngrid-1)) then
      if ((ix1 .GE. 0) .AND. (ix1 .LE. ngrid-1) .AND. (iy1 .GE. 0) .AND. (iy1 .LE. ngrid-1)) then
        q311 = iz3 + ngrid*(iy1 + ngrid*ix1) + 1
        fgrid(q311) = fgrid(q311) + wx1*wy1*wz3*fp
      end if
      if ((ix1 .GE. 0) .AND. (ix1 .LE. ngrid-1) .AND. (iy2 .GE. 0) .AND. (iy2 .LE. ngrid-1)) then
        q321 = iz3 + ngrid*(iy2 + ngrid*ix1) + 1
        fgrid(q321) = fgrid(q321) + wx1*wy2*wz3*fp
      end if
      if ((ix1 .GE. 0) .AND. (ix1 .LE. ngrid-1) .AND. (iy3 .GE. 0) .AND. (iy3 .LE. ngrid-1)) then
        q331 = iz3 + ngrid*(iy3 + ngrid*ix1) + 1
        fgrid(q331) = fgrid(q331) + wx1*wy3*wz3*fp
      end if
      if ((ix2 .GE. 0) .AND. (ix2 .LE. ngrid-1) .AND. (iy1 .GE. 0) .AND. (iy1 .LE. ngrid-1)) then
        q312 = iz3 + ngrid*(iy1 + ngrid*ix2) + 1
        fgrid(q312) = fgrid(q312) + wx2*wy1*wz3*fp
      end if
      if ((ix2 .GE. 0) .AND. (ix2 .LE. ngrid-1) .AND. (iy2 .GE. 0) .AND. (iy2 .LE. ngrid-1)) then
        q322 = iz3 + ngrid*(iy2 + ngrid*ix2) + 1
        fgrid(q322) = fgrid(q322) + wx2*wy2*wz3*fp
      end if
      if ((ix2 .GE. 0) .AND. (ix2 .LE. ngrid-1) .AND. (iy3 .GE. 0) .AND. (iy3 .LE. ngrid-1)) then
        q332 = iz3 + ngrid*(iy3 + ngrid*ix2) + 1
        fgrid(q332) = fgrid(q332) + wx2*wy3*wz3*fp
      end if
      if ((ix3 .GE. 0) .AND. (ix3 .LE. ngrid-1) .AND. (iy1 .GE. 0) .AND. (iy1 .LE. ngrid-1)) then
        q313 = iz3 + ngrid*(iy1 + ngrid*ix3) + 1
        fgrid(q313) = fgrid(q313) + wx3*wy1*wz3*fp
      end if
      if ((ix3 .GE. 0) .AND. (ix3 .LE. ngrid-1) .AND. (iy2 .GE. 0) .AND. (iy2 .LE. ngrid-1)) then
        q323 = iz3 + ngrid*(iy2 + ngrid*ix3) + 1
        fgrid(q323) = fgrid(q323) + wx3*wy2*wz3*fp
      end if
      if ((ix3 .GE. 0) .AND. (ix3 .LE. ngrid-1) .AND. (iy3 .GE. 0) .AND. (iy3 .LE. ngrid-1)) then
        q333 = iz3 + ngrid*(iy3 + ngrid*ix3) + 1
        fgrid(q333) = fgrid(q333) + wx3*wy3*wz3*fp
      end if
    end if
  end do

end subroutine p2g_tsc_3d_nonperiodic
