subroutine xgrid(index, dl, xg)

  implicit none

  integer, intent(in) :: index
  real, intent(in) :: dl
  real, intent(out) :: xg

  xg = dl/2. + real(index) * dl

end subroutine xgrid

subroutine weight_cic(xp, xg, dl, w)

  implicit none

  real, intent(in) :: xp, xg, dl
  real, intent(out) :: w

  real :: dx

  dx = abs(xp - xg)

  if (dx .LE. dl) then
    w = (1. - (dx/dl))/dl
  else
    w = 0.
  end if

end subroutine weight_cic

subroutine weight_tsc(xp, xg, dl, w)

  implicit none

  real, intent(in) :: xp, xg, dl
  real, intent(out) :: w

  real :: dx

  dx = abs(xp - xg)

  if (dx .LE. dl/2.) then
    w = 3./4. - (dx/dl)**2.
  else if (dx .LE. 3.*dl/2.) then
    w = (1./2.)*((3./2. - (dx/dl))**2.)
  else
    w = 0.
  end if

  w = w / dl

end subroutine weight_tsc

subroutine part2grid_ngp(x, y, z, f, boxsize, ngrid, npart, fgrid)
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

  ! DEFINE INPUT & OUTPUT VARIABLES

  integer, intent(in) :: ngrid, npart
  real, intent(in) :: x(npart), y(npart), z(npart), f(npart), boxsize
  real, intent(out) :: fgrid(ngrid*ngrid*ngrid)

  ! INTERNAL VARIABLES

  real :: dl, xp, yp, zp, fp
  integer :: i, ix, iy, iz, q111

  ! MAIN BODY

  do i = 1, ngrid*ngrid*ngrid
    fgrid(i) = 0.
  end do

  dl = boxsize / real(ngrid)

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

end subroutine part2grid_ngp


subroutine part2grid_cic(x, y, z, f, boxsize, ngrid, npart, fgrid)
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

  ! DEFINE INPUT & OUTPUT VARIABLES

  integer, intent(in) :: ngrid, npart
  real, intent(in) :: x(npart), y(npart), z(npart), f(npart), boxsize
  real, intent(out) :: fgrid(ngrid*ngrid*ngrid)

  ! INTERNAL VARIABLES

  real :: xp, yp, zp, fp, dl, xg, yg, zg, wx1, wx2, wy1, wy2, wz1, wz2
  integer :: ix1, ix2, iy1, iy2, iz1, iz2
  integer :: q111, q112, q121, q122, q211, q212, q221, q222
  integer :: i

  ! MAIN BODY

  do i = 1, ngrid*ngrid*ngrid
    fgrid(i) = 0.
  end do

  dl = boxsize / real(ngrid)

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

end subroutine part2grid_cic


subroutine part2grid_tsc(x, y, z, f, boxsize, ngrid, npart, fgrid)
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

  ! DEFINE INPUT & OUTPUT VARIABLES

  integer, intent(in) :: ngrid, npart
  real, intent(in) :: x(npart), y(npart), z(npart), f(npart), boxsize
  real, intent(out) :: fgrid(ngrid*ngrid*ngrid)

  ! INTERNAL VARIABLES

  real :: xp, yp, zp, fp, xg, yg, zg, dl
  real :: wx1, wx2, wx3, wy1, wy2, wy3, wz1, wz2, wz3
  integer :: i, ix1, ix2, ix3, iy1, iy2, iy3, iz1, iz2, iz3
  integer :: q111, q112, q113, q121, q122, q123, q131, q132, q133
  integer :: q211, q212, q213, q221, q222, q223, q231, q232, q233
  integer :: q311, q312, q313, q321, q322, q323, q331, q332, q333

  ! MAIN BODY

  do i = 1, ngrid*ngrid*ngrid
    fgrid(i) = 0.
  end do

  dl = boxsize / real(ngrid)

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

end subroutine part2grid_tsc


subroutine part2grid_ngp_2d(x, y, f, boxsize, ngrid, npart, fgrid)
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

  ! DEFINE INPUT & OUTPUT VARIABLES

  integer, intent(in) :: ngrid, npart
  real, intent(in) :: x(npart), y(npart), f(npart), boxsize
  real, intent(out) :: fgrid(ngrid*ngrid)

  ! INTERNAL VARIABLES

  real :: dl, xp, yp, fp
  integer :: i, ix, iy, q11

  ! MAIN BODY

  do i = 1, ngrid*ngrid
    fgrid(i) = 0.
  end do

  dl = boxsize / real(ngrid)

  do i = 1, npart
    xp = x(i)
    yp = y(i)
    fp = f(i)
    ix = int(xp / dl)
    iy = int(yp / dl)
    q11 = iy + ngrid*ix + 1
    fgrid(q11) = fgrid(q11) + fp/(dl*dl)
  end do

end subroutine part2grid_ngp_2d


subroutine part2grid_cic_2d(x, y, f, boxsize, ngrid, npart, fgrid)
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

  ! DEFINE INPUT & OUTPUT VARIABLES

  integer, intent(in) :: ngrid, npart
  real, intent(in) :: x(npart), y(npart), f(npart), boxsize
  real, intent(out) :: fgrid(ngrid*ngrid)

  ! INTERNAL VARIABLES

  real :: xp, yp, fp, dl, xg, yg, wx1, wx2, wy1, wy2
  integer :: ix1, ix2, iy1, iy2
  integer :: q11, q12, q21, q22
  integer :: i

  ! MAIN BODY

  do i = 1, ngrid*ngrid
    fgrid(i) = 0.
  end do

  dl = boxsize / real(ngrid)

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

end subroutine part2grid_cic_2d


subroutine part2grid_tsc_2d(x, y, f, boxsize, ngrid, npart, fgrid)
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

  ! DEFINE INPUT & OUTPUT VARIABLES

  integer, intent(in) :: ngrid, npart
  real, intent(in) :: x(npart), y(npart), f(npart), boxsize
  real, intent(out) :: fgrid(ngrid*ngrid)

  ! INTERNAL VARIABLES

  real :: xp, yp, fp, xg, yg, dl, wx1, wx2, wx3, wy1, wy2, wy3
  integer :: i, ix1, ix2, ix3, iy1, iy2, iy3
  integer :: q11, q12, q13, q21, q22, q23, q31, q32, q33

  ! MAIN BODY

  do i = 1, ngrid*ngrid
    fgrid(i) = 0.
  end do

  dl = boxsize / real(ngrid)

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

end subroutine part2grid_tsc_2d
