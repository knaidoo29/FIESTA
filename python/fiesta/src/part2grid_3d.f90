include "part2grid_pix.f90"
include "part2grid_wei.f90"


subroutine part2grid_ngp_3d(x, y, z, f, boxsize, npart, ngrid, fgrid)

  implicit none
  integer, parameter :: dp = kind(1.d0)

  integer, intent(in) :: npart, ngrid
  real(kind=dp), intent(in) :: x(npart), y(npart), z(npart), f(npart), boxsize
  real(kind=dp), intent(out) :: fgrid(ngrid*ngrid*ngrid)

  integer :: i, xpix, ypix, zpix, pix
  real(kind=dp) :: wngp, dx, xp, yp, zp, fp

  dx = boxsize / real(ngrid)
  wngp = 1./(dx*dx*dx)

  do i = 1, ngrid*ngrid*ngrid
    fgrid(i) = 0.
  end do

  do i = 1, npart

    xp = x(i)
    yp = y(i)
    zp = z(i)
    fp = f(i)

    call ngp_pix(xp, dx, xpix)
    call ngp_pix(yp, dx, ypix)
    call ngp_pix(zp, dx, zpix)

    if ((xpix .GE. 0) .AND. (xpix .LT. ngrid) .AND. (ypix .GE. 0) .AND. (ypix .LT. ngrid) &
    .AND. (zpix .GE. 0) .AND. (zpix .LT. ngrid)) then

      call pix1dto3d(xpix, ypix, zpix, ngrid, pix)
      fgrid(pix+1) = fgrid(pix+1) + fp*wngp

    end if

  end do

end subroutine part2grid_ngp_3d


subroutine part2grid_cic_3d(x, y, z, f, boxsize, npart, ngrid, periodic, fgrid)

  implicit none
  integer, parameter :: dp = kind(1.d0)

  integer, intent(in) :: npart, ngrid
  logical, intent(in) :: periodic
  real(kind=dp), intent(in) :: x(npart), y(npart), z(npart), f(npart), boxsize
  real(kind=dp), intent(out) :: fgrid(ngrid*ngrid*ngrid)

  integer :: i, j1, j2, j3, xpix(2), ypix(2), zpix(2), pix
  real(kind=dp) :: wx, wy, wz, dx, xp, yp, zp, xg(2), yg(2), zg(2), fp

  dx = boxsize / real(ngrid)

  do i = 1, ngrid*ngrid*ngrid
    fgrid(i) = 0.
  end do

  do i = 1, npart

    xp = x(i)
    yp = y(i)
    zp = z(i)
    fp = f(i)

    call cic_pix(xp, dx, xpix)
    call cic_pix(yp, dx, ypix)
    call cic_pix(zp, dx, zpix)
    call xgrids(xpix, 2, dx, xg)
    call xgrids(ypix, 2, dx, yg)
    call xgrids(zpix, 2, dx, zg)

    if (periodic .EQV. .TRUE.) then
      call periodic_pix(xpix, 2, ngrid)
      call periodic_pix(ypix, 2, ngrid)
      call periodic_pix(zpix, 2, ngrid)
    end if
    
    do j1 = 1, 2
      do j2 = 1, 2
        do j3 = 1, 2
          if ((xpix(j1) .GE. 0) .AND. (xpix(j1) .LT. ngrid) .AND. (ypix(j2) .GE. 0) .AND. (ypix(j2) .LT. ngrid) &
          .AND. (zpix(j3) .GE. 0) .AND. (zpix(j3) .LT. ngrid)) then

            call pix1dto3d(xpix(j1), ypix(j2), zpix(j3), ngrid, pix)
            call weight_cic(xp, xg(j1), dx, wx)
            call weight_cic(yp, yg(j2), dx, wy)
            call weight_cic(zp, zg(j3), dx, wz)

            fgrid(pix+1) = fgrid(pix+1) + fp*wx*wy*wz

          end if
        end do
      end do
    end do

  end do

end subroutine part2grid_cic_3d


subroutine part2grid_tsc_3d(x, y, z, f, boxsize, npart, ngrid, periodic, fgrid)

  implicit none
  integer, parameter :: dp = kind(1.d0)

  integer, intent(in) :: npart, ngrid
  logical, intent(in) :: periodic
  real(kind=dp), intent(in) :: x(npart), y(npart), z(npart), f(npart), boxsize
  real(kind=dp), intent(out) :: fgrid(ngrid*ngrid*ngrid)

  integer :: i, j1, j2, j3, xpix(3), ypix(3), zpix(3), pix
  real(kind=dp) :: wx, wy, wz, dx, xp, yp, zp, xg(3), yg(3), zg(3), fp

  dx = boxsize / real(ngrid)

  do i = 1, ngrid*ngrid*ngrid
    fgrid(i) = 0.
  end do

  do i = 1, npart

    xp = x(i)
    yp = y(i)
    zp = z(i)
    fp = f(i)

    call tsc_pix(xp, dx, xpix)
    call tsc_pix(yp, dx, ypix)
    call tsc_pix(zp, dx, zpix)
    call xgrids(xpix, 3, dx, xg)
    call xgrids(ypix, 3, dx, yg)
    call xgrids(zpix, 3, dx, zg)

    if (periodic .EQV. .TRUE.) then
      call periodic_pix(xpix, 2, ngrid)
      call periodic_pix(ypix, 2, ngrid)
      call periodic_pix(zpix, 2, ngrid)
    end if

    do j1 = 1, 3
      do j2 = 1, 3
        do j3 = 1, 3
          if ((xpix(j1) .GE. 0) .AND. (xpix(j1) .LT. ngrid) .AND. (ypix(j2) .GE. 0) .AND. (ypix(j2) .LT. ngrid) &
          .AND. (zpix(j3) .GE. 0) .AND. (zpix(j3) .LT. ngrid)) then

            call pix1dto3d(xpix(j1), ypix(j2), zpix(j3), ngrid, pix)
            call weight_tsc(xp, xg(j1), dx, wx)
            call weight_tsc(yp, yg(j2), dx, wy)
            call weight_tsc(zp, zg(j3), dx, wz)

            fgrid(pix+1) = fgrid(pix+1) + fp*wx*wy*wz

          end if
        end do
      end do
    end do

  end do

end subroutine part2grid_tsc_3d
