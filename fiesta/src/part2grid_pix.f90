include "grid.f90"


subroutine pix1dto2d(xpix, ypix, ngrid, pix)

  implicit none
  integer, parameter :: dp = kind(1.d0)

  integer, intent(in) :: xpix, ypix, ngrid
  integer, intent(out) :: pix

  pix = ypix + ngrid*xpix

end subroutine pix1dto2d


subroutine pix1dto3d(xpix, ypix, zpix, ngrid, pix)

  implicit none
  integer, parameter :: dp = kind(1.d0)

  integer, intent(in) :: xpix, ypix, zpix, ngrid
  integer, intent(out) :: pix

  pix = zpix + ngrid*(ypix + ngrid*xpix)

end subroutine pix1dto3d


subroutine find_pix(x, dx, pix)

  implicit none
  integer, parameter :: dp = kind(1.d0)

  real(kind=dp), intent(in) :: x, dx
  integer, intent(out) :: pix

  pix = INT(FLOOR(x/dx))

end subroutine find_pix

subroutine periodic_pix(pix, pixlen, ngrid)

  implicit none
  integer, parameter :: dp = kind(1.d0)

  integer, intent(in) :: pixlen, ngrid
  integer, intent(inout) :: pix(pixlen)

  integer :: i

  do i = 1, pixlen
    if (pix(i) .LT. 0) then
      pix(i) = pix(i) + ngrid
    else if (pix(i) .GE. ngrid) then
      pix(i) = pix(i) - ngrid
    end if
  end do

end subroutine periodic_pix

subroutine ngp_pix(x, dx, pix)

  implicit none
  integer, parameter :: dp = kind(1.d0)

  real(kind=dp), intent(in) :: x, dx
  integer, intent(out) :: pix

  call find_pix(x, dx, pix)

end subroutine ngp_pix


subroutine cic_pix(x, dx, pix)

  implicit none
  integer, parameter :: dp = kind(1.d0)

  real(kind=dp), intent(in) :: x, dx
  integer, intent(out) :: pix(2)

  integer :: xpix
  real(kind=dp) :: xg

  call find_pix(x, dx, xpix)
  call xgrid(xpix, dx, xg)

  if (x .LT. xg) then
    pix(1) = xpix - 1
    pix(2) = xpix
  else
    pix(1) = xpix
    pix(2) = xpix + 1
  end if

end subroutine cic_pix


subroutine tsc_pix(x, dx, pix)

  implicit none
  integer, parameter :: dp = kind(1.d0)

  real(kind=dp), intent(in) :: x, dx
  integer, intent(out) :: pix(3)

  integer :: xpix

  call find_pix(x, dx, xpix)

  pix(1) = xpix - 1
  pix(2) = xpix
  pix(3) = xpix + 1

end subroutine tsc_pix
