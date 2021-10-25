
subroutine xgrid(xpix, dx, xg)

  ! Returns the grid value along one axis.
  !
  ! Parameters
  ! ----------
  ! xpix : int
  !   Index along the axis.
  ! dx : float
  !   Cell width.
  !
  ! Returns
  ! -------
  ! xg : float
  !   Grid coordinate.

  implicit none
  integer, parameter :: dp = kind(1.d0)

  ! Declare variables.

  integer, intent(in) :: xpix
  real(kind=dp), intent(in) :: dx
  real(kind=dp), intent(out) :: xg

  ! Calculate the x value on the grid.

  xg = dx/2. + real(xpix) * dx

end subroutine xgrid

subroutine xgrids(xpix, pixlen, dx, xg)

  ! Returns the grid value along one axis.
  !
  ! Parameters
  ! ----------
  ! xpix : int
  !   Index along the axis.
  ! pixlen :: int
  !   Length of xpix.
  ! dx : float
  !   Cell width.
  !
  ! Returns
  ! -------
  ! xg : float
  !   Grid coordinate.

  implicit none
  integer, parameter :: dp = kind(1.d0)

  ! Declare variables.

  integer, intent(in) :: pixlen
  integer, intent(in) :: xpix(pixlen)
  real(kind=dp), intent(in) :: dx
  real(kind=dp), intent(out) :: xg(pixlen)

  integer :: i

  ! Calculate the x value on the grid.

  do i = 1, pixlen

    call xgrid(xpix(i), dx, xg(i))

  end do

end subroutine xgrids
