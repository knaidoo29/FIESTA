
subroutine xgrid(index, dl, xg)

  ! Returns the grid value along one axis.
  !
  ! Parameters
  ! ----------
  ! index : int
  !   Index along the axis.
  ! dl : float
  !   Cell width.
  !
  ! Returns
  ! -------
  ! xg : float
  !   Grid coordinate.

  implicit none

  ! Declare variables.

  integer, intent(in) :: index
  real, intent(in) :: dl
  real, intent(out) :: xg

  ! Calculate the x value on the grid.
  
  xg = dl/2. + real(index) * dl

end subroutine xgrid
