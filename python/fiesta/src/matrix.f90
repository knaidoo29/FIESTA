
subroutine inv2by2(m, invm)

  ! Invert 2 by 2 matrix.
  !
  ! Parameters
  ! ----------
  ! m : array
  !     2 by 2 matrix.
  !
  ! Returns
  ! -------
  ! invm : array
  !     2 by 2 inverse matrix.

  implicit none
  integer, parameter :: dp = kind(1.d0)

  ! Declare variables.

  real(kind=dp), intent(in) :: m(4)
  real(kind=dp), intent(out) :: invm(4)

  real(kind=dp) :: a, b, c, d, detm

  a = m(1)
  b = m(2)
  c = m(3)
  d = m(4)

  detm = a*d - b*c

  invm(1) = d/detm
  invm(2) = -b/detm
  invm(3) = -c/detm
  invm(4) = a/detm

end subroutine inv2by2


subroutine inv3by3(m, invm)

  ! Invert 3 by 3 matrix.
  !
  ! Parameters
  ! ----------
  ! m : array
  !     3 by 3 matrix.
  !
  ! Returns
  ! -------
  ! invm : array
  !     3 by 3 inverse matrix.

  implicit none
  integer, parameter :: dp = kind(1.d0)

  ! Declare variables.

  real(kind=dp), intent(in) :: m(9)
  real(kind=dp), intent(out) :: invm(9)

  real(kind=dp) :: a, b, c, d, e, f, g, h, i, detm
  real(kind=dp) :: aa, bb, cc, dd, ee, ff, gg, hh, ii

  a = m(1)
  b = m(2)
  c = m(3)
  d = m(4)
  e = m(5)
  f = m(6)
  g = m(7)
  h = m(8)
  i = m(9)

  aa = e*i - f*h
  bb = -(d*i - f*g)
  cc = d*h - e*g
  dd = -(b*i - c*h)
  ee = a*i - c*g
  ff = -(a*h - b*g)
  gg = b*f - c*e
  hh = -(a*f - c*d)
  ii = a*e - b*d

  detM = a*aa + b*bb + c*cc

  invm(1) = aa / detm
  invm(2) = dd / detm
  invm(3) = gg / detm
  invm(4) = bb / detm
  invm(5) = ee / detm
  invm(6) = hh / detm
  invm(7) = cc / detm
  invm(8) = ff / detm
  invm(9) = ii / detm

end subroutine inv3by3
