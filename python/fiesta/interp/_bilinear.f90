subroutine xgrid(index, dl, xg)

  implicit none

  integer, intent(in) :: index
  real, intent(in) :: dl
  real, intent(out) :: xg

  xg = dl/2. + real(index) * dl

end subroutine xgrid


subroutine bilinear(fgrid, x, y, boxsize, ngrid, npart, f)

  implicit none

  integer, intent(in) :: ngrid, npart
  real, intent(in) :: fgrid(ngrid*ngrid), x(npart), y(npart), boxsize
  real, intent(out) :: f(npart)

  real :: dl, xp, yp, xg1, xg2, yg1, yg2, f11, f12, f21, f22, f1, f2
  integer :: ix1, ix2, iy1, iy2, q11, q12, q21, q22
  integer :: i

  dl = boxsize / real(ngrid)

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

end subroutine bilinear
