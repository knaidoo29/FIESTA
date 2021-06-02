subroutine xgrid(index, dl, xg)

  implicit none

  integer, intent(in) :: index
  real, intent(in) :: dl
  real, intent(out) :: xg

  xg = dl/2. + real(index) * dl

end subroutine xgrid


subroutine trilinear(fgrid, x, y, z, boxsize, ngrid, npart, f)

  implicit none

  integer, intent(in) :: ngrid, npart
  real, intent(in) :: fgrid(ngrid*ngrid*ngrid)
  real, intent(in) :: x(npart), y(npart), z(npart), boxsize
  real, intent(out) :: f(npart)

  real :: dl, xp, yp, zp, xg1, xg2, yg1, yg2, zg1, zg2, xd, yd, zd
  real :: f111, f112, f121, f122, f211, f212, f221, f222
  real :: f11, f12, f21, f22, f1, f2
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

end subroutine trilinear
