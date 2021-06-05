include "polyhedron.f90"


subroutine voronoi_3d_volume(xpoints, ypoints, zpoints, xverts, yverts, zverts, &
  ridge_point1, ridge_point2, ridge_vertices, ridge_start, ridge_end, npoints, &
  nridge, nvertices, nridge_vertices, volume)

  ! Determines the volume of the voronoi cells.
  !
  ! Parameters
  ! ----------
  ! xpoints : array
  !   X-coordinates for the points.
  ! ypoints : array
  !   Y-coordinates for the points.
  ! zpoints : array
  !   Z-coordinates for the points.
  ! xverts : array
  !   X-coordinates for the voronoi vertices.
  ! yverts : array
  !   Y-coordinates for the voronoi vertices.
  ! zverts : array
  !   Z-coordinates for the voronoi vertices.
  ! ridge_points1 : array
  !   The points on either side of a voronoi ridge.
  ! ridge_points2 : array
  !   The points on either side of a voronoi ridge.
  ! ridge_vertices : array
  !   The index of the vertices on the ridge.
  ! ridge_start : array
  !   Index in ridge vertex which specifies the start of vertices that make up a single ridge.
  ! ridge_end : array
  !   Index in ridge vertex which specifies the end of vertices that make up a single ridge.
  ! npoints : int
  !   Number of points.
  ! nridge : int
  !   Number of ridges.
  ! nvertices : int
  !   Number of vertices.
  ! nridge_vertices : int
  !   Number of ridge vertices.
  !
  ! Returns
  ! -------
  ! volume : array
  !   The volume of each voronoi cell.

  implicit none
  integer, parameter :: dp = kind(1.d0)

  ! define variables

  integer, intent(in) :: npoints, nridge, nvertices, nridge_vertices
  real(kind=dp), intent(in) :: xpoints(npoints), ypoints(npoints), zpoints(npoints)
  real(kind=dp), intent(in) :: xverts(nvertices), yverts(nvertices), zverts(nvertices)
  integer, intent(in) :: ridge_point1(nridge), ridge_point2(nridge)
  integer, intent(in) :: ridge_vertices(nridge_vertices)
  integer, intent(in) :: ridge_start(nridge), ridge_end(nridge)
  real(kind=dp), intent(out) :: volume(npoints)

  integer :: i, j, check
  real(kind=dp) :: xa1, ya1, za1, xa2, ya2, za2
  real(kind=dp) :: xb, yb, zb, xc, yc, zc, xd, yd, zd, vol1, vol2

  ! computes the volume of each voronoi cell by breaking it into tetrahedron between
  ! points and ridges.

  do i = 1, npoints
    volume(i) = 0.
  end do

  do i = 1, nridge

    check = 1

    do j = ridge_start(i)+1, ridge_end(i)

      if (ridge_vertices(j) .EQ. -1) then
        check = 0
      end if

      if (ridge_vertices(j+1) .EQ. -1) then
        check = 0
      end if

    end do

    if ((check .EQ. 1)) then

      xa1 = xpoints(ridge_point1(i)+1)
      ya1 = ypoints(ridge_point1(i)+1)
      za1 = zpoints(ridge_point1(i)+1)

      xa2 = xpoints(ridge_point2(i)+1)
      ya2 = ypoints(ridge_point2(i)+1)
      za2 = zpoints(ridge_point2(i)+1)

      xb = xverts(ridge_vertices(ridge_start(i)+1)+1)
      yb = yverts(ridge_vertices(ridge_start(i)+1)+1)
      zb = zverts(ridge_vertices(ridge_start(i)+1)+1)

      do j = ridge_start(i), ridge_end(i)-2

          xc = xverts(ridge_vertices(j+2)+1)
          yc = yverts(ridge_vertices(j+2)+1)
          zc = zverts(ridge_vertices(j+2)+1)

          xd = xverts(ridge_vertices(j+3)+1)
          yd = yverts(ridge_vertices(j+3)+1)
          zd = zverts(ridge_vertices(j+3)+1)

          call tetrahedron_volume(xa1, ya1, za1, xb, yb, zb, xc, yc, zc, xd, yd, zd, vol1)
          call tetrahedron_volume(xa2, ya2, za2, xb, yb, zb, xc, yc, zc, xd, yd, zd, vol2)

          if (volume(ridge_point1(i)+1) .NE. -1.) then
            volume(ridge_point1(i)+1) = volume(ridge_point1(i)+1) + vol1
          end if

          if (volume(ridge_point2(i)+1) .NE. -1.) then
            volume(ridge_point2(i)+1) = volume(ridge_point2(i)+1) + vol2
          end if

      end do

    else

      volume(ridge_point1(i)+1) = -1.
      volume(ridge_point2(i)+1) = -1.

    end if

  end do

end subroutine voronoi_3d_volume
