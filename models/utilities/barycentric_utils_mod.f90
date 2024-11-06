!> Barycentric interpolation routines for triangular structures within grids. These were originally
!> Developed for the mpas_atm model_mod. They're being separated into their own utilities module
!> so they can be used by multiple model_mods.

module barycentric_utils_mod

   use        types_mod, only : r8, DEG2RAD, radius => earth_radius
   
   implicit none
   private
   
   public :: inside_triangle,               &
             get_3d_weights,                &
             get_barycentric_weights,       &
             latlon_to_xyz,                 &
             barycentric_average
   
   real(r8), parameter :: roundoff = 1.0e-12_r8
   
   contains
   
   subroutine inside_triangle(t1, t2, t3, r, lat, lon, inside, weights)
   
      ! given 3 corners of a triangle and an xyz point, compute whether
      ! the point is inside the triangle.  this assumes r is coplanar
      ! with the triangle - the caller must have done the lat/lon to
      ! xyz conversion with a constant radius and then this will be
      ! true (enough).  sets inside to true/false, and returns the
      ! weights if true.  weights are set to 0 if false.
      
      real(r8), intent(in)  :: t1(3), t2(3), t3(3)
      real(r8), intent(in)  :: r(3), lat, lon
      logical,  intent(out) :: inside
      real(r8), intent(out) :: weights(3)
      
      ! check for degenerate cases first - is the test point located
      ! directly on one of the vertices?  (this case may be common
      ! if we're computing on grid point locations.)
      if (all(abs(r - t1) < roundoff)) then
         inside = .true.
         weights = (/ 1.0_r8, 0.0_r8, 0.0_r8 /)
         return
      else if (all(abs(r - t2) < roundoff)) then
         inside = .true.
         weights = (/ 0.0_r8, 1.0_r8, 0.0_r8 /)
         return
      else if (all(abs(r - t3) < roundoff)) then
         inside = .true.
         weights = (/ 0.0_r8, 0.0_r8, 1.0_r8 /)
         return
      endif
      
      ! not a vertex. compute the weights.  if any are
      ! negative, the point is outside.  since these are
      ! real valued computations define a lower bound for
      ! numerical roundoff error and be sure that the
      ! weights are not just *slightly* negative.
      call get_3d_weights(r, t1, t2, t3, lat, lon, weights)
      
      if (any(weights < -roundoff)) then
         inside = .false.
         weights = 0.0_r8
         return
      endif
      
      ! truncate barely negative values to 0
      inside = .true.
      where (weights < 0.0_r8) weights = 0.0_r8
      return
      
   end subroutine inside_triangle
   
   subroutine get_3d_weights(p, v1, v2, v3, lat, lon, weights)
   
      ! Given a point p (x,y,z) inside a triangle, and the (x,y,z)
      ! coordinates of the triangle corner points (v1, v2, v3),
      ! find the weights for a barycentric interpolation.  this
      ! computation only needs two of the three coordinates, so figure
      ! out which quadrant of the sphere the triangle is in and pick
      ! the 2 axes which are the least planar:
      !  (x,y) near the poles,
      !  (y,z) near 0 and 180 longitudes near the equator,
      !  (x,z) near 90 and 270 longitude near the equator.
      ! (lat/lon are the coords of p. we could compute them here
      ! but since in all cases we already have them, pass them
      ! down for efficiency)
      
      real(r8), intent(in)  :: p(3)
      real(r8), intent(in)  :: v1(3), v2(3), v3(3)
      real(r8), intent(in)  :: lat, lon
      real(r8), intent(out) :: weights(3)
      
      real(r8) :: cxs(3), cys(3)
      
      ! above or below 45 in latitude, where -90 < lat < 90:
      if (lat >= 45.0_r8 .or. lat <= -45.0_r8) then
         cxs(1) = v1(1)
         cxs(2) = v2(1)
         cxs(3) = v3(1)
         cys(1) = v1(2)
         cys(2) = v2(2)
         cys(3) = v3(2)
         call get_barycentric_weights(p(1), p(2), cxs, cys, weights)
         return
      endif
      
      ! nearest 0 or 180 in longitude, where 0 < lon < 360:
      if ( lon <= 45.0_r8 .or. lon >= 315.0_r8 .or. &
          (lon >= 135.0_r8 .and. lon <= 225.0_r8)) then
         cxs(1) = v1(2)
         cxs(2) = v2(2)
         cxs(3) = v3(2)
         cys(1) = v1(3)
         cys(2) = v2(3)
         cys(3) = v3(3)
         call get_barycentric_weights(p(2), p(3), cxs, cys, weights)
         return
      endif
      
      ! last option, nearest 90 or 270 in lon:
      cxs(1) = v1(1)
      cxs(2) = v2(1)
      cxs(3) = v3(1)
      cys(1) = v1(3)
      cys(2) = v2(3)
      cys(3) = v3(3)
      call get_barycentric_weights(p(1), p(3), cxs, cys, weights)
      
   end subroutine get_3d_weights
   
   subroutine get_barycentric_weights(x, y, cxs, cys, weights)
   
      ! Computes the barycentric weights for a 2d interpolation point
      ! (x,y) in a 2d triangle with the given (cxs,cys) corners.
      
      real(r8), intent(in)  :: x, y, cxs(3), cys(3)
      real(r8), intent(out) :: weights(3)
      
      real(r8) :: denom
      
      ! Get denominator
      denom = (cys(2) - cys(3)) * (cxs(1) - cxs(3)) + &
         (cxs(3) - cxs(2)) * (cys(1) - cys(3))
      
      weights(1) = ((cys(2) - cys(3)) * (x - cxs(3)) + &
         (cxs(3) - cxs(2)) * (y - cys(3))) / denom
      
      weights(2) = ((cys(3) - cys(1)) * (x - cxs(3)) + &
         (cxs(1) - cxs(3)) * (y - cys(3))) / denom
      
      weights(3) = 1.0_r8 - weights(1) - weights(2)
      
      if (any(abs(weights) < roundoff)) then
         where (abs(weights) < roundoff) weights = 0.0_r8
         where (abs(1.0_r8 - abs(weights)) < roundoff) weights = 1.0_r8
      endif
      
   end subroutine get_barycentric_weights
   
   subroutine latlon_to_xyz(lat, lon, x, y, z)
   
      ! Given a lat, lon in degrees, return the cartesian x,y,z coordinate
      ! on the surface of a specified radius relative to the origin
      ! at the center of the earth.
      
      real(r8), intent(in)  :: lat, lon
      real(r8), intent(out) :: x, y, z
      
      real(r8) :: rlat, rlon
      
      rlat = lat * DEG2RAD
      rlon = lon * DEG2RAD
      
      x = radius * cos(rlon) * cos(rlat)
      y = radius * sin(rlon) * cos(rlat)
      z = radius * sin(rlat)
      
   end subroutine latlon_to_xyz
   
   function barycentric_average(nitems, weights, vertex_temp_values) result (averaged_values)
   
      integer, intent(in)                                        :: nitems
      real(r8), dimension(3), intent(in)                         :: weights
      real(r8), dimension(3, nitems), intent(in)                 :: vertex_temp_values
   
      real(r8), dimension(nitems)                                :: averaged_values
   
      integer :: iweight, iitem
   
      averaged_values(:) = 0
   
      do iitem = 1, nitems
         do iweight = 1, 3
            averaged_values(iitem) = averaged_values(iitem) + weights(iweight)*vertex_temp_values(iweight, iitem)
         end do
      end do
   
   end function barycentric_average
   
end module barycentric_utils_mod
   