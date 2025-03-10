program test_quad_bilinear_interp

  use types_mod, only : r8
  use model_mod, only : quad_bilinear_interp

  implicit none
  real(r8) :: lon, lat, expected_obs
  real(r8) :: x_corners(4), y_corners(4), p(4)
  integer :: ens_size

  ! Define the quadrilateral corners (longitude, latitude)
  x_corners = (/ 253.978338312172,        254.698325218660,        256.152459376863,        255.464594303337 /)
  y_corners = (/ 74.4775053981321,        74.2803767364551,        74.6534788353970,        74.8496146897753 /)
 


 
  ! Define the function values at the four corners
  p = (/ 1.0, 2.0, 3.0, 4.0 /)  ! Example values for testing

  ! Define the interpolation point (inside the quadrilateral)
  lon = 255.070312500000
  lat = 74.4735717773438

  ! Call the interpolation subroutine
  call quad_bilinear_interp(lon, lat, x_corners, y_corners, p, ens_size, expected_obs)

  ! Print the result
  print *, 'Interpolated value at (', lon, ',', lat, '): ', expected_obs

end program test_quad_bilinear_interp

