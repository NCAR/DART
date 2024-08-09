program test_get_close_zero_dist

use location_mod, only : get_close_init, get_close_type, location_type
use types_mod,    only : r8
use utilities_mod, only : initialize_utilities, finalize_utilities

implicit none

integer, parameter :: n = 10
type(get_close_type) :: gc
type(location_type)  ::  locs(n)
real(r8) :: max_dist = 0.0_r8

call initialize_utilities('test_get_close')

call get_close_init(gc, n, max_dist, locs)

call finalize_utilities()

end program test_get_close_zero_dist
