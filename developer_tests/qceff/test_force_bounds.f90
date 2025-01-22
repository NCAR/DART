! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program test_force_bounds

use algorithm_info_mod, only : init_algorithm_info_mod, end_algorithm_info_mod, probit_dist_info
use distribution_params_mod,   only : GAMMA_DISTRIBUTION, BETA_DISTRIBUTION, NORMAL_DISTRIBUTION
use utilities_mod,      only : initialize_utilities, finalize_utilities
use obs_kind_mod,       only : QTY_AQUIFER_WATER, QTY_AMMONIUM_SULPHATE 
use types_mod,          only : r8, MISSING_R8

use test

implicit none

logical :: is_state, is_inflation
logical :: bounded_below, bounded_above
real(r8) :: lower_bound, upper_bound
integer :: dist_type


call initialize_utilities('test_table_read')

call init_algorithm_info_mod()

! QTY1
! inflation GAMMA
call probit_dist_info(QTY_AQUIFER_WATER, .false., .true., dist_type, &
                        bounded_below, bounded_above, lower_bound, upper_bound)

call plan(30)

call ok(dist_type == GAMMA_DISTRIBUTION)
call ok(bounded_below)
call ok(.not. bounded_above)
call ok(lower_bound == 0.0_r8)
call ok(upper_bound == MISSING_R8)

! state BETA
call probit_dist_info(QTY_AQUIFER_WATER, .true., .false., dist_type, &
                        bounded_below, bounded_above, lower_bound, upper_bound)

call ok(dist_type == BETA_DISTRIBUTION)
call ok(bounded_below)
call ok(bounded_above)
call ok(lower_bound == 0.0_r8)
call ok(upper_bound == 1.0_r8)            

! extended state NORMAL
call probit_dist_info(QTY_AQUIFER_WATER, .false., .false., dist_type, &
                        bounded_below, bounded_above, lower_bound, upper_bound)

call ok(dist_type == NORMAL_DISTRIBUTION)
call ok(.not. bounded_below)
call ok(.not. bounded_above)
call ok(lower_bound == MISSING_R8)
call ok(upper_bound == MISSING_R8)


! QTY2
! inflation BETA
call probit_dist_info(QTY_AMMONIUM_SULPHATE , .false., .true., dist_type, &
                        bounded_below, bounded_above, lower_bound, upper_bound)

call ok(dist_type == BETA_DISTRIBUTION)
call ok(bounded_below)
call ok(bounded_above)
call ok(lower_bound == 0.0_r8)
call ok(upper_bound == 1.0_r8)

! state GAMMA
call probit_dist_info(QTY_AMMONIUM_SULPHATE , .true., .false., dist_type, &
                        bounded_below, bounded_above, lower_bound, upper_bound)

call ok(dist_type == GAMMA_DISTRIBUTION)
call ok(bounded_below)
call ok(.not. bounded_above)
call ok(lower_bound == 0.0_r8)
call ok(upper_bound == MISSING_R8)        

! extended state BETA
call probit_dist_info(QTY_AMMONIUM_SULPHATE , .false., .false., dist_type, &
                        bounded_below, bounded_above, lower_bound, upper_bound)

call ok(dist_type == BETA_DISTRIBUTION)
call ok(bounded_below)
call ok(bounded_above)
call ok(lower_bound == 0.0_r8)
call ok(upper_bound == 1.0_r8)


call end_algorithm_info_mod()

call finalize_utilities()

end program test_force_bounds
