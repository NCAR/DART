! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module algorithm_info_mod

use types_mod, only : r8, i8, missing_r8

use obs_def_mod, only : obs_def_type, get_obs_def_type_of_obs, get_obs_def_error_variance
use obs_kind_mod, only : get_quantity_for_type_of_obs, get_name_for_quantity

use assim_model_mod, only : get_state_meta_data
use location_mod, only    : location_type

use distribution_params_mod, only : NORMAL_DISTRIBUTION, BOUNDED_NORMAL_RH_DISTRIBUTION, &
   GAMMA_DISTRIBUTION, BETA_DISTRIBUTION, LOG_NORMAL_DISTRIBUTION, UNIFORM_DISTRIBUTION,  &
   PARTICLE_FILTER_DISTRIBUTION

implicit none
private

! Defining parameter strings for different observation space filters
! For now, retaining backwards compatibility in assim_tools_mod requires using
! these specific integer values and there is no point in using these in assim_tools.
! That will change if backwards compatibility is removed in the future.
integer, parameter :: EAKF               = 1
integer, parameter :: ENKF               = 2
integer, parameter :: UNBOUNDED_RHF      = 8
integer, parameter :: GAMMA_FILTER       = 11
integer, parameter :: BOUNDED_NORMAL_RHF = 101 

public :: obs_error_info, probit_dist_info, obs_inc_info, &
          init_algorithm_info_mod, end_algorithm_info_mod, &
          obs_error_info_type, probit_inflation_type, probit_state_type, &
          probit_extended_state_type, obs_inc_info_type, qcf_table_data_type, &
          EAKF, ENKF, BOUNDED_NORMAL_RHF, UNBOUNDED_RHF, GAMMA_FILTER

!Creates the type definitions for the QCF table
type obs_error_info_type
   logical               :: bounded_below, bounded_above
   real                  :: lower_bound,   upper_bound
end type

type probit_inflation_type
   integer               :: dist_type
   logical               :: bounded_below, bounded_above
   real                  :: lower_bound,   upper_bound
end type

type probit_state_type
   integer               :: dist_type
   logical               :: bounded_below, bounded_above
   real                  :: lower_bound,   upper_bound
end type

type probit_extended_state_type
   integer               :: dist_type
   logical               :: bounded_below, bounded_above
   real                  :: lower_bound,   upper_bound
end type

type obs_inc_info_type
   integer :: filter_kind
   logical :: rectangular_quadrature, gaussian_likelihood_tails
   logical :: sort_obs_inc, spread_restoration
   logical :: bounded_below, bounded_above
   real :: lower_bound,   upper_bound
end type

type qcf_table_data_type
   type(obs_error_info_type) :: obs_error_info
   type(probit_inflation_type) :: probit_inflation
   type(probit_state_type) :: probit_state
   type(probit_extended_state_type) :: probit_extended_state
   type(obs_inc_info_type) :: obs_inc_info
end type

type(qcf_table_data_type), allocatable :: qcf_table_data(:)
character(len=129), allocatable :: qcf_table_row_headers(:) 

! Provides routines that give information about details of algorithms for 
! observation error sampling, observation increments, and the transformations
! for regression and inflation in probit space. 
! For now, it is convenient to have these in a single module since several
! users will be developing their own problem specific versions of these
! subroutines. This will avoid constant merge conflicts as other parts of the
! assimilation code are updated.

contains

!-------------------------------------------------------------------------


subroutine obs_error_info(obs_def, error_variance, &
   bounded_below, bounded_above, lower_bound, upper_bound)

! Computes information needed to compute error sample for this observation
! This is called by perfect_model_obs when generating noisy obs
type(obs_def_type), intent(in)  :: obs_def
real(r8),           intent(out) :: error_variance
logical,            intent(out) :: bounded_below, bounded_above
real(r8),           intent(out) :: lower_bound, upper_bound

integer     :: obs_type, obs_kind
integer(i8) :: state_var_index
type(location_type) :: temp_loc

integer :: QTY_loc(1)
character(len=129) :: kind_name

! Get the kind of the observation
obs_type = get_obs_def_type_of_obs(obs_def)
! If it is negative, it is an identity obs
if(obs_type < 0) then
   state_var_index = -1 * obs_type
   call get_state_meta_data(state_var_index, temp_loc, obs_kind)
else
   obs_kind = get_quantity_for_type_of_obs(obs_type)
endif

! Get the default error variance
error_variance = get_obs_def_error_variance(obs_def)

!get actual name of QTY from integer index
kind_name = get_name_for_quantity(obs_kind)
write(*,*) 'kind_name: ', kind_name

!find location of QTY in qcf_table_data structure
QTY_loc = findloc(qcf_table_row_headers, kind_name)
write(*,*) 'findloc of kind: ', QTY_loc(1)

if (QTY_loc(1) == 0) then
   write(*,*) 'QTY not in table, using default values'

   !use default values
   bounded_below = .false.;    bounded_above = .false.
   lower_bound   = missing_r8; upper_bound   = missing_r8

   else
      bounded_below = qcf_table_data(QTY_loc(1))%obs_error_info%bounded_below
      bounded_above = qcf_table_data(QTY_loc(1))%obs_error_info%bounded_above
      lower_bound = qcf_table_data(QTY_loc(1))%obs_error_info%lower_bound
      upper_bound = qcf_table_data(QTY_loc(1))%obs_error_info%upper_bound

endif

write(*,*) 'obs_error_info: ', bounded_below, bounded_above, lower_bound, upper_bound

end subroutine obs_error_info

!-------------------------------------------------------------------------


subroutine probit_dist_info(kind, is_state, is_inflation, dist_type, &
   bounded_below, bounded_above, lower_bound, upper_bound)

! Computes the details of the probit transform for initial experiments
! with Molly 

integer,  intent(in)  :: kind
logical,  intent(in)  :: is_state      ! True for state variable, false for obs
logical,  intent(in)  :: is_inflation  ! True for inflation transform
integer,  intent(out) :: dist_type
logical,  intent(out) :: bounded_below, bounded_above
real(r8), intent(out) :: lower_bound,   upper_bound

integer :: QTY_loc(1)
character(len=129) :: kind_name

! Have input information about the kind of the state or observation being transformed
! along with additional logical info that indicates whether this is an observation
! or state variable and about whether the transformation is being done for inflation
! or for regress. 
! Need to select the appropriate transform. At present, options are NORMAL_PRIOR
! which does nothing or BOUNDED_NORMAL_RH_PRIOR. 
! If the BNRH is selected then information about the bounds must also be set.
! The two dimensional logical array 'bounded' is set to false for no bounds and true
! for bounded. the first element of the array is for the lower bound, the second for the upper.
! If bounded is chosen, the corresponding bound value(s) must be set in the two dimensional 
! real array 'bounds'.
! For example, if my_state_kind corresponds to a sea ice fraction then an appropriate choice
! would be:
! bounded_below = .true.;  bounded_above = .true.
! lower_bound  = 0.0_r8;   upper_bounds  = 1.0_r8

! In the long run, may not have to have separate controls for each of the input possibilities
! However, for now these are things that need to be explored for science understanding

!get actual name of QTY from integer index
kind_name = get_name_for_quantity(kind)
write(*,*) 'kind_name: ', kind_name

!find location of QTY in qcf_table_data structure
QTY_loc = findloc(qcf_table_row_headers, kind_name)
write(*,*) 'findloc of kind: ', QTY_loc(1)

if (QTY_loc(1) == 0) then
   write(*,*) 'QTY not in table, using default values'
   
   !use default values
   dist_type = BOUNDED_NORMAL_RH_DISTRIBUTION
   bounded_below = .false.;    bounded_above = .false.
   lower_bound   = missing_r8; upper_bound   = missing_r8

   elseif(is_inflation) then
   ! Case for inflation transformation

      dist_type = qcf_table_data(QTY_loc(1))%probit_inflation%dist_type  !dist_type has checks in transform_to_probit, transform_from_probit
      bounded_below = qcf_table_data(QTY_loc(1))%probit_inflation%bounded_below
      bounded_above = qcf_table_data(QTY_loc(1))%probit_inflation%bounded_above
      lower_bound = qcf_table_data(QTY_loc(1))%probit_inflation%lower_bound
      upper_bound = qcf_table_data(QTY_loc(1))%probit_inflation%upper_bound

   elseif(is_state) then
   ! Case for state variable priors

      dist_type = qcf_table_data(QTY_loc(1))%probit_state%dist_type
      bounded_below = qcf_table_data(QTY_loc(1))%probit_state%bounded_below
      bounded_above = qcf_table_data(QTY_loc(1))%probit_state%bounded_above
      lower_bound = qcf_table_data(QTY_loc(1))%probit_state%lower_bound
      upper_bound = qcf_table_data(QTY_loc(1))%probit_state%upper_bound

   else
   ! This case is for observation (extended state) priors

      dist_type = qcf_table_data(QTY_loc(1))%probit_extended_state%dist_type
      bounded_below = qcf_table_data(QTY_loc(1))%probit_extended_state%bounded_below
      bounded_above = qcf_table_data(QTY_loc(1))%probit_extended_state%bounded_above
      lower_bound = qcf_table_data(QTY_loc(1))%probit_extended_state%lower_bound
      upper_bound = qcf_table_data(QTY_loc(1))%probit_extended_state%upper_bound

endif

write(*,*) 'probit_dist_info: ', dist_type, bounded_below, bounded_above, lower_bound, upper_bound

end subroutine probit_dist_info

!------------------------------------------------------------------------


subroutine obs_inc_info(obs_kind, filter_kind, rectangular_quadrature, gaussian_likelihood_tails, &
   sort_obs_inc, spread_restoration, bounded_below, bounded_above, lower_bound, upper_bound)

integer,  intent(in)  :: obs_kind
integer,  intent(inout) :: filter_kind
logical,  intent(inout) :: rectangular_quadrature, gaussian_likelihood_tails
logical,  intent(inout) :: sort_obs_inc
logical,  intent(inout) :: spread_restoration
logical,  intent(inout) :: bounded_below, bounded_above
real(r8), intent(inout) :: lower_bound,  upper_bound

integer :: QTY_loc(1)
character(len=129) :: kind_name

! The information arguments are all intent (inout). This means that if they are not set
! here, they retain the default values from the assim_tools_mod namelist. Bounds don't exist 
! in that namelist, so default values are set in assim_tools_mod just before the call to here.

! Temporary approach for setting the details of how to assimilate this observation
! This example is designed to reproduce the squared forward operator results from paper

!get actual name of QTY from integer index
kind_name = get_name_for_quantity(obs_kind)
write(*,*) 'kind_name: ', kind_name

!find location of QTY in qcf_table_data structure
QTY_loc = findloc(qcf_table_row_headers, kind_name)
write(*,*) 'findloc of kind: ', QTY_loc(1)

if (QTY_loc(1) == 0) then
   write(*,*) 'QTY not in table, using default values'

   !use default values
   filter_kind = BOUNDED_NORMAL_RHF
   bounded_below = .false.;    bounded_above = .false.
   lower_bound   = missing_r8; upper_bound   = missing_r8
   sort_obs_inc = .false.; spread_restoration = .false.
   ! Default settings for now for Icepack and tracer model tests (sort_obs_inc, spread_restoration)

   else
      filter_kind = qcf_table_data(QTY_loc(1))%obs_inc_info%filter_kind !filter_kind has a check in obs_increment
      sort_obs_inc = qcf_table_data(QTY_loc(1))%obs_inc_info%sort_obs_inc
      spread_restoration = qcf_table_data(QTY_loc(1))%obs_inc_info%spread_restoration
      bounded_below = qcf_table_data(QTY_loc(1))%obs_inc_info%bounded_below
      bounded_above = qcf_table_data(QTY_loc(1))%obs_inc_info%bounded_above
      lower_bound = qcf_table_data(QTY_loc(1))%obs_inc_info%lower_bound
      upper_bound = qcf_table_data(QTY_loc(1))%obs_inc_info%upper_bound

endif

write(*,*) 'obs_inc_info: ', filter_kind, sort_obs_inc, spread_restoration, bounded_below, bounded_above, lower_bound, upper_bound

! Only need to set these two for options the original RHF implementation
!!!rectangular_quadrature = .true.
!!!gaussian_likelihood_tails = .false.

end subroutine obs_inc_info

!------------------------------------------------------------------------


subroutine init_algorithm_info_mod(qcf_table_filename)

character(len=129), intent(in) :: qcf_table_filename

integer :: numrows 
integer :: nlines
integer :: io
integer, parameter :: fileid = 10 !file identifier

open(unit=fileid, file=qcf_table_filename)
nlines = 0

do !do loop to get number of rows (or QTY's) in the table
  read(fileid,*,iostat=io)
  if(io/=0) exit
  nlines = nlines + 1
end do
close(fileid)

numrows = nlines - 2
print *, 'numrows: ', numrows

allocate(qcf_table_data(numrows))
allocate(qcf_table_row_headers(numrows))

call read_qcf_table(qcf_table_filename)

end subroutine init_algorithm_info_mod

!------------------------------------------------------------------------


subroutine read_qcf_table(qcf_table_filename)

! Reads in the QCEFF input options from tabular data file

character(len=129), intent(in) :: qcf_table_filename

integer, parameter :: fileid = 10 !file identifier
integer :: row

character(len=129), dimension(4) :: header1
character(len=129), dimension(29) :: header2

open(unit=fileid, file=qcf_table_filename)

! skip the headers 
read(fileid, *) header1
read(fileid, *) header2
write(*, *) "header1: ", header1
write(*, *) "header2: ", header2

! read in table values directly to qcf_table_data type
do row = 1, size(qcf_table_data)
   read(fileid, *) qcf_table_row_headers(row), qcf_table_data(row)%obs_error_info%bounded_below, qcf_table_data(row)%obs_error_info%bounded_above, &
                   qcf_table_data(row)%obs_error_info%lower_bound, qcf_table_data(row)%obs_error_info%upper_bound, qcf_table_data(row)%probit_inflation%dist_type, &
                   qcf_table_data(row)%probit_inflation%bounded_below, qcf_table_data(row)%probit_inflation%bounded_above, &
                   qcf_table_data(row)%probit_inflation%lower_bound, qcf_table_data(row)%probit_inflation%upper_bound, qcf_table_data(row)%probit_state%dist_type, &
                   qcf_table_data(row)%probit_state%bounded_below, qcf_table_data(row)%probit_state%bounded_above, &
                   qcf_table_data(row)%probit_state%lower_bound, qcf_table_data(row)%probit_state%upper_bound, qcf_table_data(row)%probit_extended_state%dist_type, &
                   qcf_table_data(row)%probit_extended_state%bounded_below, qcf_table_data(row)%probit_extended_state%bounded_above, &
                   qcf_table_data(row)%probit_extended_state%lower_bound, qcf_table_data(row)%probit_extended_state%upper_bound, &
                   qcf_table_data(row)%obs_inc_info%filter_kind, qcf_table_data(row)%obs_inc_info%rectangular_quadrature, &
                   qcf_table_data(row)%obs_inc_info%gaussian_likelihood_tails, qcf_table_data(row)%obs_inc_info%sort_obs_inc, &
                   qcf_table_data(row)%obs_inc_info%spread_restoration, qcf_table_data(row)%obs_inc_info%bounded_below, qcf_table_data(row)%obs_inc_info%bounded_above, &
                   qcf_table_data(row)%obs_inc_info%lower_bound, qcf_table_data(row)%obs_inc_info%upper_bound
end do

close(fileid)

call write_qcf_table()

end subroutine read_qcf_table

!------------------------------------------------------------------------


subroutine write_qcf_table()

! DRAFT SUBROUTINE 
! write to check values were correctly assigned
! testing for findloc

character(len=30), parameter :: tester_QTY = 'QTY_GPSRO'
integer :: QTY_loc(1)

character(len=30), parameter :: tester_QTY0 = 'QTY_DUMMY'
integer :: QTY_loc0(1)

integer :: row

write(*,*), 'SIZE: ', size(qcf_table_data)

do row = 1, size(qcf_table_data)
   write(*, *) "qcf_table_row_headers(", row, "): ", qcf_table_row_headers(row)
   write(*, *) "qcf_table_data(", row, "): "
   write(*, *) qcf_table_data(row)%obs_error_info%bounded_below, qcf_table_data(row)%obs_error_info%bounded_above, &
               qcf_table_data(row)%obs_error_info%lower_bound, qcf_table_data(row)%obs_error_info%upper_bound, qcf_table_data(row)%probit_inflation%dist_type, &
               qcf_table_data(row)%probit_inflation%bounded_below, qcf_table_data(row)%probit_inflation%bounded_above, &
               qcf_table_data(row)%probit_inflation%lower_bound, qcf_table_data(row)%probit_inflation%upper_bound, qcf_table_data(row)%probit_state%dist_type, &
               qcf_table_data(row)%probit_state%bounded_below, qcf_table_data(row)%probit_state%bounded_above, &
               qcf_table_data(row)%probit_state%lower_bound, qcf_table_data(row)%probit_state%upper_bound, qcf_table_data(row)%probit_extended_state%dist_type, &
               qcf_table_data(row)%probit_extended_state%bounded_below, qcf_table_data(row)%probit_extended_state%bounded_above, &
               qcf_table_data(row)%probit_extended_state%lower_bound, qcf_table_data(row)%probit_extended_state%upper_bound, &
               qcf_table_data(row)%obs_inc_info%filter_kind, qcf_table_data(row)%obs_inc_info%rectangular_quadrature, &
               qcf_table_data(row)%obs_inc_info%gaussian_likelihood_tails, qcf_table_data(row)%obs_inc_info%sort_obs_inc, &
               qcf_table_data(row)%obs_inc_info%spread_restoration, qcf_table_data(row)%obs_inc_info%bounded_below, qcf_table_data(row)%obs_inc_info%bounded_above, &
               qcf_table_data(row)%obs_inc_info%lower_bound, qcf_table_data(row)%obs_inc_info%upper_bound
end do

QTY_loc = findloc(qcf_table_row_headers, tester_QTY)
write(*, *) 'findloc of QTY_GPSRO: ', QTY_loc(1)

QTY_loc0 = findloc(qcf_table_row_headers, tester_QTY0)
write(*, *) 'findloc of invalid QTY (QTY_DUMMY): ', QTY_loc0(1)

end subroutine write_qcf_table

!------------------------------------------------------------------------


subroutine end_algorithm_info_mod()

deallocate(qcf_table_data)
deallocate(qcf_table_row_headers)

end subroutine end_algorithm_info_mod

!----------------------------------------------------------------------

end module algorithm_info_mod
