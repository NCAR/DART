! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module algorithm_info_mod

use types_mod, only : r8, i8, missing_r8

use obs_def_mod, only : obs_def_type, get_obs_def_type_of_obs, get_obs_def_error_variance
use obs_kind_mod, only : get_quantity_for_type_of_obs, get_name_for_quantity, get_index_for_quantity

use utilities_mod, only : error_handler, E_ERR, E_MSG, open_file, close_file, to_upper

use assim_model_mod, only : get_state_meta_data
use location_mod, only    : location_type

use distribution_params_mod, only : NORMAL_DISTRIBUTION, BOUNDED_NORMAL_RH_DISTRIBUTION, &
   GAMMA_DISTRIBUTION, BETA_DISTRIBUTION, LOG_NORMAL_DISTRIBUTION, UNIFORM_DISTRIBUTION,  &
   PARTICLE_FILTER_DISTRIBUTION

implicit none
private

character(len=512) :: errstring
character(len=*), parameter :: source = 'algorithm_info_mod.f90'

logical :: module_initialized = .false.
logical :: qcf_table_listed = .false.

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
          EAKF, ENKF, BOUNDED_NORMAL_RHF, UNBOUNDED_RHF, GAMMA_FILTER

!Creates the type definitions for the QCF table
type obs_error_info_type
   logical               :: bounded_below, bounded_above
   real(r8)              :: lower_bound,   upper_bound
end type

type probit_inflation_type
   character(len=129)    :: dist_type
   logical               :: bounded_below, bounded_above
   real(r8)              :: lower_bound,   upper_bound
end type

type probit_state_type
   character(len=129)    :: dist_type
   logical               :: bounded_below, bounded_above
   real(r8)              :: lower_bound,   upper_bound
end type

type probit_extended_state_type
   character(len=129)    :: dist_type
   logical               :: bounded_below, bounded_above
   real(r8)              :: lower_bound,   upper_bound
end type

type obs_inc_info_type
   character(len=129)    :: filter_kind
   logical               :: rectangular_quadrature, gaussian_likelihood_tails
   logical               :: sort_obs_inc, spread_restoration
   logical               :: bounded_below, bounded_above
   real(r8)              :: lower_bound,   upper_bound
end type

type qcf_table_data_type
   type(obs_error_info_type)        :: obs_error_info
   type(probit_inflation_type)      :: probit_inflation
   type(probit_state_type)          :: probit_state
   type(probit_extended_state_type) :: probit_extended_state
   type(obs_inc_info_type)          :: obs_inc_info
end type

character(len=129), dimension(4) :: header1
character(len=129), dimension(29) :: header2

character(len=129), allocatable :: qcf_table_row_headers(:)
type(qcf_table_data_type), allocatable :: qcf_table_data(:)

! Provides routines that give information about details of algorithms for 
! observation error sampling, observation increments, and the transformations
! for regression and inflation in probit space. 
! For now, it is convenient to have these in a single module since several
! users will be developing their own problem specific versions of these
! subroutines. This will avoid constant merge conflicts as other parts of the
! assimilation code are updated.

contains

!-------------------------------------------------------------------------


subroutine init_algorithm_info_mod(qcf_table_filename)

! Gets number of lines/QTYs in the QCF table, allocates space for the table data

character(len=129), intent(in) :: qcf_table_filename

integer :: fileid
integer :: io

integer :: numrows
integer :: nlines

if (module_initialized) return
module_initialized = .true.

if (qcf_table_filename == '') then
   write(errstring,*) 'No QCF table file listed in namelist, using default values for all QTYs'
   call error_handler(E_MSG, 'init_algorithm_info_mod:', errstring, source)
   return
endif

qcf_table_listed = .true.
fileid = open_file(trim(qcf_table_filename), 'formatted', 'read')

! Do loop to get number of rows (or QTY's) in the table
nlines = 0
do
  read(fileid,*,iostat=io)
  if(io/=0) exit
  nlines = nlines + 1
end do

call close_file(fileid)

numrows = nlines - 2

allocate(qcf_table_data(numrows))
allocate(qcf_table_row_headers(numrows))

call read_qcf_table(qcf_table_filename)
call assert_qcf_table_version()
call verify_qcf_table_data()
call log_qcf_table_data()

end subroutine init_algorithm_info_mod

!------------------------------------------------------------------------


subroutine read_qcf_table(qcf_table_filename)

! Reads in the QCEFF input options from tabular data file

character(len=129), intent(in) :: qcf_table_filename

integer :: fileid
integer :: row

if (.not. module_initialized) call init_algorithm_info_mod(qcf_table_filename)

fileid = open_file(trim(qcf_table_filename), 'formatted', 'read')

! skip the headers
read(fileid, *) header1
read(fileid, *) header2

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

call close_file(fileid)

end subroutine read_qcf_table

!------------------------------------------------------------------------


subroutine obs_error_info(obs_def, error_variance, &
   bounded_below, bounded_above, lower_bound, upper_bound)

! Computes information needed to compute error sample for this observation
! This is called by perfect_model_obs when generating noisy obs
type(obs_def_type), intent(in)  :: obs_def
real(r8),           intent(out) :: error_variance
logical,            intent(out) :: bounded_below, bounded_above
real(r8),           intent(out) :: lower_bound, upper_bound

integer     :: obs_type, obs_qty
integer(i8) :: state_var_index
type(location_type) :: temp_loc

integer :: QTY_loc(1)
character(len=129) :: qty_name

! Get the kind of the observation
obs_type = get_obs_def_type_of_obs(obs_def)
! If it is negative, it is an identity obs
if(obs_type < 0) then
   state_var_index = -1 * obs_type
   call get_state_meta_data(state_var_index, temp_loc, obs_qty)
else
   obs_qty = get_quantity_for_type_of_obs(obs_type)
endif

! Get the default error variance
error_variance = get_obs_def_error_variance(obs_def)

!use default values if qcf_table_filename is not in namelist
if (.not. qcf_table_listed) then
   bounded_below = .false.;    bounded_above = .false.
   lower_bound   = missing_r8; upper_bound   = missing_r8
   return
endif

!get actual name of QTY from integer index
qty_name = get_name_for_quantity(obs_qty)

!find location of QTY in qcf_table_data structure
QTY_loc = findloc(qcf_table_row_headers, qty_name)

if (QTY_loc(1) == 0) then
   !use default values if QTY is not in table
   bounded_below = .false.;    bounded_above = .false.
   lower_bound   = missing_r8; upper_bound   = missing_r8

   else
      bounded_below = qcf_table_data(QTY_loc(1))%obs_error_info%bounded_below
      bounded_above = qcf_table_data(QTY_loc(1))%obs_error_info%bounded_above
      lower_bound = qcf_table_data(QTY_loc(1))%obs_error_info%lower_bound
      upper_bound = qcf_table_data(QTY_loc(1))%obs_error_info%upper_bound

endif

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
character(len=129) :: qty_name

integer :: dist_type_loc(1)
character(len=129), dimension(7) :: possible_dist_types
integer, dimension(7) :: possible_dist_type_ints
character(len=129) :: dist_type_string

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

! Fill arrays with possible dist_type strings and corresponding ints
possible_dist_types(1) = 'NORMAL_DISTRIBUTION'
possible_dist_types(2) = 'BOUNDED_NORMAL_RH_DISTRIBUTION'
possible_dist_types(3) = 'GAMMA_DISTRIBUTION'
possible_dist_types(4) = 'BETA_DISTRIBUTION'
possible_dist_types(5) = 'LOG_NORMAL_DISTRIBUTION'
possible_dist_types(6) = 'UNIFORM_DISTRIBUTION'
possible_dist_types(7) = 'PARTICLE_FILTER_DISTRIBUTION'

possible_dist_type_ints(1) = 1
possible_dist_type_ints(2) = 2
possible_dist_type_ints(3) = 3
possible_dist_type_ints(4) = 4
possible_dist_type_ints(5) = 5
possible_dist_type_ints(6) = 6
possible_dist_type_ints(7) = 7

!use default values if qcf_table_filename is not in namelist
if (.not. qcf_table_listed) then
   dist_type = BOUNDED_NORMAL_RH_DISTRIBUTION
   bounded_below = .false.;    bounded_above = .false.
   lower_bound   = missing_r8; upper_bound   = missing_r8
   return
endif

!get actual name of QTY from integer index
qty_name = get_name_for_quantity(kind)

!find location of QTY in qcf_table_data structure
QTY_loc = findloc(qcf_table_row_headers, qty_name)

if (QTY_loc(1) == 0) then
   !use default values if QTY is not in table
   dist_type = BOUNDED_NORMAL_RH_DISTRIBUTION
   bounded_below = .false.;    bounded_above = .false.
   lower_bound   = missing_r8; upper_bound   = missing_r8

   elseif(is_inflation) then
   ! Case for inflation transformation

      ! Comparing the dist_type in string format to list of potential dist_types
      dist_type_string = qcf_table_data(QTY_loc(1))%probit_inflation%dist_type
      call to_upper(dist_type_string)
      dist_type_loc = findloc(possible_dist_types, trim(dist_type_string))

      if (dist_type_loc(1) == 0) then
         write(errstring, *) 'Invalid dist_type: ', trim(dist_type_string)
         call error_handler(E_ERR, 'probit_dist_info:', errstring, source)

      else
         dist_type = possible_dist_type_ints(dist_type_loc(1))
      endif

      bounded_below = qcf_table_data(QTY_loc(1))%probit_inflation%bounded_below
      bounded_above = qcf_table_data(QTY_loc(1))%probit_inflation%bounded_above
      lower_bound = qcf_table_data(QTY_loc(1))%probit_inflation%lower_bound
      upper_bound = qcf_table_data(QTY_loc(1))%probit_inflation%upper_bound

   elseif(is_state) then
   ! Case for state variable priors

      ! Comparing the dist_type in string format to list of potential dist_types
      dist_type_string = qcf_table_data(QTY_loc(1))%probit_state%dist_type
      call to_upper(dist_type_string)
      dist_type_loc = findloc(possible_dist_types, trim(dist_type_string))

      if (dist_type_loc(1) == 0) then
         write(errstring, *) 'Invalid dist_type: ', trim(dist_type_string)
         call error_handler(E_ERR, 'probit_dist_info:', errstring, source)

      else
         dist_type = possible_dist_type_ints(dist_type_loc(1))
      endif

      bounded_below = qcf_table_data(QTY_loc(1))%probit_state%bounded_below
      bounded_above = qcf_table_data(QTY_loc(1))%probit_state%bounded_above
      lower_bound = qcf_table_data(QTY_loc(1))%probit_state%lower_bound
      upper_bound = qcf_table_data(QTY_loc(1))%probit_state%upper_bound

   else
   ! This case is for observation (extended state) priors

      ! Comparing the dist_type in string format to list of potential dist_types
      dist_type_string = qcf_table_data(QTY_loc(1))%probit_extended_state%dist_type
      call to_upper(dist_type_string)
      dist_type_loc = findloc(possible_dist_types, trim(dist_type_string))

      if (dist_type_loc(1) == 0) then
         write(errstring, *) 'Invalid dist_type: ', trim(dist_type_string)
         call error_handler(E_ERR, 'probit_dist_info:', errstring, source)

      else
         dist_type = possible_dist_type_ints(dist_type_loc(1))
      endif

      bounded_below = qcf_table_data(QTY_loc(1))%probit_extended_state%bounded_below
      bounded_above = qcf_table_data(QTY_loc(1))%probit_extended_state%bounded_above
      lower_bound = qcf_table_data(QTY_loc(1))%probit_extended_state%lower_bound
      upper_bound = qcf_table_data(QTY_loc(1))%probit_extended_state%upper_bound

endif

end subroutine probit_dist_info

!------------------------------------------------------------------------


subroutine obs_inc_info(obs_qty, filter_kind, rectangular_quadrature, gaussian_likelihood_tails, &
   sort_obs_inc, spread_restoration, bounded_below, bounded_above, lower_bound, upper_bound)

integer,  intent(in)  :: obs_qty
integer,  intent(inout) :: filter_kind
logical,  intent(inout) :: rectangular_quadrature, gaussian_likelihood_tails
logical,  intent(inout) :: sort_obs_inc
logical,  intent(inout) :: spread_restoration
logical,  intent(inout) :: bounded_below, bounded_above
real(r8), intent(inout) :: lower_bound,  upper_bound

integer :: QTY_loc(1)
character(len=129) :: qty_name

integer :: filter_kind_loc(1)
character(len=129), dimension(5) :: possible_filter_kinds
integer, dimension(5) :: possible_filter_kind_ints
character(len=129) :: filter_kind_string

! The information arguments are all intent (inout). This means that if they are not set
! here, they retain the default values from the assim_tools_mod namelist. Bounds don't exist 
! in that namelist, so default values are set in assim_tools_mod just before the call to here.

! Temporary approach for setting the details of how to assimilate this observation
! This example is designed to reproduce the squared forward operator results from paper

! Fill arrays with possible filter_kind strings and corresponding ints
possible_filter_kinds(1) = 'EAKF'
possible_filter_kinds(2) = 'ENKF'
possible_filter_kinds(3) = 'UNBOUNDED_RHF'
possible_filter_kinds(4) = 'GAMMA_FILTER'
possible_filter_kinds(5) = 'BOUNDED_NORMAL_RHF'

possible_filter_kind_ints(1) = 1
possible_filter_kind_ints(2) = 2
possible_filter_kind_ints(3) = 8
possible_filter_kind_ints(4) = 11
possible_filter_kind_ints(5) = 101

!use default values if qcf_table_filename is not in namelist
if (.not. qcf_table_listed) then
   filter_kind = BOUNDED_NORMAL_RHF
   bounded_below = .false.;    bounded_above = .false.
   lower_bound   = missing_r8; upper_bound   = missing_r8
   sort_obs_inc = .false.; spread_restoration = .false.
   return
endif

!get actual name of QTY from integer index
qty_name = get_name_for_quantity(obs_qty)

!find location of QTY in qcf_table_data structure
QTY_loc = findloc(qcf_table_row_headers, qty_name)

if (QTY_loc(1) == 0) then
   !use default values if QTY is not in table
   filter_kind = BOUNDED_NORMAL_RHF
   bounded_below = .false.;    bounded_above = .false.
   lower_bound   = missing_r8; upper_bound   = missing_r8
   sort_obs_inc = .false.; spread_restoration = .false.
   ! Default settings for now for Icepack and tracer model tests (sort_obs_inc, spread_restoration)

   else

      ! Comparing the filter_kind in string format to list of potential filter_kinds
      filter_kind_string = qcf_table_data(QTY_loc(1))%obs_inc_info%filter_kind
      call to_upper(filter_kind_string)
      filter_kind_loc = findloc(possible_filter_kinds, trim(filter_kind_string))

      if (filter_kind_loc(1) == 0) then
         write(errstring, *) 'Invalid filter_kind: ', trim(filter_kind_string)
         call error_handler(E_ERR, 'obs_inc_info:', errstring, source)

      else
         filter_kind = possible_filter_kind_ints(filter_kind_loc(1))
      endif
     
      sort_obs_inc = qcf_table_data(QTY_loc(1))%obs_inc_info%sort_obs_inc
      spread_restoration = qcf_table_data(QTY_loc(1))%obs_inc_info%spread_restoration
      bounded_below = qcf_table_data(QTY_loc(1))%obs_inc_info%bounded_below
      bounded_above = qcf_table_data(QTY_loc(1))%obs_inc_info%bounded_above
      lower_bound = qcf_table_data(QTY_loc(1))%obs_inc_info%lower_bound
      upper_bound = qcf_table_data(QTY_loc(1))%obs_inc_info%upper_bound

endif

! Only need to set these two for options the original RHF implementation
!!!rectangular_quadrature = .true.
!!!gaussian_likelihood_tails = .false.

end subroutine obs_inc_info

!------------------------------------------------------------------------


subroutine assert_qcf_table_version()

! Subroutine to ensure the correct version of the QCF table is being used

if (trim(header1(4)) /= '1') then
   write(errstring,*) 'Using outdated/incorrect version of the QCF table'
   call error_handler(E_ERR, 'assert_qcf_table_version:', errstring, source)
endif

end subroutine assert_qcf_table_version

!------------------------------------------------------------------------


subroutine verify_qcf_table_data()

! Subroutine to ensure that the data in the QCF table is valid 

integer :: varid
integer :: row

if (.not. qcf_table_listed) return

!Checks that all bounds are valid; currently checks that the lower bound in less than the upper
!Here we could add more specific checks if we have known limits on the bounds
do row = 1, size(qcf_table_data)

   if (qcf_table_data(row)%obs_error_info%bounded_below .and.  qcf_table_data(row)%obs_error_info%bounded_above) then
      if(qcf_table_data(row)%obs_error_info%lower_bound > qcf_table_data(row)%obs_error_info%upper_bound) then
         write(errstring,*) 'Invalid bounds in obs_error_info'
         call error_handler(E_ERR, 'verify_qcf_table_data:', errstring, source)
      endif
   endif
   if (qcf_table_data(row)%probit_inflation%bounded_below .and. qcf_table_data(row)%probit_inflation%bounded_above) then
      if(qcf_table_data(row)%probit_inflation%lower_bound > qcf_table_data(row)%probit_inflation%upper_bound) then
         write(errstring,*) 'Invalid bounds in probit_inflation'
         call error_handler(E_ERR, 'verify_qcf_table_data:', errstring, source)
      endif
   endif
   if(qcf_table_data(row)%probit_state%bounded_below .and. qcf_table_data(row)%probit_state%bounded_above) then
      if(qcf_table_data(row)%probit_state%lower_bound > qcf_table_data(row)%probit_state%upper_bound) then
         write(errstring,*) 'Invalid bounds in probit_state'
         call error_handler(E_ERR, 'verify_qcf_table_data:', errstring, source)
      endif
   endif
   if(qcf_table_data(row)%probit_extended_state%bounded_below .and. qcf_table_data(row)%probit_extended_state%bounded_above) then
      if(qcf_table_data(row)%probit_extended_state%lower_bound > qcf_table_data(row)%probit_extended_state%upper_bound) then
         write(errstring,*) 'Invalid bounds in probit_extended_state'
         call error_handler(E_ERR, 'verify_qcf_table_data:', errstring, source)
      endif
   endif
   if(qcf_table_data(row)%obs_inc_info%bounded_below .and. qcf_table_data(row)%obs_inc_info%bounded_above) then
       if(qcf_table_data(row)%obs_inc_info%lower_bound > qcf_table_data(row)%obs_inc_info%upper_bound) then
          write(errstring,*) 'Invalid bounds in obs_inc_info'
          call error_handler(E_ERR, 'verify_qcf_table_data:', errstring, source)
       endif
    endif
end do

!Ensures that all QTYs listed in the table exist in DART
do row = 1, size(qcf_table_data)
   varid = get_index_for_quantity(trim(qcf_table_row_headers(row)))
   if(varid == -1) then
      write(errstring,*) trim(qcf_table_row_headers(row)), ' is not a valid DART QTY'
      call error_handler(E_ERR, 'verify_qcf_table_data:', errstring, source)
   endif
end do

!Ensures that there are no duplicate QTYs in the table
do row = 1, size(qcf_table_data)
   if(count(qcf_table_row_headers==trim(qcf_table_row_headers(row))) > 1) then
      write(errstring,*) trim(qcf_table_row_headers(row)), ' has multiple entries in the table'
      call error_handler(E_ERR, 'verify_qcf_table_data:', errstring, source)
   endif
end do

end subroutine verify_qcf_table_data

!------------------------------------------------------------------------


subroutine log_qcf_table_data()

! Subroutine to write the data in QCF table to dart_log
character(len=2000) :: log_msg
integer :: row

if (.not. qcf_table_listed) return

call error_handler(E_MSG, '', '', source) !Writing blank line to log
call error_handler(E_MSG, 'log_qcf_table_data:', 'Logging the data in the QCF Table', source)

! Write the table headers to the dart_log and terminal
write(log_msg, '(A4, A6, A9, A)') header1(:)
call error_handler(E_MSG, 'log_qcf_table_data:', trim(log_msg), source)

write(log_msg,'(3A14, 2A12, 3(A10, 2A14, 2A12), A12, A23, A26, A13, A19, 2A14, 2A12)') header2(:)
call error_handler(E_MSG, 'log_qcf_table_data:', trim(log_msg), source)

! Write the table data to the dart_log and terminal
do row = 1, size(qcf_table_data)
   write(log_msg, *) trim(qcf_table_row_headers(row)), qcf_table_data(row)%obs_error_info%bounded_below, qcf_table_data(row)%obs_error_info%bounded_above, &
               qcf_table_data(row)%obs_error_info%lower_bound, qcf_table_data(row)%obs_error_info%upper_bound, trim(qcf_table_data(row)%probit_inflation%dist_type), &
               qcf_table_data(row)%probit_inflation%bounded_below, qcf_table_data(row)%probit_inflation%bounded_above, &
               qcf_table_data(row)%probit_inflation%lower_bound, qcf_table_data(row)%probit_inflation%upper_bound, trim(qcf_table_data(row)%probit_state%dist_type), &
               qcf_table_data(row)%probit_state%bounded_below, qcf_table_data(row)%probit_state%bounded_above, &
               qcf_table_data(row)%probit_state%lower_bound, qcf_table_data(row)%probit_state%upper_bound, trim(qcf_table_data(row)%probit_extended_state%dist_type), &
               qcf_table_data(row)%probit_extended_state%bounded_below, qcf_table_data(row)%probit_extended_state%bounded_above, &
               qcf_table_data(row)%probit_extended_state%lower_bound, qcf_table_data(row)%probit_extended_state%upper_bound, &
               trim(qcf_table_data(row)%obs_inc_info%filter_kind), qcf_table_data(row)%obs_inc_info%rectangular_quadrature, &
               qcf_table_data(row)%obs_inc_info%gaussian_likelihood_tails, qcf_table_data(row)%obs_inc_info%sort_obs_inc, &
               qcf_table_data(row)%obs_inc_info%spread_restoration, qcf_table_data(row)%obs_inc_info%bounded_below, qcf_table_data(row)%obs_inc_info%bounded_above, &
               qcf_table_data(row)%obs_inc_info%lower_bound, qcf_table_data(row)%obs_inc_info%upper_bound
call error_handler(E_MSG, 'log_qcf_table_data:', trim(log_msg), source)
end do

call error_handler(E_MSG, '', '', source) !Writing blank line to log

end subroutine log_qcf_table_data

!------------------------------------------------------------------------


subroutine end_algorithm_info_mod()

if (.not. module_initialized) return
module_initialized = .false.

if (.not. qcf_table_listed) return

deallocate(qcf_table_data)
deallocate(qcf_table_row_headers)

end subroutine end_algorithm_info_mod

!----------------------------------------------------------------------

end module algorithm_info_mod
