! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module algorithm_info_mod

use types_mod, only : r8, i8, MISSING_R8

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
logical :: use_qty_defaults = .true.

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
   integer               :: dist_type
   logical               :: bounded_below, bounded_above
   real(r8)              :: lower_bound,   upper_bound
end type

type probit_state_type
   integer               :: dist_type
   logical               :: bounded_below, bounded_above
   real(r8)              :: lower_bound,   upper_bound
end type

type probit_extended_state_type
   integer               :: dist_type
   logical               :: bounded_below, bounded_above
   real(r8)              :: lower_bound,   upper_bound
end type

type obs_inc_info_type
   integer               :: filter_kind
   logical               :: bounded_below, bounded_above
   real(r8)              :: lower_bound,   upper_bound
end type

type algorithm_info_type
   type(obs_error_info_type)        :: obs_error_info
   type(probit_inflation_type)      :: probit_inflation
   type(probit_state_type)          :: probit_state
   type(probit_extended_state_type) :: probit_extended_state
   type(obs_inc_info_type)          :: obs_inc_info
end type

integer, parameter :: HEADER_LINES = 2
character(len=129), dimension(4) :: header1
character(len=129), dimension(25) :: header2 ! Number of table columns plus 1

character(len=129), allocatable :: specified_qtys(:)
type(algorithm_info_type), allocatable :: qcf_table_data(:)

character(len=129), allocatable :: dist_type_string_probit_inflation(:)
character(len=129), allocatable :: dist_type_string_probit_state(:)
character(len=129), allocatable :: dist_type_string_probit_extended_state(:)
character(len=129), allocatable :: filter_kind_string(:)

! Provides routines that give information about details of algorithms for 
! observation error sampling, observation increments, and the transformations
! for regression and inflation in probit space. 

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

use_qty_defaults = .false.
fileid = open_file(trim(qcf_table_filename), 'formatted', 'read')

! Do loop to get number of rows (or QTY's) in the table
nlines = 0
do
  read(fileid,*,iostat=io)
  if(io/=0) exit
  nlines = nlines + 1
end do

call close_file(fileid)

numrows = nlines - HEADER_LINES

allocate(specified_qtys(numrows))
allocate(qcf_table_data(numrows))
allocate(dist_type_string_probit_inflation(numrows))
allocate(dist_type_string_probit_state(numrows))
allocate(dist_type_string_probit_extended_state(numrows))
allocate(filter_kind_string(numrows))

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
   read(fileid, *) specified_qtys(row), qcf_table_data(row)%obs_error_info%bounded_below, qcf_table_data(row)%obs_error_info%bounded_above, &
                   qcf_table_data(row)%obs_error_info%lower_bound, qcf_table_data(row)%obs_error_info%upper_bound, dist_type_string_probit_inflation(row), &
                   qcf_table_data(row)%probit_inflation%bounded_below, qcf_table_data(row)%probit_inflation%bounded_above, &
                   qcf_table_data(row)%probit_inflation%lower_bound, qcf_table_data(row)%probit_inflation%upper_bound, dist_type_string_probit_state(row), &
                   qcf_table_data(row)%probit_state%bounded_below, qcf_table_data(row)%probit_state%bounded_above, &
                   qcf_table_data(row)%probit_state%lower_bound, qcf_table_data(row)%probit_state%upper_bound, dist_type_string_probit_extended_state(row), &
                   qcf_table_data(row)%probit_extended_state%bounded_below, qcf_table_data(row)%probit_extended_state%bounded_above, &
                   qcf_table_data(row)%probit_extended_state%lower_bound, qcf_table_data(row)%probit_extended_state%upper_bound, &
                   filter_kind_string(row), qcf_table_data(row)%obs_inc_info%bounded_below, qcf_table_data(row)%obs_inc_info%bounded_above, &
                   qcf_table_data(row)%obs_inc_info%lower_bound, qcf_table_data(row)%obs_inc_info%upper_bound

   ! Converting the distribution types (read in from table as a string) to its corresponding int value
   if (trim(dist_type_string_probit_inflation(row)) == 'NORMAL_DISTRIBUTION') then
      qcf_table_data(row)%probit_inflation%dist_type = NORMAL_DISTRIBUTION
   elseif (trim(dist_type_string_probit_inflation(row)) == 'BOUNDED_NORMAL_RH_DISTRIBUTION') then
      qcf_table_data(row)%probit_inflation%dist_type = BOUNDED_NORMAL_RH_DISTRIBUTION
   elseif (trim(dist_type_string_probit_inflation(row)) == 'GAMMA_DISTRIBUTION') then
      qcf_table_data(row)%probit_inflation%dist_type = GAMMA_DISTRIBUTION
   elseif (trim(dist_type_string_probit_inflation(row)) == 'BETA_DISTRIBUTION') then
      qcf_table_data(row)%probit_inflation%dist_type = BETA_DISTRIBUTION
   elseif (trim(dist_type_string_probit_inflation(row)) == 'LOG_NORMAL_DISTRIBUTION') then
      qcf_table_data(row)%probit_inflation%dist_type = LOG_NORMAL_DISTRIBUTION
   elseif (trim(dist_type_string_probit_inflation(row)) == 'UNIFORM_DISTRIBUTION') then
      qcf_table_data(row)%probit_inflation%dist_type = UNIFORM_DISTRIBUTION
   elseif (trim(dist_type_string_probit_inflation(row)) == 'PARTICLE_FILTER_DISTRIBUTION') then
      qcf_table_data(row)%probit_inflation%dist_type = PARTICLE_FILTER_DISTRIBUTION
   else
      write(errstring, *) 'Invalid distribution type for probit inflation: ', trim(dist_type_string_probit_inflation(row))
      call error_handler(E_ERR, 'read_qcf_table:', errstring, source)
   endif

   if (trim(dist_type_string_probit_state(row)) == 'NORMAL_DISTRIBUTION') then
      qcf_table_data(row)%probit_state%dist_type = NORMAL_DISTRIBUTION
   elseif (trim(dist_type_string_probit_state(row)) == 'BOUNDED_NORMAL_RH_DISTRIBUTION') then
      qcf_table_data(row)%probit_state%dist_type = BOUNDED_NORMAL_RH_DISTRIBUTION
   elseif (trim(dist_type_string_probit_state(row)) == 'GAMMA_DISTRIBUTION') then
      qcf_table_data(row)%probit_state%dist_type = GAMMA_DISTRIBUTION
   elseif (trim(dist_type_string_probit_state(row)) == 'BETA_DISTRIBUTION') then
      qcf_table_data(row)%probit_state%dist_type = BETA_DISTRIBUTION
   elseif (trim(dist_type_string_probit_state(row)) == 'LOG_NORMAL_DISTRIBUTION') then
      qcf_table_data(row)%probit_state%dist_type = LOG_NORMAL_DISTRIBUTION
   elseif (trim(dist_type_string_probit_state(row)) == 'UNIFORM_DISTRIBUTION') then
      qcf_table_data(row)%probit_state%dist_type = UNIFORM_DISTRIBUTION
   elseif (trim(dist_type_string_probit_state(row)) == 'PARTICLE_FILTER_DISTRIBUTION') then
      qcf_table_data(row)%probit_state%dist_type = PARTICLE_FILTER_DISTRIBUTION
   else
      write(errstring, *) 'Invalid distribution type for probit state: ', trim(dist_type_string_probit_state(row))
      call error_handler(E_ERR, 'read_qcf_table:', errstring, source)
   endif

   if (trim(dist_type_string_probit_extended_state(row)) == 'NORMAL_DISTRIBUTION') then
      qcf_table_data(row)%probit_extended_state%dist_type = NORMAL_DISTRIBUTION
   elseif (trim(dist_type_string_probit_extended_state(row)) == 'BOUNDED_NORMAL_RH_DISTRIBUTION') then
      qcf_table_data(row)%probit_extended_state%dist_type = BOUNDED_NORMAL_RH_DISTRIBUTION
   elseif (trim(dist_type_string_probit_extended_state(row)) == 'GAMMA_DISTRIBUTION') then
      qcf_table_data(row)%probit_extended_state%dist_type = GAMMA_DISTRIBUTION
   elseif (trim(dist_type_string_probit_extended_state(row)) == 'BETA_DISTRIBUTION') then
      qcf_table_data(row)%probit_extended_state%dist_type = BETA_DISTRIBUTION
   elseif (trim(dist_type_string_probit_extended_state(row)) == 'LOG_NORMAL_DISTRIBUTION') then
      qcf_table_data(row)%probit_extended_state%dist_type = LOG_NORMAL_DISTRIBUTION
   elseif (trim(dist_type_string_probit_extended_state(row)) == 'UNIFORM_DISTRIBUTION') then
      qcf_table_data(row)%probit_extended_state%dist_type = UNIFORM_DISTRIBUTION
   elseif (trim(dist_type_string_probit_extended_state(row)) == 'PARTICLE_FILTER_DISTRIBUTION') then
      qcf_table_data(row)%probit_extended_state%dist_type = PARTICLE_FILTER_DISTRIBUTION
   else
      write(errstring, *) 'Invalid distribution type for probit extended state: ', trim(dist_type_string_probit_extended_state(row))
      call error_handler(E_ERR, 'read_qcf_table:', errstring, source)
   endif

   ! Converting the filter kind (read in from table as a string) to its corresponding int value
   if (trim(filter_kind_string(row)) == 'EAKF') then
      qcf_table_data(row)%obs_inc_info%filter_kind = EAKF
   elseif (trim(filter_kind_string(row)) == 'ENKF') then
      qcf_table_data(row)%obs_inc_info%filter_kind = ENKF
   elseif (trim(filter_kind_string(row)) == 'UNBOUNDED_RHF') then
      qcf_table_data(row)%obs_inc_info%filter_kind = UNBOUNDED_RHF
   elseif (trim(filter_kind_string(row)) == 'GAMMA_FILTER') then
      qcf_table_data(row)%obs_inc_info%filter_kind = GAMMA_FILTER
   elseif (trim(filter_kind_string(row)) == 'BOUNDED_NORMAL_RHF') then
      qcf_table_data(row)%obs_inc_info%filter_kind = BOUNDED_NORMAL_RHF
   else
      write(errstring, *) 'Invalid filter kind: ', trim(filter_kind_string(row))
      call error_handler(E_ERR, 'read_qcf_table:', errstring, source)
   endif

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
if (use_qty_defaults) then
   bounded_below = .false.;    bounded_above = .false.
   lower_bound   = MISSING_R8; upper_bound   = MISSING_R8
   return
endif

!get actual name of QTY from integer index
qty_name = get_name_for_quantity(obs_qty)

!find location of QTY in qcf_table_data structure
QTY_loc = findloc(specified_qtys, qty_name)

if (QTY_loc(1) == 0) then
   !use default values if QTY is not in table
   bounded_below = .false.;    bounded_above = .false.
   lower_bound   = MISSING_R8; upper_bound   = MISSING_R8

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

! Have input information about the kind of the state or observation being transformed
! along with additional logical info that indicates whether this is an observation
! or state variable and about whether the transformation is being done for inflation
! or for regress. 
! Need to select the appropriate transform. At present, the options are NORMAL_DISTRIBUTION,
! BOUNDED_NORMAL_RH_DISTRIBUTION, GAMMA_DISTRIBUTION, BETA_DISTRIBUTION, LOG_NORMAL_DISTRIBUTION,
! UNIFORM_DISTRIBUTION, and PARTICLE_FILTER_DISTRIBUTION.
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

!use default values if qcf_table_filename is not in namelist
if (use_qty_defaults) then
   dist_type = BOUNDED_NORMAL_RH_DISTRIBUTION
   bounded_below = .false.;    bounded_above = .false.
   lower_bound   = MISSING_R8; upper_bound   = MISSING_R8
   return
endif

!get actual name of QTY from integer index
qty_name = get_name_for_quantity(kind)

!find location of QTY in qcf_table_data structure
QTY_loc = findloc(specified_qtys, qty_name)

if (QTY_loc(1) == 0) then
   !use default values if QTY is not in table
   dist_type = BOUNDED_NORMAL_RH_DISTRIBUTION
   bounded_below = .false.;    bounded_above = .false.
   lower_bound   = MISSING_R8; upper_bound   = MISSING_R8

   elseif(is_inflation) then
   ! Case for inflation transformation

      dist_type = qcf_table_data(QTY_loc(1))%probit_inflation%dist_type
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

end subroutine probit_dist_info

!------------------------------------------------------------------------


subroutine obs_inc_info(obs_qty, filter_kind, bounded_below, bounded_above, &
                        lower_bound, upper_bound)

integer,  intent(in)  :: obs_qty
integer,  intent(out) :: filter_kind
logical,  intent(out) :: bounded_below, bounded_above
real(r8), intent(out) :: lower_bound,  upper_bound

integer :: QTY_loc(1)
character(len=129) :: qty_name

!use default values if qcf_table_filename is not in namelist
if (use_qty_defaults) then
   filter_kind = BOUNDED_NORMAL_RHF
   bounded_below = .false.;    bounded_above = .false.
   lower_bound   = MISSING_R8; upper_bound   = MISSING_R8
   return
endif

!get actual name of QTY from integer index
qty_name = get_name_for_quantity(obs_qty)

!find location of QTY in qcf_table_data structure
QTY_loc = findloc(specified_qtys, qty_name)

if (QTY_loc(1) == 0) then
   !use default values if QTY is not in table
   filter_kind = BOUNDED_NORMAL_RHF
   bounded_below = .false.;    bounded_above = .false.
   lower_bound   = MISSING_R8; upper_bound   = MISSING_R8

   else

      filter_kind = qcf_table_data(QTY_loc(1))%obs_inc_info%filter_kind
      bounded_below = qcf_table_data(QTY_loc(1))%obs_inc_info%bounded_below
      bounded_above = qcf_table_data(QTY_loc(1))%obs_inc_info%bounded_above
      lower_bound = qcf_table_data(QTY_loc(1))%obs_inc_info%lower_bound
      upper_bound = qcf_table_data(QTY_loc(1))%obs_inc_info%upper_bound

endif

end subroutine obs_inc_info

!------------------------------------------------------------------------


subroutine assert_qcf_table_version()

! Subroutine to ensure the correct version of the QCF table is being used

character(1), parameter :: QCF_VERSION = '1'

if (trim(header1(4)) /= QCF_VERSION) then
   write(errstring,*) 'Using outdated/incorrect version of the QCF table'
   call error_handler(E_ERR, 'assert_qcf_table_version:', errstring, source)
endif

end subroutine assert_qcf_table_version

!------------------------------------------------------------------------


subroutine verify_qcf_table_data()

! Subroutine to ensure that the data in the QCF table is valid 

integer :: varid
integer :: row

if (use_qty_defaults) return

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
   varid = get_index_for_quantity(trim(specified_qtys(row)))
   if(varid == -1) then
      write(errstring,*) trim(specified_qtys(row)), ' is not a valid DART QTY'
      call error_handler(E_ERR, 'verify_qcf_table_data:', errstring, source)
   endif
end do

!Ensures that there are no duplicate QTYs in the table
do row = 1, size(qcf_table_data)
   if(count(specified_qtys==trim(specified_qtys(row))) > 1) then
      write(errstring,*) trim(specified_qtys(row)), ' has multiple entries in the table'
      call error_handler(E_ERR, 'verify_qcf_table_data:', errstring, source)
   endif
end do

end subroutine verify_qcf_table_data

!------------------------------------------------------------------------


subroutine log_qcf_table_data()

! Subroutine to write the data in QCF table to dart_log
character(len=2000) :: log_msg
integer :: row

if (use_qty_defaults) return

call error_handler(E_MSG, '', '', source) !Writing blank line to log
call error_handler(E_MSG, 'log_qcf_table_data:', 'Logging the data in the QCF Table', source)

! Write the table headers to the dart_log and terminal
write(log_msg, '(A4, A6, A9, A)') header1(:)
call error_handler(E_MSG, 'log_qcf_table_data:', trim(log_msg), source)

write(log_msg,'(3A14, 2A12, 3(A10, 2A14, 2A12), A12, 2A14, 2A12)') header2(:)
call error_handler(E_MSG, 'log_qcf_table_data:', trim(log_msg), source)

! Write the table data to the dart_log and terminal
do row = 1, size(qcf_table_data)
   write(log_msg, *) trim(specified_qtys(row)), qcf_table_data(row)%obs_error_info%bounded_below, qcf_table_data(row)%obs_error_info%bounded_above, &
               qcf_table_data(row)%obs_error_info%lower_bound, qcf_table_data(row)%obs_error_info%upper_bound, trim(dist_type_string_probit_inflation(row)), &
               qcf_table_data(row)%probit_inflation%bounded_below, qcf_table_data(row)%probit_inflation%bounded_above, &
               qcf_table_data(row)%probit_inflation%lower_bound, qcf_table_data(row)%probit_inflation%upper_bound, trim(dist_type_string_probit_state(row)), &
               qcf_table_data(row)%probit_state%bounded_below, qcf_table_data(row)%probit_state%bounded_above, &
               qcf_table_data(row)%probit_state%lower_bound, qcf_table_data(row)%probit_state%upper_bound, trim(dist_type_string_probit_extended_state(row)), &
               qcf_table_data(row)%probit_extended_state%bounded_below, qcf_table_data(row)%probit_extended_state%bounded_above, &
               qcf_table_data(row)%probit_extended_state%lower_bound, qcf_table_data(row)%probit_extended_state%upper_bound, &
               trim(filter_kind_string(row)), qcf_table_data(row)%obs_inc_info%bounded_below, qcf_table_data(row)%obs_inc_info%bounded_above, &
               qcf_table_data(row)%obs_inc_info%lower_bound, qcf_table_data(row)%obs_inc_info%upper_bound
call error_handler(E_MSG, 'log_qcf_table_data:', trim(log_msg), source)
end do

call error_handler(E_MSG, '', '', source) !Writing blank line to log

end subroutine log_qcf_table_data

!------------------------------------------------------------------------


subroutine end_algorithm_info_mod()

if (.not. module_initialized) return
module_initialized = .false.

if (use_qty_defaults) return

deallocate(specified_qtys)
deallocate(qcf_table_data)

end subroutine end_algorithm_info_mod

!----------------------------------------------------------------------

end module algorithm_info_mod
