! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module algorithm_info_mod

! Provides routines that give information about details of algorithms for 
! observation error sampling, observation increments, and the transformations
! for regression and inflation in probit space. 

use types_mod,       only : r8, i8, MISSING_R8, obstypelength

use obs_def_mod,     only : obs_def_type, get_obs_def_type_of_obs, get_obs_def_error_variance
use obs_kind_mod,    only : get_quantity_for_type_of_obs, get_name_for_quantity, get_index_for_quantity

use utilities_mod,   only : error_handler, E_ERR, E_MSG, open_file, close_file, to_upper, &
                            do_nml_file, do_nml_term, nmlfileunit, check_namelist_read, &
                            find_namelist_in_file

use assim_model_mod, only : get_state_meta_data
use location_mod,    only : location_type

use distribution_params_mod, only : NORMAL_DISTRIBUTION, BOUNDED_NORMAL_RH_DISTRIBUTION, &
   GAMMA_DISTRIBUTION, BETA_DISTRIBUTION, LOG_NORMAL_DISTRIBUTION, UNIFORM_DISTRIBUTION,  &
   PARTICLE_FILTER_DISTRIBUTION

implicit none
private

character(len=512) :: errstring
character(len=*), parameter :: source = 'algorithm_info_mod.f90'

logical :: module_initialized = .false.
logical :: use_qty_defaults = .true.

! Defining parameters for different observation space filters
integer, parameter :: EAKF               = 1
integer, parameter :: ENKF               = 2
integer, parameter :: KERNEL             = 3
integer, parameter :: OBS_PARTICLE       = 4
integer, parameter :: UNBOUNDED_RHF      = 8
integer, parameter :: GAMMA_FILTER       = 11
integer, parameter :: BOUNDED_NORMAL_RHF = 101 

public :: obs_error_info, probit_dist_info, obs_inc_info, &
          init_algorithm_info_mod, end_algorithm_info_mod, &
          EAKF, ENKF, BOUNDED_NORMAL_RHF, UNBOUNDED_RHF, &
          GAMMA_FILTER, KERNEL, OBS_PARTICLE

! type definitions for the QCF table
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

integer, allocatable :: specified_qtys(:)
type(algorithm_info_type), allocatable :: qceff_table_data(:)

character(len=129), allocatable :: dist_type_string_probit_inflation(:)
character(len=129), allocatable :: dist_type_string_probit_state(:)
character(len=129), allocatable :: dist_type_string_probit_extended_state(:)
character(len=129), allocatable :: filter_kind_string(:)

! namelist
character(len=129) :: qceff_table_filename = ''

namelist /algorithm_info_nml/ qceff_table_filename

contains

!-------------------------------------------------------------------------


subroutine init_algorithm_info_mod()

! Gets number of lines/QTYs in the QCF table, allocates space for the table data


integer :: fileid
integer :: io, iunit

integer :: numrows
integer :: nlines

if (module_initialized) return
module_initialized = .true.

! Read the namelist entry
call find_namelist_in_file("input.nml", "algorithm_info_nml", iunit)
read(iunit, nml = algorithm_info_nml, iostat = io)
call check_namelist_read(iunit, io, "algorithm_info_nml")

if (do_nml_file()) write(nmlfileunit, nml=algorithm_info_nml)
if (do_nml_term()) write(     *     , nml=algorithm_info_nml)


if (qceff_table_filename == '') then
   write(errstring,*) 'No QCF table file listed in namelist, using default values for all QTYs'
   call error_handler(E_MSG, 'init_algorithm_info_mod:', errstring, source)
   return
endif

use_qty_defaults = .false.
fileid = open_file(trim(qceff_table_filename), 'formatted', 'read')

! Do loop to get number of rows (or QTYs) in the table
nlines = 0
do
  read(fileid,*,iostat=io)
  if(io/=0) exit
  nlines = nlines + 1
end do

call close_file(fileid)

numrows = nlines - HEADER_LINES

allocate(specified_qtys(numrows))
allocate(qceff_table_data(numrows))
allocate(dist_type_string_probit_inflation(numrows))
allocate(dist_type_string_probit_state(numrows))
allocate(dist_type_string_probit_extended_state(numrows))
allocate(filter_kind_string(numrows))

call read_qceff_table(qceff_table_filename)
call assert_qceff_table_version()
call verify_qceff_table_data()
call log_qceff_table_data()

end subroutine init_algorithm_info_mod

!------------------------------------------------------------------------

 
subroutine read_qceff_table(qceff_table_filename)

! Reads in the QCEFF input options from tabular data file

character(len=129), intent(in) :: qceff_table_filename

integer :: fileid
integer :: row
character(len=obstypelength) :: qty_string

if (.not. module_initialized) call init_algorithm_info_mod()

fileid = open_file(trim(qceff_table_filename), 'formatted', 'read')

! skip the headers
read(fileid, *) header1
read(fileid, *) header2

! read in table values directly to qceff_table_data type
do row = 1, size(qceff_table_data)
   read(fileid, *) qty_string, qceff_table_data(row)%obs_error_info%bounded_below, qceff_table_data(row)%obs_error_info%bounded_above, &
                   qceff_table_data(row)%obs_error_info%lower_bound, qceff_table_data(row)%obs_error_info%upper_bound, dist_type_string_probit_inflation(row), &
                   qceff_table_data(row)%probit_inflation%bounded_below, qceff_table_data(row)%probit_inflation%bounded_above, &
                   qceff_table_data(row)%probit_inflation%lower_bound, qceff_table_data(row)%probit_inflation%upper_bound, dist_type_string_probit_state(row), &
                   qceff_table_data(row)%probit_state%bounded_below, qceff_table_data(row)%probit_state%bounded_above, &
                   qceff_table_data(row)%probit_state%lower_bound, qceff_table_data(row)%probit_state%upper_bound, dist_type_string_probit_extended_state(row), &
                   qceff_table_data(row)%probit_extended_state%bounded_below, qceff_table_data(row)%probit_extended_state%bounded_above, &
                   qceff_table_data(row)%probit_extended_state%lower_bound, qceff_table_data(row)%probit_extended_state%upper_bound, &
                   filter_kind_string(row), qceff_table_data(row)%obs_inc_info%bounded_below, qceff_table_data(row)%obs_inc_info%bounded_above, &
                   qceff_table_data(row)%obs_inc_info%lower_bound, qceff_table_data(row)%obs_inc_info%upper_bound
   
   call to_upper(qty_string)
   specified_qtys(row) = get_index_for_quantity(qty_string)  

   if(specified_qtys(row) == -1) then
      write(errstring,*) trim(qty_string), ' is not a valid DART QTY'
      call error_handler(E_ERR, 'read_qceff_table:', errstring, source)
   endif

   ! Converting the distribution types (read in from table as a string) to its corresponding int value
   call to_upper(dist_type_string_probit_inflation(row))

   select case (trim(dist_type_string_probit_inflation(row)))
      case ('NORMAL_DISTRIBUTION')
         qceff_table_data(row)%probit_inflation%dist_type = NORMAL_DISTRIBUTION
      case ('BOUNDED_NORMAL_RH_DISTRIBUTION')
         qceff_table_data(row)%probit_inflation%dist_type = BOUNDED_NORMAL_RH_DISTRIBUTION
      case ('GAMMA_DISTRIBUTION')
         qceff_table_data(row)%probit_inflation%dist_type = GAMMA_DISTRIBUTION
      case ('BETA_DISTRIBUTION')
         qceff_table_data(row)%probit_inflation%dist_type = BETA_DISTRIBUTION
      case ('LOG_NORMAL_DISTRIBUTION')
         qceff_table_data(row)%probit_inflation%dist_type = LOG_NORMAL_DISTRIBUTION
      case ('UNIFORM_DISTRIBUTION')
         qceff_table_data(row)%probit_inflation%dist_type = UNIFORM_DISTRIBUTION
      !!!case ('PARTICLE_FILTER_DISTRIBUTION')
         !!!qceff_table_data(row)%probit_inflation%dist_type = PARTICLE_FILTER_DISTRIBUTION
      case default
         write(errstring, *) 'Invalid distribution type for probit inflation: ', trim(dist_type_string_probit_inflation(row))
         call error_handler(E_ERR, 'read_qceff_table:', errstring, source)
   end select


   call to_upper(dist_type_string_probit_state(row))
  
   select case (trim(dist_type_string_probit_state(row)))
      case ('NORMAL_DISTRIBUTION')
         qceff_table_data(row)%probit_state%dist_type = NORMAL_DISTRIBUTION
      case ('BOUNDED_NORMAL_RH_DISTRIBUTION')
         qceff_table_data(row)%probit_state%dist_type = BOUNDED_NORMAL_RH_DISTRIBUTION
      case ('GAMMA_DISTRIBUTION')
         qceff_table_data(row)%probit_state%dist_type = GAMMA_DISTRIBUTION
      case ('BETA_DISTRIBUTION')
         qceff_table_data(row)%probit_state%dist_type = BETA_DISTRIBUTION
      case ('LOG_NORMAL_DISTRIBUTION')
         qceff_table_data(row)%probit_state%dist_type = LOG_NORMAL_DISTRIBUTION
      case ('UNIFORM_DISTRIBUTION')
         qceff_table_data(row)%probit_state%dist_type = UNIFORM_DISTRIBUTION
      !!!case ('PARTICLE_FILTER_DISTRIBUTION')
         !!!qceff_table_data(row)%probit_state%dist_type = PARTICLE_FILTER_DISTRIBUTION
      case default
         write(errstring, *) 'Invalid distribution type for probit state: ', trim(dist_type_string_probit_state(row))
         call error_handler(E_ERR, 'read_qceff_table:', errstring, source)
   end select

   call to_upper(dist_type_string_probit_extended_state(row))

   select case (trim(dist_type_string_probit_extended_state(row)))
      case ('NORMAL_DISTRIBUTION')
         qceff_table_data(row)%probit_extended_state%dist_type = NORMAL_DISTRIBUTION
      case ('BOUNDED_NORMAL_RH_DISTRIBUTION')
         qceff_table_data(row)%probit_extended_state%dist_type = BOUNDED_NORMAL_RH_DISTRIBUTION
      case ('GAMMA_DISTRIBUTION')
         qceff_table_data(row)%probit_extended_state%dist_type = GAMMA_DISTRIBUTION
      case ('BETA_DISTRIBUTION')
         qceff_table_data(row)%probit_extended_state%dist_type = BETA_DISTRIBUTION
      case ('LOG_NORMAL_DISTRIBUTION')
         qceff_table_data(row)%probit_extended_state%dist_type = LOG_NORMAL_DISTRIBUTION
      case ('UNIFORM_DISTRIBUTION')
         qceff_table_data(row)%probit_extended_state%dist_type = UNIFORM_DISTRIBUTION
      !!!case ('PARTICLE_FILTER_DISTRIBUTION')
         !!!qceff_table_data(row)%probit_extended_state%dist_type = PARTICLE_FILTER_DISTRIBUTION
      case default
         write(errstring, *) 'Invalid distribution type for probit extended state: ', trim(dist_type_string_probit_extended_state(row))
         call error_handler(E_ERR, 'read_qceff_table:', errstring, source)
   end select


   ! Converting the filter kind (read in from table as a string) to its corresponding int value
   call to_upper(filter_kind_string(row))
 
   select case (trim(filter_kind_string(row)))
      case ('EAKF')
         qceff_table_data(row)%obs_inc_info%filter_kind = EAKF
      case ('ENKF')
         qceff_table_data(row)%obs_inc_info%filter_kind = ENKF
      case ('UNBOUNDED_RHF')
         qceff_table_data(row)%obs_inc_info%filter_kind = UNBOUNDED_RHF
      case ('GAMMA_FILTER')
         qceff_table_data(row)%obs_inc_info%filter_kind = GAMMA_FILTER
      case ('BOUNDED_NORMAL_RHF')
         qceff_table_data(row)%obs_inc_info%filter_kind = BOUNDED_NORMAL_RHF
      case default
         write(errstring, *) 'Invalid filter kind: ', trim(filter_kind_string(row))
         call error_handler(E_ERR, 'read_qceff_table:', errstring, source)
   end select

end do

call close_file(fileid)

end subroutine read_qceff_table

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

if (.not. module_initialized) call init_algorithm_info_mod()

! Get the type of the observation
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

!use default values if qceff_table_filename is not in namelist
if (use_qty_defaults) then
   bounded_below = .false.;    bounded_above = .false.
   lower_bound   = MISSING_R8; upper_bound   = MISSING_R8
   return
endif

!find location of QTY in qceff_table_data structure
QTY_loc = findloc(specified_qtys, obs_qty)

if (QTY_loc(1) == 0) then
   !use default values if QTY is not in table
   bounded_below = .false.;    bounded_above = .false.
   lower_bound   = MISSING_R8; upper_bound   = MISSING_R8

   else
      bounded_below = qceff_table_data(QTY_loc(1))%obs_error_info%bounded_below
      bounded_above = qceff_table_data(QTY_loc(1))%obs_error_info%bounded_above
      lower_bound   = qceff_table_data(QTY_loc(1))%obs_error_info%lower_bound
      upper_bound   = qceff_table_data(QTY_loc(1))%obs_error_info%upper_bound

endif

end subroutine obs_error_info

!-------------------------------------------------------------------------


subroutine probit_dist_info(qty, is_state, is_inflation, dist_type, &
   bounded_below, bounded_above, lower_bound, upper_bound)

integer,  intent(in)  :: qty
logical,  intent(in)  :: is_state      ! True for state variable, false for obs
logical,  intent(in)  :: is_inflation  ! True for inflation transform
integer,  intent(out) :: dist_type
logical,  intent(out) :: bounded_below, bounded_above
real(r8), intent(out) :: lower_bound,   upper_bound

integer :: QTY_loc(1)

! Have input information about the qty of the state or observation being transformed
! along with additional logical info that indicates whether this is an observation
! or state variable and about whether the transformation is being done for inflation
! or for regress.
 
! Selects the appropriate transform, which is specified in the QCF input table per QTY. 
! At present, the options are NORMAL_DISTRIBUTION, BOUNDED_NORMAL_RH_DISTRIBUTION, 
! GAMMA_DISTRIBUTION, BETA_DISTRIBUTION, LOG_NORMAL_DISTRIBUTION,
! UNIFORM_DISTRIBUTION
! If the BNRH is selected then information about the bounds must also be set.
! For example, if qty corresponds to a sea ice fraction then an appropriate choice
! would be:
! bounded_below = .true.;  bounded_above = .true.
! lower_bound  = 0.0_r8;   upper_bounds  = 1.0_r8

! In the long run, may not have to have separate controls for each of the input possibilities
! However, for now these are things that need to be explored for science understanding

if (.not. module_initialized) call init_algorithm_info_mod()

!use default values if qceff_table_filename is not in namelist
if (use_qty_defaults) then
   dist_type = NORMAL_DISTRIBUTION
   bounded_below = .false.;    bounded_above = .false.
   lower_bound   = MISSING_R8; upper_bound   = MISSING_R8
   return
endif

QTY_loc = findloc(specified_qtys, qty)

if (QTY_loc(1) == 0) then
   !use default values if QTY is not in table
   dist_type = NORMAL_DISTRIBUTION
   bounded_below = .false.;    bounded_above = .false.
   lower_bound   = MISSING_R8; upper_bound   = MISSING_R8

   elseif(is_inflation) then
   ! Case for inflation transformation

      dist_type     = qceff_table_data(QTY_loc(1))%probit_inflation%dist_type
      bounded_below = qceff_table_data(QTY_loc(1))%probit_inflation%bounded_below
      bounded_above = qceff_table_data(QTY_loc(1))%probit_inflation%bounded_above
      lower_bound   = qceff_table_data(QTY_loc(1))%probit_inflation%lower_bound
      upper_bound   = qceff_table_data(QTY_loc(1))%probit_inflation%upper_bound

   elseif(is_state) then
   ! Case for state variable priors

      dist_type     = qceff_table_data(QTY_loc(1))%probit_state%dist_type
      bounded_below = qceff_table_data(QTY_loc(1))%probit_state%bounded_below
      bounded_above = qceff_table_data(QTY_loc(1))%probit_state%bounded_above
      lower_bound   = qceff_table_data(QTY_loc(1))%probit_state%lower_bound
      upper_bound   = qceff_table_data(QTY_loc(1))%probit_state%upper_bound

   else
   ! This case is for observation (extended state) priors

      dist_type     = qceff_table_data(QTY_loc(1))%probit_extended_state%dist_type
      bounded_below = qceff_table_data(QTY_loc(1))%probit_extended_state%bounded_below
      bounded_above = qceff_table_data(QTY_loc(1))%probit_extended_state%bounded_above
      lower_bound   = qceff_table_data(QTY_loc(1))%probit_extended_state%lower_bound
      upper_bound   = qceff_table_data(QTY_loc(1))%probit_extended_state%upper_bound

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

if (.not. module_initialized) call init_algorithm_info_mod()

!use default values if qceff_table_filename is not in namelist
if (use_qty_defaults) then
   filter_kind = EAKF
   bounded_below = .false.;    bounded_above = .false.
   lower_bound   = MISSING_R8; upper_bound   = MISSING_R8
   return
endif

!find location of QTY in qceff_table_data structure
QTY_loc = findloc(specified_qtys, obs_qty)

if (QTY_loc(1) == 0) then
   !use default values if QTY is not in table
   filter_kind = EAKF
   bounded_below = .false.;    bounded_above = .false.
   lower_bound   = MISSING_R8; upper_bound   = MISSING_R8

   else

      filter_kind   = qceff_table_data(QTY_loc(1))%obs_inc_info%filter_kind
      bounded_below = qceff_table_data(QTY_loc(1))%obs_inc_info%bounded_below
      bounded_above = qceff_table_data(QTY_loc(1))%obs_inc_info%bounded_above
      lower_bound   = qceff_table_data(QTY_loc(1))%obs_inc_info%lower_bound
      upper_bound   = qceff_table_data(QTY_loc(1))%obs_inc_info%upper_bound

endif

end subroutine obs_inc_info

!------------------------------------------------------------------------


subroutine assert_qceff_table_version()

! Subroutine to ensure the correct version of the QCF table is being used

character(1), parameter :: QCF_VERSION = '1'

if (trim(header1(4)) /= QCF_VERSION) then
   write(errstring,*) 'Using outdated/incorrect version of the QCF table'
   call error_handler(E_ERR, 'assert_qceff_table_version:', errstring, source)
endif

end subroutine assert_qceff_table_version

!------------------------------------------------------------------------


subroutine verify_qceff_table_data()

! Subroutine to ensure that the data in the QCF table is valid 

integer :: row

if (use_qty_defaults) return

!Checks that all bounds are valid; currently checks that the lower bound in less than the upper
!Here we could add more specific checks if we have known limits on the bounds
do row = 1, size(qceff_table_data)

   if (qceff_table_data(row)%obs_error_info%bounded_below .and.  qceff_table_data(row)%obs_error_info%bounded_above) then
      if(qceff_table_data(row)%obs_error_info%lower_bound > qceff_table_data(row)%obs_error_info%upper_bound) then
         write(errstring,*) 'Invalid bounds in obs_error_info'
         call error_handler(E_ERR, 'verify_qceff_table_data:', errstring, source)
      endif
   endif
   if (qceff_table_data(row)%probit_inflation%bounded_below .and. qceff_table_data(row)%probit_inflation%bounded_above) then
      if(qceff_table_data(row)%probit_inflation%lower_bound > qceff_table_data(row)%probit_inflation%upper_bound) then
         write(errstring,*) 'Invalid bounds in probit_inflation'
         call error_handler(E_ERR, 'verify_qceff_table_data:', errstring, source)
      endif
   endif
   if(qceff_table_data(row)%probit_state%bounded_below .and. qceff_table_data(row)%probit_state%bounded_above) then
      if(qceff_table_data(row)%probit_state%lower_bound > qceff_table_data(row)%probit_state%upper_bound) then
         write(errstring,*) 'Invalid bounds in probit_state'
         call error_handler(E_ERR, 'verify_qceff_table_data:', errstring, source)
      endif
   endif
   if(qceff_table_data(row)%probit_extended_state%bounded_below .and. qceff_table_data(row)%probit_extended_state%bounded_above) then
      if(qceff_table_data(row)%probit_extended_state%lower_bound > qceff_table_data(row)%probit_extended_state%upper_bound) then
         write(errstring,*) 'Invalid bounds in probit_extended_state'
         call error_handler(E_ERR, 'verify_qceff_table_data:', errstring, source)
      endif
   endif
   if(qceff_table_data(row)%obs_inc_info%bounded_below .and. qceff_table_data(row)%obs_inc_info%bounded_above) then
       if(qceff_table_data(row)%obs_inc_info%lower_bound > qceff_table_data(row)%obs_inc_info%upper_bound) then
          write(errstring,*) 'Invalid bounds in obs_inc_info'
          call error_handler(E_ERR, 'verify_qceff_table_data:', errstring, source)
       endif
    endif
end do


!Ensures that there are no duplicate QTYs in the table
do row = 1, size(qceff_table_data)
   if(count(specified_qtys==specified_qtys(row)) > 1) then
      write(errstring,*) trim(get_name_for_quantity(specified_qtys(row))), ' has multiple entries in the table'
      call error_handler(E_ERR, 'verify_qceff_table_data:', errstring, source)
   endif
end do

end subroutine verify_qceff_table_data

!------------------------------------------------------------------------


subroutine log_qceff_table_data()

! Subroutine to write the data in QCF table to dart_log
character(len=2000) :: log_msg
integer :: row

if (use_qty_defaults) return

call error_handler(E_MSG, '', '', source) !Writing blank line to log
call error_handler(E_MSG, 'log_qceff_table_data:', 'Logging the data in the QCEFF Table', source)

! Write the table headers to the dart_log and terminal
write(log_msg, '(2A6, A9, A)') header1(:)
call error_handler(E_MSG, 'log_qceff_table_data:', trim(log_msg), source)

write(log_msg,'(A10, 2A14, 2A12, 3(A10, 2A14, 2A12), A12, 2A14, 2A12)') header2(:)
call error_handler(E_MSG, 'log_qceff_table_data:', trim(log_msg), source)

! Write the table data to the dart_log and terminal
do row = 1, size(qceff_table_data)
   write(log_msg, *) trim(get_name_for_quantity(specified_qtys(row))), qceff_table_data(row)%obs_error_info%bounded_below, qceff_table_data(row)%obs_error_info%bounded_above, &
               qceff_table_data(row)%obs_error_info%lower_bound, qceff_table_data(row)%obs_error_info%upper_bound, trim(dist_type_string_probit_inflation(row)), &
               qceff_table_data(row)%probit_inflation%bounded_below, qceff_table_data(row)%probit_inflation%bounded_above, &
               qceff_table_data(row)%probit_inflation%lower_bound, qceff_table_data(row)%probit_inflation%upper_bound, trim(dist_type_string_probit_state(row)), &
               qceff_table_data(row)%probit_state%bounded_below, qceff_table_data(row)%probit_state%bounded_above, &
               qceff_table_data(row)%probit_state%lower_bound, qceff_table_data(row)%probit_state%upper_bound, trim(dist_type_string_probit_extended_state(row)), &
               qceff_table_data(row)%probit_extended_state%bounded_below, qceff_table_data(row)%probit_extended_state%bounded_above, &
               qceff_table_data(row)%probit_extended_state%lower_bound, qceff_table_data(row)%probit_extended_state%upper_bound, &
               trim(filter_kind_string(row)), qceff_table_data(row)%obs_inc_info%bounded_below, qceff_table_data(row)%obs_inc_info%bounded_above, &
               qceff_table_data(row)%obs_inc_info%lower_bound, qceff_table_data(row)%obs_inc_info%upper_bound
call error_handler(E_MSG, 'log_qceff_table_data:', trim(log_msg), source)
end do

call error_handler(E_MSG, '', '', source) !Writing blank line to log

end subroutine log_qceff_table_data

!------------------------------------------------------------------------


subroutine end_algorithm_info_mod()

if (.not. module_initialized) return
module_initialized = .false.

if (use_qty_defaults) return

deallocate(specified_qtys)
deallocate(qceff_table_data)

end subroutine end_algorithm_info_mod

!----------------------------------------------------------------------

end module algorithm_info_mod
