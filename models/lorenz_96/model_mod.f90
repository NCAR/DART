! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module model_mod

use types_mod,             only : r8, i8, i4

use time_manager_mod,      only : time_type, set_time

use location_mod,          only : location_type, set_location, get_location,  &
                                  LocationDims, LocationName, LocationLName,  &
                                  get_close_maxdist_init, get_close_obs_init, &
                                  get_close_type, loc_get_close_obs => get_close_obs

use utilities_mod,         only : register_module, do_nml_file, do_nml_term,    &
                                  nmlfileunit, find_namelist_in_file,           &
                                  check_namelist_read, nc_check

use         obs_kind_mod,  only : RAW_STATE_VARIABLE

use ensemble_manager_mod,  only : ensemble_type

use distributed_state_mod, only : get_state

use state_structure_mod,   only : add_domain

use dart_time_io_mod,      only : read_model_time, write_model_time

implicit none
private

public :: get_model_size, &
          adv_1step, &
          get_state_meta_data, &
          get_model_time_step, &
          end_model, &
          static_init_model, &
          init_time, &
          init_conditions, &
          nc_write_model_atts, &
          nc_write_model_vars, &
          get_close_maxdist_init, get_close_obs_init, get_close_obs, &
          model_interpolate, &
          query_vert_localization_coord, &
          vert_convert, &
          construct_file_name_in, &
          pert_model_copies, &
          read_model_time, &
          write_model_time

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"


!---------------------------------------------------------------

! Namelist with default values

integer(i8) :: model_size = 40
real(r8)    :: forcing    = 8.00_r8
real(r8)    :: delta_t    = 0.05_r8
integer     :: time_step_days = 0
integer     :: time_step_seconds = 3600

namelist /model_nml/ model_size, forcing, delta_t, time_step_days, time_step_seconds

!----------------------------------------------------------------

! module global variables

! Define the location of the state variables in module storage
type(location_type), allocatable :: state_loc(:)
type(time_type) :: time_step


contains

!==================================================================

! model-specific routines


!------------------------------------------------------------------

!> Computes the time tendency of the lorenz 1996 model given current state

subroutine comp_dt(x, dt)

real(r8), intent(in)  ::  x(:)
real(r8), intent(out) :: dt(:)

integer :: j, jp1, jm1, jm2

do j = 1, model_size
   jp1 = j + 1
   if(jp1 > model_size) jp1 = 1
   jm2 = j - 2
   if(jm2 < 1) jm2 = model_size + jm2
   jm1 = j - 1
   if(jm1 < 1) jm1 = model_size
   
   dt(j) = (x(jp1) - x(jm2)) * x(jm1) - x(j) + forcing
end do

end subroutine comp_dt


!------------------------------------------------------------------

!> Do single time step advance for lorenz 96 model
!> using four-step rk time step.
!>
!> the 'time' argument is unused here but used in larger models.

subroutine adv_1step(x, time)

real(r8),        intent(inout) :: x(:)
type(time_type), intent(in)    :: time

real(r8), dimension(size(x)) :: x1, x2, x3, x4, dx, inter

call comp_dt(x, dx)        !  Compute the first intermediate step
x1    = delta_t * dx
inter = x + x1 / 2.0_r8

call comp_dt(inter, dx)    !  Compute the second intermediate step
x2    = delta_t * dx
inter = x + x2 / 2.0_r8

call comp_dt(inter, dx)    !  Compute the third intermediate step
x3    = delta_t * dx
inter = x + x3

call comp_dt(inter, dx)    !  Compute fourth intermediate step
x4 = delta_t * dx

!  Compute new value for x

x = x + x1/6.0_r8 + x2/3.0_r8 + x3/3.0_r8 + x4/6.0_r8

end subroutine adv_1step

!------------------------------------------------------------------

!> This routine is called once and initializes any information needed
!> by other routines in this file.

subroutine static_init_model()

real(r8) :: x_loc
integer  :: i, dom_id

! Do any initial setup needed, including reading the namelist values
call initialize()

! Create storage for locations
allocate(state_loc(model_size))

! Define the locations of the model state variables
do i = 1, model_size
   x_loc = (i - 1.0_r8) / model_size
   state_loc(i) =  set_location(x_loc)
end do

! The time_step in terms of a time type must also be initialized. 
! (Need to determine appropriate non-dimensionalization conversion 
! for L96 from Shree Khare.)
time_step = set_time(time_step_seconds, time_step_days)

! Tell the DART I/O routines how large the model data is so they
! can read/write it.
dom_id = add_domain(model_size)

end subroutine static_init_model

!------------------------------------------------------------------

!> Supply initial conditions for lorenz 96 if not reading restart
!> data from a file.

subroutine init_conditions(x)

real(r8), intent(out) :: x(:)

x    = forcing
x(1) = 1.001_r8 * forcing

end subroutine init_conditions

!------------------------------------------------------------------

!> Sets the initial time for a state from the model if not starting
!> from a restart file.

subroutine init_time(time)

type(time_type), intent(out) :: time

! Day 0, seconds 0
time = set_time(0, 0)

end subroutine init_time

!------------------------------------------------------------------

! Returns size of model

function get_model_size()

integer :: get_model_size

get_model_size = model_size

end function get_model_size

!------------------------------------------------------------------

!> Return the minimum amount of time the model can be advanced by.
!> This is unrelated to the internal RK timesteps.

function get_model_time_step()

type(time_type) :: get_model_time_step

get_model_time_step = time_step

end function get_model_time_step

!------------------------------------------------------------------

!> Given a location, return the interpolated state value.
!> expected_obs() and istatus() are arrays.

subroutine model_interpolate(state_handle, ens_size, location, itype, expected_obs, istatus)

type(ensemble_type),  intent(in) :: state_handle
integer,              intent(in) :: ens_size
type(location_type),  intent(in) :: location
integer,              intent(in) :: itype
real(r8),            intent(out) :: expected_obs(ens_size)
integer,             intent(out) :: istatus(ens_size)

integer(i8) :: lower_index, upper_index, i
real(r8) :: lctn, lctnfrac
real(r8) :: x_lower(ens_size) !< the lower piece of state vector
real(r8) :: x_upper(ens_size) !< the upper piece of state vector

! All forward operators supported
istatus(:) = 0

! Convert location to real
lctn = get_location(location)
! Multiply by model size assuming domain is [0, 1] cyclic
lctn = model_size * lctn

lower_index = int(lctn) + 1
upper_index = lower_index + 1
if(lower_index > model_size) lower_index = lower_index - model_size
if(upper_index > model_size) upper_index = upper_index - model_size

lctnfrac = lctn - int(lctn)

! Grab the correct pieces of state vector

! Lower value
x_lower(:) = get_state(lower_index, state_handle)

! Upper value
x_upper(:) = get_state(upper_index, state_handle)

! calculate the obs value

expected_obs(:) = (1.0_r8 - lctnfrac) * x_lower(:) + lctnfrac * x_upper(:)

end subroutine model_interpolate

!------------------------------------------------------------------

!> @TODO state_handle shouldn't be needed here - IF we can prohibit
!> this routine from using the mean to do vertical conversions.

!> @todo state_loc is state vector size, do we care?
!> all tasks store locations for the entire state vector
!> not just the locations we have.  it could be computed 
!> on demand.

!> Given an integer index into the state vector structure, returns the
!> associated location and variable type.

subroutine get_state_meta_data(state_handle, index_in, location, var_type)

type(ensemble_type), intent(in)  :: state_handle !< some large models need this
integer(i8),         intent(in)  :: index_in
type(location_type), intent(out) :: location
integer,             intent(out), optional :: var_type


location = state_loc(index_in)
if (present(var_type)) var_type = RAW_STATE_VARIABLE 

end subroutine get_state_meta_data

!------------------------------------------------------------------

!> Do any initialization/setup, including reading the
!> namelist values.

subroutine initialize()

integer :: iunit, io

! Print module information
call register_module(source, revision, revdate)

! Read the namelist 
call find_namelist_in_file("input.nml", "model_nml", iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, "model_nml")

! Output the namelist values if requested
if (do_nml_file()) write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)

end subroutine initialize

!------------------------------------------------------------------

! Nothing to do in this model

subroutine end_model()

end subroutine end_model

!------------------------------------------------------------------
! Perturbs a model state copies for generating initial ensembles.
! Routine which could provide a custom perturbation routine to
! generate initial ensembles.  The default (if interface is not
! provided) is to add gaussian noise to each item in the state vector.
subroutine pert_model_copies(state_ens_handle, ens_size, pert_amp, interf_provided)

 type(ensemble_type), intent(inout) :: state_ens_handle
 integer,   intent(in) :: ens_size
 real(r8),  intent(in) :: pert_amp
 logical,  intent(out) :: interf_provided

interf_provided = .false.

end subroutine pert_model_copies

!--------------------------------------------------------------------

!> construct restart file name for reading

function construct_file_name_in(basename, domain, copy)

character(len=512), intent(in) :: basename
integer,            intent(in) :: domain
integer,            intent(in) :: copy
character(len=1024)            :: construct_file_name_in

write(construct_file_name_in, '(A, i4.4)') TRIM(basename), copy

end function construct_file_name_in

!--------------------------------------------------------------------

!> pass the vertical localization coordinate to assim_tools_mod

function query_vert_localization_coord()

integer :: query_vert_localization_coord

!> @TODO should define some parameters including something
!> like HAS_NO_VERT for this use.

query_vert_localization_coord = -1

end function query_vert_localization_coord

!--------------------------------------------------------------------

!> Unused in this model.

subroutine vert_convert(state_handle, location, obs_kind, istatus)

type(ensemble_type), intent(in)  :: state_handle
type(location_type), intent(in)  :: location
integer,             intent(in)  :: obs_kind
integer,             intent(out) :: istatus

istatus = 0

end subroutine vert_convert

!--------------------------------------------------------------------

!> Pass through to the code in the locations module

subroutine get_close_obs(gc, base_obs_loc, base_obs_kind, obs_loc, &
                         obs_kind, num_close, close_ind, dist, state_handle)

type(ensemble_type),         intent(in)     :: state_handle
type(get_close_type),        intent(in)     :: gc
type(location_type),         intent(inout)  :: base_obs_loc, obs_loc(:)
integer,                     intent(in)     :: base_obs_kind, obs_kind(:)
integer,                     intent(out)    :: num_close, close_ind(:)
real(r8),                    intent(out)    :: dist(:)


call loc_get_close_obs(gc, base_obs_loc, base_obs_kind, obs_loc, obs_kind, &
                          num_close, close_ind, dist)

end subroutine get_close_obs

!------------------------------------------------------------------

!> @TODO can we replace most of this with:
!> 
!> call put_standard_attributes()
!> call put_my_model_name("Lorenz_96")
!> call put_model_info("model_forcing", forcing)
!> call put_model_info("model_delta_t", delta_t)


! Can't write netcdf i8 - do we just have 
! a 32 bit version of the netcdf library?
! Bgrid errors out if we have an model size greater than largest i4

! Writes the model-specific attributes to a netCDF file
! TJH Jan 24 2003
!
! TJH 29 July 2003 -- for the moment, all errors are fatal, so the
! return code is always '0 == normal', since the fatal errors stop execution.
!
! For the lorenz_96 model, each state variable is at a separate location.
! that's all the model-specific attributes I can think of ...
!
! assim_model_mod:init_diag_output uses information from the location_mod
!     to define the location dimension and variable ID. All we need to do
!     is query, verify, and fill ...
!
! Typical sequence for adding new dimensions,variables,attributes:
! NF90_OPEN             ! open existing netCDF dataset
!    NF90_redef         ! put into define mode 
!    NF90_def_dim       ! define additional dimensions (if any)
!    NF90_def_var       ! define variables: from name, type, and dims
!    NF90_put_att       ! assign attribute values
! NF90_ENDDEF           ! end definitions: leave define mode
!    NF90_put_var       ! provide values for variable
! NF90_CLOSE            ! close: save updated netCDF dataset

function nc_write_model_atts( ncFileID, model_mod_writes_state_variables ) result (ierr)

use typeSizes
use netcdf

integer, intent(in)  :: ncFileID      ! netCDF file identifier
integer              :: ierr          ! return value of function
logical, intent(out) :: model_mod_writes_state_variables

!--------------------------------------------------------------------
! General netCDF variables
!--------------------------------------------------------------------

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID

!--------------------------------------------------------------------
! netCDF variables for Location
!--------------------------------------------------------------------

integer :: LocationVarID
integer :: StateVarDimID, StateVarVarID
integer :: StateVarID, MemberDimID, TimeDimID

!--------------------------------------------------------------------
! local variables
!--------------------------------------------------------------------

character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic
character(len=NF90_MAX_NAME) :: str1

integer(i8)         :: i
type(location_type) :: lctn 
character(len=128)  :: filename

type(ensemble_type) :: state_handle ! needed for compilation, not used here

ierr = 0                      ! assume normal termination
model_mod_writes_state_variables = .true. 

!--------------------------------------------------------------------
! we only have a netcdf handle here so we do not know the filename
! or the fortran unit number.  but construct a string with at least
! the netcdf handle, so in case of error we can trace back to see
! which netcdf file is involved.
!--------------------------------------------------------------------

write(filename,*) 'ncFileID', ncFileID

!--------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file 
!--------------------------------------------------------------------

call nc_check(nf90_Inquire(ncFileID, nDimensions, nVariables, nAttributes, unlimitedDimID), &
              'nc_write_model_atts', 'inquire, '//trim(filename))
call nc_check(nf90_sync(ncFileID), & ! Ensure netCDF file is current
              'nc_write_model_atts', 'sync, '//trim(filename))
call nc_check(nf90_Redef(ncFileID), &
              'nc_write_model_atts', 'redef, '//trim(filename))

!--------------------------------------------------------------------
! Determine ID's from stuff already in the netCDF file
!--------------------------------------------------------------------

! make sure time is unlimited dimid

call nc_check(nf90_inq_dimid(ncFileID,"copy",dimid=MemberDimID), &
              'nc_write_model_atts', 'inq_dimid copy, '//trim(filename))
call nc_check(nf90_inq_dimid(ncFileID,"time",dimid=TimeDimID), &
              'nc_write_model_atts', 'inq_dimid time, '//trim(filename))

!--------------------------------------------------------------------
! Write Global Attributes 
!--------------------------------------------------------------------
call DATE_AND_TIME(crdate,crtime,crzone,values)
write(str1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
                  values(1), values(2), values(3), values(5), values(6), values(7)

call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "creation_date",str1), &
              'nc_write_model_atts', 'put_att creation_date, '//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_source", source ), &
              'nc_write_model_atts', 'put_att model_source, '//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revision", revision ), &
              'nc_write_model_atts', 'put_att model_revision, '//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revdate", revdate ), &
              'nc_write_model_atts', 'put_att model_revdate, '//trim(filename))

call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model", "Lorenz_96"), &
              'nc_write_model_atts', 'put_att model, '//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_forcing", forcing ), &
              'nc_write_model_atts', 'put_att model_forcing, '//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_delta_t", delta_t ), &
              'nc_write_model_atts', 'put_att model_delta_t, '//trim(filename))

!--------------------------------------------------------------------
! Define the model size, state variable dimension ... whatever ...
!--------------------------------------------------------------------

call nc_check(nf90_def_dim(ncid=ncFileID, name="StateVariable", &
                           len=int(model_size, i4), dimid = StateVarDimID), &
                          'nc_write_model_atts', 'def_dim StateVariable, '//trim(filename))

!--------------------------------------------------------------------
! Define the Location Variable and add Attributes
! Some of the atts come from location_mod (via the USE: stmnt)
! CF standards for Locations:
! http://www.cgd.ucar.edu/cms/eaton/netcdf/CF-working.html#ctype
!--------------------------------------------------------------------

call nc_check(NF90_def_var(ncFileID, name=trim(adjustl(LocationName)), xtype=nf90_double, &
              dimids = StateVarDimID, varid=LocationVarID), &
              'nc_write_model_atts', 'check, '//trim(LocationName)//', '//trim(filename))
call nc_check(nf90_put_att(ncFileID, LocationVarID, "long_name", trim(adjustl(LocationLName))), &
              'nc_write_model_atts', 'put_att long_name, '//trim(filename))
call nc_check(nf90_put_att(ncFileID, LocationVarID, "dimension", LocationDims), &
              'nc_write_model_atts', 'put_att dimension, '//trim(filename))
call nc_check(nf90_put_att(ncFileID, LocationVarID, "units", "nondimensional"), &
              'nc_write_model_atts', 'put_att units, '//trim(filename))
call nc_check(nf90_put_att(ncFileID, LocationVarID, "valid_range", (/ 0.0_r8, 1.0_r8 /)), &
              'nc_write_model_atts', 'put_att valid_range, '//trim(filename))

!--------------------------------------------------------------------
! Define either the "state vector" variables -OR- the "prognostic" variables.
!--------------------------------------------------------------------

! Define the state vector coordinate variable
call nc_check(nf90_def_var(ncid=ncFileID,name="StateVariable", xtype=nf90_int, &
              dimids=StateVarDimID, varid=StateVarVarID), &
             'nc_write_model_atts', 'def_var StateVariable, '//trim(filename))
call nc_check(nf90_put_att(ncFileID, StateVarVarID, "long_name", "State Variable ID"), &
              'nc_write_model_atts', 'put_att long_name, '//trim(filename))
call nc_check(nf90_put_att(ncFileID, StateVarVarID, "units",     "indexical"), &
              'nc_write_model_atts', 'put_att units, '//trim(filename))
call nc_check(nf90_put_att(ncFileID, StateVarVarID, "valid_range", (/ 1, int(model_size, i4) /)), &
              'nc_write_model_atts', 'put_att valid_range, '//trim(filename)) 

! Define the actual state vector
call nc_check(nf90_def_var(ncid=ncFileID, name="state", xtype=nf90_double, &
           dimids = (/ StateVarDimID, MemberDimID, TimeDimID /), varid=StateVarID), &
           'nc_write_model_atts', 'def_var state, '//trim(filename))
call nc_check(nf90_put_att(ncFileID, StateVarID, "long_name", "model state or fcopy"), &
              'nc_write_model_atts', 'put_att long_name, '//trim(filename))

! Leave define mode so we can fill
call nc_check(nf90_enddef(ncfileID), &
              'nc_write_model_atts', 'enddef, '//trim(filename))

! Fill the state variable coordinate variable
call nc_check(nf90_put_var(ncFileID, StateVarVarID, (/ (i,i=1,int(model_size, i4)) /) ), &
              'nc_write_model_atts', 'put_var state variable coordinate, '//trim(filename))

!--------------------------------------------------------------------
! Fill the location variable
!--------------------------------------------------------------------

do i = 1,model_size
   call get_state_meta_data(state_handle, i,lctn)
   call nc_check(nf90_put_var(ncFileID, LocationVarID, get_location(lctn), (/ int(i, i4) /) ), &
              'nc_write_model_atts', 'check locationVarId, '//trim(filename))
enddo

!--------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!--------------------------------------------------------------------
call nc_check(nf90_sync(ncFileID), &
              'nc_write_model_atts', 'sync, '//trim(filename))

! write (*,*)'Model attributes written, netCDF file synched ...'

end function nc_write_model_atts

!------------------------------------------------------------------

! Writes the model-specific attributes to a netCDF file
! TJH 23 May 2003
!
! TJH 29 July 2003 -- for the moment, all errors are fatal, so the
! return code is always '0 == normal', since the fatal errors stop execution.
!
! For the lorenz_96 model, each state variable is at a separate location.
! that's all the model-specific attributes I can think of ...
!
! assim_model_mod:init_diag_output uses information from the location_mod
!     to define the location dimension and variable ID. All we need to do
!     is query, verify, and fill ...
!
! Typical sequence for adding new dimensions,variables,attributes:
! NF90_OPEN             ! open existing netCDF dataset
!    NF90_redef         ! put into define mode
!    NF90_def_dim       ! define additional dimensions (if any)
!    NF90_def_var       ! define variables: from name, type, and dims
!    NF90_put_att       ! assign attribute values
! NF90_ENDDEF           ! end definitions: leave define mode
!    NF90_put_var       ! provide values for variable
! NF90_CLOSE            ! close: save updated netCDF dataset

function nc_write_model_vars( ncFileID, statevec, copyindex, timeindex, write_my_own_variables) result (ierr)

use typeSizes
use netcdf

integer,                intent(in) :: ncFileID      ! netCDF file identifier
real(r8), dimension(:), intent(in) :: statevec
integer,                intent(in) :: copyindex
integer,                intent(in) :: timeindex
logical, optional,     intent(out) :: write_my_own_variables ! model_mod writes variables to 
                                                             ! diagnostic files
integer                            :: ierr          ! return value of function
!--------------------------------------------------------------------
! General netCDF variables
!--------------------------------------------------------------------

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID
integer :: StateVarID

!--------------------------------------------------------------------
! local variables
!--------------------------------------------------------------------

character(len=128) :: filename


if (present(write_my_own_variables)) then
   write_my_own_variables = .false.
   return
endif

!--------------------------------------------------------------------
! we only have a netcdf handle here so we do not know the filename
! or the fortran unit number.  but construct a string with at least
! the netcdf handle, so in case of error we can trace back to see
! which netcdf file is involved.
!--------------------------------------------------------------------

write(filename,*) 'ncFileID', ncFileID


ierr = 0                      ! assume normal termination

!--------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file
!--------------------------------------------------------------------

call nc_check(nf90_Inquire(ncFileID, nDimensions, nVariables, nAttributes, unlimitedDimID), &
              'nc_write_model_vars', 'inquire, '//trim(filename))

! no matter the value of "output_state_vector", we only do one thing.

call nc_check(NF90_inq_varid(ncFileID, "state", StateVarID), & 
              'nc_write_model_vars', 'inq_varid state, '//trim(filename))
call nc_check(NF90_put_var(ncFileID, StateVarID, statevec,  &
              start=(/ 1, copyindex, timeindex /)),  &
              'nc_write_model_vars', 'put_var state vector, '//trim(filename))

! write (*,*)'Finished filling variables ...'
call nc_check(nf90_sync(ncFileID), 'nc_write_model_vars', 'sync, '//trim(filename))
! write (*,*)'netCDF file is synched ...'

end function nc_write_model_vars



!===================================================================
! End of model_mod
!===================================================================
end module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
