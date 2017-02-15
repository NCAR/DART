! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module model_mod

! Implements Lorenz 84 3 variable model with intermediate level
! attractor.

use        types_mod,      only : r8, i8, i4

use time_manager_mod,      only : time_type, set_time

use     location_mod,      only : location_type, set_location, get_location, &
                                  LocationDims, LocationName, LocationLName, &
                                  get_close_maxdist_init, get_close_obs_init, &
                                  loc_get_close_obs => get_close_obs, get_close_type

use    utilities_mod,      only : register_module, error_handler, E_ERR, E_MSG, nmlfileunit, &
                                  do_output, find_namelist_in_file, check_namelist_read,     &
                                  do_nml_file, do_nml_term, nc_check

use         obs_kind_mod,  only : RAW_STATE_VARIABLE

use ensemble_manager_mod,  only : ensemble_type

use distributed_state_mod, only : get_state

use state_structure_mod,   only : add_domain, state_structure_info

use dart_time_io_mod,      only : read_model_time, write_model_time


implicit none
private

public :: get_model_size, &
          adv_1step, &
          get_state_meta_data, &
          model_interpolate, &
          get_model_time_step, &
          end_model, &
          static_init_model, &
          init_time, &
          init_conditions, &
          nc_write_model_atts, &
          nc_write_model_vars, &
          pert_model_copies, &
          get_close_maxdist_init, &
          get_close_obs_init, &
          get_close_obs, &
          vert_convert, &
          query_vert_localization_coord, &
          read_model_time, &
          write_model_time


! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

!  define model parameters

! Model size is fixed
integer(i8), parameter :: model_size = 3

!-------------------------------------------------------------
! Namelist with default values
!
real(r8) ::      a = 0.25_r8
real(r8) ::      b = 4.00_r8
real(r8) ::      f = 8.00_r8 
real(r8) ::      g = 1.25_r8
real(r8) :: deltat = 0.01_r8
integer  :: time_step_days = 0
integer  :: time_step_seconds = 3600

namelist /model_nml/ a, b, f, g, deltat, time_step_days, time_step_seconds
!---------------------------------------------------------------

! Define the location of the state variables in module storage
type(location_type) :: state_loc(model_size)
type(time_type)     :: time_step


contains

!==================================================================



subroutine static_init_model()
!------------------------------------------------------------------
! Initializes class data for L63 model and outputs I.D.
!

real(r8) :: x_loc
integer(i8) :: i
integer  :: iunit, io, dom_id

! Print module information to log file and stdout.
call register_module(source, revision, revdate)

! Read the namelist entry
call find_namelist_in_file("input.nml", "model_nml", iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, "model_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)

! Define the locations of the model state variables
do i = 1, model_size
   x_loc = (i - 1.0_r8) / model_size
   state_loc(i) =  set_location(x_loc)
end do

! The time_step in terms of a time type must also be initialized.
time_step = set_time(time_step_seconds, time_step_days)

! Tell the DART I/O routines how large the model data is so they
! can read/write it.
dom_id = add_domain(model_size)

end subroutine static_init_model



subroutine init_model_instance()
!------------------------------------------------------------------
! subroutine init_model_instance
!
! Initializes instance dependent state for model. Null for L84.

end subroutine init_model_instance



subroutine comp_dt(x, dt)
!------------------------------------------------------------------
! subroutine comp_dt(x, dt)
!
! compute time tendency of lorenz-84 3-variable model given current state

real(r8), intent( in)  ::  x(:)
real(r8), intent(out) :: dt(:)

!  compute lorenz 84 model time tendency

dt(1) = -x(2)**2 -x(3)**2 - a*x(1) + a*f
dt(2) = x(1)*x(2) - b*x(1)*x(3) - x(2) + g
dt(3) = b*x(1)*x(2) + x(1)*x(3) - x(3)

end subroutine comp_dt



subroutine init_conditions(x)
!------------------------------------------------------------------
! subroutine init_conditions(x)
!
!  defines off-attractor initial conditions for lorenz-84

real(r8), intent(out) :: x(:)

x(1) = 0.0_r8
x(2) = 1.0_r8
x(3) = 0.0_r8

end subroutine init_conditions



subroutine adv_1step(x, time)
!------------------------------------------------------------------
! subroutine adv_1step(x, time)
!
! does single time step advance 

real(r8), intent(inout) :: x(:)
type(time_type), intent(in) :: time

real(r8) :: fract

fract = 1.0_r8
call adv_single(x, fract)

return
end  subroutine adv_1step



subroutine adv_single(x, fract)
!------------------------------------------------------------------
! subroutine adv_single(x, fract)
!
!  does single time step advance for lorenz 3 variable model using two step rk


real(r8), intent(inout) :: x(:)
real(r8), intent(in)    :: fract

real(r8) :: x1(3), x2(3), dx(3)

call comp_dt(x, dx)            !  compute the first intermediate step
x1 = x + fract * deltat * dx

call comp_dt(x1, dx)           !  compute the second intermediate step
x2 = x1 + fract * deltat * dx

!  new value for x is average of original value and second intermediate

x = (x + x2) / 2.0_r8

return
end subroutine adv_single



function get_model_size()
!------------------------------------------------------------------
! function get_model_size()
!
! Returns size of model

integer(i8) :: get_model_size

get_model_size = model_size

end function get_model_size



subroutine init_time(time)
!------------------------------------------------------------------
!
! Gets the initial time for a state from the model. Where should this info
! come from in the most general case?

type(time_type), intent(out) :: time

! For now, just set to 0
time = set_time(0, 0)

end subroutine init_time



subroutine model_interpolate(state_handle, ens_size, location, itype, expected_obs, istatus)

!------------------------------------------------------------------
!
! Interpolates from state vector x to the location. It's not particularly
! happy dumping all of this straight into the model. Eventually some
! concept of a grid underlying models but above locations is going to
! be more general. May want to wait on external infrastructure projects
! for this?

! Argument itype is not used here because there is only one type of variable.
! Type is needed to allow swap consistency with more complex models.

type(ensemble_type),  intent(in) :: state_handle
integer,              intent(in) :: ens_size
type(location_type),  intent(in) :: location
integer,              intent(in) :: itype
real(r8),            intent(out) :: expected_obs(ens_size)
integer,             intent(out) :: istatus(ens_size)

integer(i8) :: lower_index, upper_index
real(r8) :: lctn, lctnfrac

! All interpolations okay for now
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
expected_obs(:) = (1.0_r8 - lctnfrac) * get_state(lower_index, state_handle) + &
                            lctnfrac  * get_state(upper_index, state_handle)

end subroutine model_interpolate



function get_model_time_step()
!------------------------------------------------------------------
! function get_model_time_step()
!
! Returns the the time step of the model. In the long run should be repalced
! by a more general routine that returns details of a general time-stepping
! capability.

type(time_type) :: get_model_time_step

get_model_time_step = time_step

end function get_model_time_step



subroutine get_state_meta_data(state_handle, index_in, location, var_type)

!------------------------------------------------------------------
!
! Given an integer index into the state vector structure, returns the
! associated location. This is not a function because the more general
! form of the call has a second intent(out) optional argument kind.
! Maybe a functional form should be added?


type(ensemble_type), intent(in)  :: state_handle !< some large models need this
integer(i8),         intent(in)  :: index_in
type(location_type), intent(out) :: location
integer,             intent(out), optional :: var_type

location = state_loc(index_in)
if (present(var_type)) var_type = RAW_STATE_VARIABLE    ! default variable type

end subroutine get_state_meta_data


!------------------------------------------------------------------
! Does any shutdown and clean-up needed for model. Nothing for L84 for now.

subroutine end_model()


end subroutine end_model



!------------------------------------------------------------------
function nc_write_model_atts( ncFileID, model_mod_writes_state_variables ) result (ierr)

use netcdf

integer, intent(in)  :: ncFileID      ! netCDF file identifier
logical, intent(out) :: model_mod_writes_state_variables
integer              :: ierr          ! return value of function

!--------------------------------------------------------------------
! General netCDF variables
!--------------------------------------------------------------------

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID

!--------------------------------------------------------------------
! netCDF variables for Location
!--------------------------------------------------------------------

integer :: LocationVarID, LocationDimID

!--------------------------------------------------------------------
! local variables
!--------------------------------------------------------------------

character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic
character(len=NF90_MAX_NAME) :: str1

integer             :: i
type(location_type) :: lctn 
ierr = 0                      ! assume normal termination

model_mod_writes_state_variables = .false.

!--------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file 
!--------------------------------------------------------------------

call check(nf90_Inquire(ncFileID, nDimensions, nVariables, nAttributes, unlimitedDimID))
call check(nf90_sync(ncFileID)) ! Ensure netCDF file is current
call check(nf90_Redef(ncFileID))

!--------------------------------------------------------------------
! Write Global Attributes 
!--------------------------------------------------------------------

call DATE_AND_TIME(crdate,crtime,crzone,values)
write(str1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
                  values(1), values(2), values(3), values(5), values(6), values(7)

call check(nf90_put_att(ncFileID, NF90_GLOBAL, "creation_date",str1))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_source", source ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revision", revision ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revdate", revdate ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model", "Lorenz_84"))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_a", a ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_b", b ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_f", f ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_g", g ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_deltat", deltat ))

!--------------------------------------------------------------------
! Define the model size, state variable dimension ... whatever ...
!--------------------------------------------------------------------

call nc_check(nf90_def_dim(ncid=ncFileID, name="location", &
                           len=int(model_size, i4), dimid = LocationDimID), &
                          'nc_write_model_atts', 'def_dim location')

!--------------------------------------------------------------------
! Define the Location Variable and add Attributes
! Some of the atts come from location_mod (via the USE: stmnt)
! CF standards for Locations:
! http://www.cgd.ucar.edu/cms/eaton/netcdf/CF-working.html#ctype
!--------------------------------------------------------------------

call check(NF90_def_var(ncFileID, name="location", xtype=nf90_double, &
              dimids = LocationDimID, varid=LocationVarID) )
call check(nf90_put_att(ncFileID, LocationVarID, "long_name", trim(adjustl(LocationLName))))
call check(nf90_put_att(ncFileID, LocationVarID, "dimension", LocationDims ))
call check(nf90_put_att(ncFileID, LocationVarID, "units", "nondimensional"))
call check(nf90_put_att(ncFileID, LocationVarID, "valid_range", (/ 0.0_r8, 1.0_r8 /)))

! Leave define mode so we can fill
call check(nf90_enddef(ncfileID))

!--------------------------------------------------------------------
! Fill the location variable
!--------------------------------------------------------------------

do i = 1,model_size
   lctn = state_loc(i)
   call nc_check(nf90_put_var(ncFileID, LocationVarID, get_location(lctn), (/ int(i,i4) /) ), 'nc_write_model_atts', 'nf90_put_var LocationVarID 2')
enddo

!--------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!--------------------------------------------------------------------
call check(nf90_sync(ncFileID))

write (*,*)'Model attributes written, netCDF file synched ...'

contains
   ! Internal subroutine - checks error status after each netcdf, prints 
   !                       text message each time an error code is returned. 
   subroutine check(istatus)
      integer, intent ( in) :: istatus
      if(istatus /= nf90_noerr) call error_handler(E_ERR,'nc_write_model_atts',&
         trim(nf90_strerror(istatus)), source, revision, revdate)
   end subroutine check
end function nc_write_model_atts

!------------------------------------------------------------------

function nc_write_model_vars( ncFileID, statevec, copyindex, timeindex ) result (ierr)         

use netcdf

integer,                intent(in) :: ncFileID      ! netCDF file identifier
real(r8), dimension(:), intent(in) :: statevec
integer,                intent(in) :: copyindex
integer,                intent(in) :: timeindex
integer                            :: ierr          ! return value of function

ierr = 0                      ! assume normal termination

end function nc_write_model_vars


  subroutine linearize(nl, l)
!-------------------------------------------------------------------
! subroutine linearize(nl, l)
!
!  compute linear operator around state nl

real(r8) :: nl(3), l(3, 3)

l(1, 1) = -1.0_r8 *     a     * deltat + 1.0_r8
l(1, 2) = -1.0_r8 * nl(2)     * deltat
l(1, 3) = -1.0_r8 * nl(3)     * deltat
l(2, 1) = (nl(2) - b * nl(3)) * deltat
l(2, 2) = (-1.0_r8 + nl(1))   * deltat + 1.0_r8
l(2, 3) = -1.0_r8 * b * nl(1) * deltat
l(3, 1) = (b * nl(2) + nl(3)) * deltat
l(3, 2) =  b * nl(1)          * deltat
l(3, 3) = (nl(1) - 1.0_r8)    * deltat + 1.0_r8

return
end subroutine linearize

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

!> Unused in this model.

subroutine vert_convert(state_handle, location, obs_kind, istatus)

type(ensemble_type), intent(in)  :: state_handle
type(location_type), intent(in)  :: location
integer,             intent(in)  :: obs_kind
integer,             intent(out) :: istatus

istatus = 0

end subroutine vert_convert

!--------------------------------------------------------------------

!> pass the vertical localization coordinate to assim_tools_mod

function query_vert_localization_coord()

integer :: query_vert_localization_coord

!> @TODO should define some parameters including something
!> like HAS_NO_VERT for this use.

query_vert_localization_coord = -1

end function query_vert_localization_coord

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

!--------------------------------------------------------------------

!===================================================================
! End of model_mod
!===================================================================
end module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
