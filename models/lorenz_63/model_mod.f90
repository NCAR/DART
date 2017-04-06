! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

!> Lorenz 63 model interfaces to DART

module model_mod

! Revised assim_model version of Lorenz-63 3-variable model

use        types_mod,      only : r8, i8, i4

use time_manager_mod,      only : time_type, set_time

use     location_mod,      only : location_type, set_location, get_location, &
                                  LocationDims, LocationName, LocationLName

use    utilities_mod,      only : register_module, error_handler, E_ERR, E_MSG, nmlfileunit, &
                                  do_output, find_namelist_in_file, check_namelist_read,     &
                                  do_nml_file, do_nml_term, nc_check

use         obs_kind_mod,  only : RAW_STATE_VARIABLE

use ensemble_manager_mod,  only : ensemble_type

use distributed_state_mod, only : get_state

use state_structure_mod,   only : add_domain

use dart_time_io_mod,      only : read_model_time, write_model_time

use default_model_mod,     only : end_model, pert_model_copies, vert_convert, &
                                  query_vert_localization_coord, get_close_obs, &
                                  get_close_maxdist_init, get_close_obs_init, &
                                  get_close_state_init, get_close_state, &
                                  get_close_type, nc_write_model_vars

implicit none
private

!>@todo list real routines first, passthrus at end,
!> or two public lists w/ !default comments
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
          get_close_state_init, &
          get_close_state, &
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

! Model size is fixed for Lorenz-63
integer(i8), parameter :: model_size = 3

!-------------------------------------------------------------
! Namelist with default values

real(r8) ::  sigma = 10.0_r8
real(r8) ::      r = 28.0_r8
real(r8) ::      b = 8.0_r8 / 3.0_r8
real(r8) :: deltat = 0.01_r8
integer  :: time_step_days = 0
integer  :: time_step_seconds = 3600

namelist /model_nml/ sigma, r, b, deltat, time_step_days, time_step_seconds

!---------------------------------------------------------------

! Define the location of the state variables in module storage
type(location_type) :: state_loc(model_size)
type(time_type)     :: time_step


contains

!==================================================================



!------------------------------------------------------------------
!> Initializes class data for L63 model 

subroutine static_init_model()

real(r8) :: x_loc
integer  :: i, iunit, io, dom_id

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

!------------------------------------------------------------------
!> Computes time tendency of the lorenz 1963 3-variable model given 
!> current state

subroutine comp_dt(x, dt)

real(r8), intent( in) ::  x(:)
real(r8), intent(out) :: dt(:)

! compute the lorenz model dt from standard equations

dt(1) = sigma * (x(2) - x(1))
dt(2) = -1.0_r8*x(1)*x(3) + r*x(1) - x(2)
dt(3) = x(1)*x(2) - b*x(3)

end subroutine comp_dt


!------------------------------------------------------------------
!>  off-attractor initial conditions for lorenz 63

subroutine init_conditions(x)

real(r8), intent(out) :: x(:)

! Initial conditions that move nicely onto attractor
x = 0.10_r8

end subroutine init_conditions


!------------------------------------------------------------------
!> does single time step advance for lorenz convective 3 variable model
!> using two step rk time step

subroutine adv_1step(x, time)

real(r8), intent(inout) :: x(:)
type(time_type), intent(in) :: time

real(r8) :: fract

fract = 1.0_r8
call adv_single(x, fract)

end  subroutine adv_1step


!------------------------------------------------------------------
!> does single time step advance for lorenz convective 3 variable model
!> using two step rk time step

subroutine adv_single(x, fract)

real(r8), intent(inout) :: x(:)
real(r8), intent(in)    :: fract

real(r8) :: x1(3), x2(3), dx(3)

call comp_dt(x, dx)            !  compute the first intermediate step
x1 = x + fract * deltat * dx

call comp_dt(x1, dx)           !  compute the second intermediate step
x2 = x1 + fract * deltat * dx

!  new value for x is average of original value and second intermediate

x = (x + x2) / 2.0_r8

end subroutine adv_single


!------------------------------------------------------------------
!> Returns size of model

function get_model_size()

integer :: get_model_size

get_model_size = model_size

end function get_model_size


!------------------------------------------------------------------
!> Sets the initial time for a state from the model.

subroutine init_time(time)

type(time_type), intent(out) :: time

! Set to 0
time = set_time(0, 0)

end subroutine init_time


!------------------------------------------------------------------
!> Interpolates from state vector x to the location.
!>
!> Argument itype is not used here because there is only one type of variable.

subroutine model_interpolate(state_handle, ens_size, location, itype, expected_obs, istatus)

type(ensemble_type),  intent(in) :: state_handle
integer,              intent(in) :: ens_size
type(location_type),  intent(in) :: location
integer,              intent(in) :: itype
real(r8),            intent(out) :: expected_obs(ens_size)
integer,             intent(out) :: istatus(ens_size)

integer(i8)  :: lower_index, upper_index
real(r8) :: lctn, lctnfrac

! All obs okay for now
istatus = 0

! Convert location to real
lctn = get_location(location)
! Multiply by model size assuming domain is [0, 1] cyclic
lctn = model_size * lctn

lower_index = int(lctn) + 1
upper_index = lower_index + 1
if(lower_index > model_size) lower_index = lower_index - model_size
if(upper_index > model_size) upper_index = upper_index - model_size

lctnfrac = lctn - int(lctn)
expected_obs = (1.0_r8 - lctnfrac) * get_state(lower_index, state_handle) + lctnfrac * get_state(upper_index, state_handle)

end subroutine model_interpolate


!------------------------------------------------------------------
!> Returns the mininum time step of the model.

function get_model_time_step()

type(time_type) :: get_model_time_step

get_model_time_step = time_step

end function get_model_time_step

!------------------------------------------------------------------
!> Given an integer index into the state vector structure, returns the
!> associated location.

subroutine get_state_meta_data(state_handle, index_in, location, var_type)

type(ensemble_type), intent(in)  :: state_handle !< some large models need this
integer(i8),         intent(in)  :: index_in
type(location_type), intent(out) :: location
integer,             intent(out), optional :: var_type

location = state_loc(index_in)
if (present(var_type)) var_type = RAW_STATE_VARIABLE    ! default variable type

end subroutine get_state_meta_data

!------------------------------------------------------------------
!> old version of linearized lorenz 63 model time tendency computation

subroutine linear_dt(x, dx, dt)
   
real(r8), intent(in)  :: x(:), dx(:)
real(r8), intent(out) :: dt(:)

!  compute linear model lorenz time tendency

dt(1) = -1.0_r8 * sigma * dx(1) + sigma*dx(2)
dt(2) = (r - x(3))*dx(1) - dx(2) - x(1)*dx(3)
dt(3) = x(2)*dx(1) + x(1)*dx(2) - b*dx(3)

end subroutine linear_dt

!------------------------------------------------------------------
!> does single time step advance for lorenz convective 3 variable model
!> using four step rk time step

subroutine adv_single_rk4(x, fract)

real(r8), intent(inout) :: x(:)
real(r8), intent(in)    :: fract

real(r8) :: x1(3), x2(3), x3(3), x4(3), dx(3), inter(3)

call comp_dt(x, dx)         !  compute the first intermediate step
x1    = fract * deltat * dx
inter = x + x1 / 2.0_r8

call comp_dt(inter, dx)     !  compute the second intermediate step
x2    = fract * deltat * dx
inter = x + x2 / 2.0_r8

call comp_dt(inter, dx)     !  compute the third intermediate step
x3    = fract * deltat * dx
inter = x + x3

call comp_dt(inter, dx)     !  compute fourth intermediate step
x4 = fract * deltat * dx

!  compute new value for x

x = x + x1/6.0_r8 + x2/3.0_r8 + x3/3.0_r8 + x4/6.0_r8

end subroutine adv_single_rk4

!------------------------------------------------------------------
!>  compute inv linear model lorenz time tendency (see notes 13mar94)
!>  for now assumes stupid leap frog, will this be sufficient?

subroutine inv_linear_dt(x, dx, px)

real(r8), intent(in)  :: x(:), dx(:)
real(r8), intent(out) :: px(3)

real(r8) a(3, 3), fact, tdx(3)

a(1, 1) = -sigma * deltat + 1.0_r8
a(1, 2) =  sigma * deltat
a(1, 3) = 0.0_r8

a(2, 1) = (r - x(3)) * deltat
a(2, 2) = -1.0_r8 * deltat + 1.0_r8
a(2, 3) = -x(1) * deltat

a(3, 1) =  x(2) * deltat
a(3, 2) =  x(1) * deltat
a(3, 3) =    -b * deltat + 1.0_r8
   
!  initialize copy of dx

  call error_handler(E_ERR,'inv_linear_dt', 'this routine is not up to date', source, revision, revdate)

!      tdx(i) = dx(i)
   tdx = dx

!  get rid of a(2, 3)

fact = a(2, 3) / a(3, 3)
   a(2, :) = a(2, :) - fact * a(3, :)
tdx(2) = tdx(2) - fact * tdx(3)

!  get rid of a(1, 2)

fact = a(1, 2) / a(2, 2)
   a(1, :) = a(1, :) - fact * a(2, :)
tdx(1) = tdx(1) - fact * tdx(2)

!  solve for the previous step linear perturbation

px(1) = tdx(1) / a(1, 1)
px(2) = (tdx(2) - a(2, 1) * px(1)) / a(2, 2)
px(3) = (tdx(3) - a(3, 1) * px(1) - a(3, 2) * px(2)) / a(3, 3)

end subroutine inv_linear_dt


!------------------------------------------------------------------
!> compute linear operator around state nl

subroutine linearize(nl, l)

real(r8), intent(in)  :: nl(3)
real(r8), intent(out) :: l(3, 3)

l(1, 1) = -1.0_r8 * sigma * deltat + 1.0_r8
l(1, 2) =           sigma * deltat
l(1, 3) =          0.0_r8 * deltat

l(2, 1) = (r - nl(3)) * deltat
l(2, 2) = -1.0_r8 * deltat + 1.0_r8
l(2, 3) = -1.0_r8 * nl(1) * deltat

l(3, 1) =       nl(2) * deltat
l(3, 2) =       nl(1) * deltat
l(3, 3) = -1.0_r8 * b * deltat + 1.0_r8

end subroutine linearize


!------------------------------------------------------------------
!> Writes the model-specific attributes to a netCDF file
!> For the lorenz_63 model, each state variable is at a separate location.

function nc_write_model_atts( ncFileID, model_mod_writes_state_variables ) result (ierr)

use typeSizes
use netcdf

integer, intent(in)  :: ncFileID      ! netCDF file identifier
logical, intent(out) :: model_mod_writes_state_variables
integer              :: ierr          ! return value of function

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID
integer :: LocationDimID, LocationVarID
integer :: StateVarID, MemberDimID, TimeDimID
integer :: i
type(location_type) :: lctn 

ierr = 0                             ! assume normal termination

! dart code will write the state into the file
! so this routine just needs to write any model-specific
! attributes it wants to record.

model_mod_writes_state_variables = .false.


! make sure ncFileID refers to an open netCDF file 

! are the inq and sync needed?  can we just redef here?
call nc_check(nf90_Inquire(ncFileID, nDimensions, nVariables, nAttributes, unlimitedDimID), &
                           'nc_write_model_atts', 'nf90_Inquire')
call nc_check(nf90_sync(ncFileID), 'nc_write_model_atts', 'nf90_sync') ! Ensure netCDF file is current
call nc_check(nf90_Redef(ncFileID), 'nc_write_model_atts', 'nf90_Redef')


! Write Global Attributes 

call put_global_creation_time(ncFileID)

call put_global_char_att(ncFileID, "model_source", source )
call put_global_char_att(ncFileID, "model_revision", revision )
call put_global_char_att(ncFileID, "model_revdate", revdate )
call put_global_char_att(ncFileID, "model", "Lorenz_63")
call put_global_real_att(ncFileID, "model_r", r )
call put_global_real_att(ncFileID, "model_b", b )
call put_global_real_att(ncFileID, "model_sigma", sigma )
call put_global_real_att(ncFileID, "model_deltat", deltat )


!--------------------------------------------------------------------
!--------------------------------------------------------------------
!>@todo
!> this same location attr and write code is in the oned locations module.
!> this code should call it to fill in the locations attrs and data.

! call nc_write_location_atts(ncFileID)
! call nc_get_location_varids(ncFileID, LocationVarID)
! do i=1, model_size
!    call nc_write_location(ncFileID, LocationVarID, loc, index)
! enddo

call nc_check(nf90_def_dim(ncid=ncFileID, name="location", len=int(model_size,i4), dimid = LocationDimID), &
              'nc_write_model_atts', 'nf90_def_dim location') 
call nc_check(NF90_def_var(ncFileID, name="location", xtype=nf90_double, dimids = LocationDimID, varid=LocationVarID) , &
             'nc_write_model_atts', 'nf90_def_var location')

call nc_check(nf90_put_att(ncFileID, LocationVarID, "short_name", trim(adjustl(LocationLName))), &
              'nc_write_model_atts', 'nf90_put_att short_name')
call nc_check(nf90_put_att(ncFileID, LocationVarID, "long_name", "location on a unit circle"), &
              'nc_write_model_atts', 'nf90_put_att long_name')
call nc_check(nf90_put_att(ncFileID, LocationVarID, "dimension", LocationDims ), &
              'nc_write_model_atts', 'nf90_put_att dimension')
call nc_check(nf90_put_att(ncFileID, LocationVarID, "valid_range", (/ 0.0_r8, 1.0_r8 /)), &
              'nc_write_model_atts', 'nf90_put_att valid_range')

! Leave define mode so we can fill
call nc_check(nf90_enddef(ncfileID), 'nc_write_model_atts', 'nf90_enddef')

! Fill the state variable coordinate variable
call nc_check(nf90_put_var(ncFileID, LocationVarID, (/ (i,i=1,int(model_size,i4)) /) ), &
              'nc_write_model_atts', 'nf90_put_var LocationVarID')

! Fill the location variable
do i = 1,model_size
   lctn = state_loc(i)
   call nc_check(nf90_put_var(ncFileID, LocationVarID, get_location(lctn), (/ i /) ), &
                 'nc_write_model_atts', 'nf90_put_var LocationVarID 2')
enddo
!--------------------------------------------------------------------
!--------------------------------------------------------------------

call nc_check(nf90_sync(ncFileID), 'nc_write_model_atts', 'nf90_sync')

end function nc_write_model_atts

!--------------------------------------------------------------------

subroutine put_global_char_att(ncid, name, val)

use netcdf

integer,          intent(in) :: ncid
character(len=*), intent(in) :: name
character(len=*), intent(in) :: val

integer :: ret

ret = nf90_put_att(ncid, NF90_GLOBAL, name, val)
call nc_check(ret, 'put_global_char_att', 'adding the global attribute: '//trim(name))

end subroutine put_global_char_att

!--------------------------------------------------------------------

subroutine put_global_real_att(ncid, name, val)

use netcdf

integer,          intent(in) :: ncid
character(len=*), intent(in) :: name
real(r8),         intent(in) :: val

integer :: ret

ret = nf90_put_att(ncid, NF90_GLOBAL, name, val)
call nc_check(ret, 'put_global_real_att', 'adding the global attribute: '//trim(name))

end subroutine put_global_real_att

!--------------------------------------------------------------------

subroutine put_global_creation_time(ncid)
integer, intent(in) :: ncid

character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic

character(len=128) :: str1

call DATE_AND_TIME(crdate,crtime,crzone,values)
write(str1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
                  values(1), values(2), values(3), values(5), values(6), values(7)

call put_global_char_att(ncid, "creation_date",str1)

end subroutine put_global_creation_time

!===================================================================
! End of model_mod
!===================================================================
end module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
