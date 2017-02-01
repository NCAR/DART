! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module model_mod

use        types_mod, only : r8
use time_manager_mod, only : time_type, set_time
use     location_mod, only : location_type, set_location, get_location, &
                             LocationDims, LocationName, LocationLName, &
                             get_close_maxdist_init, get_close_obs_init, get_close_obs

use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG, nmlfileunit, &
                             do_output, find_namelist_in_file, check_namelist_read,     &
                             do_nml_file, do_nml_term

use   random_seq_mod, only : random_seq_type, init_random_seq, random_gaussian

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
          pert_model_state, &
          get_close_maxdist_init, get_close_obs_init, get_close_obs, ens_mean_for_model

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

! Basic model parameters controlled by nameslist; have defaults

!---------------------------------------------------------------
! Namelist with default values
!
integer  :: num_state_vars = 40
real(r8) :: forcing    = 8.00_r8
real(r8) :: delta_t    = 0.05_r8
integer  :: time_step_days = 0
integer  :: time_step_seconds = 3600
logical  :: reset_forcing = .false.
real(r8) :: random_forcing_amplitude = 0.0_r8

namelist /model_nml/ num_state_vars, forcing, delta_t, time_step_days, &
   time_step_seconds, reset_forcing, random_forcing_amplitude
!----------------------------------------------------------------

! Define the location of the state variables in module storage
type(location_type), allocatable :: state_loc(:)
type(time_type) :: time_step
integer :: model_size

! Adding random noise
type(random_seq_type) :: random



contains

!==================================================================



subroutine static_init_model()
!------------------------------------------------------------------
! Initializes class data for this model. For now, simply outputs the
! identity info, sets the location of the state variables, and initializes
! the time type for the time stepping (is this general enough for time???)

real(r8) :: x_loc
integer  :: i, iunit, io

! Print module information to log file and stdout.
call register_module(source, revision, revdate)

! Read the namelist entry
call find_namelist_in_file("input.nml", "model_nml", iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, "model_nml")

! Model size is twice the number of state_vars
model_size = 2 * num_state_vars

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)

! Create storage for locations
allocate(state_loc(model_size))

! Define the locations of the model state variables
do i = 1, num_state_vars
   x_loc = (i - 1.0_r8) / num_state_vars
   state_loc(i) =  set_location(x_loc)
   ! Forcing is at same location as corresponding state variable
   state_loc(i + num_state_vars) = set_location(x_loc)
end do

! The time_step in terms of a time type must also be initialized. Need
! to determine appropriate non-dimensionalization conversion for L96 from
! Shree Khare.
time_step = set_time(time_step_seconds, time_step_days)

! Initialize the random sequence
call init_random_seq(random)

end subroutine static_init_model



subroutine comp_dt(x, dt)
!------------------------------------------------------------------
! subroutine comp_dt(x, dt)
! 
! Computes the time tendency of the lorenz 1996 model given current state

real(r8), intent( in) ::  x(:)
real(r8), intent(out) :: dt(:)

integer :: j, jp1, jm1, jm2

do j = 1, num_state_vars
   jp1 = j + 1
   if(jp1 > num_state_vars) jp1 = 1
   jm2 = j - 2
   if(jm2 < 1) jm2 = num_state_vars + jm2
   jm1 = j - 1
   if(jm1 < 1) jm1 = num_state_vars
   
   dt(j) = (x(jp1) - x(jm2)) * x(jm1) - x(j) + x(num_state_vars + j)
end do

! Time tendency for the forcing variables
dt(num_state_vars + 1 : model_size) = 0.0_r8


! Try adding in some random spread; fixed across the instances (basically a global var)
if(.not. reset_forcing) &
   dt(num_state_vars + 1 : model_size) = &
      random_gaussian(random, 0.0_r8, random_forcing_amplitude)

! Try adding in some random noise to each forcing variable, completely local
!if(.not. reset_forcing) then
!   do j = num_state_vars + 1, model_size
!      dt(j) = random_gaussian(random, 0.0_r8, random_forcing_amplitude)
!   end do
!endif 
   





end subroutine comp_dt



subroutine init_conditions(x)
!------------------------------------------------------------------
! subroutine init_conditions(x)
!
! Initial conditions for lorenz 96
! It is assumed that this is called before any other routines in this
! module. Should probably make that more formal and perhaps enforce for
! more comprehensive models.


real(r8), intent(out) :: x(:)

! Set state variable to the value of forcing with 1st element slightly perturbed
! Set forcing parameters to 8.0 if being assimilated
x    = forcing
x(1) = 1.001_r8 * forcing

end subroutine init_conditions



subroutine adv_1step(x, time)
!------------------------------------------------------------------
! subroutine adv_1step(x, time)
!
! Does single time step advance for lorenz 96 model
! using four-step rk time step
! The Time argument is needed for compatibility with more complex models
! that need to know the time to compute their time tendency and is not
! used in L96. Is there a better way to do this in F90 than to just hang
! this argument out everywhere?


real(r8), intent(inout) :: x(:)
type(time_type), intent(in) :: time

real(r8), dimension(size(x)) :: x1, x2, x3, x4, dx, inter

! If reset forcing is set, grab the value from the namelist and hold it fixed
if(reset_forcing) x(num_state_vars + 1: model_size) = forcing

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



function get_model_size()
!------------------------------------------------------------------
! function get_model_size()
!
! Returns size of model

integer :: get_model_size

! Number of state variables is really twice model size if forcing is assimed
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



subroutine model_interpolate(x, location, itype, obs_val, istatus)
!------------------------------------------------------------------
!
! Interpolates from state vector x to the location. It's not particularly
! happy dumping all of this straight into the model. Eventually some
! concept of a grid underlying models but above locations is going to
! be more general. May want to wait on external infrastructure projects
! for this?

! Argument itype is not used here because there is only one type of variable.
! Type is needed to allow swap consistency with more complex models.

! Don't allow direct interpolation of the parameter value at this point
! although this could be added in the future if useful for testing.


real(r8),            intent(in) :: x(:)
type(location_type), intent(in) :: location
integer,             intent(in) :: itype
real(r8),           intent(out) :: obs_val
integer,            intent(out) :: istatus

integer :: lower_index, upper_index, i
real(r8) :: lctn, lctnfrac

! All forward operators supported
istatus = 0

! Convert location to real
lctn = get_location(location)
! Multiply by model size assuming domain is [0, 1] cyclic
lctn = num_state_vars * lctn

lower_index = int(lctn) + 1
upper_index = lower_index + 1
if(lower_index > num_state_vars) lower_index = lower_index - num_state_vars
if(upper_index > num_state_vars) upper_index = upper_index - num_state_vars

lctnfrac = lctn - int(lctn)
obs_val = (1.0_r8 - lctnfrac) * x(lower_index) + lctnfrac * x(upper_index)

if(1 == 1) return


! All the stuff below is for strange forward operator tests; not currently used
!-----------------------------------------------
!!!obs_val = obs_val ** 2
!!!if(1 == 1) return

! Temporarily add on an observation from the other side of the domain, too
lower_index = lower_index + model_size / 2
if(lower_index > model_size) lower_index = lower_index - model_size
upper_index = upper_index + model_size / 2
if(upper_index > model_size) upper_index = upper_index - model_size
obs_val = obs_val + &
   lctnfrac * x(lower_index) + (1.0_r8 - lctnfrac) * x(upper_index)
if(1 == 1) return


! Next one does an average over a range of points
obs_val = 0.0_r8
lower_index = lower_index - 7
upper_index = upper_index - 7
if(lower_index < 1) lower_index = lower_index + model_size
if(upper_index < 1) upper_index = upper_index + model_size

do i = 1, 15
   if(lower_index > model_size) lower_index = lower_index - model_size
   if(upper_index > model_size) upper_index = upper_index - model_size
   obs_val = obs_val + &
      (1.0_r8 - lctnfrac) * x(lower_index) + lctnfrac * x(upper_index)
   lower_index = lower_index + 1
   upper_index = upper_index + 1
end do

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



subroutine get_state_meta_data(index_in, location, var_type)
!------------------------------------------------------------------
!
! Given an integer index into the state vector structure, returns the
! associated location. This is not a function because the more general
! form of the call has a second intent(out) optional argument kind.
! Maybe a functional form should be added?


integer,             intent(in)  :: index_in
type(location_type), intent(out) :: location
integer,             intent(out), optional :: var_type

location = state_loc(index_in)
if (present(var_type)) then
   if(index_in <= num_state_vars) then
      var_type = 1    ! default variable type
   else
      ! If forcing parameter is being assimilated, it has a var_type of 2
      var_type = 2
   endif 
endif

end subroutine get_state_meta_data



subroutine end_model()
!------------------------------------------------------------------
!
! Does any shutdown and clean-up needed for model. Nothing for L96 for now.


end subroutine end_model



function nc_write_model_atts( ncFileID ) result (ierr)
!------------------------------------------------------------------
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

use typeSizes
use netcdf

integer, intent(in)  :: ncFileID      ! netCDF file identifier
integer              :: ierr          ! return value of function

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

integer             :: i
type(location_type) :: lctn 
ierr = 0                      ! assume normal termination

!--------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file 
!--------------------------------------------------------------------

call check(nf90_Inquire(ncFileID, nDimensions, nVariables, nAttributes, unlimitedDimID))
call check(nf90_sync(ncFileID)) ! Ensure netCDF file is current
call check(nf90_Redef(ncFileID))

!--------------------------------------------------------------------
! Determine ID's from stuff already in the netCDF file
!--------------------------------------------------------------------

! make sure time is unlimited dimid

call check(nf90_inq_dimid(ncFileID,"copy",dimid=MemberDimID))
call check(nf90_inq_dimid(ncFileID,"time",dimid=TimeDimID))

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
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model", "Forced_Lorenz_96"))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_forcing", forcing ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_delta_t", delta_t ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_num_state_vars", num_state_vars ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_time_step_days", time_step_days ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_time_step_seconds", time_step_seconds ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_random_forcing_amplitude", random_forcing_amplitude ))
if ( reset_forcing ) then
   call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_reset_forcing", 'TRUE' ))
else
   call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_reset_forcing", 'FALSE' ))
endif

!--------------------------------------------------------------------
! Define the model size, state variable dimension ... whatever ...
!--------------------------------------------------------------------

call check(nf90_def_dim(ncid=ncFileID, name="StateVariable", &
                        len=model_size, dimid = StateVarDimID)) 

!--------------------------------------------------------------------
! Define the Location Variable and add Attributes
! Some of the atts come from location_mod (via the USE: stmnt)
! CF standards for Locations:
! http://www.cgd.ucar.edu/cms/eaton/netcdf/CF-working.html#ctype
!--------------------------------------------------------------------

call check(NF90_def_var(ncFileID, name=trim(adjustl(LocationName)), xtype=nf90_double, &
              dimids = StateVarDimID, varid=LocationVarID) )
call check(nf90_put_att(ncFileID, LocationVarID, "long_name", trim(adjustl(LocationLName))))
call check(nf90_put_att(ncFileID, LocationVarID, "dimension", LocationDims ))
call check(nf90_put_att(ncFileID, LocationVarID, "units", "nondimensional"))
call check(nf90_put_att(ncFileID, LocationVarID, "valid_range", (/ 0.0_r8, 1.0_r8 /)))

!--------------------------------------------------------------------
! Define either the "state vector" variables -OR- the "prognostic" variables.
!--------------------------------------------------------------------
!if ( output_state_vector ) then

   ! Define the state vector coordinate variable
   call check(nf90_def_var(ncid=ncFileID,name="StateVariable", xtype=nf90_int, &
              dimids=StateVarDimID, varid=StateVarVarID))
   call check(nf90_put_att(ncFileID, StateVarVarID, "long_name", "State Variable ID"))
   call check(nf90_put_att(ncFileID, StateVarVarID, "units",     "indexical") )
   call check(nf90_put_att(ncFileID, StateVarVarID, "valid_range", (/ 1, model_size /)))

   ! Define the actual state vector
   call check(nf90_def_var(ncid=ncFileID, name="state", xtype=nf90_double, &
              dimids = (/ StateVarDimID, MemberDimID, TimeDimID /), varid=StateVarID))
   call check(nf90_put_att(ncFileID, StateVarID, "long_name", "model state or fcopy"))

   ! Leave define mode so we can fill
   call check(nf90_enddef(ncfileID))

   ! Fill the state variable coordinate variable
   call check(nf90_put_var(ncFileID, StateVarVarID, (/ (i,i=1,model_size) /) ))

!This is not ready yet ... the model state is replicated. This is an attempt
! to parse it back into two pieces -- just don't know what to call them, etc. 
!else
!
!   call check(nf90_def_var(ncid=ncFileID, name="X", xtype=nf90_double, &
!        dimids = (/ , VelJDimID, levDimID, MemberDimID, unlimitedDimID /), &
!        varid  = uVarID))
!   call check(nf90_put_att(ncFileID, uVarID, "long_name", "model component"))
!
!   call check(nf90_def_var(ncid=ncFileID, name="XF", xtype=nf90_real, &
!        dimids = (/ VelIDimID, VelJDimID, levDimID, MemberDimID, unlimitedDimID /), &
!        varid  = vVarID))
!   call check(nf90_put_att(ncFileID, vVarID, "long_name", "estimate of forcing"))
!
!endif

!--------------------------------------------------------------------
! Fill the location variable
!--------------------------------------------------------------------

do i = 1,model_size
   call get_state_meta_data(i,lctn)
   call check(nf90_put_var(ncFileID, LocationVarID, get_location(lctn), (/ i /) ))
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



function nc_write_model_vars( ncFileID, statevec, copyindex, timeindex ) result (ierr)         
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

use typeSizes
use netcdf

integer,                intent(in) :: ncFileID      ! netCDF file identifier
real(r8), dimension(:), intent(in) :: statevec
integer,                intent(in) :: copyindex
integer,                intent(in) :: timeindex
integer                            :: ierr          ! return value of function

!--------------------------------------------------------------------
! General netCDF variables
!--------------------------------------------------------------------

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID
integer :: StateVarID

!--------------------------------------------------------------------
! local variables
!--------------------------------------------------------------------

ierr = 0                      ! assume normal termination

!--------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file
!--------------------------------------------------------------------

call check(nf90_Inquire(ncFileID, nDimensions, nVariables, nAttributes, unlimitedDimID))

! no matter the value of "output_state_vector", we only do one thing.
! if ( output_state_vector ) then

   call check(NF90_inq_varid(ncFileID, "state", StateVarID) )
   call check(NF90_put_var(ncFileID, StateVarID, statevec,  &
             start=(/ 1, copyindex, timeindex /)))

! else
!----------------------------------------------------------------------------
! Fill the variables
! TemperatureGrid : surface pressure  Var%ps(tis:tie, tjs:tje)
!                 : temperature       Var%t (tis:tie, tjs:tje, klb:kub)
!                 : tracers           Var%r (tis:tie, tjs:tje, klb:kub, 1:vars%ntrace)
! VelocityGrid    : u                 Var%u (vis:vie, vjs:vje, klb:kub)
!                 : v                 Var%v (vis:vie, vjs:tje, klb:kub)
!----------------------------------------------------------------------------
!
!   x = statevec ! Unfortunately, have to explicity cast it ...
!                ! the filter uses a type=double,
!                ! the vector_to_prog_var function expects a single.
!   call vector_to_prog_var(x, get_model_size(), Var)
!
!   call check(NF90_inq_varid(ncFileID, "ps", psVarID))
!   call check(nf90_put_var( ncFileID, psVarID, Var%ps(tis:tie, tjs:tje), &
!        start=(/ 1, 1, copyindex, timeindex /)))
!
!   call check(NF90_inq_varid(ncFileID, "t", tVarID))
!   call check(nf90_put_var( ncFileID, tVarID, Var%t( tis:tie, tjs:tje, klb:kub), &
!        start=(/ 1, 1, 1, copyindex, timeindex /)))
!
! endif



! write (*,*)'Finished filling variables ...'
call check(nf90_sync(ncFileID))
! write (*,*)'netCDF file is synched ...'

contains
   ! Internal subroutine - checks error status after each netcdf, prints
   !                       text message each time an error code is returned.
   subroutine check(istatus)
   integer, intent ( in) :: istatus
      if(istatus /= nf90_noerr) call error_handler(E_ERR,'nc_write_model_vars',&
         trim(nf90_strerror(istatus)), source, revision, revdate)
   end subroutine check
end function nc_write_model_vars



subroutine pert_model_state(state, pert_state, interf_provided)
!------------------------------------------------------------------
! subroutine pert_model_state(state, pert_state, interf_provided)
!
! Perturbs a model state for generating initial ensembles
! Returning interf_provided means go ahead and do this with uniform
! small independent perturbations.

real(r8), intent(in)  :: state(:)
real(r8), intent(in)  :: pert_state(:)
logical,  intent(out) :: interf_provided

interf_provided = .false.

end subroutine pert_model_state


subroutine ens_mean_for_model(ens_mean)
!------------------------------------------------------------------
! Not used in low-order models

real(r8), intent(in) :: ens_mean(:)

end subroutine ens_mean_for_model



!===================================================================
! End of model_mod
!===================================================================
end module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
