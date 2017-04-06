! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module model_mod

! This is the model module for the Lorenz 96 2-scale model, documented in
! Lorenz (1995).  It also has the option of the variant on the model 
! from Smith (2001), and it is invoked by setting local_y = .true. in the 
! namelist.  The time step, coupling, forcing, number of X variables, and the
! number of Ys per X are all specified in the namelist.  Defaults are chosen
! depending on whether the Lorenz or Smith option is specified in the namelist.
! Lorenz is the default model.
!
! May 06 2004, modification of T. Hoar's L96 1-scale model by Josh Hacker

use        types_mod, only : r8
use time_manager_mod, only : time_type, set_time
use     location_mod, only : location_type, set_location, get_location, &
                             LocationDims, LocationName, LocationLName, &
                             get_close_maxdist_init, get_close_obs_init, get_close_obs

use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG, nmlfileunit, &
                             do_output, find_namelist_in_file, check_namelist_read,     &
                             do_nml_file, do_nml_term

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
          pert_model_state, TYPE_X, TYPE_Y, &
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
integer  :: model_size_x = 36
integer  :: y_per_x    = 10
real(r8) :: forcing    = 15.0_r8
real(r8) :: delta_t    = 0.005_r8
real(r8) :: coupling_b = 10.0_r8
real(r8) :: coupling_c = 10.0_r8
real(r8) :: coupling_h = 1.0_r8
logical  :: output_state_vector = .false.
logical  :: local_y = .false.  ! default Lorenz' approach
integer  :: time_step_days = 0
integer  :: time_step_seconds = 3600

namelist /model_nml/ model_size_x, y_per_x, forcing, delta_t, &
                     coupling_b, coupling_c, coupling_h, &
                     output_state_vector, local_y, time_step_days, time_step_seconds
!----------------------------------------------------------------

! Definition of variable types
integer, parameter :: TYPE_X = 1, TYPE_Y = 2

! Define the location of the state variables in module storage
type(location_type), allocatable :: state_loc(:)
type(time_type) :: time_step

type static_data
  integer :: x_size, y_size, model_size
  real(r8) :: b,c,h,f
  real(r8), dimension(:), pointer :: x_loc,y_loc
end type static_data

type(static_data) :: l96

contains

!==================================================================



subroutine static_init_model()
!------------------------------------------------------------------
! Initializes class data for this model. For now, simply outputs the
! identity info, sets the location of the state variables, and initializes
! the time type for the time stepping (is this general enough for time???)

integer  :: i, iunit, io, kount

! Print module information to log file and stdout.
call register_module(source, revision, revdate)

! Read the namelist entry
call find_namelist_in_file("input.nml", "model_nml", iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, "model_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)

! constants
l96%h = coupling_h
l96%b = coupling_b
l96%c = coupling_c
l96%f = forcing

! size of domain and state
l96%x_size = model_size_x
l96%y_size = model_size_x*y_per_x
l96%model_size = l96%x_size + l96%y_size

! Create storage for locations
allocate(l96%x_loc(l96%x_size))
allocate(l96%y_loc(l96%y_size))
allocate(state_loc(l96%model_size))

! Define the locations of the model state variables
! Carrying all three for now - may not be necessary
kount = 1
do i = 1, l96%x_size
   l96%x_loc(i) = (i - 1.0_r8) / l96%x_size
   state_loc(kount) =  set_location(l96%x_loc(i))
   kount = kount+1
end do
do i = 1, l96%y_size
   l96%y_loc(i) = (i - 1.0_r8) / l96%y_size
   state_loc(kount) =  set_location(l96%y_loc(i))
   kount = kount+1
end do

! The time_step in terms of a time type must also be initialized. Need
! to determine appropriate non-dimensionalization conversion for L96 from
! Shree Khare.
time_step = set_time(time_step_seconds, time_step_days)

end subroutine static_init_model



subroutine init_model_instance()
!------------------------------------------------------------------
! subroutine init_model_instance
!
! Initializes instance dependent state for model. Null for L96.

end subroutine init_model_instance



subroutine comp_dt(x, dt)
!------------------------------------------------------------------
! subroutine comp_dt(x, dt) (note used for both x and y together)
! 
! Computes the time tendency of the lorenz 1996 model given current state

real(r8), intent( in) ::  x(:)
real(r8), intent(out) :: dt(:)

integer :: j, jp1, jp2, jm1, jm2 
integer :: k, kp1,      km1, km2 
integer :: jk
integer :: xs, xe, ys, ye, js, je
real(r8) :: fast_sum, c1, c2, c3
real(r8), dimension(y_per_x) :: tmp_Y, tmp_dt

c1 = l96%c * l96%b
c2 = l96%c
c3 = l96%h * l96%c / l96%b

ys = l96%x_size + 1
ye = l96%model_size
! first small-scale, then large

if ( local_y ) then

! Smith's version treats each group of Y variables independently (could be
! a mistake in the paper)

   do jk = ys, ye, y_per_x
      tmp_Y = x(jk:jk+y_per_x-1)
      do j = 1, y_per_x
         jp1 = j + 1
         if(jp1 > y_per_x) jp1 = 1
         jp2 = j + 2
         if(jp2 > y_per_x) jp2 = jp2 - y_per_x 
         jm2 = j - 2
         if(jm2 < 1) jm2 = y_per_x + jm2
         jm1 = j - 1
         if(jm1 < 1) jm1 = y_per_x
         tmp_dt(j) = c1 * tmp_Y(jp1) * (tmp_Y(jm1)-tmp_Y(jp2)) - c2 * tmp_Y(j) &
             + c3 * x(idint(dble(jk+j-1-ys)/dble(y_per_x))+1)
      enddo
      dt(jk:jk+y_per_x-1) = tmp_dt
   enddo
else
! Lorenz' version treats all the Y variables together so that they
! are aware of each other across the Xs

   do j = ys, ye
      jp1 = j + 1
      if(jp1 > ye) jp1 = ys
      jp2 = j + 2
      if(jp2 > ye) jp2 = ys + jp2 - ye - 1
      jm2 = j - 2
      if(jm2 < ys) jm2 = ye + jm2 - ys + 1
      jm1 = j - 1
      if(jm1 < ys) jm1 = ye
      dt(j) = c1 * x(jp1) * (x(jm1)-x(jp2)) - c2 * x(j) &
           + c3 * x(idint(dble(j-ys)/dble(y_per_x))+1)
   enddo
endif

! Now large scale
xs = 1
xe = l96%x_size
do k = xs, xe
   js = (k-1)*y_per_x + l96%x_size + 1
   je = k*y_per_x + l96%x_size 
   fast_sum = sum(x(js:je))
   kp1 = k + 1
   if(kp1 > xe) kp1 = 1
   km2 = k - 2
   if(km2 < 1) km2 = xe + km2
   km1 = k - 1
   if(km1 < 1) km1 = xe
   dt(k) = (x(kp1) - x(km2)) * x(km1) - x(k) + l96%f &
           - c3 * fast_sum
end do
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

! Need noise in both the large and small scales.  The small-scale noise
! should be present in each group of Ys in case the local_y option is selected

x(1:l96%x_size) = l96%f
x(l96%x_size+1:l96%model_size) = 0.0_r8
x(1) = 1.001_r8 * l96%f

x(l96%x_size+1:) = 0.01_r8 * x(2)
x(l96%x_size+1:l96%model_size:y_per_x) = 0.011_r8 * x(2) 

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

get_model_size = l96%model_size
!print*, 'model size is ',l96%model_size

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

! Argument itype IS used here 


real(r8),            intent(in) :: x(:)
real(r8)                        :: obs_val
type(location_type), intent(in) :: location
integer,             intent(in) :: itype
integer,            intent(out) :: istatus


integer :: lower_index, upper_index,    base_index, top_index
real(r8) :: lctn, lctnfrac

! All interpolations okay for now
istatus = 0

! Convert location to real
lctn = get_location(location)
! Multiply by model size assuming domain is [0, 1] cyclic
if ( itype == 1 ) then
   lctn = l96%x_size * lctn
   base_index = 1
   top_index = l96%x_size
elseif ( itype == 2 ) then
   lctn = l96%y_size * lctn
   base_index = l96%x_size + 1
   top_index = l96%model_size
else
   call error_handler(E_ERR,'model_interpolate', 'cannot handle this type', &
                     source, revision, revdate)
endif

lower_index = int(lctn) + base_index
upper_index = lower_index + 1
if(lower_index > top_index) lower_index = lower_index - top_index
if(upper_index > top_index) upper_index = upper_index - top_index

lctnfrac = lctn - int(lctn)
obs_val = (1.0_r8 - lctnfrac) * x(lower_index) + lctnfrac * x(upper_index)

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
  var_type = 1    
  if ( index_in > model_size_x ) var_type = 2
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
! JPH 6 May 2004 -- additional attributes for the Y variables if requested
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

!-----------------------------------------------------------------------------------------
! General netCDF variables
!-----------------------------------------------------------------------------------------

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID

!-----------------------------------------------------------------------------------------
! netCDF variables for Location
!-----------------------------------------------------------------------------------------

integer :: LocationVarID, XLocationVarID, YLocationVarID
integer :: StateVarDimID, StateVarVarID
integer :: StateVarID, MemberDimID, TimeDimID
integer :: XDimID, XVarID
integer :: YDimID, YVarID
integer :: XID, YID

!--------------------------------------------------------------------
!Bounds of X and Y in vector
!--------------------------------------------------------------------

integer              :: xs, xe, ys, ye 

!-----------------------------------------------------------------------------------------
! local variables
!-----------------------------------------------------------------------------------------

character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic
character(len=NF90_MAX_NAME) :: str1

integer             :: i
type(location_type) :: lctn 
ierr = 0                      ! assume normal termination

!-------------------------------------------------------------------------------
! Get bounds
!-------------------------------------------------------------------------------

xs = 1
xe = l96%x_size
ys = xe + 1
ye = xe + l96%y_size

!-------------------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file 
!-------------------------------------------------------------------------------

call check(nf90_Inquire(ncFileID, nDimensions, nVariables, nAttributes, unlimitedDimID))
call check(nf90_sync(ncFileID)) ! Ensure netCDF file is current
call check(nf90_Redef(ncFileID))

!-------------------------------------------------------------------------------
! Determine ID's from stuff already in the netCDF file
!-------------------------------------------------------------------------------

! make sure time is unlimited dimid

call check(nf90_inq_dimid(ncFileID,"copy",dimid=MemberDimID))
call check(nf90_inq_dimid(ncFileID,"time",dimid=TimeDimID))

!-------------------------------------------------------------------------------
! Write Global Attributes 
!-------------------------------------------------------------------------------

call DATE_AND_TIME(crdate,crtime,crzone,values)
write(str1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
                  values(1), values(2), values(3), values(5), values(6), values(7)

call check(nf90_put_att(ncFileID, NF90_GLOBAL, "creation_date",str1))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_source", source ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revision", revision ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revdate", revdate ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model", "Lorenz_96_2scale"))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_delta_t",    delta_t ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_coupling_b", l96%b ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_coupling_c", l96%c ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_coupling_h", l96%h ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_forcing",    l96%f ))

!-------------------------------------------------------------------------------
! Define the model size, state variable dimension ... whatever ...
!-------------------------------------------------------------------------------

call check(nf90_def_dim(ncid=ncFileID, name="StateVariable", &
                        len=l96%model_size, dimid = StateVarDimID)) 

!-------------------------------------------------------------------------------
! Define the dimensions IDs for X and Y dimensions
!-------------------------------------------------------------------------------
call check(nf90_def_dim(ncid=ncFileID, name="Xdim", &
                        len=l96%x_size, dimid = XDimID)) 
call check(nf90_def_dim(ncid=ncFileID, name="Ydim", &
                        len=l96%y_size, dimid = YDimID)) 

!-------------------------------------------------------------------------------
! Define the Location Variable and add Attributes
! Some of the atts come from location_mod (via the USE: stmnt)
! CF standards for Locations:
! http://www.cgd.ucar.edu/cms/eaton/netcdf/CF-working.html#ctype
!-------------------------------------------------------------------------------

!--------------------------------------------------------------------
! The location info depends on output_state_vector (don't need both)
!--------------------------------------------------------------------

if ( output_state_vector) then

   call check(NF90_def_var(ncFileID, name=trim(adjustl(LocationName)), xtype=nf90_double, &
              dimids = StateVarDimID, varid=LocationVarID) )
   call check(nf90_put_att(ncFileID, LocationVarID, "long_name", trim(adjustl(LocationLName))))
   call check(nf90_put_att(ncFileID, LocationVarID, "dimension", LocationDims ))
   call check(nf90_put_att(ncFileID, LocationVarID, "valid_range", (/ 0.0_r8, 1.0_r8 /)))

else

   call check(NF90_def_var(ncFileID, name="XLocation", xtype=nf90_double, &
              dimids = XDimID, varid=XLocationVarID) )
   call check(nf90_put_att(ncFileID, XLocationVarID, "long_name", trim(adjustl(LocationLName))))
   call check(nf90_put_att(ncFileID, XLocationVarID, "dimension", LocationDims ))
   call check(nf90_put_att(ncFileID, XLocationVarID, "valid_range", (/ 0.0_r8, 1.0_r8 /)))

   call check(NF90_def_var(ncFileID, name="YLocation", xtype=nf90_double, &
              dimids = YDimID, varid=YLocationVarID) )
   call check(nf90_put_att(ncFileID, YLocationVarID, "long_name", trim(adjustl(LocationLName))))
   call check(nf90_put_att(ncFileID, YLocationVarID, "dimension", LocationDims ))
   call check(nf90_put_att(ncFileID, YLocationVarID, "valid_range", (/ 0.0_r8, 1.0_r8 /)))

endif

!-------------------------------------------------------------------------------
! Define either the "state vector" variables -OR- the "prognostic" variables.
!-------------------------------------------------------------------------------

if ( output_state_vector ) then

! Define the state vector coordinate variable
   call check(nf90_def_var(ncid=ncFileID,name="StateVariable", xtype=nf90_int, &
           dimids=StateVarDimID, varid=StateVarVarID))
   call check(nf90_put_att(ncFileID, StateVarVarID, "long_name", "State Variable ID"))
   call check(nf90_put_att(ncFileID, StateVarVarID, "units",     "indexical") )
   call check(nf90_put_att(ncFileID, StateVarVarID, "valid_range", (/ 1, l96%model_size /)))

! Define the actual state vector
   call check(nf90_def_var(ncid=ncFileID, name="state", xtype=nf90_double, &
           dimids = (/ StateVarDimID, MemberDimID, TimeDimID /), varid=StateVarID))
   call check(nf90_put_att(ncFileID, StateVarID, "long_name", "model state or fcopy"))

! Leave define mode so we can fill
   call check(nf90_enddef(ncfileID))

! Fill the state variable coordinate variable
   call check(nf90_put_var(ncFileID, StateVarVarID, (/ (i,i=1,l96%model_size) /) ))

!-------------------------------------------------------------------------------
! Fill the location variable
!-------------------------------------------------------------------------------

   do i = 1,l96%model_size
      call get_state_meta_data(i,lctn)
      call check(nf90_put_var(ncFileID, LocationVarID, get_location(lctn), (/ i /) ))
   enddo

else ! output prognostic variables

! Define the coordinate variables
   call check(nf90_def_var(ncid=ncFileID,name="Xdim", xtype=nf90_int, &
           dimids=XDimID, varid=XVarID))
   call check(nf90_put_att(ncFileID, XVarID, "long_name", "X ID"))
   call check(nf90_put_att(ncFileID, XVarID, "units",     "indexical") )
   call check(nf90_put_att(ncFileID, XVarID, "valid_range", (/ 1, l96%x_size /)))

   call check(nf90_def_var(ncid=ncFileID,name="Ydim", xtype=nf90_int, &
           dimids=YDimID, varid=YVarID))
   call check(nf90_put_att(ncFileID, YVarID, "long_name", "Y ID"))
   call check(nf90_put_att(ncFileID, YVarID, "units",     "indexical") )
   call check(nf90_put_att(ncFileID, YVarID, "valid_range", (/ 1, l96%y_size /)))

! Define the actual state variables
   call check(nf90_def_var(ncid=ncFileID, name="X", xtype=nf90_double, &
           dimids = (/ XDimID, MemberDimID, TimeDimID /), varid=XID))
   call check(nf90_put_att(ncFileID, XID, "long_name", "slow variables X"))

   call check(nf90_def_var(ncid=ncFileID, name="Y", xtype=nf90_double, &
           dimids = (/ YDimID, MemberDimID, TimeDimID /), varid=YID))
   call check(nf90_put_att(ncFileID, YID, "long_name", "fast variables Y"))

! Leave define mode so we can fill
   call check(nf90_enddef(ncfileID))

! Fill the state variable coordinate variables
   call check(nf90_put_var(ncFileID, XVarID, (/ (i,i=1,l96%x_size) /) ))

   call check(nf90_put_var(ncFileID, YVarID, (/ (i,i=1,l96%y_size) /) ))

!-------------------------------------------------------------------------------
! Fill the location variables
!-------------------------------------------------------------------------------

   do i = xs, xe
      call get_state_meta_data(i,lctn)
      call check(nf90_put_var(ncFileID, XLocationVarID, get_location(lctn), (/ i /) ))
   enddo

   do i = ys, ye
      call get_state_meta_data(i,lctn)
      call check(nf90_put_var(ncFileID, YLocationVarID, get_location(lctn), (/ i - ys + 1 /) ))
   enddo

endif

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------
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

!-------------------------------------------------------------------------------
! General netCDF variables
!-------------------------------------------------------------------------------

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID
integer :: StateVarID, XID, YID

!--------------------------------------------------------------------
!Bounds of X and Y in vector
!--------------------------------------------------------------------

integer              :: xs, xe, ys, ye 

!-------------------------------------------------------------------------------
! local variables
!-------------------------------------------------------------------------------

ierr = 0                      ! assume normal termination

!-------------------------------------------------------------------------------
! Get bounds
!-------------------------------------------------------------------------------

xs = 1
xe = l96%x_size
ys = xe + 1
ye = xe + l96%y_size

!-------------------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file
!-------------------------------------------------------------------------------

call check(nf90_Inquire(ncFileID, nDimensions, nVariables, nAttributes, unlimitedDimID))

!------------------------------------------------------------------------
! Branch for state vector OR prognostic variables
!------------------------------------------------------------------------

if ( output_state_vector ) then

   call check(NF90_inq_varid(ncFileID, "state", StateVarID) )
   call check(NF90_put_var(ncFileID, StateVarID, statevec,  &
             start=(/ 1, copyindex, timeindex /)))

else

   call check(NF90_inq_varid(ncFileID, "X", XID) )
   call check(NF90_put_var(ncFileID, XID, statevec(xs:xe),  &
             start=(/ 1, copyindex, timeindex /)))

   call check(NF90_inq_varid(ncFileID, "Y", YID) )
   call check(NF90_put_var(ncFileID, YID, statevec(ys:ye),  &
             start=(/ 1, copyindex, timeindex /)))

endif

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
