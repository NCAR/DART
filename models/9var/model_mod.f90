! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module model_mod

use        types_mod,      only : r8, i8, i4

use     location_mod,      only : location_type, set_location, get_location, &
                                  LocationDims, LocationName, LocationLName, &
                                  get_close_maxdist_init, get_close_obs_init, &
                                  loc_get_close_obs => get_close_obs, get_close_type

use    utilities_mod,      only : register_module, error_handler, E_ERR, E_MSG, do_output, &
                                  nmlfileunit, find_namelist_in_file, check_namelist_read, &
                                  do_nml_file, do_nml_term

use ensemble_manager_mod,  only : ensemble_type

use distributed_state_mod, only : get_state

use state_structure_mod,   only : add_domain

use dart_time_io_mod,      only : read_model_time, write_model_time


! All random_seq_mod calls were suppressed because a) they are not being used,
! and b) they make the pg5.02 compiler complain about a gap in the common block.
! TJH 29 April 2004
! TJH use   random_seq_mod, only : random_seq_type, random_gaussian, &
! TJH                             init_random_seq, several_random_gaussians
use time_manager_mod, only : time_type, set_time

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

integer(i8), parameter :: model_size = 9

!  define model parameters
! c is sqrt(0.75)

real(r8), private, parameter :: a(3) = (/  1.0_r8,  1.0_r8, 3.0_r8 /), &
                                b(3) = (/ -1.5_r8, -1.5_r8, 0.5_r8 /), &
                                f(3) = (/ 0.10_r8,  0.0_r8, 0.0_r8 /), &
                                h(3) = (/ -1.0_r8,  0.0_r8, 0.0_r8 /), &
                                nu = 1.0_r8 / 48.0_r8, kappa = nu, c = 0.8660254_r8

!-------------------------------------------------------------
! Namelist with default values
!
real(r8) :: g = 8.0_r8     ! lorenz default
!real(r8) :: g = 9.90_r8    ! higher dimension attractor

real(r8) :: deltat = 1.0_r8 / 12.0_r8     ! model time step
integer :: time_step_days = 0
integer :: time_step_seconds = 3600

namelist /model_nml/ g, deltat, time_step_days, time_step_seconds
!---------------------------------------------------------------


! Define the location of the state variables in module storage
! This is used for distance dependence stuff in more general models but is 
! currently just set to give 0 separation distance for all 9 vars here.

type(location_type) :: state_loc(model_size)
type(time_type)     :: time_step

! Need reproducible sequence of noise added so that different runs
! can be cleanly compared

!! used in code that is currently commented out.
!!logical :: first_ran_call = .true.
!!type(random_seq_type) :: ran_seq



contains

!======================================================================


subroutine static_init_model()
!----------------------------------------------------------------------
! subroutine static_init_model()
!
! Initializes class data for this model. For now, simply outputs the
! identity info, sets the location of the state variables, and initializes
! the time type for the time stepping (is this general enough for time???)

implicit none
real(r8) :: x_loc
integer :: i, iunit, io, dom_id

! Register the module into the logfile
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

! The time_step in terms of a time type must also be initialized. Need
! to determine appropriate non-dimensionalization conversion for L96 from
! Shree Khare.
time_step = set_time(time_step_seconds, time_step_days)

! Tell the DART I/O routines how large the model data is so they
! can read/write it.
dom_id = add_domain(model_size)

end subroutine static_init_model




subroutine comp_dt(xxx, dxxx)
!----------------------------------------------------------------------
! subroutine comp_dt(xxx, dxxx)
!
! computes time tendency for the 9 variable lorenz model given the
! current values of the 9 variables (x, y, z each with three elements)
! and the values of a number of parameters.

implicit none

real(r8), intent(in)  :: xxx(:)
real(r8), intent(out) :: dxxx(:)

real(r8) :: x(3), y(3), z(3), dx(3), dy(3), dz(3)
!!! real(r8) :: rnum(9)
integer  :: i, j, k

!  unpack the 9-vectors into the x, y and z 3-vectors

call unpack9var(xxx, x, y, z)

!  equations 33-35 from lorenz, 1980, jas, p.1688
!  equations are defined with cyclic indices

do i = 1, 3
   J = MOD(I, 3) + 1
   K = MOD(I + 1, 3) + 1

   !  BEGIN WITH THE FIRST EQUATION (33)

   dx(i) = (a(i)*b(i)*x(j)*x(k) - c*(a(i) - a(k))*x(j)*y(k) + &
           c*(a(i) - a(j))*y(j)*x(k) - 2*c**2*y(j)*y(k) - &
           nu*a(i)**2*x(i) + a(i)*y(i) - a(i)*z(i)) / a(i)

   !  equation (34)

   dy(i) = (-1.*a(k)*b(k)*x(j)*y(k) - a(j)*b(j)*y(j)*x(k) + &
           c*(a(k) - a(j))*y(j)*y(k) - a(i)*x(i) - nu*a(i)**2*y(i)) / a(i)

   !  equation (35)

   dz(i) = -1.*b(k)*x(j)*(z(k) - h(k)) - &
           b(j)*(z(j) - h(j))*x(k) + c*y(j)*(z(k) -h(k)) - &
           c*(z(j) - h(j))*y(k) + g*a(i)*x(i) - kappa*a(i)*z(i) + f(i)
end do   

call pack9var(dx, dy, dz, dxxx)     !  pack the results into 9 vector

! OPTIONAL ADDItion OF NOISE
! ADDITION OF SOME NOISE AT 1/10 the amplitude of DT
! We need to initialize the repeatable random gen
!!!if(first_ran_call) then
!!!   first_ran_call = .false.
!!!   call init_random_seq(ran_seq)
!!!end if

!!!call several_random_gaussians(ran_seq, dble(0.0), dxxx(i) / 10.0, 9, rnum)
!!!do i = 1, 9
!!!   dxxx(i) = dxxx(i) + rnum(i)
!!!end do

end subroutine comp_dt




subroutine pack9var(x, y, z, pert)
!---------------------------------------------------------------------------
! subroutine pack9var(x, y, z, pert)
!
! set of routines used to switch between 3 3-variable sets and 9-variable
! set for lorenz 9 variable pe model

implicit none

!  pack9var and unpack9var convert from x, y, z to full 9 vector format

real(r8), intent(in)  :: x(3), y(3), z(3)
real(r8), intent(out) :: pert(9)

pert(1:3) = x
pert(4:6) = y
pert(7:9) = z

end subroutine pack9var



subroutine unpack9var(pert, x, y, z)
!---------------------------------------------------------------------------
! subroutine unpack9var(pert, x, y, z)
!
! inverse of pack9var above

implicit none

real(r8), intent(in)  :: pert(9)
real(r8), intent(out) :: x(3), y(3), z(3)

x = pert(1:3)
y = pert(4:6)
z = pert(7:9)

end subroutine unpack9var



subroutine advance(x, num, xnew, time)
!-----------------------------------------------------------------------
! subroutine advance(x, num, xnew, time)
!
! advance advances the 9 variable model by a given number of steps

implicit none

real(r8),        intent(in)  :: x(9)
integer,         intent(in)  :: num
real(r8),        intent(out) :: xnew(9)
type(time_type), intent(in)  :: time

integer :: i

xnew = x                 !  copy initial conditions to avoid overwrite

do i = 1, num            !  advance the appropriate number of steps
   call adv_1step(xnew, time)
end do

end subroutine advance





subroutine adv_1step(x, time)
!-------------------------------------------------------------------------
! subroutine adv_1step(x, time)
!
! does one time step advance for 9 variable model using two-step rk.
! The Time argument is needed for compatibility with more complex models
! that need to know the time to compute their time tendency and is not
! used in L96. Is there a better way to do this in F90 than to just hang
! this argument out everywhere?

implicit none

real(r8),        intent(inout) :: x(:)
type(time_type), intent(in)    :: time

real(r8) :: fract = 1.0_r8

call adv_single(x, fract)

end subroutine adv_1step




subroutine adv_single(x, fract)
!-------------------------------------------------------------------------
! subroutine adv_single(x, fract)
!
! does one time step advance for 9 variable model using two-step rk

implicit none

real(r8), intent(inout) :: x(:)
real(r8), intent(in)    :: fract

real(r8) :: x1(9), x2(9), dx(9)

!  compute the first intermediate step

call comp_dt(x, dx)
x1 = x + fract * deltat * dx

!  compute the second intermediate step

call comp_dt(x1, dx)
x2 = x1 + fract * deltat * dx

!  new value for x is average of original value and second intermediate

x = (x + x2) / 2.0_r8

end subroutine adv_single



subroutine init_conditions(x)
!----------------------------------------------------------------------
! subroutine init_conditions(x)
!
! generates set of random off-attractor initial conditions for 9 variable

implicit none

real(r8), intent(out) ::  x(:)

x = 0.10_r8

end subroutine init_conditions



subroutine linearize(nl, l)
!----------------------------------------------------------------------
! subroutine linearize(nl, l)
!

implicit none

real(r8) :: nl(:), l(:, :)

!  no-op subroutine header for linking in standard packages

end subroutine linearize



subroutine balance_init(xxx, init_xxx)
!---------------------------------------------------------------------------
! subroutine balance_init(xxx, init_xxx)
!
! application of balance eqtn as initialization method
! to the 9 variable lorenz model given the
! current values of the 9 variables (x, y, z each with three elements)
! and the values of a number of parameters.

implicit none

real(r8), intent(in)  :: xxx(9)
real(r8), intent(out) :: init_xxx(9)

real(r8) :: x(3), y(3), z(3)
!!!real(r8) :: wkspce(3)
real(r8) :: lhs(3,3),rhs(3)
integer  :: i, j, k
!!!integer  :: ifail

!  unpack the 9-vectors into the x, y and z 3-vectors

call unpack9var(xxx, x, y, z)

! application of balance eqtn as initialization method, curry ET.AL, tellus,
!    1995, p. 153-154, section 3.3
! equation 11 from gent & mcwilliams, 1982, jas, p.4

DO i = 1, 3
   j    = mod(i, 3) + 1
   k    = mod(i + 1, 3) + 1
   z(i) = (a(i)*y(i)-2*c**2*y(j)*y(k))/a(i)
end do

!  equation 29/30 from gent & mcwilliams, 1982, jas, p.6

do i = 1, 3
   j = mod(i, 3) + 1
   k = mod(i + 1, 3) + 1
   lhs(i,i) = ( a(i)*a(j)*a(k)*(1+g*a(i))-2*c**2*(a(j)**2* &
      b(j)*y(j)**2+a(k)**2*b(k)*y(k)**2) )
   lhs(j,i) = -( a(j)*a(k)*( y(k)*(2*c**2-a(k)*b(k))+a(i)* &
      b(k)*(z(k)-h(k)) )+2*c**2*a(i)*a(j)*b(i)*y(i) *y(j)  )
   lhs(k,i) = -( a(j)*a(k)*( y(j)*(2*c**2-a(j)*b(j))+a(i)* &
      b(j)*(z(j)-h(j)) )+2*c**2*a(i)*a(k)*b(i)*y(i) *y(k)  )
   rhs(i) = a(j)*a(k)*(  c*(a(k)-a(j))*y(j)*y(k)+c*a(i)*( &
      (z(j)-h(j))*y(k)-y(j)*(z(k)-h(k)) )+a(i)*( nu* &
      a(i)*(z(i)-y(i))-f(i) )  )-2*c**2*(  c*a(j)* &
      (a(j)-a(i))*y(i)*y(j)**2+c*a(k)*(a(i)-a(k))*y(i) &
      *y(k)**2-nu*a(j)*a(k)*(a(j)+a(k))*y(j)*y(k)  )
end do

!!!ifail=0
!  n.a.g. ROUTINE TO CALCULATE APPROXIMATE SOLUTION X TO aX=B (I.E. lhs*X=rhs)
!!!call f04arf_wrap(lhs,3,rhs,3,x,wkspce,ifail)
!        print *,'XINIT=',X

!  pack the results into 9 vector

call pack9var(x, y, z, init_xxx)

end subroutine balance_init



subroutine get_close_pts(list, num)
!-------------------------------------------------------------------------
! subroutine get_close_pts(list, num)
!

implicit none

integer, intent(in)    :: num
integer, intent(inout) :: list(model_size, num)

integer :: i, offset, indx, temp

do i = 1, model_size

   do offset = -num/2, -num/2 + num - 1
      indx = i + offset
      if(indx > model_size) indx = indx - model_size
      if(indx < 1         ) indx = model_size + indx
      list(i, offset + num/2 + 1) = indx
   end do

   ! Always need the actual point first in list

   temp = list(i, 1)
   list(i, 1) =  list(i, num / 2 + 1)
   list(i, num / 2 + 1) = temp

end do

end subroutine get_close_pts



function get_model_size()
!-------------------------------------------------------------------------
! function get_model_size()
!
! Returns size of model

integer :: get_model_size

get_model_size = model_size

end function get_model_size




function get_model_time_step()
!------------------------------------------------------------------------
! function get_model_time_step()
!
! Returns the the time step of the model. In the long run should be repalced
! by a more general routine that returns details of a general time-stepping
! capability.

type(time_type) :: get_model_time_step

get_model_time_step = time_step

end function get_model_time_step



subroutine model_interpolate(state_handle, ens_size, location, itype, expected_obs, istatus)
!---------------------------------------------------------------------
!
! Interpolates from state vector x to the location. It's not particularly
! happy dumping all of this straight into the model. Eventually some
! concept of a grid underlying models but above locations is going to
! be more general. May want to wait on external infrastructure projects
! for this?

! Argument itype is not used here because there is only one type of variable.
! itype is needed to allow swap consistency with more complex models.

implicit none

type(ensemble_type),  intent(in) :: state_handle
integer,              intent(in) :: ens_size
type(location_type),  intent(in) :: location
integer,              intent(in) :: itype
real(r8),            intent(out) :: expected_obs(ens_size)
integer,             intent(out) :: istatus(ens_size)

integer(i8)  :: lower_index, upper_index
real(r8) :: lctn, lctnfrac

! All interps okay for now
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
expected_obs = (1.0_r8 - lctnfrac) * get_state(lower_index, state_handle) + &
                         lctnfrac  * get_state(upper_index, state_handle)

end subroutine model_interpolate



subroutine get_state_meta_data(state_handle, index_in, location, var_type)
!---------------------------------------------------------------------
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
if (present(var_type)) var_type = 1    ! default variable type

end subroutine get_state_meta_data



  subroutine init_model()
!-------------------------------------------------------------------------
! subroutine init_model()
!
! Stub for model initialization, not needed for 9var

end subroutine init_model




subroutine init_time(time)
!----------------------------------------------------------------------
!
! Gets the initial time for a state from the model. Where should this info
! come from in the most general case?

implicit none

type(time_type), intent(out) :: time

! For now, just set to 0
time = set_time(0, 0)

end subroutine init_time



subroutine end_model()
!------------------------------------------------------------------------
!
! Does any shutdown and clean-up needed for model. Nothing for L96 for now.


end subroutine end_model






function nc_write_model_atts( ncFileID, model_mod_writes_state_variables ) result (ierr)
!--------------------------------------------------------------------
! Writes the model-specific attributes to a netCDF file
! TJH Jan 24 2003; added by JLA 18 June, 2003
!
! TJH 29 July 2003 -- for the moment, all errors are fatal, so the
! return code is always '0 == normal', since the fatal errors stop execution.
!
! For the 9var model, each state variable is at a separate location.
! that's all the model-specific attributes I can think of ...
! In addition, there are technically three kinds of variables; this
! level of detail will have to be added by TJH at a later date.
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
!

use typesizes           ! comes from F90 netCDF interface
use netcdf              ! comes from F90 netCDF interface
implicit none

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
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model", "9var" ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_g", g ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_deltat", deltat ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_a", a ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_b", b ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_f", f ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_h", h ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_nu", nu ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_kappa", kappa ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_c", c ))

!--------------------------------------------------------------------
! Define the model size, state variable dimension ... whatever ...
!--------------------------------------------------------------------
call check(nf90_def_dim(ncid=ncFileID, name="location", &
                        len=int(model_size,i4), dimid = LocationDimID))

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
   call check(nf90_put_var(ncFileID, LocationVarID, get_location(state_loc(i)), (/ int(i, i4) /) ))
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

    if(istatus /= nf90_noerr) call error_handler(E_ERR, 'nc_write_model_atts', &
      trim(nf90_strerror(istatus)), source, revision, revdate)

  end subroutine check

end function nc_write_model_atts



function nc_write_model_vars( ncFileID, statevec, copyindex, timeindex ) result (ierr)
!--------------------------------------------------------------------
! Writes the model-specific attributes to a netCDF file
! TJH 25 June 2003
!
! TJH 29 July 2003 -- for the moment, all errors are fatal, so the
! return code is always '0 == normal', since the fatal errors stop execution.
!
! For the 9var model, there are actually 3 triplets of variables, 
! which generate a state vector of 9 variables. Completing the "prognostic"
! netCDF write will occurr at a later date.
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
!

use typeSizes      ! comes from F90 netCDF interface
use netcdf         ! comes from F90 netCDF interface
implicit none

integer,                intent(in) :: ncFileID      ! netCDF file identifier
real(r8), dimension(:), intent(in) :: statevec
integer,                intent(in) :: copyindex
integer,                intent(in) :: timeindex
integer                            :: ierr          ! return value of function

ierr = 0                      ! assume normal termination

end function nc_write_model_vars


subroutine pert_model_copies(state_ens_handle, ens_size, pert_amp, interf_provided)
!------------------------------------------------------------------
! Perturbs a model state copies for generating initial ensembles.
! Routine which could provide a custom perturbation routine to
! generate initial ensembles.  The default (if interface is not
! provided) is to add gaussian noise to each item in the state vector.

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
! End of 9var model_mod 
!===================================================================
end module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
