! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module model_mod

! Assimilation interface for Held-Suarez Bgrid

!-----------------------------------------------------------------------
!
!     interface for B-grid dynamics using the Held-Suarez forcing
!
!-----------------------------------------------------------------------
!---------------- m o d u l e   i n f o r m a t i o n ------------------

use bgrid_core_driver_mod, only: bgrid_dynam_type,       &
                                 bgrid_core_driver_init, &
                                 bgrid_core_driver,      &
                                 bgrid_core_time_diff,   &
                                 bgrid_core_driver_end,  &
                                 get_bottom_data,        &
                                 put_bottom_data

!use diag_manager_mod, only: diag_manager_init, get_base_date


use   bgrid_prog_var_mod, only: prog_var_type, var_init, prog_var_init, &
                                open_prog_var_file, read_prog_var

use      bgrid_horiz_mod, only: get_horiz_grid_size,        &
                                get_horiz_grid_bound, TGRID, VGRID

use time_manager_mod, only: time_type, set_time, get_time,  &
                            operator(+), operator (<), operator (>), &
                            operator (/=), operator (/), operator (*)

use              fms_mod, only: file_exist, open_namelist_file, &
                                error_mesg, FATAL,              &
                                stdlog,        &
                                write_version_number,           &
                                close_file,         &
                                fms_init,        &
                                open_restart_file


! routines used by subroutine bgrid_physics
use bgrid_change_grid_mod, only: vel_to_mass, mass_to_vel
use       bgrid_horiz_mod, only: horiz_grid_type
use        bgrid_vert_mod, only: vert_grid_type, &
                                 compute_pres_full, compute_pres_half
use        bgrid_halo_mod, only: update_halo, UWND, VWND, TEMP, &
                                 NORTH, EAST, WEST, SOUTH
use        hs_forcing_mod, only: hs_forcing_init, hs_forcing

! DART routines 
use          location_mod, only: location_type, get_location, set_location, &
                                 vert_is_level, query_location, &
                                 LocationDims, LocationName, LocationLName, &
                                 vert_is_pressure, vert_is_surface, &
                                 get_close_maxdist_init, get_close_obs_init, get_close_obs


use        random_seq_mod, only: random_seq_type, init_random_seq, random_gaussian
use             types_mod, only: r8, pi, MISSING_R8
! combined duplicate use lines for utilities_mod; the intel 8.x compiler
! was unhappy about the repetition.  nsc 11apr06
use        utilities_mod, only : open_file, error_handler, E_ERR, E_MSG, &
                                 nmlfileunit, register_module, &
                                 find_namelist_in_file, check_namelist_read, &
                                 do_nml_file, do_nml_term, do_output

use          obs_kind_mod, only: KIND_U_WIND_COMPONENT, KIND_V_WIND_COMPONENT, &
                                 KIND_SURFACE_PRESSURE, KIND_TEMPERATURE

!-----------------------------------------------------------------------

implicit none
private

public  get_model_size, adv_1step, get_state_meta_data, model_interpolate, &
        get_model_time_step, end_model, static_init_model, init_time, &
        init_conditions, TYPE_PS, TYPE_T, TYPE_U, TYPE_V, TYPE_TRACER, &
        nc_write_model_atts, nc_write_model_vars, &
        pert_model_state, &
        get_close_maxdist_init, get_close_obs_init, get_close_obs, ens_mean_for_model


!-----------------------------------------------------------------------
! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"
character(len=32 ), parameter :: version = "$Revision$"
character(len=128), parameter :: tag = "$Id$"

logical, save :: module_initialized = .false.

!-----------------------------------------------------------------------
! bgrid_prog_var_mod:prog_var_type
! data structure that contains all prognostic fields and tracers
! type prog_var_type
!      integer       :: nlon, nlat, nlev, ntrace
!      integer       :: ilb, iub, jlb, jub, klb, kub
!      real, pointer :: ps(:,:), pssl(:,:)
!      real, pointer :: u(:,:,:), v(:,:,:), t(:,:,:), r(:,:,:,:)
! end type prog_var_type
! It would be great if the data structure could have character names and some metadata
! for each variable ... units, etc. 

! integer, parameter :: N_prog_Var = 6
! character(len=NF90_MAX_NAME), dimension(N_prog_Var) :: Prog_Var_Names = (/ "ps", "pssl", "u","v","t","r" /)
! character(len=NF90_MAX_NAME), dimension(N_prog_Var) :: Prog_Var_Units = (/ "hPa", "hPa", "m/s","m/s","degrees Kelvin","depends" /)


!-----------------------------------------------------------------------
!---- namelist (saved in file input.nml) ----
!
! physics_window  The number of "i" by "j" rows processed each time
!                 the modular physics is called. To process the entire
!                 domain use physics_window = (/0,0/).
!                   [integer, default: physics_window = 0,1]

   integer, dimension(2) :: physics_window = (/0,0/)

   namelist /atmosphere_nml/ physics_window

!-----------------------------------------------------------------------
! Additional stuff currently in main_nml moved down from atmos_solo_driver
! This time stuff may not really belong here???

   integer, dimension(4) :: current_time = (/ 0, 0, 0, 0 /)
   logical  :: override = .false.  ! override restart values for date
   integer  :: days=10, hours=0, minutes=0, seconds=0
   integer  :: dt_atmos = 3600
   real(r8) :: noise_sd = 0.0_r8
   integer  :: dt_bias  = -1
   logical  :: output_state_vector = .false.  ! output prognostic variables

   namelist /model_nml/ current_time, override, dt_atmos, &
                       days, hours, minutes, seconds, noise_sd, &
                       dt_bias, output_state_vector 

!-----------------------------------------------------------------------
! More stuff from atmos_solo driver
! ----- model time -----

   type (time_type) :: mTime, Time_init, Time_end, Time_step_atmos

! ----- coupled model initial date -----

   integer :: date_init(6)

!-----------------------------------------------------------------------
! Public definition of variable types
integer, parameter :: TYPE_PS = 0, TYPE_T = 1, TYPE_U = 2, TYPE_V = 3, TYPE_TRACER = 4

!-----------------------------------------------------------------------
!---- private data ----

type (bgrid_dynam_type) :: Dynam
! added global_Var here and removed several local instances of Var
! which had space allocated to the data arrays inside Var but were
! then never deallocated (and possibly reallocated each time).
! (Var_dt was already here as a module global variable.)  nsc 31mar06
type    (prog_var_type) :: global_Var, Var_dt

integer                            :: model_size
!real                               :: dt_atmos
real(r8),    dimension(:,:,:), pointer :: omega
integer, dimension(4)              :: atmos_axes
integer                            :: num_levels
integer                            :: ntracers 
! Havana no longer distinguishes dynamic tracers
!!!integer                            :: nprognostic

real(r8), dimension(:), pointer        :: v_lons, v_lats, t_lons, t_lats

!------------------------------------------------------------
! ----- timing flags -----

integer :: id_init, id_loop, id_end
integer, parameter :: timing_level = 1


!-----------------------------------------------------------------------

! Stuff to allow addition of random 'sub-grid scale' noise
!!!!!!type(random_seq_type) :: randseed
logical :: first_call = .true.


contains


!!!#######################################################################
!! this subroutine appears to be unused (was originally from the atmos_solo
!!  driver program).  nsc 31mar06
!!
!! subroutine atmosphere (Var, mTime)
!!
!!
!! type (time_type), intent(in) :: mTime
!! type (prog_var_type), intent(inout) :: Var
!!
!!  type(time_type) :: Time_prev, Time_next
!!!-----------------------------------------------------------------------
!!
!!   Time_prev = mTime                       ! two time-level scheme
!!   Time_next = mTime + Time_step_atmos
!!
!!!---- dynamics -----
!!
!!   call bgrid_core_driver ( Time_next, Var, Var_dt, Dynam, omega )
!!
!!!---- call physics -----
!!
!!   call bgrid_physics ( physics_window, 1.0_r8 * dt_atmos, Time_next,  &
!!                        Dynam%Hgrid,   Dynam%Vgrid,    Dynam,             &
!!                        Var,     Var_dt                       )
!!
!!!---- time differencing and diagnostics -----
!!
!!   call bgrid_core_time_diff ( omega, Time_next, Dynam, Var, Var_dt )
!!
!!!-----------------------------------------------------------------------
!!
!! end subroutine atmosphere
!!
!!!#######################################################################

subroutine adv_1step(x, mTime)

! Does single time-step advance for B-grid model with vector state as
! input and output. This is a modified version of subroutine atmosphere
! in original bgrid_solo_atmosphere driver.


real(r8), intent(inout) :: x(:)
! Time is needed for more general models like this; need to add in to 
! low-order models
type(time_type), intent(in) :: mTime
! removed local instance of Var (since it contains pointers which are
! then left hanging if Var goes out of scope), and replaced all references
! to Var with global_Var below.    nsc 31mar06
!type(prog_var_type), intent(inout) :: Var

type(time_type) :: Time_next

!!!integer :: i, j, k

if ( .not. module_initialized ) call static_init_model

! Convert the vector to a B-grid native state representation
call vector_to_prog_var(x, get_model_size(), global_Var)

! Compute the end of the time interval
Time_next = mTime + Time_step_atmos

! Do dynamics
! Dynam, Var_dt and omega currently in static storage, is that where they should stay?
call bgrid_core_driver(Time_next, global_Var, Var_dt, Dynam, omega)

! Call physics; physics_window is also in global static storage
! dt_atmos is seconds for time step from namelist for now
call bgrid_physics(physics_window, 1.0_r8 * dt_atmos, Time_next, &
   Dynam%Hgrid, Dynam%Vgrid, Dynam, global_Var, Var_dt)

!!!write(*, *) 'max tendency ', maxval(abs(Var_dt%t))
!!!write(*, *) 'mean tendency ' , sum(abs(Var_dt%t)) / (size(Var_dt%t))

! First pass at creating some random 'sub-grid' noise 
! Need to initialize random sequence on first call
if(first_call) then
   !write(*, *) 'NOISE_SD is ', noise_sd
   !!!!!!call init_random_seq(randseed)
   first_call = .false.
endif

!write(*, *) 'in model mod temp_t loop ', size(Var_dt%t, 1), size(Var_dt%t, 2),size(Var_dt%t, 3)
! Just modify T for now by multiplying by 1 + r(0.0, noise_sd) * dt
!do i = 1, size(Var_dt%t, 1)
!   do j = 1, size(Var_dt%t, 2)
!      do k = 1, size(Var_dt%t, 3)
!         temp_t = Var_dt%t(i, j, k)
         !!!!!!Var_dt%t(i, j, k) = &
            !!!!!!(1.0_r8 + random_gaussian(randseed, 0.0_r8, dble(noise_sd))) * Var_dt%t(i, j, k)
!      end do
!   end do
!end do

!!!write(*, *) 'max tendency with noise ', maxval(abs(Var_dt%t))
!!!write(*, *) 'mean tendency with noise ' , sum(abs(Var_dt%t)) / (size(Var_dt%t))

! Time differencing and diagnostics
call bgrid_core_time_diff(omega, Time_next, Dynam, global_Var, Var_dt)

! Convert back to vector form
call prog_var_to_vector(global_Var, x, get_model_size())

end subroutine adv_1step

!#######################################################################

subroutine static_init_model()

! INitializes class data for a B-grid model (all the stuff that needs to
! be done once.

if ( module_initialized ) return ! only need to do this once.

module_initialized = .true.

call fms_init()
call atmos_model_init()

end subroutine static_init_model

!#######################################################################

subroutine init_model_instance(Var)

! Initializes an instance (Var) of a B-grid model state variable

type(prog_var_type), intent(out) :: Var

call prog_var_init(Dynam%Hgrid, num_levels, ntracers, Var)
! Havana no longer distinguishes prognostic tracers
!!!call prog_var_init(Dynam%Hgrid, num_levels, ntracers, nprognostic, Var)

end subroutine init_model_instance

!#######################################################################

subroutine init_conditions(x)


! Reads in restart initial conditions from B-grid and converts to vector

! Following changed to intent(inout) for ifc compiler;should be like this
real(r8), intent(inout) :: x(:)

! removed local instance of Var (since it contains pointers which are
! then left hanging if Var goes out of scope), and replaced all references
! to Var with global_Var below.    nsc 31mar06
!type(prog_var_type) :: Var
real(r8), dimension(Dynam%Hgrid%ilb:Dynam%Hgrid%iub, &
                    Dynam%Hgrid%jlb:Dynam%Hgrid%jub) :: fis, res
real(r8), allocatable, dimension(:) :: eta, peta
integer :: ix, jx, kx
! Havana no longer distinguishes prognostic tracers
!!!integer :: ix, jx, kx, nt, ntp

if ( .not. module_initialized ) call static_init_model

! Need to initialize var???
call init_model_instance(global_Var)

! FOR NOW, TRY TO READ IN CURRENT ICS via read_prog_var
call open_prog_var_file(ix, jx, kx)
! Havana no longer distinguishes prognostic tracers
! Havana no longer returns the number of tracers on this call
!!!call open_prog_var_file(ix, jx, kx, nt, ntp)

allocate(eta(kx + 1), peta(kx + 1))

!!! WARNING; MAY BE DANGEROUS TO USE Dynam%Hgrid here, it could get changed
!!! inappropriately???
call read_prog_var(Dynam%Hgrid, global_Var, eta, peta, fis, res)

deallocate(eta, peta)

call prog_var_to_vector(global_Var, x, get_model_size())

! Probably need to release allocated storage from Var, too???

end subroutine init_conditions

!#######################################################################

! THIS SUBROUTINE WAS ORIGINALLY IN ATMOS_SOLO DRIVER

   subroutine atmos_model_init()

!-----------------------------------------------------------------------

   integer :: iunit, io
   integer :: idate(6)
   type (time_type) :: Run_length
   logical :: use_namelist

!-----------------------------------------------------------------------

   !----- initialization timing identifiers ----

   call register_module(source,revision,revdate)

   !----- read namelist -------
   ! Read the namelist entry
   call find_namelist_in_file("input.nml", "model_nml", iunit)
   read(iunit, nml = model_nml, iostat = io)
   call check_namelist_read(iunit, io, "model_nml")

   !----- write namelist to logfile -----

   call write_version_number (version,tag)
   if (do_nml_file()) write (nmlfileunit, nml=model_nml)
   if (do_nml_term()) write (stdlog(),    nml=model_nml)

   if(dt_atmos == 0) then
     call error_mesg ('program atmos_model', 'dt_atmos has not been specified', FATAL)
   endif

   !----- read restart file -----

   if (file_exist('INPUT/atmos_model.res')) then
       iunit = open_restart_file ('INPUT/atmos_model.res', 'read')
       read  (iunit) idate
       close (iunit)
       use_namelist = .false.
   else
       use_namelist = .true.
   endif

!----- override date with namelist values ------
!----- (either no restart or override flag on) ---

   if ( use_namelist .or. override ) then
      idate(1:2) = 0
      idate(3:6) = current_time
   endif

!----- write current/initial date actually used to logfile file -----

   write (stdlog(),16) idate(3:6)

 16 format ('  current time used = day',i5,' hour',i3,2(':',i2.2))

!-----------------------------------------------------------------------
!------ initialize diagnostics manager ------

!    call diag_manager_init

!----- always override initial/base date with diag_manager value -----

!    call get_base_date ( date_init(1), date_init(2), date_init(3), &
!                         date_init(4), date_init(5), date_init(6)  )

   if ( date_init(1)+date_init(2) /= 0 ) then
        call error_mesg('program atmos_model', 'invalid base base - &
                          &must have year = month = 0', FATAL)
   endif

!----- set initial and current time types ------
!----- set run length and compute ending time -----

    Time_init  = set_time_increment (date_init(3), date_init(4), date_init(5), date_init(6))
    mTime      = set_time_increment (idate    (3), idate    (4), idate    (5), idate    (6))
    Run_length = set_time_increment (days        , hours       , minutes     , seconds     )
    Time_end   = mTime + Run_length

!-----------------------------------------------------------------------
!----- write time stamps (for start time and end time) ------

      iunit = open_file ('time_stamp.out', 'FORMATTED')

       write (iunit,20) idate

!     compute ending time in days,hours,minutes,seconds
      call get_time_increment (Time_end, idate(3), idate(4), idate(5), idate(6))

       write (iunit,20) idate

      close (iunit)

  20  format (6i7,2x,'day')   ! can handle day <= 999999

!-----------------------------------------------------------------------
!----- compute the time steps ------
!----- determine maximum number of iterations per loop ------

      Time_step_atmos = set_time (dt_atmos,0)

!-----------------------------------------------------------------------
!----- initial (base) time must not be greater than current time -----

   if ( Time_init > mTime ) call error_mesg ('program atmos_model',  &
                   'initial time is greater than current time', FATAL)

!-----------------------------------------------------------------------
!------ initialize atmospheric model ------
      call atmosphere_init (Time_init, mTime, Time_step_atmos)

!-----------------------------------------------------------------------
!---- open and close output restart to make sure directory is there ----
! Don't want this test in BGRID if possible
      !iunit = open_restart_file ('RESTART/atmos_model.res', 'write')
      !call close (iunit)


!  ---- terminate timing ----

!-----------------------------------------------------------------------

   end subroutine atmos_model_init



!#######################################################################

 subroutine atmosphere_init (Time_init, mTime, Time_step )

 type (time_type),     intent(in)    :: Time_init, mTime, Time_step

! removed local instance of Var (since it contains pointers which are
! then left hanging if Var goes out of scope), and replaced all references
! to Var with global_Var below.    nsc 31mar06
!type(prog_var_type) :: Var

 integer :: i, iunit, sec, io
 integer :: tnlon, tnlat, vnlon, vnlat
 integer :: t_horiz_size, v_horiz_size

 real(r8), allocatable :: t_lat_bnds(:), t_lon_bnds(:), v_lat_bnds(:), v_lon_bnds(:)

!-----------------------------------------------------------------------
! Read the namelist entry
call find_namelist_in_file("input.nml", "atmosphere_nml", iunit)
read(iunit, nml = atmosphere_nml, iostat = io)
call check_namelist_read(iunit, io, "atmosphere_nml")

!----- write version and namelist to log file -----

   call write_version_number ( version, tag )
   if (do_nml_file()) write (stdlog(),    nml=atmosphere_nml)
   if (do_nml_term()) write (nmlfileunit, nml=atmosphere_nml)

!---- compute physics/atmos time step in seconds ----

   Time_step_atmos = Time_step
   call get_time (Time_step_atmos, sec)
!   dt_atmos = real(sec) LEFT OVER FROM COMBINE OF ATMOS_SOLO DRIVER

!----- initialize dynamical core -----
   call bgrid_core_driver_init ( Time_init, mTime, Time_step,    &
                                 global_Var, Var_dt, Dynam, atmos_axes )

!----- initialize storage needed for vert motion ----

   omega => var_init (Dynam%Hgrid, Dynam%Vgrid%nlev)

!----- initialize physics interface -----

   call hs_forcing_init ( atmos_axes, mTime )

!   ----- use entire grid as window ? -----

   if (physics_window(1) <= 0) physics_window(1) = Dynam%Hgrid%Tmp%ie-Dynam%Hgrid%Tmp%is+1
   if (physics_window(2) <= 0) physics_window(2) = Dynam%Hgrid%Tmp%je-Dynam%Hgrid%Tmp%js+1

!-----------------------------------------------------------------------

! Initialize model_size variables
! Fixed to use only global grid for now (need to modify for mpp)

call get_horiz_grid_size (Dynam % Hgrid, TGRID, tnlon, tnlat, .true.)
call get_horiz_grid_size (Dynam % Hgrid, VGRID, vnlon, vnlat, .true.)
t_horiz_size = tnlon * tnlat
v_horiz_size = vnlon * vnlat
! Num_levels is in global storage
num_levels = Dynam%Vgrid%nlev

! U and V for size
model_size = 2 * num_levels * v_horiz_size
! T and PS for size
model_size = model_size + (1 + num_levels) * t_horiz_size
! Tracers for size
model_size = model_size + global_Var%ntrace * num_levels * t_horiz_size
if (do_output()) write(*, *) 'model_size ', model_size

! Also static store the number of levels, ntracers, and prognostic tracers
ntracers = global_Var%ntrace
! Havana no longer distinguishes prognostic tracers
!!!nprognostic = global_Var%ntprog

! Get the lats and lons of the actual grid points; keep in static store
allocate(t_lons(tnlon), t_lats(tnlat), v_lons(vnlon), v_lats(vnlat))
allocate(t_lon_bnds(tnlon + 1), t_lat_bnds(tnlat + 1), &
         v_lon_bnds(vnlon + 1), v_lat_bnds(vnlat + 1))
call get_horiz_grid_bound(Dynam%Hgrid, TGRID, t_lon_bnds, t_lat_bnds)
call get_horiz_grid_bound(Dynam%Hgrid, VGRID, v_lon_bnds, v_lat_bnds)

! Grid points are at the mid points of the bounds
do i = 1, tnlon
   t_lons(i) = (t_lon_bnds(i) + t_lon_bnds(i + 1)) / 2.0 * (180.0 / 3.14159265)
end do
do i = 1, tnlat
   t_lats(i) = (t_lat_bnds(i) + t_lat_bnds(i + 1)) / 2.0 * (180.0 / 3.14159265)
end do
do i = 1, vnlon
   v_lons(i) = (v_lon_bnds(i) + v_lon_bnds(i + 1)) / 2.0  * (180.0 / 3.14159265)
end do
do i = 1, vnlat
   v_lats(i) = (v_lat_bnds(i) + v_lat_bnds(i + 1)) / 2.0  * (180.0 / 3.14159265)
end do

deallocate(t_lon_bnds, t_lat_bnds, v_lon_bnds, v_lat_bnds)

 end subroutine atmosphere_init

!!!#######################################################################
!! the following subroutines appears to be unused (were originally from 
!!  the atmos_solo driver program).  nsc 31mar06
!!
!! subroutine atmosphere_end(Var)
!!
!! type(prog_var_type), intent(in) :: Var
!!
!!    call bgrid_core_driver_end ( Var, Dynam )
!!
!! end subroutine atmosphere_end
!!
!!!#######################################################################
!!!    returns the number of longitude and latitude grid points
!!!    for either the local PEs grid (default) or the global grid
!!
!! subroutine atmosphere_resolution (nlon, nlat, global)
!!
!!  integer, intent(out)          :: nlon, nlat
!!  logical, intent(in), optional :: global
!!
!!!---- return the size of the grid used for physics computations ----
!!
!!    call get_horiz_grid_size (Dynam % Hgrid, TGRID, nlon, nlat, global)
!!
!! end subroutine atmosphere_resolution
!!
!!!#######################################################################
!!!    returns the longitude and latitude grid box edges
!!!    for either the local PEs grid (default) or the global grid
!!
!! subroutine atmosphere_boundary (blon, blat, global)
!!
!!    real(r8),    intent(out)          :: blon(:), blat(:)
!!    logical, intent(in), optional :: global
!!
!!!----- return the longitudinal and latitudinal grid box edges ----------
!!
!!    call get_horiz_grid_bound (Dynam % Hgrid, TGRID, blon, blat, global)
!!
!! end subroutine atmosphere_boundary
!!
!!!#######################################################################
!!!    returns the axis indices associated with the coupling grid
!!
!! subroutine get_atmosphere_axes ( axes )
!!
!!   integer, intent(out) :: axes (:)
!!
!!!----- returns the axis indices for the atmospheric (mass) grid -----
!!
!!     if ( size(axes) < 0 .or. size(axes) > 4 ) call error_mesg (    &
!!                    'get_atmosphere_axes in atmosphere_mod', &
!!                           'size of argument is incorrect', FATAL   )
!!
!!     axes (1:size(axes)) = atmos_axes (1:size(axes))
!!
!! end subroutine get_atmosphere_axes
!!
!#######################################################################

subroutine bgrid_physics ( window, dt_phys, mTime, Hgrid, Vgrid, &
                           Dynam, Var, Var_dt )

!-----------------------------------------------------------------------
!
!   Time      =  current time (time_type, see time manager)
!
!-----------------------------------------------------------------------
  integer, intent(in)                :: window(2)
  real(r8),    intent(in)                :: dt_phys
       type(time_type),intent(in)    :: mTime
type (horiz_grid_type),intent(inout) :: Hgrid
type  (vert_grid_type),intent(in)    :: Vgrid
type(bgrid_dynam_type),intent(in)    :: Dynam
type   (prog_var_type),intent(in)    :: Var
type   (prog_var_type),intent(inout) :: Var_dt

!-----------------------------------------------------------------------
  integer :: is, ie, js, je, i1, i2, j1, j2, nt
  integer :: ix, jx, dimi, jdim
!-----------------------------------------------------------------------

   real(r8), dimension(window(1),window(2),Vgrid%nlev) :: p_full, u_dt, v_dt

   real(r8), dimension(window(1),window(2),Vgrid%nlev+1) :: p_half

   real(r8), dimension(Hgrid%ilb:Hgrid%iub, &
                   Hgrid%jlb:Hgrid%jub, &
                   Vgrid%nlev) :: uh, vh, uh_dt, vh_dt

   real(r8), dimension(window(1),window(2)) :: pssl_new
!-----------------------------------------------------------------------
!---------------------------- do physics -------------------------------

    dimi = window(1)
    jdim = window(2)

!   --- momentum and momentum tendency on mass grid ---

    call update_halo (Hgrid, UWND, Var_dt%u, flags=SOUTH+WEST)
    call update_halo (Hgrid, VWND, Var_dt%v, flags=SOUTH+WEST)

    call vel_to_mass (Hgrid, Var%u, Var%v,   &
                      uh, vh, Dynam%Masks%Vel%mask)
    call vel_to_mass (Hgrid, Var_dt %u, Var_dt %v, &
                      uh_dt, vh_dt, Dynam%Masks%Vel%mask)

!   --- loop through physics windows ---

    nt = Var%ntrace
    js = Hgrid%Tmp%js

    do while ( js <= Hgrid%Tmp%je )

       je = min ( js+jdim-1, Hgrid%Tmp%je )
       jx = je-js+1
       is = Hgrid%Tmp%is

    do while ( is <= Hgrid%Tmp%ie )

       ie = min ( is+dimi-1, Hgrid%Tmp%ie )
       ix = ie-is+1

!      ---- pass updated surface pressure ----
       pssl_new(1:ix,1:jx) = Var%pssl(is:ie,js:je) + &
                             Var_dt%pssl(is:ie,js:je) * dt_phys

       call compute_pres_full (Vgrid, pssl_new(1:ix,1:jx), &
                                        p_full(1:ix,1:jx,:))
       call compute_pres_half (Vgrid, pssl_new(1:ix,1:jx), &
                                        p_half(1:ix,1:jx,:))


       u_dt(1:ix,1:jx,:) = uh_dt(is:ie,js:je,:)
       v_dt(1:ix,1:jx,:) = vh_dt(is:ie,js:je,:)


!      ---- j-axis indices in the global physics grid ----

       j1 = js-Hgrid%Tmp%js+1; j2 = j1+(je-js)
       i1 = is-Hgrid%Tmp%is+1; i2 = i1+(ie-is)

!-----------------------------------------------------------------------
!-------------------------- call physics -------------------------------
!------------ (need to add leap-frog option for uh,vh) -----------------
!-----------------------------------------------------------------------
  if (.not.Dynam%Masks%sigma) then
!------------ eta coordinate -------------------------------------------

      call hs_forcing ( i1, i2, j1, j2, dt_phys, mTime ,&
                            Hgrid%Tmp%aph(is:ie,js:je)    ,&
                            p_half   ( 1:ix, 1:jx,:)      ,&
                            p_full   ( 1:ix, 1:jx,:)      ,&
                            uh       (is:ie,js:je,:)      ,&
                            vh       (is:ie,js:je,:)      ,&
                            Var%t    (is:ie,js:je,:)      ,&
                            Var%r    (is:ie,js:je,:,:)    ,&
                            uh       (is:ie,js:je,:)      ,&
                            vh       (is:ie,js:je,:)      ,&
                            Var%t    (is:ie,js:je,:)      ,&
                            Var%r    (is:ie,js:je,:,:)    ,&
                            u_dt     ( 1:ix, 1:jx,:)      ,&
                            v_dt     ( 1:ix, 1:jx,:)      ,&
                            Var_dt%t (is:ie,js:je,:)      ,&
                            Var_dt%r (is:ie,js:je,:,:)    ,&
                            mask=Dynam%Masks%Tmp%mask(is:ie,js:je,:) ,&
                            kbot=Dynam%Masks%Tmp%kbot(is:ie,js:je)    )

  else
!------------- sigma coordinate ----------------------------------------

      call hs_forcing ( i1, i2, j1, j2, dt_phys, mTime ,&
                            Hgrid%Tmp%aph(is:ie,js:je)    ,&
                            p_half   ( 1:ix, 1:jx,:)      ,&
                            p_full   ( 1:ix, 1:jx,:)      ,&
                            uh       (is:ie,js:je,:)      ,&
                            vh       (is:ie,js:je,:)      ,&
                            Var%t    (is:ie,js:je,:)      ,&
                            Var%r    (is:ie,js:je,:,:)    ,&
                            uh       (is:ie,js:je,:)      ,&
                            vh       (is:ie,js:je,:)      ,&
                            Var%t    (is:ie,js:je,:)      ,&
                            Var%r    (is:ie,js:je,:,:)    ,&
                            u_dt     ( 1:ix, 1:jx,:)      ,&
                            v_dt     ( 1:ix, 1:jx,:)      ,&
                            Var_dt%t (is:ie,js:je,:)      ,&
                            Var_dt%r (is:ie,js:je,:,:)     )

  endif

!       ---- physics tendency on mass grid ----
        uh(is:ie,js:je,:) = u_dt(1:ix,1:jx,:) - uh_dt(is:ie,js:je,:)
        vh(is:ie,js:je,:) = v_dt(1:ix,1:jx,:) - vh_dt(is:ie,js:je,:)

        is = is + dimi

     enddo

        js = js + jdim

     enddo

!-----------------------------------------------------------------------
!---- move momentum tendencies from mass to momentum points ----
!---- zero out unused polar row, no harm if not polar row ----

     call update_halo (Hgrid, TEMP, uh, flags=NORTH+EAST)
     call update_halo (Hgrid, TEMP, vh, flags=NORTH+EAST)

     call mass_to_vel (Hgrid, uh, uh)
     call mass_to_vel (Hgrid, vh, vh)

     uh(:,Hgrid%jub,:) = 0.0_r8
     vh(:,Hgrid%jub,:) = 0.0_r8

!---- update momentum tendencies ----

     Var_dt%u = Var_dt%u + uh * Dynam%Masks%Vel%mask
     Var_dt%v = Var_dt%v + vh * Dynam%Masks%Vel%mask

!---- update halo rows ----

     call update_halo (Hgrid, TEMP, Var_dt%t)
     call update_halo (Hgrid, TEMP, Var_dt%r)
     call update_halo (Hgrid, UWND, Var_dt%u)
     call update_halo (Hgrid, VWND, Var_dt%v)

!-----------------------------------------------------------------------

end subroutine bgrid_physics

!#######################################################################

function get_model_size()

integer :: get_model_size

if ( .not. module_initialized ) call static_init_model

get_model_size = model_size

end function get_model_size

!#######################################################################

subroutine prog_var_to_vector(vars, x, isize)


integer, intent(in) :: isize
type(prog_var_type), intent(in) :: vars
real(r8), intent(out) :: x(isize)

integer :: i, j, k, nt, indx
integer :: num_levs
integer :: tis, tie, tjs, tje, vis, vie, vjs, vje
character(len=129) :: errstring

! Get the bounds for storage on Temp and Velocity grids
tis = Dynam%Hgrid%Tmp%is; tie = Dynam%Hgrid%Tmp%ie
tjs = Dynam%Hgrid%Tmp%js; tje = Dynam%Hgrid%Tmp%je
vis = Dynam%Hgrid%Vel%is; vie = Dynam%Hgrid%Vel%ie
vjs = Dynam%Hgrid%Vel%js; vje = Dynam%Hgrid%Vel%je
num_levs = vars%kub - vars%klb + 1

! Start copying fields to straight vector
! Do everything on the T grid first
indx = 0
do i = tis, tie
   do j = tjs, tje
! Surface pressure is first
      indx = indx + 1
      x(indx) = vars%ps(i, j)
! Now do t and tracers at successive levels
      do k = vars%klb, vars%kub
         indx = indx + 1
         x(indx) = vars%t(i, j, k)
         do nt = 1, vars%ntrace
            indx = indx + 1
            x(indx) = vars%r(i, j, k, nt)
         end do
      end do
   end do
end do

! Now do the velocity grid, u and v
do i = vis, vie
   do j = vjs, vje
      do k = vars%klb, vars%kub
         indx = indx + 1
         x(indx) = vars%u(i, j, k)
         indx = indx + 1
         x(indx) = vars%v(i, j, k)
      end do
   end do
end do
   
! Temporary check
if(indx /= isize) then
   write(errstring, *) 'indx (',indx,') should be equal to (',isize,') isize'
   call error_handler(E_ERR,'prog_var_to_vector', errstring, source, revision, revdate)
endif

end subroutine prog_var_to_vector

!#######################################################################

subroutine vector_to_prog_var(x, isize, vars)


integer, intent(in) :: isize
! Following line fouls up the ibm/mac compiler; doesn't like to get
! A fixed size argument from an assumed shape.
!real(r8), intent(in) :: x(isize)
real(r8), intent(in) :: x(:)
type(prog_var_type), intent(inout) :: vars

integer :: i, j, k, nt, indx
integer :: num_levs
integer :: tis, tie, tjs, tje, vis, vie, vjs, vje

character(len=129) :: errstring

! Initialize the static parts of the prog var type
! Don't want to initialize prog_var_type vars if already done (memory)
! Modified on 11 Dec. 2002 to slow B-grid memory leak
if(.not. associated(vars%ps)) &
   call prog_var_init(Dynam%Hgrid, num_levels, ntracers, vars)

! Havana no longer distinguishes prognostic tracers
!!!call prog_var_init(Dynam%Hgrid, num_levels, ntracers, nprognostic, vars)

! Get the bounds for storage on Temp and Velocity grids
tis = Dynam%Hgrid%Tmp%is; tie = Dynam%Hgrid%Tmp%ie
tjs = Dynam%Hgrid%Tmp%js; tje = Dynam%Hgrid%Tmp%je
vis = Dynam%Hgrid%Vel%is; vie = Dynam%Hgrid%Vel%ie
vjs = Dynam%Hgrid%Vel%js; vje = Dynam%Hgrid%Vel%je
num_levs = vars%kub - vars%klb + 1

! Start copying fields from straight vector
! Everything on T grid first
indx = 0
do i = tis, tie
   do j = tjs, tje
! Surface pressure is first
      indx = indx + 1
      vars%ps(i, j) = x(indx)
! For non-eta models, pssl is same as ps??? Need to change?
      vars%pssl(i, j) = vars%ps(i, j)
! Now do t and tracers at successive levels
      do k = vars%klb, vars%kub
         indx = indx + 1
         vars%t(i, j, k) = x(indx)
         do nt = 1, vars%ntrace
            indx = indx + 1
            vars%r(i, j, k, nt) = x(indx)
         end do
      end do
   end do
end do

! Now do the velocity grid, u and v
do i = vis, vie
   do j = vjs, vje
      do k = vars%klb, vars%kub
         indx = indx + 1
         vars%u(i, j, k) = x(indx)
         indx = indx + 1
         vars%v(i, j, k) = x(indx)
      end do
   end do
end do

! Need to do halo updates to fill in the rest
!!!call halo_update()
call update_halo (Dynam%Hgrid, UWND, vars%u)
call update_halo (Dynam%Hgrid, VWND, vars%v)
call update_halo (Dynam%Hgrid, TEMP, vars%t)
call update_halo (Dynam%Hgrid, TEMP, vars%r)
call update_halo (Dynam%Hgrid, TEMP, vars%ps)
call update_halo (Dynam%Hgrid, TEMP, vars%pssl)


! Temporary check
if(indx /= isize) then
   write(errstring, *) 'indx (',indx,') should be equal to (',isize,') isize'
   call error_handler(E_ERR,'vector_to_prog_var', errstring, source, revision, revdate)
endif


end subroutine vector_to_prog_var

!#######################################################################

function get_model_time_step()
!------------------------------------------------------------------------
! function get_model_time_step()
!
! Returns the the time step of the model. In the long run should be repalced
! by a more general routine that returns details of a general time-stepping
! capability.

type(time_type) :: get_model_time_step

if ( .not. module_initialized ) call static_init_model

! Time_step_atmos is global static storage
get_model_time_step =  Time_step_atmos


! CODE ADDED TO SIMULATE MODEL ERROR BY MAKING MODEL THINK IT IS
! ADVANCING AT A DIFFERENT RATE THAN IT IS
! THIS auxiliary timestep is only used if namelist parameter
! dt_bias is set to a non-zero value 

if(dt_bias > 0) get_model_time_step = set_time(dt_bias, 0)

end function get_model_time_step



  subroutine get_state_meta_data(index_in, location, var_type)
!---------------------------------------------------------------------
! subroutine get_state_meta_data(index_in, location, var_type)
!
! Given an integer index into the state vector structure, returns the
! associated location. This is not a function because the more general
! form of the call has a second intent(out) optional argument kind.
! Maybe a functional form should be added?
! Types for this bgrid model are, TYPE_PS, TYPE_T, TYPE_U, TYPE_V, TYPE_TRACER

integer,             intent(in)  :: index_in
type(location_type), intent(out) :: location
integer, optional,   intent(out) :: var_type

integer :: indx, local_var_type, var_type_temp
integer :: tis, tie, tjs, tje, vis, vie, vjs, vje
integer :: t_size, v_size, t_grid_size, v_grid_size, t_per_col, v_per_col
integer :: num_t_lons, num_t_lats, num_v_lons, num_v_lats
integer :: col_num, col_elem, v_index
integer :: lat_index, lon_index
real(r8) :: lon, lat, lev

if ( .not. module_initialized ) call static_init_model

! Get the bounds for storage on Temp and Velocity grids
tis = Dynam%Hgrid%Tmp%is; tie = Dynam%Hgrid%Tmp%ie
num_t_lons = tie - tis + 1
tjs = Dynam%Hgrid%Tmp%js; tje = Dynam%Hgrid%Tmp%je
num_t_lats = tje - tjs + 1
t_size = num_t_lons *num_t_lats
vis = Dynam%Hgrid%Vel%is; vie = Dynam%Hgrid%Vel%ie
num_v_lons = vie - vis + 1
vjs = Dynam%Hgrid%Vel%js; vje = Dynam%Hgrid%Vel%je
num_v_lats = vje - vjs + 1
v_size = num_v_lons * num_v_lats


! Compute size of t_grid storage
t_per_col = 1 + (1 + ntracers) * num_levels
t_grid_size = t_per_col * t_size
v_per_col = 2 * num_levels
v_grid_size = v_per_col * v_size

! Easier to compute with a 0 to size - 1 index
indx = index_in - 1

! Is this point in the t_grid
if(indx < t_grid_size) then
   col_num = indx / t_per_col
   col_elem = indx - col_num * t_per_col
   !   write(*, *) 't_grid col and element ', col_num, col_elem
   lon_index = col_num / num_t_lats
   lat_index = col_num - lon_index * num_t_lats
   !   write(*, *) 'lon and lat index ', lon_index, lat_index

   ! Note, array runs from 1, index runs from 0
   lon = t_lons(lon_index + 1)
   lat = t_lats(lat_index + 1)

   if(col_elem == 0) then ! First variable on mass grid is PS
      lev = -1
      local_var_type = TYPE_PS 
   else ! Rest of variables are temperature and tracers

      lev = int(col_elem / (1 + ntracers))
      var_type_temp = mod(col_elem - 1, 1 + ntracers)

      ! First element on each level is T, 
      ! remainder are tracers from 1 to ntracers
      if(var_type_temp == 0) then
         local_var_type = TYPE_T
      else
         local_var_type = TYPE_TRACER + var_type_temp - 1
      endif
   endif

else

   ! It's in the v_grid
   v_index = indx - t_grid_size
   col_num = v_index / v_per_col
   col_elem = v_index - col_num * v_per_col

   !   write(*, *) 'v_grid col and element ', col_num, col_elem
   lon_index = col_num / num_v_lats
   lat_index = col_num - lon_index * num_v_lats

   !   write(*, *) 'lon and lat index ', lon_index, lat_index

   ! Note, array runs from 1, index runs from 0
   lon = v_lons(lon_index + 1)

   ! Problems with round-off over 360.0
   if(abs(lon - 360.0_r8) < 0.00001_r8) lon = 360.0_r8

   lat = v_lats(lat_index + 1)
   lev = int((col_elem + 2) / 2)

   ! Compute u or v, u is even, v is odd
   if(col_elem / 2 * 2 == col_elem) then
      local_var_type = TYPE_U
   else
      local_var_type = TYPE_V
   endif
endif

!write(*, *) 'lon, lat, and lev ', lon, lat, lev

location = set_location(lon, lat, lev, 1) ! 1 == level (indexical)

! If the type is wanted, return it
if(present(var_type)) var_type = local_var_type


end subroutine get_state_meta_data



recursive subroutine model_interpolate(x, location, itype, obs_val, istatus)

real(r8),            intent(in) :: x(:)
type(location_type), intent(in) :: location
integer,             intent(in) :: itype
real(r8),            intent(out):: obs_val
integer,             intent(out):: istatus

integer :: num_lons, num_lats, lon_below, lon_above, lat_below, lat_above, i
real(r8) :: bot_lon, top_lon, delta_lon, bot_lat, top_lat, delta_lat
real(r8) :: lon_fract, lat_fract, val(2, 2), temp_lon, a(2)
real(r8) :: lon, lat, level, lon_lat_lev(3), pressure

if ( .not. module_initialized ) call static_init_model

! All interps okay for now
istatus = 0

! Would it be better to pass state as prog_var_type (model state type) to here?
! As opposed to the stripped state vector. YES. This would give time interp.
! capability here; do we really want this or should it be pushed up?

! Get the position, determine if it is model level or pressure in vertical
lon_lat_lev = get_location(location)
lon = lon_lat_lev(1); lat = lon_lat_lev(2); 
if(vert_is_level(location)) then 
   level = lon_lat_lev(3)
else if(vert_is_pressure(location)) then
   pressure = lon_lat_lev(3)
else if(vert_is_surface(location)) then
   ! level is not used for surface pressure observations
   level = -1
else
   call error_handler(E_ERR,'model_interpolate', &
      'Bgrid can only handle pressure or model level for obs vertical coordinate', &
      source, revision, revdate)
endif
   
! Depending on itype, get appropriate lon and lat grid specs
! Types temporarily defined as 1=u, 2=v, 3=ps, 4=t, n=tracer number n-4
if(itype == KIND_U_WIND_COMPONENT .or. itype == KIND_V_WIND_COMPONENT) then
   num_lons = size(v_lons)
   num_lats = size(v_lats)
   bot_lon = v_lons(1)
   top_lon = v_lons(num_lons)
   delta_lon = v_lons(2) - v_lons(1)
   bot_lat = v_lats(1)
   top_lat = v_lats(num_lats)
   delta_lat = v_lats(2) - v_lats(1)
else
   num_lons = size(t_lons)
   num_lats = size(t_lats)
   bot_lon = t_lons(1)
   top_lon = t_lons(num_lons)
   delta_lon = t_lons(2) - t_lons(1)
   bot_lat = t_lats(1)
   top_lat = t_lats(num_lats)
   delta_lat = t_lats(2) - t_lats(1)
endif

! Compute bracketing lon indices
if(lon >= bot_lon .and. lon <= top_lon) then
   lon_below = int((lon - bot_lon) / delta_lon) + 1
   lon_above = lon_below + 1
!   write(*, *) 'lon, delta_lon, bot_lon', lon, delta_lon, bot_lon
!   write(*, *) 'prod ', ((lon_below - 1) * delta_lon + bot_lon)
   lon_fract = (lon - ((lon_below - 1) * delta_lon + bot_lon)) / delta_lon
else
! At wraparound point
   lon_below = num_lons
   lon_above = 1
   if(lon < bot_lon) then 
      temp_lon = lon + 360.0_r8
   else
      temp_lon = lon
   endif
   lon_fract = (temp_lon - top_lon) / delta_lon
endif
   

! Next, compute neighboring lat rows
! NEED TO BE VERY CAREFUL ABOUT POLES; WHAT'S BEING DONE MAY BE WRONG
if(lat >= bot_lat .and. lat <= top_lat) then
   lat_below = int((lat - bot_lat) / delta_lat) + 1
   lat_above = lat_below + 1
   lat_fract = (lat - ((lat_below - 1) * delta_lat + bot_lat)) / delta_lat
else if(lat <= bot_lat) then
! South of bottom lat NEED TO DO BETTER: NOT REALLY REGULAR
   lat_below = 1
   lat_above = 2
   lat_fract = 0.0_r8
else
! North of top lat NEED TO DO BETTER: NOT REALLY REGULAR
   lat_below = num_lats - 1
   lat_above = num_lats
   lat_fract = 1.0_r8
endif

! Case 1: model level specified in vertical
if(vert_is_level(location) .or. vert_is_surface(location)) then
! Now, need to find the values for the four corners
   val(1, 1) =  get_val(x, lon_below, lat_below, nint(level), itype)
   val(1, 2) =  get_val(x, lon_below, lat_above, nint(level), itype)
   val(2, 1) =  get_val(x, lon_above, lat_below, nint(level), itype)
   val(2, 2) =  get_val(x, lon_above, lat_above, nint(level), itype)

else
! Case of pressure specified in vertical
   val(1, 1) =  get_val_pressure(x, lon_below, lat_below, pressure, itype, istatus)
   val(1, 2) =  get_val_pressure(x, lon_below, lat_above, pressure, itype, istatus)
   val(2, 1) =  get_val_pressure(x, lon_above, lat_below, pressure, itype, istatus)
   val(2, 2) =  get_val_pressure(x, lon_above, lat_above, pressure, itype, istatus)
endif

! Do the weighted average for interpolation
!write(*, *) 'fracts ', lon_fract, lat_fract
do i = 1, 2
   a(i) = lon_fract * val(2, i) + (1.0 - lon_fract) * val(1, i)
end do

obs_val = lat_fract * a(2) + (1.0 - lat_fract) * a(1)

! Since there is no possibility of having obs_val == missing_r ...
! the return codes are always "good" i.e. zero
!
! normally set istatus here 

end subroutine model_interpolate

!#######################################################################

function get_val(x, lon_index, lat_index, level, itype)


real(r8) :: get_val
real(r8), intent(in) :: x(:)
integer, intent(in) :: lon_index, lat_index, level, itype

character(len = 129) :: msg_string
integer :: model_type

! Need to change the obs kind defined itype to the appropriate model type
! The itype_in variable uses types defined in the kinds module. The whole bgrid 
! model_mod should be modified to use this correctly. However, as a fast patch
! for the initial I-release, will just map the u, v, t, and ps kinds to the
! default expectations of the bgrid which are hard coded as u=1, v=2, ps =3,
! t = 4, and tracers for numbers larger than 4. For now, the tracer observations
! are not implemented.
if(itype == KIND_U_WIND_COMPONENT) then
   model_type = 1
else if(itype == KIND_V_WIND_COMPONENT) then
   model_type = 2
else if(itype == KIND_SURFACE_PRESSURE) then
   model_type = 3
else if(itype == KIND_TEMPERATURE) then
   model_type = 4
else
   ! Error for higher or lower for now
   write(msg_string, *) 'Only know how to do u, v, ps, t observations'
   call error_handler(E_ERR, 'get_val', msg_string, &
      source, revision, revdate)
endif 

! Find the index into state array and return this value
get_val = x(get_state_index(lon_index, lat_index, level, model_type))

end function get_val



  recursive function get_val_pressure(x, lon_index, lat_index, pressure, itype, istatus)
!================================================================================
! function get_val_pressure(x, lon_index, lat_index, pressure, itype, istatus)
!
! Gets the vertically interpolated value on pressure for variable type
! at lon_index, lat_index horizontal grid point

real(r8) :: get_val_pressure
real(r8), intent(in) :: x(:), pressure
integer,  intent(in) :: lon_index, lat_index, itype
integer, intent(out) :: istatus

type(location_type) :: ps_location
real(r8) :: ps(1, 1), pfull(1, 1, Dynam%Vgrid%nlev), rfrac
integer  :: top_lev, bot_lev, i
real(r8) :: bot_val, top_val, ps_lon

! Need to get the surface pressure at this point.
! For t or tracers (on mass grid with ps) this is trivial
! For u or v (on velocity grid)

if(itype == KIND_TEMPERATURE .or. itype == KIND_SURFACE_PRESSURE) then

   ps = get_val(x, lon_index, lat_index, -1, KIND_SURFACE_PRESSURE)

else

   ! Bgrid model has nasty habit of having longitudes truncated to > 360.0
   ! Need to fix this here to avoid death in location module
   ps_lon = v_lons(lon_index)
   if(ps_lon > 360.00_r8 .and. ps_lon < 360.00001_r8) ps_lon = 360.0_r8

   ! The vertical is not important for this interpolation -- still --
   ! mark it as missing (-1.0) but give it some type information (2==pressure)
   ps_location = set_location(ps_lon, v_lats(lat_index), -1.0_r8, 2 )
   call model_interpolate(x, ps_location, KIND_SURFACE_PRESSURE, ps(1,1), istatus)

endif

! Next, get the values on the levels for this ps
call compute_pres_full(Dynam%Vgrid, ps, pfull)

! Interpolate in vertical to get two bounding levels

!write(*, *) 'itype is ', itype
!write(*, *) 'pfull values for ps = ', ps(1, 1)
!write(*, *) pfull(1, 1, :)

! What to do about pressures above top??? Just use the top for now.
! Could extrapolate, but that would be very tricky. Might need to
! reject somehow.
!write(*, *) 'pressure , ps ', pressure, ps(1, 1)
if(pressure < pfull(1, 1, 1)) then
   top_lev = 1
   bot_lev = 2
   rfrac = 1.0_r8
   ! Actually, just fail using istatus
   istatus = 1
else if(pressure > pfull(1, 1, Dynam%Vgrid%nlev)) then
! Same for bottom
   bot_lev = Dynam%Vgrid%nlev
   top_lev = bot_lev - 1
   rfrac = 0.0_r8
   ! Actually, just fail using istatus
   istatus = 1
else

! Search down through pressures
   do i = 2, Dynam%Vgrid%nlev
      if(pressure < pfull(1, 1, i)) then
         top_lev = i -1
         bot_lev = i
         rfrac = (pfull(1, 1, i) - pressure) / &
            (pfull(1, 1, i) - pfull(1, 1, i - 1))
         goto 21
      endif
   end do
end if

! Get the value at these two points
21 bot_val = get_val(x, lon_index, lat_index, bot_lev, itype)
top_val = get_val(x, lon_index, lat_index, top_lev, itype)
!write(*, *) 'bot_lev, top_lev, fraction', bot_lev, top_lev, rfrac
   
get_val_pressure = (1.0_r8 - rfrac) * bot_val + rfrac * top_val
!write(*, *) 'bot_val, top_val, val', bot_val, top_val, get_val_pressure

end function get_val_pressure

!#######################################################################

function get_state_index(lon_index, lat_index, level, itype)


integer :: get_state_index
integer, intent(in) :: lon_index, lat_index, level, itype

! Returns the index in the state vector for variable of itype at
! the lon_index, lat_index, level position on the grid
! Types are currently hard-coded as u=1, v=2, ps=3, t=4, tracers = n-4
integer :: tis, tie, num_t_lons, tjs, tje, num_t_lats, t_size
integer :: vis, vie, num_v_lons, vjs, vje, num_v_lats, v_size
integer :: t_per_col, t_grid_size, v_per_col, v_grid_size

! Get the bounds for storage on Temp and Velocity grids
tis = Dynam%Hgrid%Tmp%is; tie = Dynam%Hgrid%Tmp%ie
num_t_lons = tie - tis + 1
tjs = Dynam%Hgrid%Tmp%js; tje = Dynam%Hgrid%Tmp%je
num_t_lats = tje - tjs + 1
t_size = num_t_lons *num_t_lats
vis = Dynam%Hgrid%Vel%is; vie = Dynam%Hgrid%Vel%ie
num_v_lons = vie - vis + 1
vjs = Dynam%Hgrid%Vel%js; vje = Dynam%Hgrid%Vel%je
num_v_lats = vje - vjs + 1
v_size = num_v_lons * num_v_lats


! Compute size of t_grid storage
t_per_col = 1 + (1 + ntracers) * num_levels
t_grid_size = t_per_col * t_size
v_per_col = 2 * num_levels
v_grid_size = v_per_col * v_size

! Is this point in the t_grid; ps, t or tracer
if(itype > 2) then
   get_state_index = t_per_col * (lat_index - 1 + (lon_index - 1) * num_t_lats)
   if(itype == 3) then
      get_state_index = get_state_index + 1
   else 
      get_state_index = get_state_index + &
         1 + (level - 1) * (1 + ntracers) + (itype - 3)
   endif

! Type 1 or 2 is u or v
else
! It's in the v_grid
   get_state_index = t_grid_size + &
      v_per_col * (lat_index - 1 + (lon_index - 1) * num_v_lats)
   get_state_index = get_state_index + (level - 1) * 2 + itype
endif

end function get_state_index

!#######################################################################

subroutine end_model()


! At some point, this stub should coordinate with atmosphere_end but
! that requires an instance variable.

end subroutine end_model

!#######################################################################
! routines to set/get date when no calendar is set (i.e., yr=0 and mo=0)
!#######################################################################
! return the time increment for the given
! number of days, hours, minutes, and seconds

 function Set_time_increment ( d, h, m, s )
 integer, intent(in) :: d, h, m, s
 type(time_type) :: Set_time_increment

   Set_time_increment = set_time ( h*3600+m*60+s, d )

 end function Set_time_increment

!#######################################################################
! compute time in days, hours, minutes, seconds ----

 subroutine get_time_increment ( T, d, h, m, s )
 type(time_type), intent(in)  :: T
 integer,         intent(out) :: d, h, m, s

   call get_time ( T, s, d )

 ! compute hours and minutes
   h = s/3600 ;   s = s - h*3600
   m = s/60   ;   s = s - m*60

 end subroutine get_time_increment

!#######################################################################

subroutine init_time(i_time)


! For now returns value of Time_init which is set in initialization routines.

type(time_type), intent(out) :: i_time

if ( .not. module_initialized ) call static_init_model

i_time = mTime

end subroutine init_time



!--------------------------------------------------------------------




function nc_write_model_atts( ncFileID ) result (ierr)
!-----------------------------------------------------------------------------------------
! Writes the model-specific attributes to a netCDF file
! TJH Dec 5 2002
!
! TJH 29 July 2003 -- for the moment, all errors are fatal, so the
! return code is always '0 == normal', since the fatal errors stop execution.
!                                                                                
! There are two different (staggered) 3D grids being used simultaneously here. 
! The routine "prog_var_to_vector" packs the prognostic variables into
! the requisite array for the data assimilation routines. That routine
! is the basis for the information stored in the netCDF files.
!
! TemperatureGrid : surface pressure  vars%ps(tis:tie, tjs:tje) 
!                 : temperature       vars%t (tis:tie, tjs:tje, klb:kup)
!                 : tracers           vars%r (tis:tie, tjs:tje, klb:kub, 1:vars%ntrace)
! VelocityGrid    : u                 vars%u (vis:vie, vjs:vje, klb:kub) 
!                 : v                 vars%v (vis:vie, vjs:tje, klb:kup)
!
! So there are six different dimensions and five different variables as long as
! simply lump "tracers" into one. 
!
! TJH 22 May 2003 -- It is now possible to output the prognostic variables 
! _or_ the state vector. If a "state" variable exists, we do nothing. 
! If it does not exist, we can assume we need to define the prognostic variables.

use typeSizes
use netcdf

integer, intent(in)  :: ncFileID        ! netCDF file identifier
integer              :: ierr            ! return value of function

!-----------------------------------------------------------------------------------------

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID
integer :: TmpIDimID, TmpJDimID, levDimID, tracerDimID, VelIDimID, VelJDimID, MemberDimID
integer :: TmpIVarID, TmpJVarID, levVarID, tracerVarID, VelIVarID, VelJVarID, StateVarID
integer :: StateVarDimID, StateVarVarID, TimeDimID
integer :: psVarID, tVarID, rVarID, uVarID, vVarID
integer :: tis, tie, tjs, tje       ! temperature grid start/stop
integer :: vis, vie, vjs, vje       ! velocity    grid start/stop
integer :: klb, kub
integer :: nTmpI, nTmpJ, nVelI, nVelJ, nlev, ntracer, i

character(len=129) :: errstring

character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic
character(len=NF90_MAX_NAME) :: str1
!-----------------------------------------------------------------------------------------

if ( .not. module_initialized ) call static_init_model

ierr = 0     ! assume normal termination

!-------------------------------------------------------------------------------
! Get the bounds for storage on Temp and Velocity grids
! I hate using these in this manner.
!-------------------------------------------------------------------------------

tis = Dynam%Hgrid%Tmp%is; tie = Dynam%Hgrid%Tmp%ie
tjs = Dynam%Hgrid%Tmp%js; tje = Dynam%Hgrid%Tmp%je
vis = Dynam%Hgrid%Vel%is; vie = Dynam%Hgrid%Vel%ie
vjs = Dynam%Hgrid%Vel%js; vje = Dynam%Hgrid%Vel%je

!write(*,*)'klb = ',klb,'  kub = ',kub

nTmpI   = tie - tis + 1
nTmpJ   = tje - tjs + 1
nlev    = Var_dt%kub - Var_dt%klb + 1
ntracer = Var_dt%ntrace 
nVelI   = vie - vis + 1
nVelJ   = vje - vjs + 1

!-------------------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file, 
! and then put into define mode.
!-------------------------------------------------------------------------------

call check(nf90_Inquire(ncFileID, nDimensions, nVariables, nAttributes, unlimitedDimID), &
           "inquire")
call check(nf90_Redef(ncFileID),"redef")

!-------------------------------------------------------------------------------
! We need the dimension ID for the number of copies 
!-------------------------------------------------------------------------------

call check(nf90_inq_dimid(ncid=ncFileID, name="copy", dimid=MemberDimID),"copy dimid")
call check(nf90_inq_dimid(ncid=ncFileID, name="time", dimid=  TimeDimID),"time dimid")

if ( TimeDimID /= unlimitedDimId ) then
   write(errstring,*)"Time Dimension ID ",TimeDimID, &
                     " should equal Unlimited Dimension ID",unlimitedDimID
   call error_handler(E_ERR,"nc_write_model_atts", errstring, source, revision, revdate)
endif

!-------------------------------------------------------------------------------
! Define the model size, state variable dimension ... whatever ...
!-------------------------------------------------------------------------------
call check(nf90_def_dim(ncid=ncFileID, name="StateVariable", &
                        len=model_size, dimid = StateVarDimID),"state def_dim")

!-------------------------------------------------------------------------------
! Write Global Attributes 
!-------------------------------------------------------------------------------

call DATE_AND_TIME(crdate,crtime,crzone,values)
write(str1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
                  values(1), values(2), values(3), values(5), values(6), values(7)

call check(nf90_put_att(ncFileID, NF90_GLOBAL, "creation_date" ,str1    ),"creation put")
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_source"  ,source  ),"source put")
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revision",revision),"revision put")
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revdate" ,revdate ),"revdate put")
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model","FMS_Bgrid"      ),"model put")

! how about namelist input? might be nice to save ...

!-------------------------------------------------------------------------------
! Define the new dimensions IDs
!-------------------------------------------------------------------------------

call check(nf90_def_dim(ncid=ncFileID, name="TmpI", &
                                      len = nTmpI, dimid = TmpIDimID),"TmpI def_dim") 
call check(nf90_def_dim(ncid=ncFileID, name="TmpJ", &
                                      len = nTmpJ, dimid = TmpJDimID),"TmpJ def_dim") 
call check(nf90_def_dim(ncid=ncFileID, name="lev",  &
                                      len = nlev,  dimid = levDimID), "lev  def_dim") 
call check(nf90_def_dim(ncid=ncFileID, name="VelI", &
                                      len = nVelI, dimid = VelIDimID),"VelI def_dim") 
call check(nf90_def_dim(ncid=ncFileID, name="VelJ", &
                                      len = nVelJ, dimid = VelJDimID),"VelJ def_dim") 
if ( ntracer > 0 ) then
   call check(nf90_def_dim(ncid=ncFileID, name="tracers", &
                                         len = ntracer, dimid = tracerDimID),"tracer def_dim") 
endif

! should implement "trajectory-like" coordinate defn ... a'la section 5.4, 5.5 of CF standard
! call check(nf90_def_dim(ncid=ncFileID, name="locationrank", &
!   len = LocationDims, dimid = LocationDimID))

!-------------------------------------------------------------------------------
! Create the (empty) Variables and the Attributes
!-------------------------------------------------------------------------------

! Temperature Grid Longitudes
call check(nf90_def_var(ncFileID, name="TmpI", &
              xtype=nf90_double, dimids=TmpIDimID, varid=TmpIVarID),   "TmpI def_var")
call check(nf90_put_att(ncFileID, TmpIVarID, "long_name", "longitude"),"TmpI long_name")
call check(nf90_put_att(ncFileID, TmpIVarID, "cartesian_axis", "X"),   "TmpI cartesian_axis")
call check(nf90_put_att(ncFileID, TmpIVarID, "units", "degrees_east"), "TmpI units")
call check(nf90_put_att(ncFileID, TmpIVarID, "valid_range", (/ 0.0_r8, 360.0_r8 /)), &
                                 "TmpI valid_range")

! Temperature Grid Latitudes
call check(nf90_def_var(ncFileID, name="TmpJ", &
              xtype=nf90_double, dimids=TmpJDimID, varid=TmpJVarID),   "TmpJ def_var" )
call check(nf90_put_att(ncFileID, TmpJVarID, "long_name", "latitude"), "TmpJ long_name")
call check(nf90_put_att(ncFileID, TmpJVarID, "cartesian_axis", "Y"),   "TmpJ cartesian_axis")
call check(nf90_put_att(ncFileID, TmpJVarID, "units", "degrees_north"),"TmpJ units")
call check(nf90_put_att(ncFileID, TmpJVarID, "valid_range", (/ -90.0_r8, 90.0_r8 /)), &
                                 "TmpJ valid_range")

! (Common) grid levels
call check(nf90_def_var(ncFileID, name="lev", &
              xtype=nf90_int, dimids=levDimID, varid=levVarID) ,    "lev def_var")
call check(nf90_put_att(ncFileID, levVarID, "long_name", "level"),  "lev long_name")
call check(nf90_put_att(ncFileID, levVarID, "cartesian_axis", "Z"), "lev cartesian_axis")
call check(nf90_put_att(ncFileID, levVarID, "units", "hPa"),        "lev units")
call check(nf90_put_att(ncFileID, levVarID, "positive", "down"),    "lev positive")

! Velocity Grid Longitudes
call check(nf90_def_var(ncFileID, name="VelI", &
              xtype=nf90_double, dimids=VelIDimID, varid=VelIVarID) ,  "VelI def_var")
call check(nf90_put_att(ncFileID, VelIVarID, "long_name", "longitude"),"VelI long_name")
call check(nf90_put_att(ncFileID, VelIVarID, "cartesian_axis", "X"),   "VelI cartesian_axis")
call check(nf90_put_att(ncFileID, VelIVarID, "units", "degrees_east"), "VelI units")
call check(nf90_put_att(ncFileID, VelIVarID, "valid_range", (/ 0.0_r8, 360.1_r8 /)), &
                                 "VelI valid_range")

! Velocity Grid Latitudes
call check(nf90_def_var(ncFileID, name="VelJ", &
              xtype=nf90_double, dimids=VelJDimID, varid=VelJVarID) ,   "VelJ def_var")
call check(nf90_put_att(ncFileID, VelJVarID, "long_name", "latitude"),  "VelJ long_name")
call check(nf90_put_att(ncFileID, VelJVarID, "cartesian_axis", "Y"),    "VelJ cartesian_axis")
call check(nf90_put_att(ncFileID, VelJVarID, "units", "degrees_north"), "VelJ units")
call check(nf90_put_att(ncFileID, VelJVarID, "valid_range", (/ -90.0_r8, 90.0_r8 /)), &
                                 "VelJ valid_range")

! Number of Tracers
if ( ntracer > 0 ) then
   call check(nf90_def_var(ncFileID, name="tracers", &
                 xtype=nf90_int, dimids=tracerDimID, varid=tracerVarID) , "tracers def_var")
   call check(nf90_put_att(ncFileID, tracerVarID, "long_name", "tracer identifier"), &
                                    "tracers long_name")
endif

if ( output_state_vector ) then

   !----------------------------------------------------------------------------
   ! Create attributes for the state vector 
   !----------------------------------------------------------------------------

  ! Define the state vector coordinate variable
   call check(nf90_def_var(ncid=ncFileID,name="StateVariable", xtype=nf90_int, &
              dimids=StateVarDimID, varid=StateVarVarID), "statevariable def_var")
   call check(nf90_put_att(ncFileID, StateVarVarID, "long_name", "State Variable ID"), &
                                    "statevariable long_name")
   call check(nf90_put_att(ncFileID, StateVarVarID, "units",     "indexical"), &
                                    "statevariable units")
   call check(nf90_put_att(ncFileID, StateVarVarID, "valid_range", (/ 1, model_size /)), &
                                    "statevariable valid_range")

   ! Define the actual state vector
   call check(nf90_def_var(ncid=ncFileID, name="state", xtype=nf90_real, &
              dimids = (/ StateVarDimID, MemberDimID, unlimitedDimID /), &
              varid=StateVarID), "state def_var")
   call check(nf90_put_att(ncFileID, StateVarID, "long_name", "model state or fcopy"), &
                                           "state long_name")
   call check(nf90_put_att(ncFileID, StateVarId, "vector_to_prog_var","FMS-Bgrid"), &
                                           "state vector_to_prog_var")
   call check(nf90_put_att(ncFileID, StateVarId, "temperature_units","degrees Kelvin"), &
                                           "state temperature_units")
   call check(nf90_put_att(ncFileID, StateVarId, "pressure_units","Pa"), &
                                           "state pressure_units")
   call check(nf90_put_att(ncFileID, StateVarId, "U_units","m/s"),"state U_units")
   call check(nf90_put_att(ncFileID, StateVarId, "V_units","m/s"),"state V_units")

   ! Leave define mode so we can fill
   call check(nf90_enddef(ncfileID),"state enddef")

   ! Fill the state variable coordinate variable
   call check(nf90_put_var(ncFileID, StateVarVarID, (/ (i,i=1,model_size) /) ), &
                                    "state put_var")

else

   !----------------------------------------------------------------------------
   ! We need to process the prognostic variables.
   !----------------------------------------------------------------------------
   ! TemperatureGrid : surface pressure  vars%ps(tis:tie, tjs:tje) 
   !                 : temperature       vars%t (tis:tie, tjs:tje, klb:kub)
   !                 : tracers           vars%r (tis:tie, tjs:tje, klb:kub, 1:vars%ntrace)
   ! VelocityGrid    : u                 vars%u (vis:vie, vjs:vje, klb:kub) 
   !                 : v                 vars%v (vis:vie, vjs:tje, klb:kub)
   !----------------------------------------------------------------------------
   ! Create the (empty) Variables and the Attributes
   !----------------------------------------------------------------------------
 
   call check(nf90_def_var(ncid=ncFileID, name="ps", xtype=nf90_real, &
         dimids = (/ TmpIDimID, TmpJDimID, MemberDimID, unlimitedDimID /), &
         varid  = psVarID), "ps def_var")
   call check(nf90_put_att(ncFileID, psVarID, "long_name", "surface pressure"), &
                                           "ps long_name")
   call check(nf90_put_att(ncFileID, psVarID, "units", "Pa"), &
                                           "ps units")
   call check(nf90_put_att(ncFileID, psVarID, "units_long_name", "pascals"), &
                                           "ps units_long_name")


   call check(nf90_def_var(ncid=ncFileID, name="t", xtype=nf90_real, &
         dimids = (/ TmpIDimID, TmpJDimID, levDimID, MemberDimID, unlimitedDimID /), &
         varid  = tVarID), "t def_var")
   call check(nf90_put_att(ncFileID, tVarID, "long_name", "temperature"), "t long_name")
   call check(nf90_put_att(ncFileID, tVarID, "units", "degrees Kelvin"), "t units")


   call check(nf90_def_var(ncid=ncFileID, name="u", xtype=nf90_real, &
         dimids = (/ VelIDimID, VelJDimID, levDimID, MemberDimID, unlimitedDimID /), &
         varid  = uVarID), "u def_var")
   call check(nf90_put_att(ncFileID, uVarID, "long_name", "zonal wind component"), &
                                           "u long_name")
   call check(nf90_put_att(ncFileID, uVarID, "units", "m/s"), "u units")


   call check(nf90_def_var(ncid=ncFileID, name="v", xtype=nf90_real, &
         dimids = (/ VelIDimID, VelJDimID, levDimID, MemberDimID, unlimitedDimID /), &
         varid  = vVarID), "v def_var")
   call check(nf90_put_att(ncFileID, vVarID, "long_name", "meridional wind component"), &
                                           "v long_name")
   call check(nf90_put_att(ncFileID, vVarID, "units", "m/s"), "v units")

   if ( ntracer > 0 ) then
      call check(nf90_def_var(ncid=ncFileID, name="r", xtype=nf90_real, &
      dimids = (/TmpIDimID, TmpJDimID, levDimID, tracerDimID, MemberDimID, unlimitedDimID/),&
         varid  = rVarID), "r def_var")
      call check(nf90_put_att(ncFileID, rVarID, "long_name", "various tracers"), "r long_name")
   endif

   call check(nf90_enddef(ncfileID), "prognostic enddef")

endif

!-------------------------------------------------------------------------------
! Fill the variables
!-------------------------------------------------------------------------------

call check(nf90_put_var(ncFileID,   TmpIVarID, t_lons(tis:tie) ), "TmpI put_var")
call check(nf90_put_var(ncFileID,   TmpJVarID, t_lats(tjs:tje) ), "TmpJ put_var")
call check(nf90_put_var(ncFileID,   VelIVarID, v_lons(vis:vie) ), "VelI put_var")
call check(nf90_put_var(ncFileID,   VelJVarID, v_lats(vjs:vje) ), "VelJ put_var")
call check(nf90_put_var(ncFileID,    levVarID, (/ (i,i=1,   nlev) /) ), "lev put_var")
if ( ntracer > 0 ) then
   call check(nf90_put_var(ncFileID, tracerVarID, (/ (i,i=1,ntracer) /) ), "tracer put_var")
endif

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------
call check(nf90_sync(ncFileID),"atts sync")

write (*,*)'nc_write_model_atts: netCDF file ',ncFileID,' is synched ...'

contains

  ! Internal subroutine - checks error status after each netcdf, prints 
  !                       text message each time an error code is returned. 
  subroutine check(istatus, string1)
    integer, intent ( in) :: istatus
    character(len=*), intent(in), optional :: string1

    character(len=20)  :: myname = 'nc_write_model_atts '
    character(len=129) :: mystring
    integer            :: indexN

    if( istatus /= nf90_noerr) then

       if (present(string1) ) then
          if ((len_trim(string1)+len(myname)) <= len(mystring) ) then
             mystring = myname // trim(adjustl(string1))
          else
             indexN = len(mystring) - len(myname)
             mystring = myname // string1(1:indexN)
          endif 
       else
          mystring = myname
       endif

       call error_mesg(mystring, trim(nf90_strerror(istatus)), FATAL)

    endif

  end subroutine check

end function nc_write_model_atts



function nc_write_model_vars( ncFileID, statevec, copyindex, timeindex ) result (ierr)
!-----------------------------------------------------------------------------------------
! Writes the model-specific variables to a netCDF file
! TJH 25 June 2003
!
! TJH 29 July 2003 -- for the moment, all errors are fatal, so the
! return code is always '0 == normal', since the fatal errors stop execution.
!                                                                                
! There are two different (staggered) 3D grids being used simultaneously here. 
! The routines "prog_var_to_vector" and "vector_to_prog_var", 
! packs the prognostic variables into
! the requisite array for the data assimilation routines. That routine
! is the basis for the information stored in the netCDF files.
!
! TemperatureGrid : surface pressure  vars%ps(tis:tie, tjs:tje) 
!                 : temperature       vars%t (tis:tie, tjs:tje, klb:kup)
!                 : tracers           vars%r (tis:tie, tjs:tje, klb:kub, 1:vars%ntrace)
! VelocityGrid    : u                 vars%u (vis:vie, vjs:vje, klb:kub) 
!                 : v                 vars%v (vis:vie, vjs:tje, klb:kup)
!
! So there are six different dimensions and five different variables as long as
! simply lump "tracers" into one. 

use typeSizes
use netcdf

integer,                intent(in) :: ncFileID      ! netCDF file identifier
real(r8), dimension(:), intent(in) :: statevec
integer,                intent(in) :: copyindex
integer,                intent(in) :: timeindex
integer                            :: ierr          ! return value of function

!-------------------------------------------------------------------------------
real(r8), dimension(SIZE(statevec)) :: x
! removed local instance of Var (since it contains pointers which are
! then left hanging if Var goes out of scope), and replaced all references
! to Var with global_Var below.    nsc 31mar06
!type(prog_var_type) :: Var

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID
integer :: StateVarID, psVarID, tVarID, rVarID, uVarID, vVarID
integer :: tis, tie, tjs, tje       ! temperature grid start/stop
integer :: vis, vie, vjs, vje       ! velocity    grid start/stop
integer :: kub, klb
integer :: nTmpI, nTmpJ, nVelI, nVelJ, nlev, ntracer

if ( .not. module_initialized ) call static_init_model

ierr = 0     ! assume normal termination

!-------------------------------------------------------------------------------
! Get the bounds for storage on Temp and Velocity grids
! ?JEFF why can't I use the components of the prog_var_type?  
! More precisely, why doesn't prog_var_type drag around the necessary
! indices instead of just the extents?
!-------------------------------------------------------------------------------

tis = Dynam%Hgrid%Tmp%is; tie = Dynam%Hgrid%Tmp%ie
tjs = Dynam%Hgrid%Tmp%js; tje = Dynam%Hgrid%Tmp%je
vis = Dynam%Hgrid%Vel%is; vie = Dynam%Hgrid%Vel%ie
vjs = Dynam%Hgrid%Vel%js; vje = Dynam%Hgrid%Vel%je
kub = Var_dt%kub
klb = Var_dt%klb

nTmpI   = tie - tis + 1
nTmpJ   = tje - tjs + 1
nlev    = Var_dt%kub - Var_dt%klb + 1
ntracer = Var_dt%ntrace 
nVelI   = vie - vis + 1
nVelJ   = vje - vjs + 1

!-------------------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file, 
! then get all the Variable ID's we need.
!-------------------------------------------------------------------------------

call check(nf90_Inquire(ncFileID, nDimensions, nVariables, nAttributes, unlimitedDimID), &
        "inquire")

if ( output_state_vector ) then

   call check(NF90_inq_varid(ncFileID, "state", StateVarID), "state inq_varid" )
   call check(NF90_put_var(ncFileID, StateVarID, statevec,  &
                start=(/ 1, copyindex, timeindex /)), "state put_var")                               

else
   
   !----------------------------------------------------------------------------
   ! Fill the variables
   ! TemperatureGrid : surface pressure  Var%ps(tis:tie, tjs:tje) 
   !                 : temperature       Var%t (tis:tie, tjs:tje, klb:kub)
   !                 : tracers           Var%r (tis:tie, tjs:tje, klb:kub, 1:vars%ntrace)
   ! VelocityGrid    : u                 Var%u (vis:vie, vjs:vje, klb:kub) 
   !                 : v                 Var%v (vis:vie, vjs:tje, klb:kub)
   !----------------------------------------------------------------------------

   x = statevec ! Unfortunately, have to explicity cast it ...
                ! the filter uses a type=double,
                ! the vector_to_prog_var function expects a single.
   call vector_to_prog_var(x, get_model_size(), global_Var)
   
   
   call check(NF90_inq_varid(ncFileID, "ps", psVarID), "ps inq_varid")
   call check(nf90_put_var( ncFileID, psVarID, global_Var%ps(tis:tie, tjs:tje), &
                            start=(/ 1, 1, copyindex, timeindex /) ), "ps put_var")

   call check(NF90_inq_varid(ncFileID,  "t",  tVarID), "t inq_varid")
   call check(nf90_put_var( ncFileID,  tVarID, global_Var%t( tis:tie, tjs:tje, klb:kub ), &
                            start=(/ 1, 1, 1, copyindex, timeindex /) ), "t put_var")

   call check(NF90_inq_varid(ncFileID,  "u",  uVarID), "u inq_varid")
   call check(nf90_put_var( ncFileID,  uVarId, global_Var%u( vis:vie, vjs:vje, klb:kub ), &
                            start=(/ 1, 1, 1, copyindex, timeindex /) ), "u put_var")

   call check(NF90_inq_varid(ncFileID,  "v",  vVarID), "v inq_varid")
   call check(nf90_put_var( ncFileID,  vVarId, global_Var%v( vis:vie, vjs:vje, klb:kub ), &
                            start=(/ 1, 1, 1, copyindex, timeindex /) ), "v put_var")

   if ( ntracer > 0 ) then
      call check(NF90_inq_varid(ncFileID,  "r",  rVarID), "r inq_varid")
      call check(nf90_put_var( ncFileID,  rVarID, &
                    global_Var%r( tis:tie, tjs:tje, klb:kub, 1:ntracer ), & 
                   start=(/  1,   1,   1,   1, copyindex, timeindex /) ), "r put_var")
   endif
endif

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------

! write (*,*)'Finished filling variables ...'
call check(nf90_sync(ncFileID), "sync")
! write (*,*)'netCDF file is synched ...'

contains

  ! Internal subroutine - checks error status after each netcdf, prints 
  !                       text message each time an error code is returned. 
  subroutine check(istatus, string1)
    integer, intent ( in) :: istatus
    character(len=*), intent(in), optional :: string1

    character(len=20)  :: myname = 'nc_write_model_vars '
    character(len=129) :: mystring
    integer            :: indexN

    if( istatus /= nf90_noerr) then

       if (present(string1) ) then
          if ((len_trim(string1)+len(myname)) <= len(mystring) ) then
             mystring = myname // trim(adjustl(string1))
          else
             indexN = len(mystring) - len(myname)
             mystring = myname // string1(1:indexN)
          endif 
       else
          mystring = myname
       endif

       call error_mesg(mystring, trim(nf90_strerror(istatus)), FATAL)

    endif

  end subroutine check

end function nc_write_model_vars




subroutine pert_model_state(state, pert_state, interf_provided)


! Perturbs a model state for generating initial ensembles
! Returning interf_provided means go ahead and do this with uniform
! small independent perturbations.

real(r8), intent(in)  :: state(:)
real(r8), intent(out) :: pert_state(:)
logical,  intent(out) :: interf_provided

if ( .not. module_initialized ) call static_init_model

interf_provided = .false.

! you *cannot* set this to junk. in at least one
! location the caller is passing the same
! array into the first arg as the second.  doing this
! corrupts the state vector completely. set it to the
! incoming state if you must do something.
pert_state = state ! Just to satisfy INTENT(OUT)

end subroutine pert_model_state





subroutine ens_mean_for_model(ens_mean)
!------------------------------------------------------------------
! Not used in this model.

real(r8), intent(in) :: ens_mean(:)

if ( .not. module_initialized ) call static_init_model

end subroutine ens_mean_for_model



!#######################################################################
end module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
