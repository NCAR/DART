!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                       !!
!!                   GNU General Public License                          !!
!!                                                                       !!
!! This file is part of the Flexible Modeling System (FMS).              !!
!!                                                                       !!
!! FMS is free software; you can redistribute it and/or modify           !!
!! it and are expected to follow the terms of the GNU General Public     !!
!! License as published by the Free Software Foundation.                 !!
!!                                                                       !!
!! FMS is distributed in the hope that it will be useful,                !!
!! but WITHOUT ANY WARRANTY; without even the implied warranty of        !!
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         !!
!! GNU General Public License for more details.                          !!
!!                                                                       !!
!! You should have received a copy of the GNU General Public License     !!
!! along with FMS; if not, write to:                                     !!
!!          Free Software Foundation, Inc.                               !!
!!          59 Temple Place, Suite 330                                   !!
!!          Boston, MA  02111-1307  USA                                  !!
!! or see:                                                               !!
!!          http://www.gnu.org/licenses/gpl.txt                          !!
!!                                                                       !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

use diag_manager_mod, only: diag_manager_init, diag_manager_end, get_base_date


use   bgrid_prog_var_mod, only: prog_var_type, var_init, prog_var_init, &
                                open_prog_var_file, read_prog_var

use      bgrid_horiz_mod, only: get_horiz_grid_size,        &
                                get_horiz_grid_bound, TGRID, VGRID

use time_manager_mod, only: time_type, set_time, get_time,  &
                            operator(+), operator (<), operator (>), &
                            operator (/=), operator (/), operator (*)

use              fms_mod, only: file_exist, open_namelist_file, &
                                error_mesg, FATAL,              &
                                check_nml_error, stdlog,        &
                                write_version_number,           &
                                mpp_pe, mpp_root_pe,            &
                                close_file, set_domain,         &
                                fms_init, mpp_clock_init,       &
                                MPP_CLOCK_SYNC,                 &
                                open_restart_file, mpp_clock_end, &
                                mpp_clock_begin

use       mpp_io_mod, only: mpp_open, mpp_close, MPP_ASCII, MPP_OVERWR, &
                            MPP_SEQUENTIAL, MPP_SINGLE, MPP_DELETE


! routines used by subroutine bgrid_physics
use bgrid_change_grid_mod, only: vel_to_mass, mass_to_vel
use bgrid_horiz_mod      , only: horiz_grid_type
use bgrid_vert_mod       , only: vert_grid_type, &
                                 compute_pres_full, compute_pres_half
use bgrid_halo_mod       , only: update_halo, UWND, VWND, TEMP, &
                                 NORTH, EAST, WEST, SOUTH
use hs_forcing_mod       , only: hs_forcing_init, hs_forcing

use location_mod         , only: location_type, get_location, set_location, get_dist, &
                                 nc_write_location

use types_mod

!-----------------------------------------------------------------------

implicit none
private

public  atmosphere,      &
        atmosphere_init, &
        atmosphere_end,  &
        atmosphere_resolution,  &
        atmosphere_boundary,    &
        get_atmosphere_axes, &
        prog_var_type, &
        get_model_size, &
        prog_var_to_vector, &
        vector_to_prog_var, &
        adv_1step, &
        get_state_meta_data, &
        model_interpolate, &
        get_val, &              ! THIS SHOULD ONLY BE FOR TEMP TEST
        model_interpolate, &
        get_model_time_step, &
        end_model, &
        static_init_model,  &
        init_model_instance, &
        init_time, &
        init_conditions

!-----------------------------------------------------------------------

character(len=128) :: version = '$Id$'
character(len=128) :: tag = '$Name$'

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
      logical :: override = .false.  ! override restart values for date
      integer :: days=0, hours=0, minutes=0, seconds=0
      integer :: dt_atmos = 0

      namelist /main_nml/ current_time, override, dt_atmos, &
                          days, hours, minutes, seconds

!-----------------------------------------------------------------------
! More stuff from atmos_solo driver
! ----- model time -----

   type (time_type) :: Time, Time_init, Time_end, Time_step_atmos

! ----- coupled model initial date -----

   integer :: date_init(6)



!-----------------------------------------------------------------------
!---- private data ----

type (bgrid_dynam_type) :: Dynam
type    (prog_var_type) :: Var_dt

integer                            :: model_size
!real                               :: dt_atmos
real,    dimension(:,:,:), pointer :: omega
integer, dimension(4)              :: atmos_axes
integer                            :: num_levels
integer                            :: ntracers 
integer                            :: nprognostic

real, dimension(:), pointer        :: v_lons, v_lats, t_lons, t_lats

!------------------------------------------------------------
! ----- timing flags -----

   integer :: id_init, id_loop, id_end
   integer, parameter :: timing_level = 1



!-----------------------------------------------------------------------

contains

!#######################################################################

 subroutine atmosphere (Var, Time)

implicit none

 type (time_type), intent(in) :: Time
 type (prog_var_type), intent(inout) :: Var

  type(time_type) :: Time_prev, Time_next
!-----------------------------------------------------------------------

   Time_prev = Time                       ! two time-level scheme
   Time_next = Time + get_model_time_step()

!---- dynamics -----

   call bgrid_core_driver ( Time_next, Var, Var_dt, Dynam, omega )

!---- call physics -----

   call bgrid_physics ( physics_window, real(dt_atmos), Time_next,  &
                        Dynam%Hgrid,   Dynam%Vgrid,    Dynam,             &
                        Var,     Var_dt                       )

!---- time differencing and diagnostics -----

   call bgrid_core_time_diff ( omega, Time_next, Dynam, Var, Var_dt )

!-----------------------------------------------------------------------

 end subroutine atmosphere

!#######################################################################

subroutine adv_1step(x, Time)

! Does single time-step advance for B-grid model with vector state as
! input and output. This is a modified version of subroutine atmosphere
! in original bgrid_solo_atmosphere driver.

implicit none

real, intent(inout) :: x(:)
! Time is needed for more general models like this; need to add in to 
! low-order models
type(time_type), intent(in) :: Time
!!!type(prog_var_type), intent(inout) :: Var

type(prog_var_type) :: Var
type(time_type) :: Time_next

! Convert the vector to a B-grid native state representation
call vector_to_prog_var(x, get_model_size(), Var)

! Compute the end of the time interval
Time_next = Time + get_model_time_step()

! Do dynamics
! Dynam, Var_dt and omega currently in static storage, is that where they should stay?
call bgrid_core_driver(Time_next, Var, Var_dt, Dynam, omega)

! Call physics; physics_window is also in global static storage
! dt_atmos is seconds for time step from namelist for now
call bgrid_physics(physics_window, real(dt_atmos), Time_next, &
   Dynam%Hgrid, Dynam%Vgrid, Dynam, Var, Var_dt)

! Time differencing and diagnostics
call bgrid_core_time_diff(omega, Time_next, Dynam, Var, Var_dt)

! Convert back to vector form
call prog_var_to_vector(Var, x, get_model_size())

end subroutine adv_1step

!#######################################################################

subroutine static_init_model()

! INitializes class data for a B-grid model (all the stuff that needs to
! be done once.

implicit none

call fms_init()
call atmos_model_init()

end subroutine static_init_model

!#######################################################################

subroutine init_model_instance(Var)

! Initializes an instance (Var) of a B-grid model state variable

implicit none

type(prog_var_type), intent(out) :: Var

call prog_var_init(Dynam%Hgrid, num_levels, ntracers, nprognostic, Var)

end subroutine init_model_instance

!#######################################################################

subroutine init_conditions(x)

implicit none

! Reads in restart initial conditions from B-grid and converts to vector

real, intent(in) :: x(:)

type(prog_var_type) :: Var
real, dimension(Dynam%Hgrid%ilb:Dynam%Hgrid%iub, Dynam%Hgrid%jlb:Dynam%Hgrid%jub) :: fis, res
real, allocatable, dimension(:) :: eta, peta
integer :: ix, jx, kx, nt, ntp

! Need to initialize var???
call init_model_instance(Var)

! FOR NOW, TRY TO READ IN CURRENT ICS via read_prog_var
call open_prog_var_file(ix, jx, kx, nt, ntp)

allocate(eta(kx + 1), peta(kx + 1))

!!! WARNING; MAY BE DANGEROUS TO USE Dynam%Hgrid here, it could get changed
!!! inappropriately???
call read_prog_var(Dynam%Hgrid, Var, eta, peta, fis, res)

deallocate(eta, peta)

call prog_var_to_vector(Var, x, get_model_size())

! Probably need to release allocated storage from Var, too???

end subroutine init_conditions

!#######################################################################

! THIS SUBROUTINE WAS ORIGINALLY IN ATMOS_SOLO DRIVER

   subroutine atmos_model_init()

!-----------------------------------------------------------------------
    integer :: total_days, total_seconds, unit, ierr, io, id, jd, kd
    integer :: date(6)
    type (time_type) :: Run_length
    logical :: use_namelist

    integer :: num_atmos_calls
!-----------------------------------------------------------------------
!----- initialization timing identifiers ----

 id_init = mpp_clock_init ('MAIN: initialization', timing_level, flags=MPP_CLOCK_SYNC)
 id_loop = mpp_clock_init ('MAIN: time loop'     , timing_level, flags=MPP_CLOCK_SYNC)
 id_end  = mpp_clock_init ('MAIN: termination'   , timing_level, flags=MPP_CLOCK_SYNC)

 call mpp_clock_begin (id_init)

!----- read namelist -------

   unit = open_namelist_file ( )
   ierr=1; do while (ierr /= 0)
          read  (unit, nml=main_nml, iostat=io, end=10)
          ierr = check_nml_error (io, 'main_nml')
   enddo
10 call mpp_close (unit)

!----- write namelist to logfile -----

   call write_version_number (version,tag)
   if ( mpp_pe() == mpp_root_pe() ) write (stdlog(), nml=main_nml)

   if(dt_atmos == 0) then
     call error_mesg ('program atmos_model', 'dt_atmos has not been specified', FATAL)
   endif

!----- read restart file -----

   if (file_exist('INPUT/atmos_model.res')) then
       unit = open_restart_file ('INPUT/atmos_model.res', 'read')
       read  (unit) date
       call mpp_close (unit)
       use_namelist = .false.
   else
       use_namelist = .true.
   endif

!----- override date with namelist values ------
!----- (either no restart or override flag on) ---

 if ( use_namelist .or. override ) then
      date(1:2) = 0
      date(3:6) = current_time
 endif

!----- write current/initial date actually used to logfile file -----

    if ( mpp_pe() == mpp_root_pe() ) then
      write (stdlog(),16) date(3:6)
    endif

 16 format ('  current time used = day',i5,' hour',i3,2(':',i2.2))

!-----------------------------------------------------------------------
!------ initialize diagnostics manager ------

    call diag_manager_init

!----- always override initial/base date with diag_manager value -----

    call get_base_date ( date_init(1), date_init(2), date_init(3), &
                         date_init(4), date_init(5), date_init(6)  )

    if ( date_init(1)+date_init(2) /= 0 ) then
         call error_mesg ('program atmos_model', 'invalid base base - &
                          &must have year = month = 0', FATAL)
    endif

!----- set initial and current time types ------
!----- set run length and compute ending time -----

    Time_init  = set_time_increment (date_init(3), date_init(4), date_init(5), date_init(6))
    Time       = set_time_increment (date     (3), date     (4), date     (5), date     (6))
    Run_length = set_time_increment (days        , hours       , minutes     , seconds     )
    Time_end   = Time + Run_length

!-----------------------------------------------------------------------
!----- write time stamps (for start time and end time) ------

      call mpp_open (unit, 'time_stamp.out', form=MPP_ASCII, action=MPP_OVERWR, &
                     access=MPP_SEQUENTIAL, threading=MPP_SINGLE, nohdrs=.true. )

      if ( mpp_pe() == mpp_root_pe() ) write (unit,20) date

!     compute ending time in days,hours,minutes,seconds
      call get_time_increment (Time_end, date(3), date(4), date(5), date(6))

      if ( mpp_pe() == mpp_root_pe() ) write (unit,20) date

      call mpp_close (unit)

  20  format (6i7,2x,'day')   ! can handle day <= 999999

!-----------------------------------------------------------------------
!----- compute the time steps ------
!----- determine maximum number of iterations per loop ------

      Time_step_atmos = set_time (dt_atmos,0)

!-----------------------------------------------------------------------
!----- initial (base) time must not be greater than current time -----

   if ( Time_init > Time ) call error_mesg ('program atmos_model',  &
                   'initial time is greater than current time', FATAL)

!-----------------------------------------------------------------------
!------ initialize atmospheric model ------

      call atmosphere_init (Time_init, Time, Time_step_atmos)

!-----------------------------------------------------------------------
!---- open and close output restart to make sure directory is there ----

      unit = open_restart_file ('RESTART/atmos_model.res', 'write')
      call mpp_close (unit, action=MPP_DELETE)


!  ---- terminate timing ----
   call mpp_clock_end (id_init)

!-----------------------------------------------------------------------

   end subroutine atmos_model_init



!#######################################################################

 subroutine atmosphere_init (Time_init, Time, Time_step )

 type (time_type),     intent(in)    :: Time_init, Time, Time_step


!!! WARNING: This PROG_VAR_TYPE MAY HOG STORAGE FOREVER, NEED TO DEALLOCATE
type(prog_var_type) :: Var

  integer :: unit, sec, ierr, io
  integer :: tnlon, tnlat, vnlon, vnlat
  integer :: t_horiz_size, v_horiz_size

real, allocatable :: t_lat_bnds(:), t_lon_bnds(:), v_lat_bnds(:), v_lon_bnds(:)

integer :: i


!-----------------------------------------------------------------------
!----- read namelist -----

    if (file_exist('input.nml')) then
        unit = open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
           read (unit, nml=atmosphere_nml, iostat=io, end=10)
           ierr = check_nml_error (io, 'atmosphere_nml')
        enddo
 10     call close_file (unit)
    endif

!----- write version and namelist to log file -----

    call write_version_number ( version, tag )
    if ( mpp_pe() == mpp_root_pe() ) write (stdlog(), nml=atmosphere_nml)

!---- compute physics/atmos time step in seconds ----

   Time_step_atmos = Time_step
   call get_time (Time_step_atmos, sec)
!   dt_atmos = real(sec) LEFT OVER FROM COMBINE OF ATMOS_SOLO DRIVER

!----- initialize dynamical core -----

   call bgrid_core_driver_init ( Time_init, Time, Time_step,    &
                                 Var, Var_dt, Dynam, atmos_axes )

!----- initialize storage needed for vert motion ----

    omega => var_init (Dynam%Hgrid, Dynam%Vgrid%nlev)

!----- initialize physics interface -----

    call hs_forcing_init ( atmos_axes, Time )

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
model_size = model_size + Var%ntrace * num_levels * t_horiz_size
write(*, *) 'model_size ', model_size

! Also static store the number of levels, ntracers, and prognostic tracers
ntracers = Var%ntrace
nprognostic = Var%ntprog

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

!#######################################################################

 subroutine atmosphere_end(Var)

 type(prog_var_type), intent(in) :: Var

 integer :: unit

    call bgrid_core_driver_end ( Var, Dynam )

 end subroutine atmosphere_end

!#######################################################################
!    returns the number of longitude and latitude grid points
!    for either the local PEs grid (default) or the global grid

 subroutine atmosphere_resolution (nlon, nlat, global)

  integer, intent(out)          :: nlon, nlat
  logical, intent(in), optional :: global

!---- return the size of the grid used for physics computations ----

    call get_horiz_grid_size (Dynam % Hgrid, TGRID, nlon, nlat, global)

 end subroutine atmosphere_resolution

!#######################################################################
!    returns the longitude and latitude grid box edges
!    for either the local PEs grid (default) or the global grid

 subroutine atmosphere_boundary (blon, blat, global)

    real,    intent(out)          :: blon(:), blat(:)
    logical, intent(in), optional :: global

!----- return the longitudinal and latitudinal grid box edges ----------

    call get_horiz_grid_bound (Dynam % Hgrid, TGRID, blon, blat, global)

 end subroutine atmosphere_boundary

!#######################################################################
!    returns the axis indices associated with the coupling grid

 subroutine get_atmosphere_axes ( axes )

   integer, intent(out) :: axes (:)

!----- returns the axis indices for the atmospheric (mass) grid -----

     if ( size(axes) < 0 .or. size(axes) > 4 ) call error_mesg (    &
                    'get_atmosphere_axes in atmosphere_mod', &
                           'size of argument is incorrect', FATAL   )

     axes (1:size(axes)) = atmos_axes (1:size(axes))

 end subroutine get_atmosphere_axes

!#######################################################################

subroutine bgrid_physics ( window, dt_phys, Time, Hgrid, Vgrid, &
                           Dynam, Var, Var_dt )

!-----------------------------------------------------------------------
!
!   Time      =  current time (time_type, see time manager)
!
!-----------------------------------------------------------------------
  integer, intent(in)                :: window(2)
  real,    intent(in)                :: dt_phys
       type(time_type),intent(in)    :: Time
type (horiz_grid_type),intent(inout) :: Hgrid
type  (vert_grid_type),intent(in)    :: Vgrid
type(bgrid_dynam_type),intent(in)    :: Dynam
type   (prog_var_type),intent(in)    :: Var
type   (prog_var_type),intent(inout) :: Var_dt

!-----------------------------------------------------------------------
  integer :: j, k, n, is, ie, js, je, i1, i2, j1, j2, nt
  integer :: ix, jx, idim, jdim
!-----------------------------------------------------------------------

   real, dimension(window(1),window(2),Vgrid%nlev) :: p_full, u_dt, v_dt

   real, dimension(window(1),window(2),Vgrid%nlev+1) :: p_half

   real, dimension(Hgrid%ilb:Hgrid%iub, &
                   Hgrid%jlb:Hgrid%jub, &
                   Vgrid%nlev) :: uh, vh, uh_dt, vh_dt

   real, dimension(window(1),window(2)) :: pssl_new
!-----------------------------------------------------------------------
!---------------------------- do physics -------------------------------

    idim = window(1)
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

       ie = min ( is+idim-1, Hgrid%Tmp%ie )
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

      call hs_forcing ( i1, i2, j1, j2, dt_phys, Time ,&
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

      call hs_forcing ( i1, i2, j1, j2, dt_phys, Time ,&
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

        is = is + idim

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

     uh(:,Hgrid%jub,:) = 0.0
     vh(:,Hgrid%jub,:) = 0.0

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

implicit none

integer :: get_model_size

get_model_size = model_size

end function get_model_size

!#######################################################################

subroutine prog_var_to_vector(vars, x, size)

implicit none

integer, intent(in) :: size
type(prog_var_type), intent(in) :: vars
real, intent(out) :: x(size)

integer :: i, j, k, nt, index
integer :: num_levs
integer :: tis, tie, tjs, tje, vis, vie, vjs, vje

! Get the bounds for storage on Temp and Velocity grids
tis = Dynam%Hgrid%Tmp%is; tie = Dynam%Hgrid%Tmp%ie
tjs = Dynam%Hgrid%Tmp%js; tje = Dynam%Hgrid%Tmp%je
vis = Dynam%Hgrid%Vel%is; vie = Dynam%Hgrid%Vel%ie
vjs = Dynam%Hgrid%Vel%js; vje = Dynam%Hgrid%Vel%je
num_levs = vars%kub - vars%klb + 1

! Start copying fields to straight vector
! Do everything on the T grid first
index = 0
do i = tis, tie
   do j = tjs, tje
! Surface pressure is first
      index = index + 1
      x(index) = vars%ps(i, j)
! Now do t and tracers at successive levels
      do k = vars%klb, vars%kub
         index = index + 1
         x(index) = vars%t(i, j, k)
         do nt = 1, vars%ntrace
            index = index + 1
            x(index) = vars%r(i, j, k, nt)
         end do
      end do
   end do
end do

! Now do the velocity grid, u and v
do i = vis, vie
   do j = vjs, vje
      do k = vars%klb, vars%kub
         index = index + 1
         x(index) = vars%u(i, j, k)
         index = index + 1
         x(index) = vars%v(i, j, k)
      end do
   end do
end do
   



! Temporary check
if(index /= size) then
   write(*, *) 'prog_var_to_vector bad index sum '
   write(*, *) 'index, size ', index, size
   stop
endif

end subroutine prog_var_to_vector

!#######################################################################

subroutine vector_to_prog_var(x, size, vars)

implicit none

integer, intent(in) :: size
real, intent(in) :: x(size)
type(prog_var_type), intent(inout) :: vars

integer :: i, j, k, nt, index
integer :: num_levs
integer :: tis, tie, tjs, tje, vis, vie, vjs, vje

! Initialize the static parts of the prog var type
call prog_var_init(Dynam%Hgrid, num_levels, ntracers, nprognostic, vars)

! Get the bounds for storage on Temp and Velocity grids
tis = Dynam%Hgrid%Tmp%is; tie = Dynam%Hgrid%Tmp%ie
tjs = Dynam%Hgrid%Tmp%js; tje = Dynam%Hgrid%Tmp%je
vis = Dynam%Hgrid%Vel%is; vie = Dynam%Hgrid%Vel%ie
vjs = Dynam%Hgrid%Vel%js; vje = Dynam%Hgrid%Vel%je
num_levs = vars%kub - vars%klb + 1

! Start copying fields from straight vector
! Everything on T grid first
index = 0
do i = tis, tie
   do j = tjs, tje
! Surface pressure is first
      index = index + 1
      vars%ps(i, j) = x(index)
! For non-eta models, pssl is same as ps??? Need to change?
      vars%pssl(i, j) = vars%ps(i, j)
! Now do t and tracers at successive levels
      do k = vars%klb, vars%kub
         index = index + 1
         vars%t(i, j, k) = x(index)
         do nt = 1, vars%ntrace
            index = index + 1
            vars%r(i, j, k, nt) = x(index)
         end do
      end do
   end do
end do

! Now do the velocity grid, u and v
do i = vis, vie
   do j = vjs, vje
      do k = vars%klb, vars%kub
         index = index + 1
         vars%u(i, j, k) = x(index)
         index = index + 1
         vars%v(i, j, k) = x(index)
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
if(index /= size) then
   write(*, *) 'prog_var_to_vector bad index sum '
   write(*, *) 'index, size ', index, size
   stop
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

! Time_step_atmos is global static storage
get_model_time_step =  Time_step_atmos

end function get_model_time_step

!#######################################################################

subroutine get_state_meta_data(index_in, location)
!---------------------------------------------------------------------
!
! Given an integer index into the state vector structure, returns the
! associated location. This is not a function because the more general
! form of the call has a second intent(out) optional argument kind.
! Maybe a functional form should be added?
!  SHOULD THIS ALSO RETURN THE TYPE OF THIS VARIABLE???

implicit none

integer, intent(in) :: index_in
! Temporary kluge of location type
!integer, intent(in) :: location
type(location_type), intent(out) :: location

integer :: i, j, k, nt, index
integer :: tis, tie, tjs, tje, vis, vie, vjs, vje
integer :: t_size, v_size, t_grid_size, v_grid_size, t_per_col, v_per_col
integer :: num_t_lons, num_t_lats, num_v_lons, num_v_lats
integer :: col_num, col_elem, v_index
integer :: lat_index, lon_index
real :: lon, lat, lev

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
index = index_in - 1

! Is this point in the t_grid
if(index < t_grid_size) then
   col_num = index / t_per_col
   col_elem = index - col_num * t_per_col
!   write(*, *) 't_grid col and element ', col_num, col_elem
   lon_index = col_num / num_t_lats
   lat_index = col_num - lon_index * num_t_lats
!   write(*, *) 'lon and lat index ', lon_index, lat_index
! Note, array runs from 1, index runs from 0
   lon = t_lons(lon_index + 1)
   lat = t_lats(lat_index + 1)
   if(col_elem == 0) then
      lev = -1
   else
      lev = int((col_elem  + 1)/ 2)
   endif

else
! It's in the v_grid
   v_index = index - t_grid_size
   col_num = v_index / v_per_col
   col_elem = v_index - col_num * v_per_col
!   write(*, *) 'v_grid col and element ', col_num, col_elem
   lon_index = col_num / num_v_lats
   lat_index = col_num - lon_index * num_v_lats
!   write(*, *) 'lon and lat index ', lon_index, lat_index
! Note, array runs from 1, index runs from 0
   lon = v_lons(lon_index + 1)
! Problems with round-off over 360.0
   if(abs(lon - 360.0) < 0.00001) lon = 360.0
   lat = v_lats(lat_index + 1)
   lev = int((col_elem + 2) / 2)
endif

!write(*, *) 'lon, lat, and lev ', lon, lat, lev
location = set_location(lon, lat, lev)


end subroutine get_state_meta_data

!#######################################################################

function model_interpolate(x, location, type)
!!!function model_interpolate(x, lon, lat, level, type)

implicit none

!!! EVENTUALLY NEED TO DO SOMETHING WITH TYPE; BUT FOR NOW THIS JUST FINDS
!!! THE HORIZONTAL PART OF THE INTPERPOLATION??? SHOULD ARGUMENT BE A HYBRID
!!! LOCATION TYPE VARIABLE HERE???

real :: model_interpolate
real, intent(in) :: x(:)
type(location_type), intent(in) :: location
integer, intent(in) :: type

integer :: num_lons, num_lats, lon_below, lon_above, lat_below, lat_above, i
real :: bot_lon, top_lon, delta_lon, bot_lat, top_lat, delta_lat
real :: lon_fract, lat_fract, val(2, 2), temp_lon, a(2)
real :: lon, lat, level, lon_lat_lev(3)

! Would it be better to pass state as prog_var_type (model state type) to here?

! First version only allows identity obs (level specifies verical???)
lon_lat_lev = get_location(location)
lon = lon_lat_lev(1); lat = lon_lat_lev(2); level = lon_lat_lev(3)

! Depending on type, get appropriate lon and lat grid specs
! Types temporarily defined as 1=u, 2=v, 3=ps, 4=t, n=tracer number n-4
if(type == 1 .or. type == 2) then
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
   if(lon < 180.0) temp_lon = lon + 360.0
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
   lat_fract = 0.0
else
! North of top lat NEED TO DO BETTER: NOT REALLY REGULAR
   lat_below = num_lats - 1
   lat_above = num_lats
   lat_fract = 1.0
endif

! Level is obvious for now

! Now, need to find the values for the four corners
val(1, 1) =  get_val(x, lon_below, lat_below, int(level), type)
val(1, 2) =  get_val(x, lon_below, lat_above, int(level), type)
val(2, 1) =  get_val(x, lon_above, lat_below, int(level), type)
val(2, 2) =  get_val(x, lon_above, lat_above, int(level), type)

! Do the weighted average for interpolation
!write(*, *) 'fracts ', lon_fract, lat_fract
do i = 1, 2
   a(i) = lon_fract * val(2, i) + (1.0 - lon_fract) * val(1, i)
end do

model_interpolate = lat_fract * a(2) + (1.0 - lat_fract) * a(1)


end function model_interpolate

!#######################################################################

function get_val(x, lon_index, lat_index, level, type)

implicit none

real :: get_val
real, intent(in) :: x(:)
integer, intent(in) :: lon_index, lat_index, level, type

integer :: tis, tie, num_t_lons, tjs, tje, num_t_lats, t_size
integer :: vis, vie, num_v_lons, vjs, vje, num_v_lats, v_size
integer :: t_per_col, t_grid_size, v_per_col, v_grid_size, index

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
if(type > 2) then
   index = t_per_col * (lat_index - 1 + (lon_index - 1) * num_t_lats)
   if(type == 3) then
      index = index + 1
   else 
      index = index + 1 + (level - 1) * (1 + ntracers) + (type - 3)
   endif

! Type 1 or 2 is u or v
else
! It's in the v_grid
   index = t_grid_size + v_per_col * (lat_index - 1 + (lon_index - 1) * num_v_lats)
   index = index + (level - 1) * 2 + type
endif

get_val = x(index)

end function get_val

!#######################################################################

subroutine end_model()

implicit none

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
   h = h/3600 ;   s = s - h*3600
   m = s/60   ;   s = s - m*60

 end subroutine get_time_increment

!#######################################################################

subroutine init_time(i_time)

implicit none

! For now returns value of Time_init which is set in initialization routines.

type(time_type), intent(out) :: i_time

i_time = Time

end subroutine init_time



!#######################################################################
end module model_mod
