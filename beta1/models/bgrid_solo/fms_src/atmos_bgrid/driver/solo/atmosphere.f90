!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                   !!
!!                   GNU General Public License                      !!
!!                                                                   !!
!! This file is part of the Flexible Modeling System (FMS).          !!
!!                                                                   !!
!! FMS is free software; you can redistribute it and/or modify       !!
!! it and are expected to follow the terms of the GNU General Public !!
!! License as published by the Free Software Foundation.             !!
!!                                                                   !!
!! FMS is distributed in the hope that it will be useful,            !!
!! but WITHOUT ANY WARRANTY; without even the implied warranty of    !!
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     !!
!! GNU General Public License for more details.                      !!
!!                                                                   !!
!! You should have received a copy of the GNU General Public License !!
!! along with FMS; if not, write to:                                 !!
!!          Free Software Foundation, Inc.                           !!
!!          59 Temple Place, Suite 330                               !!
!!          Boston, MA  02111-1307  USA                              !!
!! or see:                                                           !!
!!          http://www.gnu.org/licenses/gpl.txt                      !!
!!                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module atmosphere_mod

!-----------------------------------------------------------------------
!
!     interface for B-grid dynamics using the Held-Suarez forcing
!
!-----------------------------------------------------------------------
!---------------- m o d u l e   i n f o r m a t i o n ------------------
use types_mod, only : r8
use bgrid_core_driver_mod, only: bgrid_dynam_type,       &
                                 bgrid_core_driver_init, &
                                 bgrid_core_driver,      &
                                 bgrid_core_time_diff,   &
                                 bgrid_core_driver_end,  &
                                 get_bottom_data,        &
                                 put_bottom_data

use   bgrid_prog_var_mod, only: prog_var_type, var_init

use      bgrid_horiz_mod, only: get_horiz_grid_size,        &
                                get_horiz_grid_bound, TGRID

use     time_manager_mod, only: time_type, get_time, operator(+)

use utilities_mod, only : find_namelist_in_file, check_namelist_read
use              fms_mod, only: file_exist, open_namelist_file, &
                                error_mesg, FATAL,              &
                                stdlog,        &
                                write_version_number,           &
                                close_file

! routines used by subroutine bgrid_physics
use bgrid_change_grid_mod, only: vel_to_mass, mass_to_vel
use bgrid_horiz_mod      , only: horiz_grid_type
use bgrid_vert_mod       , only: vert_grid_type, &
                                 compute_pres_full, compute_pres_half
use bgrid_halo_mod       , only: update_halo, UWND, VWND, TEMP, &
                                 NORTH, EAST, WEST, SOUTH
use hs_forcing_mod       , only: hs_forcing_init, hs_forcing

!-----------------------------------------------------------------------

implicit none
private

public  atmosphere,      &
        atmosphere_init, &
        atmosphere_end,  &
        atmosphere_resolution,  &
        atmosphere_boundary,    &
        get_atmosphere_axes

!-----------------------------------------------------------------------

character(len=128) :: version = '$Revision$'
character(len=128) :: tag = '$Id$'

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
!---- private data ----

type (bgrid_dynam_type) :: Dynam
type    (prog_var_type) :: Var, Var_dt
type        (time_type) :: Time_step_atmos

real(r8)                               :: dt_atmos
real(r8),    dimension(:,:,:), pointer :: omega
integer, dimension(4)              :: atmos_axes

!-----------------------------------------------------------------------

contains

!#######################################################################

 subroutine atmosphere (Time)
 type (time_type), intent(in) :: Time

  type(time_type) :: Time_prev, Time_next
!-----------------------------------------------------------------------

   Time_prev = Time                       ! two time-level scheme
   Time_next = Time + Time_step_atmos

!---- dynamics -----

   call bgrid_core_driver ( Time_next, Var, Var_dt, Dynam, omega )

!---- call physics -----

   call bgrid_physics ( physics_window, dt_atmos, Time_next,  &
                        Dynam%Hgrid,   Dynam%Vgrid,    Dynam,             &
                        Var,     Var_dt                       )

!---- time differencing and diagnostics -----

   call bgrid_core_time_diff ( omega, Time_next, Dynam, Var, Var_dt )

!-----------------------------------------------------------------------

 end subroutine atmosphere

!#######################################################################

 subroutine atmosphere_init ( Time_init, Time, Time_step )

 type (time_type),     intent(in)    :: Time_init, Time, Time_step

  integer :: unit, sec, ierr, io, iunit

!-----------------------------------------------------------------------
!----- read namelist -----
call find_namelist_in_file("input.nml", "atmosphere_nml", iunit)
read(iunit, nml = atmosphere_nml, iostat = io)
call check_namelist_read(iunit, io, "atmosphere_nml")


!!!    if (file_exist('input.nml')) then
!!!        unit = open_namelist_file ( )
!!!        ierr=1; do while (ierr /= 0)
!!!           read (unit, nml=atmosphere_nml, iostat=io, end=10)
!!!           ierr = check_nml_error (io, 'atmosphere_nml')
!!!        enddo
!!! 10     call close_file (unit)
!!!    endif

!----- write version and namelist to log file -----

    call write_version_number ( version, tag )
    write (stdlog(), nml=atmosphere_nml)

!---- compute physics/atmos time step in seconds ----

   Time_step_atmos = Time_step
   call get_time (Time_step_atmos, sec)
   dt_atmos = sec

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

 end subroutine atmosphere_init

!#######################################################################

 subroutine atmosphere_end

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

    real(r8),    intent(out)          :: blon(:), blat(:)
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
  real(r8),    intent(in)                :: dt_phys
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

   real(r8), dimension(window(1),window(2),Vgrid%nlev) :: p_full, u_dt, v_dt

   real(r8), dimension(window(1),window(2),Vgrid%nlev+1) :: p_half

   real(r8), dimension(Hgrid%ilb:Hgrid%iub, &
                   Hgrid%jlb:Hgrid%jub, &
                   Vgrid%nlev) :: uh, vh, uh_dt, vh_dt

   real(r8), dimension(window(1),window(2)) :: pssl_new
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

end module atmosphere_mod

