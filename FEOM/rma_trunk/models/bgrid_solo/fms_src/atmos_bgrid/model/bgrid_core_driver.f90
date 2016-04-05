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

module bgrid_core_driver_mod

!-----------------------------------------------------------------------
!
!      wrapper module for initializing the b-grid dynamics core
!
!         reads namelist
!         initializes bgrid_core_mod
!
!-----------------------------------------------------------------------
!---------------- m o d u l e   i n f o r m a t i o n ------------------

use types_mod, only : r8
use bgrid_core_mod           , only: bgrid_dynam_type,  &
                                     bgrid_core_init, update_bgrid_core,  &
                                     bgrid_core_end
use bgrid_horiz_mod          , only: horiz_grid_type, horiz_grid_init
use bgrid_vert_mod           , only: vert_grid_type, vert_grid_init
use bgrid_prog_var_mod       , only: prog_var_type, prog_var_init, &
                                     prog_var_time_diff, var_init, &
                                     open_prog_var_file,           &
                                     read_prog_var, write_prog_var
use bgrid_halo_mod           , only: horiz_grid_boundary
use bgrid_diagnostics_mod    , only: bgrid_diagnostics,      &
                                     bgrid_diagnostics_tend, &
                                     bgrid_diagnostics_init
use bgrid_integrals_mod      , only: bgrid_integrals, bgrid_integrals_init, &
                                     bgrid_integrals_end
use bgrid_conserve_energy_mod, only: bgrid_conserve_energy_init, &
                                     bgrid_conserve_energy

use  field_manager_mod, only: MODEL_ATMOS
use tracer_manager_mod, only: get_number_tracers
use   time_manager_mod, only: time_type, get_time
use utilities_mod, only : find_namelist_in_file, check_namelist_read
use            fms_mod, only: error_mesg, FATAL, file_exist, open_namelist_file,  &
                              write_version_number,     &
                              close_file, stdlog

!-----------------------------------------------------------------------

implicit none
private

public  bgrid_dynam_type, bgrid_core_driver_init, &
        bgrid_core_driver, bgrid_core_driver_end, &
        bgrid_core_time_diff, get_bottom_data, put_bottom_data

!-----------------------------------------------------------------------
character(len=128) :: version =  '$Revision$'
character(len=128) :: tag =  '$Id$'
!-----------------------------------------------------------------------
!
!             NAMELIST INPUT: bgrid_core_driver_nml
!
!           This namelist is read from file input.nml.
!           See the on-line documentation for more details.
!
!-----------------------------------------------------------------------
!
!  damp_scheme          Determines how horizontal damping coefficients
!                       vary with latitude.
!
!  damp_order_wind      The horizontal damping order for momentum,
!  damp_order_temp      temperature, and default order for all 
!  damp_order_tracer    prognostic tracers.
!
!  damp_coeff_wind      The horizontal damping coefficients for
!  damp_coeff_temp      momentum, temperature, and default value for
!  damp_coeff_tracer    all prognostic tracers.
!
!  slope_corr_wind      The topography slope correction for horizontal
!  slope_corr_temp      damping of momentum and temperature (including
!                       all prognostic tracers).
!
!  advec_order_wind     The advection order for momentum, temperature,
!  advec_order_temp     and default order for all prognostic tracers.
!  advec_order_tracer
!
!  advec_coeff_wind     Coefficients for modified Euler-backward advection
!  advec_coeff_temp     scheme for momentum, temperature, and all
!  advec_coeff_tracer   prognostic tracers.
!
!  num_fill_pass        The number of successive passes applied in the tracer
!                       borowing/filling scheme.
!
!  filter_option        Determines how polar filtering is performed.
! 
!  filter_weight       Weight applied to the polar filter that will
!                      increase (or decrease) the strength of the standard
!                      polar filter response function.
!
!  ref_lat_filter      The reference latitude at which polar filtering
!                      (in each hemisphere) will begin to be applied.
!
!  num_sponge_levels   Number of uppermost model level where a band-pass
!                      filter is applied to damp undesireable waves.
!
!  sponge_coeff_wind     Damping coefficients for the sponge layer in
!  sponge_coeff_temp     the uppermost model level.
!  sponge_coeff_tracer
!
!  halo                The number of halo rows along all (NEWS) boundaries.
!                      There is currently no namelist option that allows unequal
!                      halo boundary.
!
!  num_adjust_dt       The number of adjustment time steps for each advection
!                      time step, where num_adjust_dt >= 1.
!
!  num_advec_dt        The number of advection/dynamics time steps for each
!                      atmospheric/physics time step, where num_advec_dt >= 1.
!
!  decomp              The domain decomposition, where decomp(1) = x-axis
!                      decomposition, decomp(2) = y-axis decomposition.
!                      * If decomp(1)*decomp(2) does not equal the number
!                        of processors the model will fail.
!                      * If decomp(1)=decomp(2)=0 then default rules apply.
!                      * By default, one-dimensional decomposition (in Y) is used.
!                        When there is fewer than 2 points per processor, then 2-D
!                        decomposition is used.
!
!  do_conserve_energy  If TRUE the temperature tendency will be updated to
!                      guarantee that the dynamical core conserves total energy.
!                      The correction is applied to a uniform global value.
!
!  grid_sep_coeff      Coefficient to suppress grid-separation problem
!                      associated with the B-grid. Currently, this option 
!                      has been disabled within the model, so that this
!                      coefficient does nothing.
!
!  verbose             Flag that control additional printed output
!                      Currently, this option is not being used.
!


   integer   ::     damp_scheme       = 1
   integer   ::     damp_order_wind   = 4
   integer   ::     damp_order_temp   = 4
   integer   ::     damp_order_tracer = 4
   real(r8)      ::     damp_coeff_wind   = 0.50
   real(r8)      ::     damp_coeff_temp   = 0.50
   real(r8)      ::     damp_coeff_tracer = 0.50

   integer   ::     advec_order_wind   = 2
   integer   ::     advec_order_temp   = 2
   integer   ::     advec_order_tracer = 2
   real(r8)      ::     advec_weight_wind   = 0.7
   real(r8)      ::     advec_weight_temp   = 0.7
   real(r8)      ::     advec_weight_tracer = 0.7

   integer   ::     num_fill_pass = 1
   integer   ::     filter_option = 2
   integer   ::     filter_weight = 1
   real(r8)      ::     ref_lat_filter = 60.
   integer   ::     verbose = -1
   real(r8)      ::     grid_sep_coeff = 0.00

   integer   ::     num_sponge_levels = 0
   real(r8)      ::     sponge_coeff_wind   = 0.0
   real(r8)      ::     sponge_coeff_temp   = 0.0
   real(r8)      ::     sponge_coeff_tracer = 0.0

   real(r8), dimension(4) :: slope_corr_wind = (/0.,0.,0.,0./)
   real(r8), dimension(4) :: slope_corr_temp = (/0.,0.,0.,0./)

   logical   :: do_conserve_energy = .false.
   integer   :: num_adjust_dt = 3
   integer   :: num_advec_dt = 1

   integer   :: halo = 1
   integer, dimension(2) :: decomp = (/0,0/)


   namelist /bgrid_core_driver_nml/ damp_order_wind,                 &
                                    damp_order_temp,                 &
                                    damp_order_tracer,               &
                                    damp_coeff_wind,                 &
                                    damp_coeff_temp,                 &
                                    damp_coeff_tracer,               &
                                    damp_scheme,                     &
                                    advec_order_wind,                &
                                    advec_order_temp,                &
                                    advec_order_tracer,              &
                                    advec_weight_wind,               &
                                    advec_weight_temp,               &
                                    advec_weight_tracer,             &
                                    num_fill_pass, grid_sep_coeff,   &
                                    filter_option, filter_weight,    &
                                    ref_lat_filter,                  &
                                    num_sponge_levels,               &
                                    sponge_coeff_wind,               &
                                    sponge_coeff_temp,               &
                                    sponge_coeff_tracer,             &
                                    slope_corr_wind,                 &
                                    slope_corr_temp,                 &
                                    verbose,                         &
                                    do_conserve_energy,              &
                                    num_adjust_dt, num_advec_dt,     &
                                    halo, decomp

!-----------------------------------------------------------------------
!------ private data ------

real(r8), dimension(:,:),   pointer :: fis, res
real(r8), dimension(:), allocatable :: eta, peta

type  (horiz_grid_type), target      :: Hgrid
type   (vert_grid_type), target      :: Vgrid

integer, dimension(4) :: mass_axes, vel_axes

real(r8) :: dt_atmos

logical :: channel_model = .false.
real(r8)    :: tph0d_in = 0.0, tlm0d_in = 0.0

!-----------------------------------------------------------------------

contains

!#######################################################################

 subroutine bgrid_core_driver_init ( Time_init, Time, Time_step, &
                                     Var, Var_dt, Dynam, phys_axes )

 type       (time_type), intent(in)    :: Time_init, Time, Time_step
 type   (prog_var_type), intent(inout) :: Var, Var_dt
 type(bgrid_dynam_type), intent(inout) :: Dynam
 integer,                intent(out)   :: phys_axes(4)

 integer :: unit, io, ierr, iunit
 integer :: ix, jx, kx
 integer :: sec, ntrace, ntprog, ntdiag
!-----------------------------------------------------------------------
!  ----- read namelist -----
! fms replaced with dart 8 June, 2006
call find_namelist_in_file("input.nml", "bgrid_core_driver_nml", iunit)
read(iunit, nml = bgrid_core_driver_nml, iostat = io)
call check_namelist_read(iunit, io, "bgrid_core_driver_nml")

!!!    if (file_exist('input.nml')) then
!!!        unit = open_namelist_file ( )
!!!        ierr=1; do while (ierr /= 0)
!!!           read (unit, nml=bgrid_core_driver_nml, iostat=io, end=10)
!!!           ierr = check_nml_error (io, 'bgrid_core_driver_nml')
!!!        enddo
!!! 10     call close_file (unit)
!!!    endif

!-----------------------------------------------------------------------
!  ----- read restart header records and set up grid resolution -----
   call open_prog_var_file (ix, jx, kx)

!  ---- horizontal grid initialization ----
   Hgrid = horiz_grid_init (ix, jx, ihalo=halo, jhalo=halo,       &
                                    decomp=decomp,                &
                                    channel=channel_model,        &
                                    tph0d=tph0d_in, tlm0d=tlm0d_in)
   call horiz_grid_boundary (Hgrid)

! how many tracers have been registered?
   call get_number_tracers ( MODEL_ATMOS, num_tracers=ntrace, num_prog=ntprog, num_diag=ntdiag )
!  ----- write version, namelist and tracer info to log file -----

    call write_version_number (version, tag)
        write (stdlog(), nml=bgrid_core_driver_nml)
        write (stdlog(), '(a,i3)') 'Number of tracers =', ntrace
        write (stdlog(), '(a,i3)') 'Number of prognostic tracers =', ntprog
        write (stdlog(), '(a,i3)') 'Number of diagnostic tracers =', ntdiag

!  ---- prognostic variable initialization -----
   call prog_var_init (Hgrid, kx, ntrace, Var)     ! prognostic+diagnostic tracers
   call prog_var_init (Hgrid, kx, ntprog, Var_dt)

!----- read data -----

   fis   => var_init (Hgrid)
   res   => var_init (Hgrid)

   allocate (eta(kx+1), peta(kx+1))

   call read_prog_var (Hgrid, Var, eta, peta, fis, res)


!---- vertical grid initialization ----

   Vgrid = vert_grid_init (eta, peta)

   deallocate (eta, peta)

!---- compute time step in seconds ----

   call get_time (Time_step, sec)
   dt_atmos = sec

!-----------------------------------------------------------------------
!  ----- initialize dynamical core -----
   Dynam = bgrid_core_init (Hgrid, Vgrid, fis, res, dt_atmos,          &
                          num_adjust_dt, num_advec_dt,                 &
                          damp_order_wind, damp_order_temp,            &
                          damp_order_tracer,                           &
                          damp_coeff_wind, damp_coeff_temp,            &
                          damp_coeff_tracer, damp_scheme,              &
                          slope_corr_wind, slope_corr_temp,            &
                          num_fill_pass, num_fill_pass,                &
                          advec_order_wind, advec_order_temp,          &
                          advec_order_tracer,  advec_weight_wind,      &
                          advec_weight_temp, advec_weight_tracer,      &
                          grid_sep_coeff,                              &
                          filter_option, filter_weight, ref_lat_filter,&
                          num_sponge_levels, sponge_coeff_wind,        &
                          sponge_coeff_temp, sponge_coeff_tracer,      &
                          verbose )

!-----------------------------------------------------------------------
!---- initialize history (netcdf) file and integrals -------

   call bgrid_diagnostics_init ( Time, Hgrid, Vgrid, Var,      &
                                 fis, res, mass_axes, vel_axes )
   phys_axes = mass_axes

   call bgrid_integrals_init (Time_init, Time)

!---- initialize integrals ----

   call bgrid_integrals (Time, Hgrid, Vgrid, Var, Dynam%Masks)

!---- initialize energy conservation module ----

   if (do_conserve_energy) call bgrid_conserve_energy_init (Time, mass_axes)

!-----------------------------------------------------------------------

 end subroutine bgrid_core_driver_init

!#######################################################################

 subroutine bgrid_core_driver (Time_diag, Var, Var_dt, Dynam, omega)

!-----------------------------------------------------------------------
   type       (time_type), intent(in)    :: Time_diag
   type   (prog_var_type), intent(in)    :: Var
   type   (prog_var_type), intent(inout) :: Var_dt
   type(bgrid_dynam_type), intent(inout) :: Dynam
   real(r8),                   intent(out)   :: omega(:,:,:)
!-----------------------------------------------------------------------
!  dynamics
   call update_bgrid_core (Var, Var_dt, Dynam, omega)

!  energy conservation
   if (do_conserve_energy) then
       call bgrid_conserve_energy ( dt_atmos, Time_diag, Hgrid, Vgrid, &
                                    Dynam%Masks, Var, Var_dt )
   endif

!  diagnostics for dynamics tendencies
   call bgrid_diagnostics_tend ( Hgrid, Var_dt, Dynam%Masks, Time_diag )

!-----------------------------------------------------------------------

 end subroutine bgrid_core_driver

!#######################################################################

 subroutine bgrid_core_time_diff ( omega, Time_diag, Dynam, Var, Var_dt )

!-----------------------------------------------------------------------
   real(r8),                   intent(in)    :: omega(:,:,:)
   type       (time_type), intent(in)    :: Time_diag
   type(bgrid_dynam_type), intent(in)    :: Dynam
   type   (prog_var_type), intent(inout) :: Var
   type   (prog_var_type), intent(inout) :: Var_dt
!-----------------------------------------------------------------------

!  time differencing

   call prog_var_time_diff ( dt_atmos, Var_dt, Var )

!  global integrals and diagnostics

   call bgrid_integrals ( Time_diag, Hgrid, Vgrid, Var, Dynam%Masks )

   call bgrid_diagnostics ( Hgrid, Vgrid, Var, Dynam%Masks, &
                            Time_diag, omega ) 

!-----------------------------------------------------------------------

 end subroutine bgrid_core_time_diff

!#######################################################################

 subroutine bgrid_core_driver_end ( Var, Dynam )

!-----------------------------------------------------------------------
   type   (prog_var_type), intent(in)    :: Var
   type(bgrid_dynam_type), intent(inout) :: Dynam
!-----------------------------------------------------------------------
!  terminate dynamics

   call bgrid_core_end ( Dynam )

!  write restart for prognostic variables

   call write_prog_var ( Var, Hgrid, Vgrid, fis, res )

!  terminate integrals
   call bgrid_integrals_end

!-----------------------------------------------------------------------

 end subroutine bgrid_core_driver_end

!#######################################################################
!    The following routines do not really belong in this module
!    but since they are not used within the core itself they
!    will reside here temporarily.
!#######################################################################

 subroutine get_bottom_data ( a, b, a_bot, b_bot, k_bot ) 

  real(r8)   , intent(in) , dimension(:,:,:) :: a    , b
  real(r8)   , intent(out), dimension(:,:)   :: a_bot, b_bot
  integer, intent(in) , dimension(:,:), optional :: k_bot

! returns the lowest level data (a_bot,b_bot) from 3d fields (a,b)

  integer :: i, j, kb, kb_min
    
                      kb_min = size(a,3)
  if (present(k_bot)) kb_min = minval (k_bot) 

  if ( kb_min == size(a,3) ) then
          a_bot = a(:,:,kb_min)
          b_bot = b(:,:,kb_min)
  else
       do j = 1, size(a,2)
       do i = 1, size(a,1)
          kb = k_bot(i,j)
          a_bot(i,j) = a(i,j,kb)
          b_bot(i,j) = b(i,j,kb)
       enddo   
       enddo   
  endif

 end subroutine get_bottom_data

!#######################################################################

 subroutine put_bottom_data ( a_bot, b_bot, a, b, k_bot ) 

  real(r8)   , intent(in)   , dimension(:,:)   :: a_bot, b_bot
  real(r8)   , intent(inout), dimension(:,:,:) :: a    , b
  integer, intent(in)   , dimension(:,:), optional :: k_bot

! inserts the lowest level data (a_bot,b_bot) into 3d fields (a,b)

  integer :: i, j, kb, kb_min

                      kb_min = size(a,3)
  if (present(k_bot)) kb_min = minval (k_bot)

  if ( kb_min == size(a,3) ) then
          a(:,:,kb_min) = a_bot
          b(:,:,kb_min) = b_bot
  else
       do j = 1, size(a,2)
       do i = 1, size(a,1)
          kb = k_bot(i,j)
          a(i,j,kb) = a_bot(i,j)
          b(i,j,kb) = b_bot(i,j)
       enddo
       enddo
  endif

 end subroutine put_bottom_data

!#######################################################################

end module bgrid_core_driver_mod

