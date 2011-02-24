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

module bgrid_core_mod

!-----------------------------------------------------------------------
!
!              gfdl global b-grid dynamical core
!
!           (with eta and hybrid pressure coordinate)
!
!-----------------------------------------------------------------------
!--------------------------- modules -----------------------------------
!-----------------------------------------------------------------------

use types_mod, only : r8
use bgrid_prog_var_mod    , only: prog_var_type, var_init,   &
                                  prog_var_times_scalar
use bgrid_horiz_mod       , only: horiz_grid_type
use bgrid_vert_mod        , only: vert_grid_type, compute_pres_depth, &
                                  compute_pressures
use bgrid_masks_mod       , only: grid_mask_type, grid_masks_init
use bgrid_advection_mod   , only: advec_control_type, &
                                  advection_init, advection
use bgrid_horiz_diff_mod  , only: hdiff_control_type, &
                                  horiz_diff_init, horiz_diff
use bgrid_horiz_adjust_mod, only: horiz_adjust_vel, horiz_adjust_mass,  &
                                  press_grad, compute_grad_pres
use bgrid_vert_adjust_mod , only: vert_adjust
use bgrid_polar_filter_mod, only: pfilt_control_type, polar_filter_init, &
                                  polar_filter, polar_filter_wind
use bgrid_halo_mod        , only: update_halo, TEMP, UWND, VWND
use bgrid_sponge_mod      , only: sponge_driver, sponge_init

use                fms_mod, only: error_mesg, FATAL, write_version_number
use          constants_mod, only:  CP, RDGAS, RVGAS

use  field_manager_mod, only: MODEL_ATMOS
use tracer_manager_mod, only: get_tracer_index

!-----------------------------------------------------------------------

implicit none
private

public  update_bgrid_core, bgrid_core_init, bgrid_core_end

public  bgrid_dynam_type

!-----------------------------------------------------------------------
!  ----- defined data types -----
!    Hgrid = horizontal grid constants
!    Vgrid = vertical grid constants
!    Masks = eta coordinate topography masks and indices
!    Pfilt = polar filter constants
!    Advec = advection constants
!    Hdiff = horizontal diffusion constants
!
!  ----- 2-dimensional (nlon,nlat) fields -----
!    fis  = geopotential height of the surface
!    fisl = geopotential height at eta=1. (for eta coord = 0.0,
!    res  = reciprical of eta at the surface
!
!  ----- time step terms ----
!    nt_adv = no. of advection time steps per atmosphere step (integer)
!    nt_adj = no. of adjustment time steps per advection step (integer)
!    dt_adj = adjustment time step in seconds (real(r8))
!
!    verbose   = verbose flag [integer]
!
!  ----- miscellaneous ----
!    wcorr     = coefficient for the grid separation smoothing operator
!    fopt      = filtering option [integer]
!-----------------------------------------------------------------------

type bgrid_dynam_type
   type(horiz_grid_type), pointer  :: Hgrid
   type (vert_grid_type), pointer  :: Vgrid
   type (grid_mask_type)           :: Masks
   type (pfilt_control_type)       :: Pfilt
   type (advec_control_type)       :: Advec
   type (hdiff_control_type)       :: Hdiff
   real(r8), pointer, dimension(:,:)   :: fis, fisl, res
   integer                         :: nt_adv, nt_adj
   real(r8)                            :: dt_adj
   real(r8)                            :: wcorr
   integer                         :: fopt, verbose
end type bgrid_dynam_type

!-----------------------------------------------------------------------
!------- internal options ---------!  alpha_implicit determines how the
                                   !  coriolis and press grad force
   real(r8)  :: alpha_implicit = 0.5   !  terms are solved
                                   !    = 0.5  trapezoidal implicit
                                   !    = 1.0        fully implicit

!---- un-used leapfrog options ----

   real(r8)  :: smooth = 0.2475
   real(r8)  :: robert = 0.0

!---- internal parameters ----

   real(r8), parameter :: d608 = (RVGAS-RDGAS)/RDGAS

!---- version number ----

   character(len=128) :: version='$Revision$'
   character(len=128) :: tagname='$Id$'

!---- timing data ----

   integer, parameter :: time_level = 4
   integer :: id_
   logical :: do_clock_init = .true.

!-----------------------------------------------------------------------

contains

!#######################################################################

 function bgrid_core_init (Hgrid, Vgrid, fis, res, dt, ntadj, ntadv,   &
              damp_order_vel,  damp_order_tmp,  damp_order_trs,        &
              damp_coeff_vel,  damp_coeff_tmp,  damp_coeff_trs,        &
              damp_scheme, damp_slope_coeff_vel, damp_slope_coeff_tmp, &
              num_horiz_fill, num_vert_fill,                           &
              advec_order_vel, advec_order_tmp, advec_order_trs,       &
              advec_coeff_vel, advec_coeff_tmp, advec_coeff_trs,       &
              grid_sep_coeff, filter_option, filter_weight,            &
              ref_lat_filter, num_sponge_levels,                       &
              sponge_coeff_vel, sponge_coeff_tmp, sponge_coeff_trs,    &
              verbose) result (Dynam)

   type(horiz_grid_type), intent(in), target :: Hgrid
   type (vert_grid_type), intent(in), target :: Vgrid
   real(r8), intent(in), dimension(:,:),  target :: fis, res
   real(r8),            intent(in)               :: dt
   integer,         intent(in)               :: ntadj, ntadv

!          ---- optional arguments ----

   integer,         intent(in), optional :: damp_order_vel, &
                                            damp_order_tmp, &
                                            damp_order_trs
   real(r8),            intent(in), optional :: damp_coeff_vel, &
                                            damp_coeff_tmp, &
                                            damp_coeff_trs, &
                                            damp_slope_coeff_vel(4), &
                                            damp_slope_coeff_tmp(4)
   integer,         intent(in), optional :: damp_scheme,    &
                                            num_horiz_fill, &
                                            num_vert_fill,  &
                                            advec_order_vel,&
                                            advec_order_tmp,&
                                            advec_order_trs
   real(r8),            intent(in), optional :: advec_coeff_vel, &
                                            advec_coeff_tmp, &
                                            advec_coeff_trs, &
                                            grid_sep_coeff
   integer,         intent(in), optional :: filter_option,   &
                                            filter_weight,   &
                                            verbose,         &
                                            num_sponge_levels
   real(r8),            intent(in), optional :: ref_lat_filter,   &
                                            sponge_coeff_vel, &
                                            sponge_coeff_tmp, &
                                            sponge_coeff_trs

   type(bgrid_dynam_type) :: Dynam

!-----------------------------------------------------------------------
!
!    performs initialization for b-grid dynamics type
!
! input:  Hgrid      horizontal grid constants
!         Vgrid      vertical grid constants (see vert_grid_mod)
!         fis        geopotential height of the surface,
!                    dimensioned by (nlon,nlat)
!         res        reciprocal of eta at the surface,
!                    dimensioned by (nlon,nlat)
!
!         dt         time step in seconds for each update call
!         ntadj      number of adjustment time steps for each
!                    advective time step [integer]
!         ntadv      number of advection time steps for each
!                    update call [integer]
!
!         IMPORTANT:  The input arguments (Hgrid, Vgrid, fis, res)
!                     must have permanent storage for the length
!                     of the model integration at some higher level
!
! input (optional):
!
!  damp_order_vel       The horizontal damping order for momentum,
!  damp_order_tmp       temperature, and default order for all 
!  damp_order_trs       prognostic tracers.
!
!  damp_coeff_vel       The horizontal damping coefficients for
!  damp_coeff_tmp       momentum, temperature, and default value for
!  damp_coeff_trs       all prognostic tracers.
!
!  damp_scheme          Determines how horizontal damping coefficients
!                       vary with latitude.
!
!  damp_slope_corr_vel  The topography slope correction for horizontal
!  damp_slope_corr_tmp  damping of momentum and temperature (including
!                       all prognostic tracers).
!
!  advec_order_vel     The advection order for momentum, temperature,
!  advec_order_tmp     and default order for all prognostic tracers.
!  advec_order_trs   
!
!  advec_coeff_vel     Coefficients for modified Euler-backward advection
!  advec_coeff_tmp     scheme for momentum, temperature, and all
!  advec_coeff_trs     prognostic tracers.
!
!  num_horiz_fill      The number of successive horizontal and vertical passes
!  num_vert_fill       applied in the tracer borrowing/filling scheme.
!  
!  grid_sep_coeff      Coefficient to suppress grid-separation problem
!                      associated with the B-grid. Currently, this option 
!                      has been disabled within the model, so that this
!                      coefficient does nothing.
!
!  filter_option       Determines how polar filtering is performed.
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
!  sponge_coeff_vel    Damping coefficients for the sponge layer in
!  sponge_coeff_tmp    the uppermost model level.
!  sponge_coeff_trs    
!
!  verbose             Flag that control additional printed output
!                      Currently, this option is not being used.
!
! NOTE: also see bgrid_core_driver for description of optional arguments
!
!-----------------------------------------------------------------------
!  ---- required time step arguments -----

   Dynam % nt_adj = ntadj
   if (Dynam % nt_adj <= 0) call error_mesg ('bgrid_core_init',  &
                            'input argument ntadj must be >= 1', FATAL)

   Dynam % nt_adv = ntadv
   if (Dynam % nt_adv <= 0) call error_mesg ('bgrid_core_init',  &
                            'input argument ntadv must be >= 1', FATAL)

   Dynam % dt_adj = dt / float(Dynam%nt_adv*Dynam%nt_adj)
   if (Dynam % dt_adj <= 0.0) call error_mesg ('bgrid_core_init',  &
                             'input argument dtadv must be > 0.', FATAL)

!  ---- optional arguments ----

   Dynam % wcorr = 0.25
   if (present(grid_sep_coeff)) Dynam % wcorr = grid_sep_coeff

   Dynam % fopt = 2
   if (present(filter_option)) Dynam % fopt = filter_option

   Dynam % verbose = 0
   if (present(verbose)) Dynam % verbose = max(0, verbose)

!-----------------------------------------------------------------------
!     ------ pointers ------

   Dynam % Hgrid => Hgrid
   Dynam % Vgrid => Vgrid
   Dynam % fis   => fis
   Dynam % res   => res

!-----------------------------------------------------------------------
!     ------- allocate space for other variables -------

   Dynam%fisl => var_init (Dynam%Hgrid)

!-----------------------------------------------------------------------
!  ---- eta coordinate masks ----

   Dynam%Masks = grid_masks_init (Dynam%Hgrid, Dynam%Vgrid, Dynam%res)

!  ---- define sea level geop height ----
   if (Dynam%Masks%sigma) then
        Dynam%fisl = Dynam%fis
   else
        Dynam%fisl = 0.0
   endif

!-----------------------------------------------------------------------
!       ---- initialize polar filtering routines -----

   Dynam%Pfilt = polar_filter_init (Dynam%Hgrid, reflat=ref_lat_filter, &
                         weight=filter_weight, sigma=Dynam%Masks%sigma, &
                         verbose=Dynam%verbose )

!-----------------------------------------------------------------------
!------- initialize calls to bgrid_core routines ----------

     Dynam%Advec = advection_init (Dynam%Hgrid,                 &
               advec_order_vel, advec_order_tmp, advec_order_trs, &
               advec_coeff_vel, advec_coeff_tmp, advec_coeff_trs, &
               num_horiz_fill, num_vert_fill)

     Dynam%Hdiff = horiz_diff_init (Hgrid, Dynam%nt_adj, Vgrid%nplev,  &
                   damp_order_vel, damp_order_tmp, damp_order_trs, &
                   damp_coeff_vel, damp_coeff_tmp, damp_coeff_trs, &
                   damp_slope_coeff_vel, damp_slope_coeff_tmp,     &
                   fis=Dynam%fisl, Vgrid=Dynam%Vgrid,                  &
                   damping_scheme=damp_scheme                          )

     call sponge_init ( sponge_coeff_vel, sponge_coeff_tmp, sponge_coeff_trs, &
                        num_sponge_levels )

!-----------------------------------------------------------------------

     call write_version_number (version,tagname)

!-----------------------------------------------------------------------

 end function bgrid_core_init

!#######################################################################

 subroutine update_bgrid_core (Var, Var_dt, Dynam, omega)

!-----------------------------------------------------------------------
   type (prog_var_type), intent(in)    :: Var
   type (prog_var_type), intent(inout) :: Var_dt
   type(bgrid_dynam_type), intent(inout) :: Dynam
   real(r8), intent(out), dimension(Dynam%Hgrid%ilb:Dynam%Hgrid%iub, &
                                Dynam%Hgrid%jlb:Dynam%Hgrid%jub, &
                                Dynam%Vgrid%nlev) :: omega
!-----------------------------------------------------------------------

real(r8), dimension(Dynam%Hgrid%ilb:Dynam%Hgrid%iub, &
                Dynam%Hgrid%jlb:Dynam%Hgrid%jub) :: psdt

real(r8), dimension(Dynam%Hgrid%ilb:Dynam%Hgrid%iub, &
                Dynam%Hgrid%jlb:Dynam%Hgrid%jub) :: pssl

real(r8), dimension(Dynam%Hgrid%ilb:Dynam%Hgrid%iub, &
                Dynam%Hgrid%jlb:Dynam%Hgrid%jub, Dynam%Vgrid%nlev) :: &
                dpde, few, fns, div, pgfew, pgfns, dpde_old,          &
                pfull, wta, wtb, cew, cns, u, v, tq, uo, vo, to

real(r8), dimension(Dynam%Hgrid%ilb:Dynam%Hgrid%iub, &
                Dynam%Hgrid%jlb:Dynam%Hgrid%jub, Dynam%Vgrid%nlev+1)  &
                :: phalf, etadot

   real(r8)    :: fadv, tdt_adj, rscale
   integer :: i, k, n, m, nt
!-----------------------------------------------------------------------
!    ---- definition of local variables ----
!
!    ps    = surface pressure
!    psdt  = surface pressure tendency
!    pssl  = surface pressure adjusted to eta=1.
!    tq    = (virtual) temperature
!    div   = mass divergence
!    omega = thermodynamic (omega) term in some form
!    pgfew = zonal pressure gradient force (at v pts)
!    pgfns = meridional pressure gradient force (at v pts)
!    dpde  = pressure thickness of model layers
!    pfull = pressure at full model levels
!    phalf = pressure at half model levels (between full levels)
!    wta,
!      wtb = weights for obtaining the pressure at full levels
!    cew,
!      cns = grad(p)/p term (zonal and meridional components)
!    few,
!      fns = mass fluxes time summed over advection interval
!            (zonal and meridional components)
!   etadot = vertical mass flux summed over advection interval
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!------- number of tracers and variable (time) levels ---------

    nt = Var%ntrace

!   ---- set-up time step ----

    tdt_adj = Dynam % dt_adj

!   ------ pressure variables ------

    pssl = Var % pssl + tdt_adj * Var_dt % pssl

    call compute_pressures ( Dynam%Vgrid, pssl,           &
                             phalf, pfull, dpde, wta, wtb )

    call compute_grad_pres ( Dynam%Hgrid, Dynam%Vgrid%nplev, &
                             phalf, dpde, wta, wtb, cew, cns )

!   ------ zero fluxes -----

    few = 0.0;  fns = 0.0;  etadot = 0.0

!-----------------------------------------------------------------------
!************** start of bgrid_core time step loop *********************
 do m = 1, Dynam % nt_adv

!   ------ save this time level for advection ------

    dpde_old = dpde

    to = Var % t + tdt_adj * Var_dt % t
    uo = Var % u + tdt_adj * Var_dt % u
    vo = Var % v + tdt_adj * Var_dt % v

!   ------ variables at current time level ------

    tq = to
    call compute_virt_temp ( tq, Var%r, Var_dt%r, tdt_adj )
    u  = uo
    v  = vo

!-----------------------------------------------------------------------
!*************** start of adjustment time step loop ********************

 do n = 1, Dynam % nt_adj

!-----------------------------------------------------------------------
!------------------------- compute divergence --------------------------
!--------------------- horizontal alpha-omega term ---------------------

   call horiz_adjust_mass ( Dynam%Vgrid%nplev, Dynam%Hgrid, &
                            Dynam%Masks,       u,           &
                            v,                 dpde,        &
                            cew,               cns,         &
                            few,               fns,         &
                            div,               omega        )

!    ---- polar filtering -----

     if (Dynam%fopt >= 1) then
        call polar_filter (Dynam%Pfilt, div, omega, 1, Dynam%Masks%Tmp%mask)
     endif

     call update_halo (Dynam%Hgrid, TEMP, div  )
     call update_halo (Dynam%Hgrid, TEMP, omega)

!-----------------------------------------------------------------------
!------- compute sfc pres tendency, vert. alpha-omega, & vert vel.------
    call vert_adjust ( Dynam%res, div, wta, wtb,          &
                       Dynam%Masks%Tmp%mask, Dynam%Vgrid, &
                       omega, etadot, psdt                )     

!-----------------------------------------------------------------------
!-------- update surface pressure tendency -----------------------------

    Var_dt % ps   = Var_dt % ps   + psdt
    Var_dt % pssl = Var_dt % pssl + psdt * Dynam % res

!-------- recompute pressure variables at next time level -----------

    pssl = Var % pssl + tdt_adj * Var_dt % pssl
    call compute_pressures ( Dynam%Vgrid, pssl, phalf, pfull, dpde, &
                             wta, wtb )

!-------- update thermodynamic tendency ----------

    omega = omega/dpde
    Var_dt % t = Var_dt % t + omega * tq * RDGAS/CP

!-----------------------------------------------------------------------
!------------ compute geopotential height and rt/p ---------------------
!     ----- (use smoothed value for pssl if leapfrog) ----

    tq = Var % t + tdt_adj * Var_dt % t
    call compute_virt_temp ( tq, Var%r, Var_dt%r, tdt_adj )

    call compute_grad_pres ( Dynam%Hgrid, Dynam%Vgrid%nplev, &
                             phalf, dpde, wta, wtb, cew, cns )

    call press_grad ( Dynam%Hgrid, Dynam%Vgrid, Dynam%Masks,    &
                      Dynam%fisl, tq, dpde, wta, wtb, cew, cns, &
                      pgfew, pgfns )

!------------------- adjustment of wind components ---------------------
    call horiz_adjust_vel ( Dynam%Hgrid, Dynam%Masks,                &
                                        tdt_adj, pgfew, pgfns,       &
                            u, v, Var_dt%u, Var_dt%v, alpha_implicit )

  ! polar filtering of momentum (new scheme)
    if (Dynam%fopt == 2)  &
    call prog_var_filter_vel ( Dynam%Pfilt, tdt_adj, Dynam%Masks%Vel%mask, &
                               Var%u, Var%v, Var_dt%u, Var_dt%v            )

!-----------------------------------------------------------------------
!----------------------- advection -------------------------------------
    if (n ==  Dynam%nt_adj) then
         call advection ( Dynam%Advec, Dynam%Pfilt,              &
                          Dynam%Hgrid, Dynam%Vgrid, Dynam%Masks, &
                          Dynam%fopt,  tdt_adj,                  &
                          dpde_old, dpde, few, fns, etadot,      &
                          uo, vo, to,   Var, Var_dt              )
    endif

!-----------------------------------------------------------------------

  ! polar filtering of momentum (old scheme)
    if (Dynam%fopt == 1)  &
    call prog_var_filter_vel ( Dynam%Pfilt, tdt_adj, Dynam%Masks%Vel%mask, &
                               Var%u, Var%v, Var_dt%u, Var_dt%v            )
!   ---- recompute momentum at next time level ----

    call update_halo (Dynam%Hgrid, UWND, Var_dt%u )
    call update_halo (Dynam%Hgrid, VWND, Var_dt%v )

    u = Var % u + tdt_adj * Var_dt % u
    v = Var % v + tdt_adj * Var_dt % v

!-----------------------------------------------------------------------

 enddo

!-------- polar filtering of temp and tracers --------------------------

!!  if (Dynam%fopt == 1) then
!!      call prog_var_filter_mass &
!!                      ( Dynam%Pfilt, tdt_adj, dpde, Dynam%Masks%Tmp%mask, &
!!                        Var%t, Var%r, Var_dt%t, Var_dt%r                  )
!!      call update_halo (Dynam%Hgrid, TEMP, Var_dt%t )
!!      call update_halo (Dynam%Hgrid, TEMP, Var_dt%r )
!!  endif

!-----------------------------------------------------------------------
!----------------- horizontal diffusion --------------------------------
    call horiz_diff ( Dynam%Hdiff, Dynam%Masks, tdt_adj, dpde, pfull, &
                      Var, Var_dt )

!-----------------------------------------------------------------------

    call sponge_driver ( Dynam%Hgrid, Dynam%Vgrid%nplev, &
                         tdt_adj, dpde, Var, Var_dt      )

!    ---------------- halo updates --------------------
!    skip on the last pass will update last when needed

    if (m == Dynam % nt_adv) cycle
    call update_halo (Dynam%Hgrid, TEMP, Var_dt%t )
    call update_halo (Dynam%Hgrid, TEMP, Var_dt%r )
    call update_halo (Dynam%Hgrid, UWND, Var_dt%u )
    call update_halo (Dynam%Hgrid, VWND, Var_dt%v )

!-----------------------------------------------------------------------

 enddo

!-----------------------------------------------------------------------
!  ---- compute omega diagnostic ----

        omega = omega * pfull

!  ---- scale tendencies before physics for time step difference ----

   rscale = 1.0 / float(Dynam%nt_adj*Dynam%nt_adv)
   call prog_var_times_scalar (Var_dt, rscale)

!-----------------------------------------------------------------------

 end subroutine update_bgrid_core
        
!#######################################################################

 subroutine bgrid_core_end (Dynam)

!-----------------------------------------------------------------------
   type(bgrid_dynam_type), intent(inout) :: Dynam
!-----------------------------------------------------------------------

!  **** may want to deallocate arrays ? ****

!-----------------------------------------------------------------------

 end subroutine bgrid_core_end

!#######################################################################

 subroutine compute_virt_temp ( tq, r, rdt, dt )

  real(r8), intent(inout) :: tq(:,:,:)
  real(r8), intent(in)    :: r (:,:,:,:), rdt(:,:,:,:), dt

  real(r8), dimension(size(r,1),size(r,2),size(r,3)) :: q
  integer :: sphum

!  given the temp and tracers (where tracer 1 is assumed to be
!  specific humidity) the virtual temperature is computed

    if ( size(r,4) == 0 ) return
    sphum = get_tracer_index ( MODEL_ATMOS, 'sphum' )
    if ( sphum <= 0 ) return

    q  = r(:,:,:,sphum) + dt * rdt(:,:,:,sphum)
    tq = tq * ( 1. + d608 * q )

 end subroutine compute_virt_temp

!#######################################################################
!#### various polar filter routines (some of which may be used) ########
!#######################################################################

 subroutine prog_var_filter_mass ( Pfilt, dt, dpde, mask,  &
                                   t, r, tdt, rdt )

 type(pfilt_control_type), intent(in)    :: Pfilt
 real(r8),                     intent(in)    :: dt
 real(r8), dimension(:,:,:),   intent(in)    :: dpde, mask, t
 real(r8), dimension(:,:,:,:), intent(in)    :: r
 real(r8), dimension(:,:,:),   intent(inout) :: tdt
 real(r8), dimension(:,:,:,:), intent(inout) :: rdt

   real(r8), dimension(size(t,1),size(t,2),size(t,3)) :: work
   integer :: n

!  ---- temperature ----
     work = (t + dt*tdt) * dpde 
     call polar_filter ( Pfilt, work, 1, mask )
     tdt = (work/dpde-t)/dt

!  ---- tracers ----
   do n = 1, size(r,4)
     work = (r(:,:,:,n) + dt*rdt(:,:,:,n)) * dpde
     call polar_filter ( Pfilt, work, 1, mask )
     rdt(:,:,:,n) = (work/dpde-r(:,:,:,n))/dt
   enddo

 end subroutine prog_var_filter_mass

!---------------------------------------------------------------------

 subroutine prog_var_filter_vel ( Pfilt, dt, mask, u, v, udt, vdt )

 type(pfilt_control_type), intent(in)    :: Pfilt
 real(r8),                     intent(in)    :: dt
 real(r8), dimension(:,:,:),   intent(in)    :: mask
 real(r8), dimension(:,:,:),   intent(in)    :: u, v
 real(r8), dimension(:,:,:),   intent(inout) :: udt, vdt

   real(r8), dimension(size(u,1),size(u,2),size(u,3)) :: ut, vt

!  ---- momentum ----
     ut = u + dt*udt
     vt = v + dt*vdt
     call polar_filter_wind ( Pfilt, ut, vt, mask )
     udt = (ut-u)/dt
     vdt = (vt-v)/dt

 end subroutine prog_var_filter_vel

!---------------------------------------------------------------------

 subroutine prog_tend_filter_mass ( Pfilt, mask, dpde, tdt, rdt,  &
                                    dpdt, t, r )

 type(pfilt_control_type), intent(in)    :: Pfilt
 real(r8), dimension(:,:,:),   intent(in)    :: mask, dpde
 real(r8), dimension(:,:,:),   intent(inout) :: tdt
 real(r8), dimension(:,:,:,:), intent(inout) :: rdt

 real(r8), dimension(:,:,:),   intent(in), optional :: dpdt, t
 real(r8), dimension(:,:,:,:), intent(in), optional :: r

   real(r8), dimension(size(tdt,1),size(tdt,2),size(tdt,3)) :: work
   integer :: n

!  ---- temperature ----
     work = tdt * dpde
     if (present(dpdt)) work = work - t*dpdt
     call polar_filter ( Pfilt, work, 1, mask )
     tdt = work/dpde

!  ---- tracers ----
   do n = 1, size(r,4)
     work = rdt(:,:,:,n) * dpde
     if (present(dpdt)) work = work - r(:,:,:,n)*dpdt
     call polar_filter ( Pfilt, work, 1, mask )
     rdt(:,:,:,n) = work/dpde
   enddo

 end subroutine prog_tend_filter_mass

!#######################################################################
!#######################################################################

end module bgrid_core_mod

