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

module bgrid_diagnostics_mod

!-----------------------------------------------------------------------

use types_mod, only : r8
use       bgrid_horiz_mod, only: horiz_grid_type
use        bgrid_vert_mod, only: vert_grid_type, compute_pres_full,  &
                                 compute_pres_half, compute_pres_depth
use       bgrid_masks_mod, only: grid_mask_type
use    bgrid_prog_var_mod, only: prog_var_type
use bgrid_change_grid_mod, only: mass_to_vel

!use      diag_manager_mod, only: diag_axis_init, register_diag_field, &
!                                 register_static_field, send_data
use      time_manager_mod, only: time_type
use            fms_mod, only: file_exist, open_namelist_file,    &
                              error_mesg, NOTE, &
                              stdlog,       &
                              close_file, write_version_number
use      constants_mod, only: GRAV
use  field_manager_mod, only: MODEL_ATMOS
use tracer_manager_mod, only: get_tracer_names, get_number_tracers



implicit none
private

public :: bgrid_diagnostics_init, &
          bgrid_diagnostics,      &
          bgrid_diagnostics_tend

!-----------------------------------------------------------------------
!------------------------- axis names ----------------------------------

character(len=8) :: axiset = 'dynamics'
character(len=8) :: mod_name = 'dynamics'

!-----------------------------------------------------------------------

   integer, parameter :: mxtr = 10
   real(r8),    parameter :: ginv = 1./GRAV

!-----------------------------------------------------------------------

integer :: id_hlonb, id_hlon , id_hlatb, id_hlat , &
           id_vlonb, id_vlon , id_vlatb, id_vlat , &
           id_phalf, id_pfull, id_hlat_wgt, id_vlat_wgt

integer :: id_bk   , id_pk   , id_zsurf, id_res  , id_wspd,          &
           id_ps   , id_ucomp, id_vcomp, id_temp ,                   &
           id_omega, id_div  , id_vor  , id_pgfx , id_pgfy,          &
           id_udt  , id_vdt  , id_tdt
integer, allocatable :: id_tracer(:), id_tracer_tend(:)

integer :: id_ucomp_sq, id_vcomp_sq, id_temp_sq, id_omega_sq, &
           id_ucomp_vcomp, id_omega_temp

!-----------------------------------------------------------------------

 character(len=128) :: version = '$Revision$'
 character(len=128) :: tag = '$Id$'

!-----------------------------------------------------------------------

contains

!#######################################################################

subroutine bgrid_diagnostics_init ( Time, Hgrid, Vgrid, Var, &
                                    fis, res,                &
                                    mass_axes, vel_axes      )

   type(time_type),       intent(in)  :: Time
   type(horiz_grid_type), intent(in)  :: Hgrid
   type (vert_grid_type), intent(in)  :: Vgrid
   type  (prog_var_type), intent(in)  :: Var
   real(r8), intent(in), dimension(Hgrid%ilb:,Hgrid%jlb:)  :: fis, res
   integer, dimension(4), intent(out) :: mass_axes, vel_axes

!-----------------------------------------------------------------------
!           subroutine to set up netcdf info for this run
!-----------------------------------------------------------------------
      real(r8), dimension(Hgrid%Tmp%isg:Hgrid%Tmp%ieg+1) :: hlonb
      real(r8), dimension(Hgrid%Vel%isg:Hgrid%Vel%ieg+1) :: vlonb
      real(r8), dimension(Hgrid%Tmp%jsg:Hgrid%Tmp%jeg+1) :: hlatb
      real(r8), dimension(Hgrid%Vel%jsg:Hgrid%Vel%jeg+1) :: vlatb

      real(r8), dimension(Hgrid%Tmp%isg:Hgrid%Tmp%ieg) :: hlon
      real(r8), dimension(Hgrid%Vel%isg:Hgrid%Vel%ieg) :: vlon
      real(r8), dimension(Hgrid%Tmp%jsg:Hgrid%Tmp%jeg) :: hlat
      real(r8), dimension(Hgrid%Vel%jsg:Hgrid%Vel%jeg) :: vlat

      real(r8), dimension(Hgrid%Tmp%js:Hgrid%Tmp%je) :: hlat_wgt
      real(r8), dimension(Hgrid%Vel%js:Hgrid%Vel%je) :: vlat_wgt

      real(r8), dimension(1,1)              :: psurf
      real(r8), dimension(1,1,Vgrid%nlev)   :: pfull
      real(r8), dimension(1,1,Vgrid%nlev+1) :: phalf


      real(r8)    :: vrange(2), trange(2)
      real(r8)    :: rad2deg
      integer :: i, j, n, unit, io, ierr, ntprog
      integer :: isg, ieg, hsg, heg, vsg, veg
      integer :: is, ie, hs, he, vs, ve
      logical :: used
      character(len=128) :: tname
      character(len=256) :: longname, units

!-----------------------------------------------------------------------
!      ----- write version (to log file) -----

    call write_version_number (version,tag)

!--------------------------- set up axes -------------------------------

      is = Hgrid % Tmp % is;  ie = Hgrid % Tmp % ie
      hs = Hgrid % Tmp % js;  he = Hgrid % Tmp % je
      vs = Hgrid % Vel % js;  ve = Hgrid % Vel % je

      isg = Hgrid % Tmp % isg;  ieg = Hgrid % Tmp % ieg
      hsg = Hgrid % Tmp % jsg;  heg = Hgrid % Tmp % jeg
      vsg = Hgrid % Vel % jsg;  veg = Hgrid % Vel % jeg

      rad2deg = 90./acos(0.0)

      hlonb(isg:ieg+1) = Hgrid % Tmp % blong(isg:ieg+1) * rad2deg
      hlatb(hsg:heg+1) = Hgrid % Tmp % blatg(hsg:heg+1) * rad2deg
      vlonb(isg:ieg+1) = Hgrid % Vel % blong(isg:ieg+1) * rad2deg
      vlatb(vsg:veg+1) = Hgrid % Vel % blatg(vsg:veg+1) * rad2deg

   do i = isg, ieg
      hlon(i) = 0.5*(hlonb(i)+hlonb(i+1))
      vlon(i) = 0.5*(vlonb(i)+vlonb(i+1))
   enddo

   do j = hsg, heg
      hlat(j) = 0.5*(hlatb(j)+hlatb(j+1))
   enddo

   do j = vsg, veg
      vlat(j) = 0.5*(vlatb(j)+vlatb(j+1))
   enddo

!---------- reference profile -----------

    psurf = reshape ( (/ 100000. /), (/ 1, 1 /) )
    call compute_pres_full (Vgrid, psurf, pfull)
    call compute_pres_half (Vgrid, psurf, phalf)
! --- in units of hPa ---
    pfull = pfull*0.01
    phalf = phalf*0.01


!----- mass axes ------

! id_hlonb = diag_axis_init ( 'lonb', hlonb, 'degrees_E', 'x',     &
!                             'longitude edges', set_name='atmos', &
!                             Domain2=Hgrid%Tmp%Domain_nohalo      )

! id_hlon  = diag_axis_init ( 'lon', hlon, 'degrees_E', 'x',       &
!                             'longitude', set_name='atmos',       &
!                             edges=id_hlonb,                      &
!                             Domain2=Hgrid%Tmp%Domain_nohalo      )

! id_hlatb = diag_axis_init ( 'latb', hlatb, 'degrees_N', 'y',     &
!                             'latitude edges', set_name='atmos',  &
!                             Domain2=Hgrid%Tmp%Domain_nohalo      )

! id_hlat  = diag_axis_init ( 'lat', hlat, 'degrees_N', 'y',       &
!                             'latitude', set_name='atmos',        &
!                             edges=id_hlatb,                      &
!                             Domain2=Hgrid%Tmp%Domain_nohalo      )

!----- velocity axes ------

! id_vlonb = diag_axis_init ( 'vlonb', vlonb, 'degrees_E', 'x',    &
!                             'longitude edges', set_name='atmos', &
!                             Domain2=Hgrid%Vel%Domain_nohalo      )

! id_vlon  = diag_axis_init ( 'vlon', vlon, 'degrees_E', 'x',      &
!                             'longitude', set_name='atmos',       &
!                             edges=id_vlonb,                      &
!                             Domain2=Hgrid%Vel%Domain_nohalo      )

! id_vlatb = diag_axis_init ( 'vlatb', vlatb, 'degrees_N', 'y',    &
!                             'latitude edges', set_name='atmos',  &
!                             Domain2=Hgrid%Vel%Domain_nohalo      )

! id_vlat  = diag_axis_init ( 'vlat', vlat, 'degrees_N', 'y',      &
!                             'latitude', set_name='atmos',        &
!                             edges=id_vlatb,                      &
!                             Domain2=Hgrid%Vel%Domain_nohalo      )

!----- vertical axes -----

! id_phalf = diag_axis_init ( 'phalf', phalf(1,1,:), 'hPa', 'z', &
!                             'approx half pressure level',      &
!                             direction=-1, set_name='atmos'     )

! id_pfull = diag_axis_init ( 'pfull', pfull(1,1,:), 'hPa', 'z', &
!                             'approx full pressure level',      &
!                             direction=-1, edges=id_phalf,      &
!                             set_name='atmos'                   )

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-------- initialize and output variables with no time axis ------------


    mass_axes = (/ id_hlon, id_hlat, id_pfull, id_phalf /)
     vel_axes = (/ id_vlon, id_vlat, id_pfull, id_phalf /)

     vrange = (/ -400., +400. /)
     trange = (/  100.,  400. /)

!-----------------------------------------------------------------------
!---- register static fields -------

!  id_bk    = register_static_field ( mod_name, 'bk', (/id_phalf/), &
!                        'vertical coordinate sigma value', 'none' )

!  id_pk    = register_static_field ( mod_name, 'pk', (/id_phalf/), &
!   'vertical coordinate reference pressure value (ak*pref)', 'pascals' )

!  id_zsurf = register_static_field ( mod_name, 'zsurf', mass_axes(1:2),&
!                                       'surface height', 'm' )

!  id_res   = register_static_field ( mod_name, 'res', mass_axes(1:2), &
!                     'reciprocal of sigma/eta at the surface', 'none' )

! these changes cannot be implemented until changes to diag_manager
!
! id_hlat_wgt = register_static_field ( mod_name, 'lat_wgt',  &
!                 (/id_hlat/), 'latitude weight for mass grid', 'none' )
!
! id_vlat_wgt = register_static_field ( mod_name, 'vlat_wgt',  &
!             (/id_vlat/), 'latitude weight for momentum grid', 'none' )
                    

!      if ( id_bk > 0 ) &
!      used = send_data ( id_bk, Vgrid%eta, Time )

!      if ( id_pk > 0 ) &
!      used = send_data ( id_pk, Vgrid%peta, Time )

!      if ( id_zsurf > 0 ) &
!      used = send_data ( id_zsurf, fis(is:ie,hs:he)*ginv, Time )

!      if ( id_res > 0 ) &
!      used = send_data ( id_res, res(is:ie,hs:he), Time )

!     if ( id_hlat_wgt > 0 ) then
!        hlat_wgt = sin(Hgrid%Tmp%blatg(hs+1:he+1))-sin(Hgrid%Tmp%blatg(hs:he))
!        used = send_data ( id_hlat_wgt, hlat_wgt, Time )
!     endif
!
!     if ( id_vlat_wgt > 0 ) then
!        vlat_wgt = sin(Hgrid%Vel%blatg(vs+1:ve+1))-sin(Hgrid%Vel%blatg(vs:ve))
!        used = send_data ( id_vlat_wgt, vlat_wgt, Time )
!     endif

!---- register non-static fields -------

!   id_ps   = register_diag_field ( mod_name, 'ps', mass_axes(1:2), &
!                               Time, 'surface pressure', 'pascals' )

!   id_ucomp = register_diag_field ( mod_name, 'ucomp', vel_axes(1:3), &
!                           Time, 'zonal wind component', 'm/sec',     &
!                           missing_value=vrange(1), range=vrange      )

!   id_vcomp = register_diag_field ( mod_name, 'vcomp', vel_axes(1:3), &
!                        Time, 'meridional wind component', 'm/sec',   &
!                        missing_value=vrange(1), range=vrange         )

!   id_temp = register_diag_field ( mod_name, 'temp', mass_axes(1:3), &
!                             Time, 'temperature', 'deg_k',           &
!                             missing_value=trange(1), range=trange   )

 ! diagnostics for all tracers
   allocate (id_tracer(Var%ntrace))
   write(stdlog(),100) trim(mod_name)
   do n = 1, Var%ntrace
     call get_tracer_names ( MODEL_ATMOS, n, tname, longname, units )
     write(stdlog(),110) trim(tname),trim(longname),trim(units)
!     id_tracer(n) = register_diag_field ( mod_name, trim(tname),  &
!                            mass_axes(1:3), Time, trim(longname), &
!                            trim(units), missing_value=-999.      )
   enddo
100 format ('Diagnostics for the following tracer fields are available for module name = ',a)
110 format (3x,a,' (',a,'; ',a,')')

!   id_omega = register_diag_field ( mod_name, 'omega', mass_axes(1:3),&
!                                 Time, 'omega vertical velocity',     &
!                                 'pascals/sec',                       &
!                                 missing_value=-999.                  )

!-------- register second-moment quantities -------

!   id_ucomp_sq = register_diag_field ( mod_name, 'ucomp_sq', vel_axes(1:3), &
!                           Time, 'zonal wind component squared', 'm2/s2',   &
!                           missing_value=-1., range=(/0.,vrange(2)**2/)   )

!   id_vcomp_sq = register_diag_field ( mod_name, 'vcomp_sq', vel_axes(1:3), &
!                      Time, 'meridional wind component squared', 'm2/s2',   &
!                       missing_value=-1., range=(/0.,vrange(2)**2/)   )

!   id_temp_sq = register_diag_field ( mod_name, 'temp_sq', mass_axes(1:3), &
!                             Time, 'temperature squared', 'deg_K**2',      &
!                             missing_value=-1., range=(/0.,trange(2)**2/)  )

!   id_omega_sq = register_diag_field ( mod_name, 'omega_sq', mass_axes(1:3),&
!                                 Time, 'omega vertical velocity squared',   &
!                                 'Pa**2/s**2', missing_value=-999.          )

!   id_ucomp_vcomp = register_diag_field ( mod_name, 'ucomp_vcomp', vel_axes(1:3),&
!                       Time, 'zonal times meridional wind components', 'm2/s2',  &
!                       missing_value=-1. )

!   id_omega_temp = register_diag_field ( mod_name, 'omega_temp', mass_axes(1:3),&
!                               Time, 'omega vertical velocity time temperature',&
!                                 'Pascals*deg_K/sec', missing_value=-999.       )

!-------- wind speed, divergence, vorticity ----------------------------

!   id_wspd = register_diag_field ( mod_name, 'wspd', vel_axes(1:3),   &
!                       Time, 'wind speed', 'm/s', missing_value=-999.,&
!                       range=(/0.,vrange(2)/) )

!   id_div  = register_diag_field ( mod_name, 'div', mass_axes(1:3),   &
!                       Time, 'divergence', '1/s', missing_value=-999. )

!   id_vor  = register_diag_field ( mod_name, 'vor', mass_axes(1:3),   &
!               Time, 'relative vorticity', '1/s', missing_value=-999. )

!--------------- pressure gradient components --------------------------

!  id_pgfx = register_diag_field ( mod_name, 'pgfx', vel_axes(1:3), &
!                            Time, 'zonal pressure gradient force', &
!                            'm/s2', missing_value=-999.            )

!  id_pgfy = register_diag_field ( mod_name, 'pgfy', vel_axes(1:3), &
!                       Time, 'meridional pressure gradient force', &
!                       'm/s2', missing_value=-999.                 )

!-----------------------------------------------------------------------
!         -------- tendencies ---------

!   id_udt = register_diag_field ( mod_name, 'udt_dyn', vel_axes(1:3),  &
!                            Time, 'zonal wind tendency for dynamics', &
!                            'm/s2', missing_value=-999. )

!   id_vdt = register_diag_field ( mod_name, 'vdt_dyn', vel_axes(1:3),  &
!                        Time, 'meridional wind tendency for dynamics', &
!                        'm/s2', missing_value=-999. )

!   id_tdt = register_diag_field ( mod_name, 'tdt_dyn', mass_axes(1:3), &
!                            Time, 'temperature tendency for dynamics', &
!                            'deg_k/sec', missing_value=-999. )

 ! need to do this for just prognostic variables
   call get_number_tracers ( MODEL_ATMOS, num_prog=ntprog )
   allocate (id_tracer_tend(ntprog))
   do n = 1, ntprog
     call get_tracer_names ( MODEL_ATMOS, n, tname, longname, units )
     tname    = trim(tname)   //'_dt_dyn'
     longname = trim(longname)//' tendency for dynamics'
     units    = trim(units)   //'/s'
     if (units == 'none') units = '1/sec'
     write(stdlog(),110) trim(tname),trim(longname),trim(units)
!     id_tracer_tend(n) = register_diag_field ( mod_name, trim(tname),  &
!                                 mass_axes(1:3), Time, trim(longname), &
!                                 trim(units), missing_value=-999.      )
   enddo

!-----------------------------------------------------------------------

   end subroutine bgrid_diagnostics_init

!#######################################################################

 subroutine bgrid_diagnostics ( Hgrid, Vgrid, Var, Masks, Time, omega )

!-----------------------------------------------------------------------
type(horiz_grid_type), intent(in) :: Hgrid
type (vert_grid_type), intent(in) :: Vgrid
type (prog_var_type),  intent(in) :: Var
type(grid_mask_type),  intent(in) :: Masks
type(time_type),       intent(in) :: Time
   real(r8), intent(in), dimension(Hgrid%ilb:,Hgrid%jlb:,:) :: omega

!  real(r8), intent(in), dimension(Hgrid%ilb:,Hgrid%jlb:,:), optional ::&
!                                               div, pgfx, pgfy
                                                       
   real(r8), dimension(Hgrid%ilb:Hgrid%iub,Hgrid%jlb:Hgrid%jub, &
                   Vgrid%nlev) :: wspd, vor, div, adp, udp, vdp
!-----------------------------------------------------------------------
   integer :: is, ie, hs, he, vs, ve, n, k
   logical :: used
!-----------------------------------------------------------------------

      is = Hgrid % Tmp % is;  ie = Hgrid % Tmp % ie
      hs = Hgrid % Tmp % js;  he = Hgrid % Tmp % je
      vs = Hgrid % Vel % js;  ve = Hgrid % Vel % je

!-----------------------------------------------------------------------
!---------------- surface fields ---------------------------------------

!      if ( id_ps > 0 ) &
!      used = send_data ( id_ps , Var%ps(is:ie,hs:he), Time )

!---------------- 3d momentum fields (u & v) ---------------------------

!      if ( id_ucomp > 0 ) &
!      used = send_data ( id_ucomp, Var%u(is:ie,vs:ve,:), Time,    &
!                         mask=Masks%Vel%mask(is:ie,vs:ve,:) > 0.5 )

!      if ( id_vcomp > 0 ) &
!      used = send_data ( id_vcomp, Var%v(is:ie,vs:ve,:), Time,    &
!                         mask=Masks%Vel%mask(is:ie,vs:ve,:) > 0.5 )

!      if ( id_temp > 0 ) &
!      used = send_data ( id_temp, Var%t(is:ie,hs:he,:), Time,     &
!                         mask=Masks%Tmp%mask(is:ie,hs:he,:) > 0.5 )

!      do n = 1, Var%ntrace
!        if ( id_tracer(n) > 0 ) &
!        used = send_data ( id_tracer(n), Var%r(is:ie,hs:he,:,n), Time, &
!                           mask=Masks%Tmp%mask(is:ie,hs:he,:) > 0.5    )
!      enddo

!      if ( id_omega > 0 ) &
!      used = send_data ( id_omega, omega(is:ie,hs:he,:), Time,    &
!                         mask=Masks%Tmp%mask(is:ie,hs:he,:) > 0.5 )

!--------- second moment quantities ----------

!      if ( id_ucomp_sq > 0 ) &
!      used = send_data ( id_ucomp_sq, Var%u(is:ie,vs:ve,:)**2, Time, &
!                         mask=Masks%Vel%mask(is:ie,vs:ve,:) > 0.5 )

!      if ( id_vcomp_sq > 0 ) &
!      used = send_data ( id_vcomp_sq, Var%v(is:ie,vs:ve,:)**2, Time, &
!                         mask=Masks%Vel%mask(is:ie,vs:ve,:) > 0.5 )

!      if ( id_temp_sq > 0 ) &
!      used = send_data ( id_temp_sq, Var%t(is:ie,hs:he,:)**2, Time, &
!                         mask=Masks%Tmp%mask(is:ie,hs:he,:) > 0.5 )

!      if ( id_omega_sq > 0 ) &
!      used = send_data ( id_omega_sq, omega(is:ie,hs:he,:)**2, Time, &
!                         mask=Masks%Tmp%mask(is:ie,hs:he,:) > 0.5 )

!      if ( id_ucomp_vcomp > 0 ) used = send_data ( id_ucomp_vcomp, &
!                  Var%u(is:ie,vs:ve,:)*Var%v(is:ie,vs:ve,:), Time, &
!                         mask=Masks%Vel%mask(is:ie,vs:ve,:) > 0.5 )

!      if ( id_omega_temp > 0 ) used = send_data ( id_omega_temp, &
!                omega(is:ie,hs:he,:)*Var%t(is:ie,hs:he,:), Time, &
!                        mask=Masks%Tmp%mask(is:ie,hs:he,:) > 0.5 )

!------------ wind speed, divergence, vorticity ------------------------

      if ( id_wspd > 0 ) then
          wspd(is:ie,vs:ve,:) = sqrt &
                      ( Var%u(is:ie,vs:ve,:)*Var%u(is:ie,vs:ve,:) + &
                        Var%v(is:ie,vs:ve,:)*Var%v(is:ie,vs:ve,:) )
!          used = send_data ( id_wspd, wspd(is:ie,vs:ve,:), Time,      &
!                             mask=Masks%Vel%mask(is:ie,vs:ve,:) > 0.5 )
      endif

      if ( id_vor > 0 .or. id_div > 0 ) then
!     --- precompute quantities common to both vor and div ---
         call compute_pres_depth (Vgrid, Var%pssl, adp)
         call mass_to_vel (Hgrid, adp, udp)
         vdp = Var%v * udp
         udp = Var%u * udp
         do k=1,size(adp,3)
           adp(:,:,k) = 2.0 * Hgrid%Tmp%area * adp(:,:,k)
         enddo
         if ( id_vor > 0 ) then
             call compute_vorticity (Hgrid, adp, udp, vdp, vor )
!             used = send_data ( id_vor, vor(is:ie,hs:he,:), Time,      &
!                              mask=Masks%Tmp%mask(is:ie,hs:he,:) > 0.5 )
         endif
         if ( id_div > 0 ) then
             call compute_divergence (Hgrid, adp, udp, vdp, div )
!             used = send_data ( id_div, div(is:ie,hs:he,:), Time,      &
!                              mask=Masks%Tmp%mask(is:ie,hs:he,:) > 0.5 )
         endif
      endif

!--------------- optional arguments ------------------------------------
!--------------- pressure gradient components --------------------------

!     if ( id_div > 0 .and. present(div) )  &
!     used = send_data ( id_div, div(is:ie,hs:he,:), Time,        &
!                        mask=Masks%Tmp%mask(is:ie,hs:he,:) > 0.5 )

!     if ( id_pgfx > 0 .and. present(pgfx) )  &
!     used = send_data ( id_pgfx, pgfx(is:ie,vs:ve,:), Time,      &
!                        mask=Masks%Vel%mask(is:ie,vs:ve,:) > 0.5 )

!     if ( id_pgfy > 0 .and. present(pgfy) )  &
!     used = send_data ( id_pgfy, pgfy(is:ie,vs:ve,:), Time,      &
!                        mask=Masks%Vel%mask(is:ie,vs:ve,:) > 0.5 )

!-----------------------------------------------------------------------

 end subroutine bgrid_diagnostics

!#######################################################################

 subroutine bgrid_diagnostics_tend ( Hgrid, Var_dt, Masks, Time )

!-----------------------------------------------------------------------
type(horiz_grid_type), intent(in) :: Hgrid
type (prog_var_type),  intent(in) :: Var_dt
type(grid_mask_type),  intent(in) :: Masks
type(time_type),       intent(in) :: Time
!-----------------------------------------------------------------------
   integer :: is, ie, hs, he, vs, ve, n
   logical :: used
!-----------------------------------------------------------------------

      is = Hgrid % Tmp % is;  ie = Hgrid % Tmp % ie
      hs = Hgrid % Tmp % js;  he = Hgrid % Tmp % je
      vs = Hgrid % Vel % js;  ve = Hgrid % Vel % je

!-----------------------------------------------------------------------
!---------------- 3d prognostic fields ---------------------------

!      if ( id_udt > 0 ) &
!      used = send_data ( id_udt, Var_dt%u(is:ie,vs:ve,:), Time,    &
!                          mask=Masks%Vel%mask(is:ie,vs:ve,:) > 0.5 )

!      if ( id_vdt > 0 ) &
!      used = send_data ( id_vdt, Var_dt%v(is:ie,vs:ve,:), Time,   &
!                         mask=Masks%Vel%mask(is:ie,vs:ve,:) > 0.5 )

!      if ( id_tdt > 0 ) &
!      used = send_data ( id_tdt, Var_dt%t(is:ie,hs:he,:), Time,   &
!                         mask=Masks%Tmp%mask(is:ie,hs:he,:) > 0.5 )

!      do n = 1, Var_dt%ntrace
!       if ( id_tracer_tend(n) > 0 ) &
!       used = send_data ( id_tracer_tend(n), Var_dt%r(is:ie,hs:he,:,n),  &
!                          Time, mask=Masks%Tmp%mask(is:ie,hs:he,:) > 0.5 )
!      enddo

!-----------------------------------------------------------------------

 end subroutine bgrid_diagnostics_tend

!#######################################################################

 subroutine compute_vorticity ( Hgrid, adp, udp, vdp, vor )

!-----------------------------------------------------------------------
type(horiz_grid_type), intent(in) :: Hgrid
real(r8), intent(in),  dimension(Hgrid%ilb:,Hgrid%jlb:,:) :: adp, udp, vdp
real(r8), intent(out), dimension(Hgrid%ilb:,Hgrid%jlb:,:) :: vor

real(r8),dimension(Hgrid%ilb:Hgrid%iub,Hgrid%jlb:Hgrid%jub) :: &
                          vdy,udx,few,fns
integer :: i,j,k,is,ie,js,je

   is = Hgrid % Tmp % is;  ie = Hgrid % Tmp % ie
   js = Hgrid % Tmp % js;  je = Hgrid % Tmp % je

   do k = 1, size(adp,3)
      vdy(:,:) = vdp(:,:,k)*Hgrid%Vel%dy
      udx(:,:) = udp(:,:,k)*Hgrid%Vel%dx(:,:)

      do j = js,   je
      do i = is-1, ie
         fns(i,j) = vdy(i,j-1)+vdy(i,j)
      enddo
      enddo
      do j = js-1, je
      do i = is,   ie
         few(i,j) = udx(i-1,j)+udx(i,j)
      enddo
      enddo

!  ------ vorticity ------
      do j = js, je
      do i = is, ie
         vor(i,j,k)=((fns(i,j)-fns(i-1,j))-(few(i,j)-few(i,j-1))) &
                    /adp(i,j,k)
      enddo
      enddo
   enddo

 end subroutine compute_vorticity

!#######################################################################

 subroutine compute_divergence ( Hgrid, adp, udp, vdp, div )

!-----------------------------------------------------------------------
type(horiz_grid_type), intent(in) :: Hgrid
real(r8), intent(in),  dimension(Hgrid%ilb:,Hgrid%jlb:,:) :: adp, udp, vdp
real(r8), intent(out), dimension(Hgrid%ilb:,Hgrid%jlb:,:) :: div

real(r8),dimension(Hgrid%ilb:Hgrid%iub,Hgrid%jlb:Hgrid%jub) ::  &
                          udy,vdx,few,fns
integer :: i,j,k,is,ie,js,je

   is = Hgrid % Tmp % is;  ie = Hgrid % Tmp % ie
   js = Hgrid % Tmp % js;  je = Hgrid % Tmp % je

   do k = 1, size(adp,3)
      udy(:,:) = udp(:,:,k)*Hgrid%Vel%dy
      vdx(:,:) = vdp(:,:,k)*Hgrid%Vel%dx(:,:)

      do j = js,   je
      do i = is-1, ie
         few(i,j) = udy(i,j-1)+udy(i,j)
      enddo
      enddo
      do j = js-1, je
      do i = is,   ie
         fns(i,j) = vdx(i-1,j)+vdx(i,j)
      enddo
      enddo

!  ------ divergence ------
      do j = js, je
      do i = is, ie
         div(i,j,k)=((few(i,j)+fns(i,j))-(few(i-1,j)+fns(i,j-1))) &
                    /adp(i,j,k)
      enddo
      enddo
   enddo

 end subroutine compute_divergence

!#######################################################################

end module bgrid_diagnostics_mod

