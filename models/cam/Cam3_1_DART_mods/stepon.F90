! Data Assimilation Research Testbed -- DART
! Copyright 2004-2007, Data Assimilation Research Section
! University Corporation for Atmospheric Research

#include <misc.h>
#include <params.h>
#if ( defined SCAM )
#include <max.h>
#endif

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

subroutine stepon
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Loop over time, calling driving routines for physics, dynamics, 
! transport
! 
! Method: 
! 
! Author: 
! Original version:  CCM1
! Standardized:      J. Rosinski, June 1992
! Reviewed:          B. Boville, D. Williamson, August 1992
! Reviewed:          B. Boville, D. Williamson, April 1996
! Restructured:      J. Truesdale, May 1999
! Augmented:         K. Raeder, June 2004 (kdr)
!
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use history, only: wshist, wrapup
   use pmgrid,  only: plev, plat, plevp, plon, beglat, endlat, &
                      masterproc
   use scanslt, only: advection_state, scanslt_initial, scanslt_final, qfcst
   use rgrid,   only: nlon
   use prognostics, only: ps, u3, v3, t3, q3, qminus, div, &
                          dpsl, dpsm, omga, phis, n3, n3m2, n3m1
   use prognostics, only: pdeld
   use restart, only: write_restart
   use constituents, only: pcnst
#if (defined COUP_CSM)
   use ccsm_msg, only: csmstop, ccsmfin
#endif

   use ppgrid,         only: begchunk, endchunk
   use physics_types,  only: physics_state, physics_tend
   use phys_buffer,    only: pbuf
   use dp_coupling,    only: d_p_coupling, p_d_coupling
   use commap,         only: clat
   use physconst, only: gravit
   use time_manager, only: advance_timestep, get_step_size, get_nstep, &
                           is_first_step, is_first_restart_step, &
                           is_last_step, is_end_curr_day, get_curr_date, &
! kdr
                           nelapse, nestep, stop_ymd, stop_tod
#if ( defined BFB_CAM_SCAM_IOP )
   use iop, only: init_iop_fields,fixmassav,betasav,alphasav,divt3dsav,divq3dsav,dqfx3savm1
#endif
#if ( defined SCAM )
   use scamMod, only :use_iop,use_pert_frc,doiopupdate,switch
   use iop, only: init_iop_fields,fixmassav,betasav,alphasav,divt3dsav,divq3dsav,dqfx3savm1
#endif
   implicit none

!-----------------------------------------------------------------------
#include <comctl.h>
!-----------------------------------------------------------------------
#include <comhyb.h>
!-----------------------------------------------------------------------
#include <comlun.h>
!-----------------------------------------------------------------------
#if ( defined SCAM )
#include <comfrc.h>
!-----------------------------------------------------------------------
#endif
!
   type(physics_state), allocatable :: phys_state(:)
   type(physics_tend ), allocatable :: phys_tend(:)

   real(r8) gw(plat)                ! Gaussian weights
   real(r8) detam(plev)             ! intervals between vert full levs.
   real(r8) cwava(plat)             ! weight applied to global integrals
   real(r8) etamid(plev)            ! vertical coords at midpoints 
   real(r8), allocatable :: t2(:,:,:) ! temp tendency
   real(r8), allocatable :: fu(:,:,:) ! u wind tendency
   real(r8), allocatable :: fv(:,:,:) ! v wind tendency
   real(r8) flx_net(plon,beglat:endlat)       ! net flux from physics
   real(r8) coslat(plon)
   real(r8) rcoslat(plon)
   real(r8) rpmid(plon,plev)
   real(r8) pdel(plon,plev)
   real(r8) pint(plon,plevp)
   real(r8) pmid(plon,plev)
   real(r8) dtime               ! timestep size
   real(r8) ztodt               ! twice time step unless nstep=0
   real(r8) :: wcstart, wcend   ! wallclock timestamp at start, end of timestep
   real(r8) :: usrstart, usrend ! user timestamp at start, end of timestep
   real(r8) :: sysstart, sysend ! sys timestamp at start, end of timestep
   type(advection_state) :: adv_state ! Advection state data
!
   integer i,k,lat,j            ! longitude,level,latitude indices
   integer jd1,jd2              ! latitude index for pdeld
   integer iter
! kdr
   integer nstep,  ymd2
   real(r8) days
! kdr end
   integer :: yr, mon, day    ! year, month, and day components of a date
   integer :: ncsec           ! current time of day [seconds]

#if ( defined SCAM )
   integer dummy
   real(r8) omegapert
   real(r8) pertval              ! perturbation valie
   real(r8) drand48
   external drand48
#endif
!
! Externals
!
   logical, external :: rstwr  ! whether or not to write restart files
!
!-----------------------------------------------------------------------
   call t_startf ('stepon_startup')
   dtime = get_step_size()
!
! Define eta coordinates: Used for calculation etadot vertical velocity 
! for slt.
!
   do k=1,plev
      etamid(k) = hyam(k) + hybm(k)
   end do

   call scanslt_initial( adv_state, etamid, gravit, gw, detam, cwava )
!
! Initial guess for trajectory midpoints in spherical coords.
! nstep = 0:  use arrival points as initial guess for trajectory midpoints.
! nstep > 0:  use calculated trajectory midpoints from previous time 
! step as first guess.
! NOTE:  reduce number of iters necessary for convergence after nstep = 1.
!
   if (is_first_step()) then
      do lat=beglat,endlat
#if ( !defined SCAM )
         do i=1,nlon(lat)
            coslat(i) = cos(clat(lat))
            rcoslat(i) = 1./coslat(i)
         end do
#endif
!     
! Set current time pressure arrays for model levels etc.
!
         call plevs0(nlon(lat), plon, plev, ps(1,lat,n3), pint, pmid, pdel)
!
         do k=1,plev
            do i=1,nlon(lat)
               rpmid(i,k) = 1./pmid(i,k)
            end do
         end do

#if ( ! defined SCAM )
!
! Calculate vertical motion field
!
         call omcalc (rcoslat, div(1,1,lat,n3), u3(1,1,lat,n3), v3(1,1,lat,n3), dpsl(1,lat), &
                      dpsm(1,lat), pmid, pdel, rpmid   ,pint(1,plevp), &
                      omga(1,1,lat), nlon(lat))
#else
         
         omga(1,:,lat)=wfld(:)
#endif
      end do
   end if

   allocate(phys_state(begchunk:endchunk))
   allocate(phys_tend(begchunk:endchunk))
   allocate(t2(plon,plev,beglat:endlat))
   allocate(fu(plon,plev,beglat:endlat))
   allocate(fv(plon,plev,beglat:endlat))
!
! Beginning of basic time step loop
!
   call t_stopf ('stepon_startup')


#if ( defined BFB_CAM_SCAM_IOP )
   if (is_first_step()) then
      call init_iop_fields(ps(1,beglat,n3m2), t3(1,1,beglat,n3m2), &
        u3(1,1,beglat,n3m2), v3(1,1,beglat,n3m2), q3(1,1,1,beglat,n3m2))
      fixmassav(:)=0.
      betasav(:)=0.
      alphasav(:,:)=0.
      divt3dsav(:,:,:)=0.
      divq3dsav(:,:,:,:)=0.
      dqfx3savm1(:,:,:,:)=0.
   endif
#endif

! Begin time loop.

   do
      call t_startf('stepon_st')
      if (masterproc .and. print_step_cost) then
         call t_stampf (wcstart, usrstart, sysstart)
      end if

      ztodt = 2.0*dtime
!
! If initial time step adjust dt
!
      if (is_first_step()) ztodt = dtime
!
! adjust hydrostatic matrices if the time step has changed.  This only
! happens on transition from time 0 to time 1. 
      if (get_nstep() == 1) then
         call settau(dtime)
      end if
!
! Dump state variables to IC file
!
      call diag_dynvar_ic (phis, ps(1,beglat,n3m1), t3(1,1,beglat,n3m1), u3(1,1,beglat,n3m1), &
                           v3(1,1,beglat,n3m1), q3(1,1,1,beglat,n3m1) )
!
!----------------------------------------------------------
! PHYSPKG  Call the Physics package
!----------------------------------------------------------
!
      call t_stopf('stepon_st')
      call t_startf('d_p_coupling')
      call d_p_coupling (ps(1,beglat,n3m2), t3(1,1,beglat,n3m2), u3(1,1,beglat,n3m2), &
                         v3(1,1,beglat,n3m2), q3(1,1,1,beglat,n3m2), &
                         omga, phis, phys_state, phys_tend, pbuf, pdeld(:,:,:,n3m2))
      call t_stopf('d_p_coupling')

      call t_startf('phys_driver')
      if (ideal_phys) then
         call phys_idealized(phys_state, phys_tend, ztodt, etamid)
      else if (adiabatic) then
         call phys_adiabatic(phys_state, phys_tend)
      else
! kdr land model called here; BEFORE timestep is updated;
!     write of initial file must use different criterion; 1 step "earlier"
         call physpkg(phys_state, gw, ztodt, phys_tend, pbuf)
      end if
      call t_stopf('phys_driver')
   
      call t_startf('p_d_coupling')
      call p_d_coupling (phys_state, phys_tend, t2, fu, fv, flx_net, &
                         qminus(1,1,1,beglat),  q3(1,1,1,beglat,n3))
      call t_stopf('p_d_coupling')
#if (defined SCAM )
!
! In order to provide forcasted values for u,v, and ps from a ccm derived
! iop dataset we must advance timestep and perform iop read here.  Dynamics
! residuals are also read here to allow forcasting of q and t.
 
     call advance_timestep()
     call get_curr_date(yr, mon, day, ncsec)
!
!=====================================================================
!     Determine whether it is time for an IOP update;
!     doiopupdate set to true if model time step > next available IOP
!     dataset time step
!
     if (use_iop) then
        call setiopupdate
     end if
!
!     Update IOP properties e.g. omega, divT, divQ
!
     if (doiopupdate) then
        call readiopdata(dummy) ! not checking error code
!
!  Randomly perturb omega if requested
!         
        omegapert = .30             ! +/- 30 %
!         
        if (use_pert_frc .eqv. .true. ) then
           do k=1,plev
              pertval = omegapert * (2.*(0.5 - drand48()))
              wfld(k) = wfld(k) * (1. + pertval )
           end do
        end if
!
     end if                    ! (do_iopupdate)
#endif
!----------------------------------------------------------
! DYNPKG Call the Dynamics Package
!----------------------------------------------------------

      call t_startf('dynpkg')
#if ( defined SCAM )
      ! No dynamics for the column radiation model
      if(.not.switch(CRM_SW+1)) &
#endif
! kdr 'NSTEP =' line comes from courlim call down in this dynpkg call, I think.
      call dynpkg (adv_state, t2      ,fu      ,fv      ,etamid  ,  &
                   cwava   ,detam   ,flx_net ,ztodt   )
      call t_stopf('dynpkg')

      call t_startf('stepon_st')
      if (is_first_step() .or. is_first_restart_step()) then
         call print_memusage ('stepon after dynpkg')
      end if

! Set end of run flag.

#if ( ! defined COUP_CSM )
      if (is_last_step()) nlend = .true.
#else
      if (csmstop) then
         if ( masterproc ) write(6,*)'atm: Stopping at the end of this day'
         if (is_end_curr_day()) nlend = .true.
      end if
#endif

#if ( ! defined SCAM )
!
!----------------------------------------------------------
! History and restart logic: Write and/or dispose history tapes if required
!----------------------------------------------------------
!
      call t_startf ('wshist')
      call wshist ()
      call t_stopf ('wshist')
!
! Write restart file
!
      if (rstwr() .and. nrefrq /= 0) then
         call t_startf ('write_restart')
         call write_restart
         call t_stopf ('write_restart')
      end if
!
! Dispose necessary files
!
      call t_startf ('wrapup')
      call wrapup
      call t_stopf ('wrapup')

      if (masterproc .and. print_step_cost) then
         call t_stampf (wcend, usrend, sysend)
         write(6,'(a,3f8.3,a)')'Prv timestep wallclock, usr, sys=', &
                               wcend-wcstart, usrend-usrstart, sysend-sysstart, ' seconds'
      end if
!
! Advance timestep before returning to top of loop
!
      call advance_timestep()
      call get_curr_date(yr, mon, day, ncsec)
      
#else                 
! else if running SCAM then end timestep and dealloc
      nlend = .true.
#endif
      call t_stopf('stepon_st')
!
! Check for end of run
!
      if (nlend) then
         call scanslt_final( adv_state )
         call print_memusage ('End stepon')
         deallocate(phys_state)
         deallocate(phys_tend)
         deallocate(t2)
         deallocate(fu)
         deallocate(fv)
#ifdef COUP_CSM
         call ccsmfin
#endif
         return
      end if

   end do  ! End of timestep loop

end subroutine stepon
