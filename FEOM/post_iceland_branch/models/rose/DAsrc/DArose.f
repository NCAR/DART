! Data Assimilation Research Testbed -- DART
! Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

!-------------------------------------------------------------------
!
!
!             ......................................
!             :                                    :
!             :                                    :
!             :              R O S E               :
!             :                                    :
!             :      Research on Ozone in the      :
!             :   Stratosphere and its Evolution   :
!             :                                    :
!             :                                    :
!             :       3-D Mechanistic Model        :
!             :      of the Middle atmosphere      :
!             :       with coupled chemistry       :
!             :                                    :
!             :                                    :
!             :                NCAR                :
!             :            P.O Box 3000            : 
!             :         Boulder, CO  80307         :
!             :                                    :
!             :  Major Contributors:               :
!             :                                    :
!             :  G. Brasseur (brasseur@dkrz.de)    :
!             :  A. Smith (aksmith@ucar.edu)       :
!             :  S. Walters (stacy@ucar.edu)       :  
!             :  D. Marsh (marsh@ucar.edu)         :
!             :....................................:
!
!
!
!-------------------------------------------------------------------

      program  mlt3d

! <next four lines automatically updated by CVS, do not edit>
! $Source$ 
! $Revision$ 
! $Date$ 
! $Author$ 

!===================================================================
!    3d model with annual cycle
!    This version has: 
!           chemistry for stratosphere and mesosphere
!           Hines gravity wave drag
!           polar ozone hole chemistry omitted
!===================================================================
!
!    Modified from a3d.chem.f by T.M.
!
!===================================================================

      use params 
      use chem
      use dynam
      use phys
      use output
      use diagnostic

      use utilities_mod, only : open_file, close_file, 
     $                          error_handler, E_ERR, E_MSG, E_WARN
      use DArose_mod, only : msetvar, output_prog_diag

      implicit none

!-------------------------------------------------------------------
!   local variables
!-------------------------------------------------------------------

      integer :: doy             ! day of year
      integer :: utsec           ! universal time of model in seconds
      real    :: modday          ! number of days model has run
      integer :: iyear           ! number of years model has run

      real    :: gmt_frac        ! fraction of day since 0Z
      integer :: lengdy          ! length of a day in steps

      integer :: iunit, iunit_tidal
      integer :: i, k
      integer :: dummy(3)
      integer :: daynum
      integer :: ierr, iftype
      character(len=129) :: errstring

      real, dimension (nz,nx,ny,nvars) :: vars

      real, dimension(nz,nx,ny) :: tfull, fxt
      
      real, parameter :: gcp = 9.8e-03

!-------------------------------------------------------------------
! CVS Generated file description for error handling, do not edit
      character(len=128) :: source = " "
      character(len=128) :: revision = " "
      character(len=128) :: revdate = " "
!-------------------------------------------------------------------
!  set variable program parameters including program control
!  parameters and i/o file names
!-------------------------------------------------------------------

      print *, 'msetvar'
      call msetvar()
      lengdy = ntime*24

!-------------------------------------------------------------------
!  set geometric parameters and arrays
!-------------------------------------------------------------------

      print *, 'msetfix'
      call msetfix

!-------------------------------------------------------------------
!   read in initial dynamical and chemical fields
!-------------------------------------------------------------------

      print *, 'startup: ', namf17

      iunit = open_file(namf17,form='unformatted',action='read')
      read(iunit, iostat=ierr) iftype
      if ( ierr /= 0 ) then
      write(errstring,*)'restart file is ',trim(adjustl(namf17))
      call error_handler(E_MSG,'rose main: ',
     $     errstring,source,revision,revdate)
      call error_handler(E_ERR,'rose main: ',
     $     'unknown restart file type ',source,revision,revdate)
      endif

!     open (unit   = 17, 
!    $         file   = namf17,
!    $         form   = 'unformatted',
!    $         status = 'old')

      if (iftype == 1) then !! old_restart == .true.

        read(iunit,iostat=ierr) dummy, gmt_frac, daynum, tref, treflb,
     $              trefub, un1, vn1, tn1, un0, vn0, tn0, 
     $              qn1(:,:,:,1:nbcon-1), q_o2, q_n2
        if ( ierr /= 0 ) then
          write(errstring,*)'restart file is ',trim(adjustl(namf17))
          call error_handler(E_MSG,'rose main: ',
     $         errstring,source,revision,revdate)
          call error_handler(E_ERR,'rose main: ',
     $         'read error ',source,revision,revdate)
        endif

        day0 = dummy(3)
        print *, 'day0:', day0, 'gmt_frac:', gmt_frac
        ut0 = int(gmt_frac * 24.0 * 3600.)
        iyear = (day0 - 1)/365
        doy = mod( day0, 365)

        nstart = 0

      else !! old_restart == .false.

        read(iunit,iostat=ierr) iyear,doy,utsec,year0,day0,ut0,tref, 
     $              treflb, trefub, un1, vn1, tn1, un0, vn0, tn0, 
     $              qn1, q_o2, q_n2
        if ( ierr /= 0 ) then
          write(errstring,*)'restart file is ',trim(adjustl(namf17))
          call error_handler(E_MSG,'rose main: ',
     $         errstring,source,revision,revdate)
          call error_handler(E_ERR,'rose main: ',
     $         'read error ',source,revision,revdate)
        endif

        ut0 = utsec
        gmt_frac = float(utsec)/(3600.0*24.0)
        year0 = iyear
        iyear = 0
        day0 = doy
        print *, iyear, doy, utsec, year0, day0, ut0

        nstart = 1

      endif 

      call close_file(iunit) 

!-------------------------------------------------------------------
!   other initial fields
!-------------------------------------------------------------------

!... define the reference static stability

      do k=2,nz-1
        ssref(k) = (tref(k+1) - tref(k-1))/(2.*dz) + gcp
      end do
      ssref(1) = (tref(2) - treflb)/(2.*dz) + gcp
      ssref(nz) = (trefub - tref(nz-1))/(2.*dz) + gcp

!... full temperature

      do k=1,nz
        tfull(k,:,:) = tn1(k,:,:) + tref(k)
      end do

!... if not restart run initialize O2

      if (nstart.eq.0) then
        do i=1,nx
          qn1(:,i,:,26) = q_o2(:,:)
        end do
      end if

!... calculate K for molecular diffusion of heat from reference temperature

      print *, 'dmolec'
      call dmolec

      print *, 'moldiff'
      call moldiff

!... initial w calculated from continuity equation

      call wcont

!... initial climatology fields

      call tinterp_dyn (day0)
      call tinterp_rad (day0)


!-------------------------------------------------------------------
!   open files for results
!-------------------------------------------------------------------
      
      if (output_prog_diag) then

      call main_init( namf10 )
      call diag_init( namf41 )

      vars(:,:,:,1) = un1
      vars(:,:,:,2) = vn1
      vars(:,:,:,3) = ww 
      vars(:,:,:,4) = tfull
      vars(:,:,:,5:ndyn+nbcon) = qn1

      call main_update( vars, iyear, day0, ut0 ) 

!-------------------------------------------------------------------
!   open file iunit_tidal:  fields needed for tidal analysis 
!-------------------------------------------------------------------

      iunit_tidal = open_file(namf13,form='unformatted',action='write')

!     open (unit   = 13, 
!    $      file   = namf13,
!    $      form   = 'unformatted',
!    $      status = 'unknown')

      endif

!==================================================================
!    time integration loop
!==================================================================

      do nstep = 1,nend
 
!-------------------------------------------------------------------
!  define time & date variables
!-------------------------------------------------------------------

         modday = float(nstep)/float(lengdy)

         iyear =  (      day0  + modday + float(ut0)/86400.- 1)/365
         doy   = mod(float(day0) + modday + float(ut0)/86400., 365.0)
         if(doy.eq.0) doy = 365
         utsec = mod( nstep*3600/ntime + ut0, 3600*24)
         gmt_frac = float(utsec)/(3600.0*24.0)
  
         print 101, nstep, iyear, doy, utsec 
 101     format('nstep:', i12, 'year:', i8, 'doy', i6, 'utsec', i8)

!-------------------------------------------------------------------
!  chemical boundary conditions
!-------------------------------------------------------------------

         call chembc_intrp (doy, gmt_frac)

!-------------------------------------------------------------------
!  TRANSPORT
!-------------------------------------------------------------------
!   compute semi lagrangian transport (every ntrans timestep)
!   and updated concentration fields
!-------------------------------------------------------------------
         if (ntrans.eq.1 .or. mod(nstep,ntrans).eq.1) then
            print *, 'transport'
            call transport()
         end if

         where( qn1 < 0) 
            qn1 = 0             ! enforce non-negative mixing ratios
         endwhere

!-------------------------------------------------------------------
!   vertical diffusion of chemical species
!-------------------------------------------------------------------
         print *, 'vertdiff_c'
         call vertdiff_c()

         where( qn1 < 0) 
            qn1 = 0             ! enforce non-negative mixing ratios
         endwhere

!-------------------------------------------------------------------
!  CHEMISTRY
!-------------------------------------------------------------------
!   compute chemical tendency and updated concentration fields
!-------------------------------------------------------------------

!-------------------------------------------------------------------
!     calculation of number/cm**3 from ideal gas law
!-------------------------------------------------------------------
         print *, 'density'
         call density (tfull)

!-------------------------------------------------------------------
!     calculation of 3-d column ozone for photolysis calc.
!-------------------------------------------------------------------
         print *, 'total'
         call total

!-------------------------------------------------------------------
!     calculation of solar zenith angle
!-------------------------------------------------------------------
         print *, 'zenith'
         call zenith (doy, gmt_frac)

!-------------------------------------------------------------------
!     calculation of photodissociation frequencies
!-------------------------------------------------------------------
         print *, 'phot'
         call phot(doy)

!-------------------------------------------------------------------
!     calculate rate coefficients
!-------------------------------------------------------------------
         if (nrates.eq.1 .or. mod(nstep,nrates).eq.1) then
            print *, 'rates'
            call rates (tfull)
         end if

!-------------------------------------------------------------------
!     calculate rate coefficients for heterogeneous reactions on
!     aerosols
!-------------------------------------------------------------------
!!!         if(particle) call aerosols

!-------------------------------------------------------------------
!     chemistry program to update long-lived and short-lived species
!-------------------------------------------------------------------
         print *, 'chem3d'
         call chem3d 

         where( qn1 < 0) 
            qn1 = 0             ! enforce non-negative mixing ratios
         endwhere

!-------------------------------------------------------------------
!   correction for mass conservation
!-------------------------------------------------------------------
         print *, 'no conserv' 
         ! call conserv

!-------------------------------------------------------------------
!  DYNAMICS
!-------------------------------------------------------------------
!   dynamical lower boundary conditions, including planetary waves and tides
!-------------------------------------------------------------------
         if (lengdy.eq.1 .or. mod(nstep,lengdy).eq.1) then

            print *, 'tinterp'
            call tinterp_dyn (day0)
            call tinterp_rad (day0)

            call nmc_lbc (doy,iyear)

            call tidalinp (doy)

         end if

         print *, 'dynamlbc'
         call dynamlbc( gmt_frac )

!-------------------------------------------------------------------
!   integrate vertically for fi
!-------------------------------------------------------------------
         print *, 'hydrost'
         call hydrost()

!-------------------------------------------------------------------
!   w from continuity
!-------------------------------------------------------------------
         print *, 'wcont'
         call wcont()

!-------------------------------------------------------------------
!   gravity wave parameterization
!-------------------------------------------------------------------
         print *, 'gravity'
         call hines_driver
         call gravity
         fxt = fgr*un1 - fcgr + falph

!-------------------------------------------------------------------
!   calculate heating (solar and chemical) and cooling rates
!-------------------------------------------------------------------
!   NOTE:  For now, use climatological ozone, water, & CO2 in
!   ir radiative transfer code
!-------------------------------------------------------------------

         if (ndiabat.eq.1 .or. mod(nstep,ndiabat).eq.1) then

            print *, 'radiation' 
            call radiation (doy, gmt_frac, tfull)

         endif

!-------------------------------------------------------------------
!   compute dynamical tendency and calculate updated fields
!-------------------------------------------------------------------
         print *, 'tendency'
         call tendency

!-------------------------------------------------------------------
!   vertical diffusion of potential temperature
!-------------------------------------------------------------------
         print *, 'vertdiff_t'
         call vertdiff_t

!-------------------------------------------------------------------
!   filter for dynamical variables
!-------------------------------------------------------------------
         call filter                            ! for high latitudes
         call shapirx                           ! shapiro filters
         call shapiry

!-------------------------------------------------------------------
!   write standard and diagnostic variables 
!-------------------------------------------------------------------

        if (output_prog_diag) then

         if (nout.eq.1 .or. mod(nstep,nout).eq.1) then

           print *, 'saving model data -- step=', nstep

           do k=1,nz
             tfull(k,:,:) = tn1(k,:,:) + tref(k)
           end do

           vars(:,:,:,1) = un1
           vars(:,:,:,2) = vn1
           vars(:,:,:,3) = ww 
           vars(:,:,:,4) = tfull
           vars(:,:,:,5:ndyn+nbcon) = qn1
               
           call main_update( vars, iyear, doy, utsec ) 

           call diag_update( iyear, doy, utsec )

         endif

!-------------------------------------------------------------------
!   write fields needed for tidal analysis to file unit 13
!   (every 6 timesteps on certain days)
!-------------------------------------------------------------------

            if(mod(doy,nouttid).eq.0)then
               if(mod(nstep,6).eq.0)then
                  write(13) nstep, ntime, un1, vn1, ww, tfull
               endif
            endif
 
         endif

!-------------------------------------------------------------------
!   shift variables by one timestep
!-------------------------------------------------------------------

         un0 = un1
         vn0 = vn1
         tn0 = tn1
         un1 = un2
         vn1 = vn2
         tn1 = tn2
         do k=1,nz
            tfull(k,:,:) = tn1(k,:,:) + tref(k)
         end do

!-------------------------------------------------------------------
!  save fields for restarting
!-------------------------------------------------------------------

         if((nstep.gt.2.and.mod(nstep,nsave).eq.1).or.nstep.eq.nend)then
  
            print *, 'writing restart file for step:', nstep
    
            iunit = open_file(namf14,form='unformatted',action='write')
!           open (unit = 14, file = namf14, form = 'unformatted')
   
            iftype = 2 ! ALWAYS want to write the 'new-style' restart file.
            write(iunit) iftype !  old_restart == .true.
            write(iunit) iyear, doy, utsec, year0, day0, ut0, 
     $                   tref, treflb, trefub, 
     $                   un1, vn1, tn1, un0, vn0, tn0, 
     $                   qn1, q_o2, q_n2

            call close_file(iunit)

       end if

!==================================================================
      end do                            ! end time integration loop
!==================================================================

!-------------------------------------------------------------------
!... close output files
!-------------------------------------------------------------------

      if (output_prog_diag) then 
        call main_close
        call diag_close
        call close_file(iunit_tidal)
      endif

      stop
      end

