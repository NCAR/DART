! Data Assimilation Research Testbed -- DART
! Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module DArose_mod

! <next five lines automatically updated by CVS, do not edit>
! $Source$ 
! $Revision$ 
! $Date$ 
! $Author$
! $Name$

use params, only : ntime, nout, noutdiag, nouttid, nend, nsave
use dynam,  only : deltat, dtleap, ntrans, ndiabat, ninterp, &
                   namf10, namf12, namf13, namf14, namf17, & 
                   namf20, namf23, namf26, namf27, namf28, &
                   namf30, namf40, namf41, namf52, namf53       
use chem,   only : nrates
use chem0d, only : niter, dtchem

use     types_mod, only: r8, pi
use utilities_mod, only: open_file, close_file, file_exist, &
                            register_module, error_handler, E_ERR, E_MSG, logfileunit

implicit none
private

public :: msetvar, output_prog_diag, h_tune, z_tune

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

   logical :: output_prog_diag = .false.
   character (len=50) :: input_dir = '../DAinput/'
   character (len=50) :: out_dir   = '../DAoutput/'
   character (len=30) :: ncep_file = 'nmc_lbc.02.nc'
   character (len=30) :: restart_file = 'NMC_SOC.day151_2002.dat'  
   real(kind=r8) :: h_tune = pi
   real(kind=r8) :: z_tune = 1.0
   real(kind=r8) :: target_time = 168.0 ! 7 days * 24 [hr]  
   integer       :: lengdy 


   namelist /rose_nml/ target_time, &
                       input_dir, out_dir, ncep_file, restart_file,&
                       output_prog_diag, &
                       h_tune, z_tune, &
                       ntime
   contains

subroutine msetvar

   implicit none

   integer :: iunit, io
   character(len=129) :: err_string, nml_string

   ! Print module information to log file and stdout.
   call register_module(source, revision, revdate)

   ntime = 8              ! number of steps per hour

   ! Begin by reading the namelist input
   if(file_exist('rose.nml')) then
      iunit = open_file('rose.nml', action = 'read')
      read(iunit, nml = rose_nml, iostat = io)
      if(io /= 0) then
         ! A non-zero return means a bad entry was found for this namelist
         ! Reread the line into a string and print out a fatal error message.
         BACKSPACE iunit
         read(iunit, '(A)') nml_string
         write(err_string, *) 'INVALID NAMELIST ENTRY: ', trim(adjustl(nml_string))
         call error_handler(E_ERR, 'msetvar:&rose_nml problem', &
                            err_string, source, revision, revdate)
      endif
      call close_file(iunit)
   endif

   ! Record the namelist values used for the run ...
   call error_handler(E_MSG,'msetvar','rose_nml values are',' ',' ',' ')
   write(logfileunit, nml=rose_nml)
   write(     *     , nml=rose_nml)

   !  write(*,*) 'rose_nml values  --- may come from rose.nml'
   !  write(*,*) 'target_time ', target_time
   !  write(*,*) 'input_dir ', input_dir 
   !  write(*,*) 'out_dir ', out_dir
   !  write(*,*) 'ncep_file ', ncep_file
   !  write(*,*) 'restart_file ', restart_file
   !  write(*,*) 'output_prog_diag ', output_prog_diag
   !  write(*,*) 'h_tune ', h_tune
   !  write(*,*) 'z_tune ', z_tune
   !  write(*,*) 'ntime ', ntime

!-------------------------------------------------------------------
!... set program control parameters 
!-------------------------------------------------------------------

      lengdy = ntime*24

      deltat = 3600.0_r8/float(ntime)
      dtleap = deltat*2.0_r8
      dtchem = deltat

!... length of integration

      nend = floor(target_time*real(ntime))

!... frequency at which output is saved

    if (output_prog_diag) then

!     nout = 1                      ! save every step
!     nout = 6                      ! every 45 minutes
!     nout = ntime                  ! save every hour 
      nout = ntime * 12             ! save every 12 hours 
!     nout = ntime *24              ! save every day
!     nout = lengdy * 2             ! or less frequently
      
      noutdiag = nout               ! freq to save diagnostics (typically nout)
      nouttid = 10                  ! frequency IN DAYS for higher 
                                    ! freq output for tides
    endif 

!-------------------------------------------------------------------
!  other model control parameters
!  variable - update for each run   ! need to be added to namelist
!-------------------------------------------------------------------

!      if (nstart == 0) then  
!        nseg = 1       ! segment number (=1 for initial run)
!      else             ! used for tracking NetCDF output files
!        nseg = 2       ! USED WHEN output_prog_diag = .true. 
!      endif

!-------------------------------------------------------------------
!  fixed for all portions of run
!-------------------------------------------------------------------

      ntrans = 1          ! frequency (in timestep) of transport
      ndiabat = 1         ! frequency of diabatic heating/cooling calc
      nrates = 1          ! frequency of rate constant update
      ninterp = lengdy    ! frequency of climatology b.c. update

      nsave = 15*lengdy   ! save for restarting

!-------------------------------------------------------------------
!  for chemistry
!-------------------------------------------------------------------

      niter = 32 

!-------------------------------------------------------------------
!  i/o files
!-------------------------------------------------------------------

!-------------------------------------------------------------------
!... fixed  (DEFAULT)
!-------------------------------------------------------------------

!   unit 12 zonal wind climatology for sponge layer
      namf12 = trim(input_dir)//'init.mlt.data'

!   unit 20 dynamical troposphere and lower boundary fields
!     namf20 = trim(rose_data)//'nmc_lbc.02.nc'
      namf20 = trim(input_dir)//trim(ncep_file)

!   unit 23 climatology zonal mean dynamical ubc
      namf23 = trim(input_dir)//'mlt.ubc.data'

!   unit 26 chemical upper boundary from SOCRATES
      namf26 = trim(input_dir)//'chembc.112km.data'

!   unit 27 chemical lower boundary from SOCRATES
      namf27 = trim(input_dir)//'chembc.15km.data'

!   unit 28 water vapor needed for radiation 
      namf28 = trim(input_dir)//'rad_gases.dat'

!   unit 30 table for photolysis (phot)
      namf30 = trim(input_dir)//'jtable2.tuv4_2.koppers.dat' 

!   unit 40 Curtis matrices
      namf40 = trim(input_dir)//'curtis.dat'

!   unit 52 GSWM diurnal tide at 16 km - 4 seasons
      namf52 = trim(input_dir)//'diurnal.16km.data'

!   unit 53 GSWM semidiurnal tide at 16 km - 4 seasons
      namf53 = trim(input_dir)//'semi.16km.data'

!-------------------------------------------------------------------
!... input
!-------------------------------------------------------------------

!   unit 17 - startup fields

!      if (nstart.eq.0) then
!         ! start from initial conditions
!         ! namf17 = rose_data//'NMC_SOC.day001_2002.dat'
!         ! namf17 = rose_data//'NMC_SOC.day151_2002.dat'
!         write(namf17,'(a,a)') trim(input_dir),trim(restart_file)
!      else
!         ! restart fields from previous run
!         write(namf17,'(a,"rose_restart.dat")') trim(input_dir) 
!      end if
       namf17 = trim(restart_file)

!-------------------------------------------------------------------
!... output
!-------------------------------------------------------------------

    if (output_prog_diag) then
!   unit 10 output for plotting, diagnostics, etc. - STANDARD 
      write(namf10,120) trim(out_dir)
 120  format(a, 'day151.nc')
      print *, namf10

!   unit 41 output for plotting - ADDITIONAL DIAGNOSTICS
      write(namf41,125) trim(out_dir)
 125  format(a, 'day151.diag.nc')
      print *, namf41

!   unit 13 - 32 times per day output for tidal analysis
      write(namf13,140) trim(out_dir)
 140  format(a, 'day151.tides.nc')
    endif

!   unit 14 data for restarting integration
!     write(namf14,'(a,"rose_restart.dat")') trim(out_dir)
      write(namf14,'("rose_restart.dat")') 

end subroutine msetvar

end module DArose_mod
