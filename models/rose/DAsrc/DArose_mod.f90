module DArose_mod

   use params, only : ntime, nout, noutdiag, nouttid, nend, nseg,&
       & nstart, nsave
   use dynam,  only : deltat, dtleap, ntrans, ndiabat, ninterp, &
       & namf10, namf12, namf13, namf14, namf17, & 
       & namf20, namf23, namf26, namf27, namf28, &
       & namf30, namf40, namf41, namf52, namf53       
   use chem,   only : nrates
   use chem0d, only : niter, dtchem

   use types_mod, only: r8, pi
   use utilities_mod, only: open_file, close_file, &
                            check_nml_error, file_exist

   implicit none
   private
   public :: msetvar, old_restart, output_prog_diag, & 
             h_tune, z_tune

   logical :: old_restart = .true.
   logical :: output_prog_diag = .false.
   character (len=30) :: input_dir = '../DAinput/'
   character (len=30) :: out_dir   = '../DAoutput/'
   character (len=30) :: ncep_file = 'nmc_lbc.02.nc'
   character (len=30) :: restart_file = 'NMC_SOC.day151_2002.dat'  
   real(kind=r8) :: h_tune = pi
   real(kind=r8) :: z_tune = 1.0
   real(kind=r8) :: target_time = 168.0 ! 7 days * 24 [hr]  
   integer       :: lengdy 


   namelist /rose_nml/ old_restart, nstart, target_time, &
                       input_dir, out_dir, ncep_file, restart_file,&
                       output_prog_diag, &
                       h_tune, z_tune, &
                       ntime
   contains

subroutine msetvar

      implicit none

      integer :: iunit, io, ierr
 
      !data nstart /0/                     ! =0 initial run; =1 restart
      !data ntime  /8/                     ! number of steps per hour

      if(file_exist('rose.nml')) then
         iunit = open_file('rose.nml', action = 'read')
         read(iunit, nml = rose_nml, iostat = io)
         ierr = check_nml_error(io, 'gswm_nml')
         call close_file(iunit)
      endif

      write(*,*) 'rose_nml values  --- may come from rose.nml'
      write(*,*) 'old_restart ', old_restart
      write(*,*) 'nstart ', nstart
      write(*,*) 'target_time ', target_time
      write(*,*) 'input_dir ', input_dir 
      write(*,*) 'out_dir ', out_dir
      write(*,*) 'ncep_file ', ncep_file
      write(*,*) 'restart_file ', restart_file
      write(*,*) 'output_prog_diag ', output_prog_diag
      write(*,*) 'h_tune ', h_tune
      write(*,*) 'z_tune ', z_tune
      write(*,*) 'ntime ', ntime

!-------------------------------------------------------------------
!... set program control parameters 
!-------------------------------------------------------------------

      lengdy = ntime*24

      deltat = 3600./float(ntime)
      dtleap = deltat*2.
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

      if (nstart == 0) then  
        nseg = 1       ! segment number (=1 for initial run)
      else             ! used for tracking NetCDF output files
        nseg = 2       ! USED WHEN output_prog_diag = .true. 
      endif

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

      if (nstart.eq.0) then
         ! start from initial conditions
         ! namf17 = rose_data//'NMC_SOC.day001_2002.dat'
         ! namf17 = rose_data//'NMC_SOC.day151_2002.dat'
         write(namf17,'(a,a)') trim(input_dir),trim(restart_file)
      else
         ! restart fields from previous run
         write(namf17,'(a,"rose_restart.dat")') trim(input_dir) 
      end if

!-------------------------------------------------------------------
!... output
!-------------------------------------------------------------------

    if (output_prog_diag) then
!   unit 10 output for plotting, diagnostics, etc. - STANDARD 
      write(namf10,120) out_dir, nseg
 120  format(a, 'day151.nc', i1)
      print *, namf10

!   unit 41 output for plotting - ADDITIONAL DIAGNOSTICS
      write(namf41,125) out_dir, nseg
 125  format(a, 'day151.diag.nc',i1)
      print *, namf41

!   unit 13 - 32 times per day output for tidal analysis
      write(namf13,140) out_dir, nseg
 140  format(a, 'day151.tides.nc',i1)
    endif

!   unit 14 data for restarting integration
      write(namf14,'(a,"rose_restart.dat")') trim(out_dir)

end subroutine msetvar

end module DArose_mod
