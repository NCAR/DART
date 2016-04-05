!-------------------------------------------------------------------
      subroutine msetvar
!-------------------------------------------------------------------
!  specify variable parameters for 3-d model run
!-------------------------------------------------------------------

      use params
      use dynam
      use chem,   only : nrates
      use chem0d, only : niter, dtchem

      implicit none

!... local variables

      integer :: lengdy, maxdays

      character*50 namf11            ! run log file

      character (len=11)  :: rose_data = '../DAinput/'
      character (len=12) :: out_dir  = '../DAoutput/'

!-------------------------------------------------------------------
!... set program control parameters 
!-------------------------------------------------------------------

      ntime = 8                      ! number of steps per hour

      deltat = 3600./float(ntime)
      dtleap = deltat*2.
      dtchem = deltat

      lengdy = ntime*24

!... length of integration

      maxdays = 7 
      nend = maxdays * lengdy

!... frequency at which output is saved

!      nout = 1                      ! save every step
!      nout = 6                       ! every 45 minutes
!      nout = ntime                  ! save every hour 
      nout = ntime * 12             ! save every 12 hours 
!      nout = ntime *24              ! save every day
!      nout = lengdy * 2             ! or less frequently
      
      noutdiag = nout                ! freq to save diagnostics 
                                     ! (typically nout)

      nouttid = 10                   ! frequency IN DAYS for higher 
                                     ! freq output for tides

!-------------------------------------------------------------------
!  other model control parameters
!  variable - update for each run
!-------------------------------------------------------------------

      nstart = 0                     ! =0 initial run; =1 restart
      nseg = 1                       ! segment number (=1 for initial run)

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
!... fixed  
!-------------------------------------------------------------------

!   unit 12 zonal wind climatology for sponge layer
      namf12 = rose_data//'init.mlt.data'

!   unit 20 dynamical troposphere and lower boundary fields
      namf20 = rose_data//'nmc_lbc.02.nc'

!   unit 23 climatology zonal mean dynamical ubc
      namf23 = rose_data//'mlt.ubc.data'

!   unit 26 chemical upper boundary from SOCRATES
      namf26 = rose_data//'chembc.112km.data'

!   unit 27 chemical lower boundary from SOCRATES
      namf27 = rose_data//'chembc.15km.data'

!   unit 28 water vapor needed for radiation 
      namf28 = rose_data//'rad_gases.dat'

!   unit 30 table for photolysis (phot)
      namf30 = rose_data//'jtable2.tuv4_2.koppers.dat' 

!   unit 40 Curtis matrices
      namf40 = rose_data//'curtis.dat'

!   unit 52 GSWM diurnal tide at 16 km - 4 seasons
      namf52 = rose_data//'diurnal.16km.data'

!   unit 53 GSWM semidiurnal tide at 16 km - 4 seasons
      namf53 = rose_data//'semi.16km.data'

!-------------------------------------------------------------------
!... variable
!-------------------------------------------------------------------

!-------------------------------------------------------------------
!... input
!-------------------------------------------------------------------

!   unit 17 - startup fields

      if (nstart.eq.0) then
         ! start from initial conditions
         ! namf17 = rose_data//'NMC_SOC.day001_2002.dat'
         namf17 = rose_data//'NMC_SOC.day151_2002.dat'
      else
         ! restart fields from previous run
         write(namf17,160) out_dir, nseg-1
      end if

!-------------------------------------------------------------------
!... output
!-------------------------------------------------------------------

!   unit 10 output for plotting, diagnostics, etc. - STANDARD 
      write(namf10,120) out_dir, nseg
 120  format(a, 'day151.nc',i1)
      print *, namf10

!   unit 41 output for plotting - ADDITIONAL DIAGNOSTICS
      write(namf41,125) out_dir, nseg
 125  format(a, 'day151.diag.nc',i1)
      print *, namf41

!   unit 13 - 32 times per day output for tidal analysis
      write(namf13,140) out_dir, nseg
 140  format(a, 'day151.tides.nc',i1)

!   unit 14 data for restarting integration
      write(namf14,160) out_dir, nseg
 160  format(a, 'restart.day151.seg',i1)

!-------------------------------------------------------------------
      end subroutine msetvar
!-------------------------------------------------------------------

