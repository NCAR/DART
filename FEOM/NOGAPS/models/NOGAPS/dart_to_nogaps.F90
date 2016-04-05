! This code may (or may not) be part of the NOGAPS distribution,
! So it is not protected by the DART copyright agreement.
!
! DART $Id$

program dart_to_nogaps

!----------------------------------------------------------------------
! purpose: interface between nogaps and DART
!
! method: Read dart state vector with model state
!         Reforms fields into nogaps arrays.
!         Transform into spectral space.
!         Write out NOGAPS shist000000 file.
! 
! USAGE:  The filenames is read from the dart_to_nogaps namelist
!         <edit dart_to_nogaps_input_file in input.nml:dart_to_nogaps_nml>
!         dart_to_nogaps
!
!----------------------------------------------------------------------

use        types_mod, only  : r4, r8
use    utilities_mod, only  : E_ERR, E_WARN, E_MSG, error_handler, logfileunit, &
                              initialize_utilities, timestamp, &
                              find_namelist_in_file, check_namelist_read, &
                              open_file, close_file
use        model_mod, only  : static_init_model, get_model_size, get_gridsize
use  assim_model_mod, only  : aread_state_restart, open_restart_read, close_restart
use time_manager_mod, only  : time_type, print_time, print_date, get_date, &
                              set_calendar_type, GREGORIAN, operator(-), get_time
use mpi_utilities_mod, only : initialize_mpi_utilities, finalize_mpi_utilities, &
                              iam_task0

use gcvars
use times

use netcdf
implicit none

! This is from the original assembla server we used during collaboration.
! character(len=128), parameter :: &
!    source   = "$orgURL: https://svn2.assembla.com/svn/ngdart/dart_to_nogaps.F90 $", &
!    revision = "$orgRevision: 112 $", &
!    revdate  = "$orgDate: 2010-06-11 13:02:17 -0600 (Fri, 11 Jun 2010) $"

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

!-----------------------------------------------------------------------
! namelist parameters with default values.
!-----------------------------------------------------------------------

character (len = 128) :: dart_to_nogaps_input_file   = 'dart_old_vector'
character (len = 128) :: nogaps_restart_filename     = 'specfiles/shist000000'
character (len = 128) :: dart_data_time_filename     = 'dart_data.time'
logical               :: advance_time_present        = .TRUE.

namelist /dart_to_nogaps_nml/ nogaps_restart_filename, dart_to_nogaps_input_file , &
                              dart_data_time_filename, advance_time_present

!----------------------------------------------------------------------
! global storage
!----------------------------------------------------------------------

integer               :: io, iunit, x_size
type(time_type)       :: model_time, advance_time, forecast_length
real(r8), allocatable :: statevector(:)
logical               :: verbose = .FALSE.
integer               :: numlons, numlats, numlevels
character(len=10)     :: vector_time
integer               :: year, month, day, hour, minute, second
character(len=128)    :: msgstring

!----------------------------------------------------------------------

      include 'imp.h'
      include 'param.h'
 
      integer ioc
      integer istat

      logical lfcst
      logical lgtrdy
      logical lnmode
      logical loutp
      logical ptime

      data lgtrdy/.true./, lnmode/.true./, lfcst/.true./, loutp/.true./, ptime/.false./

      namelist /rdilst/ lgtrdy, lnmode, lfcst, loutp, ptime, reprod, &
                        al2al, lmpi2, nproc, jsplit, iproc, jproc

!----------------------------------------------------------------------
! This needs to be the first thing that executes, because it starts
! MPI and does all the setup for DART.  Don't try to print or write
! before executing this.
!----------------------------------------------------------------------

call initialize_mpi_utilities(progname='dart_to_nogaps')

!----------------------------------------------------------------------
! Read the namelist to get the input and output filenames.
!----------------------------------------------------------------------

call find_namelist_in_file("input.nml", "dart_to_nogaps_nml", iunit)
read(iunit, nml = dart_to_nogaps_nml, iostat = io)
call check_namelist_read(iunit, io, "dart_to_nogaps_nml") ! closes, too.

if (iam_task0()) then
   write(*,*)
   write(*,'(''dart_to_nogaps:converting DART file '',A,   &
                               '' into nogaps file '',A)') &
   trim(dart_to_nogaps_input_file), trim(nogaps_restart_filename)
endif

!----------------------------------------------------------------------
! Call model_mod:static_init_model(), which reads the namelists
! to set calendar type, starting date, deltaT, etc.
!----------------------------------------------------------------------

call set_calendar_type(GREGORIAN)

call static_init_model()

x_size = get_model_size()
allocate(statevector(x_size))

call get_gridsize(numlons, numlats, numlevels)

!----------------------------------------------------------------------
! Read in the state vector from a dart file that has time information
! in the header. Minimally, state vectors have a time tag with the valid
! time of the ensuing state. Sometimes (under namelist control) the state
! also has an 'advance_to' time. If this is present, we can calculate
! a forecast length ... write these times to a text file that may be
! used by the scripts. (may obviate need for 'trans_time' - if we can
! resolve the CRDATE.dat conundrum)
!----------------------------------------------------------------------

if (iam_task0()) then

   iunit = open_restart_read(dart_to_nogaps_input_file)

   if ( advance_time_present ) then
      call aread_state_restart(model_time, statevector, iunit, advance_time)
   else
      call aread_state_restart(model_time, statevector, iunit )
   endif
   call close_restart(iunit)

   
   iunit = open_file(dart_data_time_filename, 'formatted', 'write')
   
   call get_date(model_time, year, month, day, hour, minute, second)
   write(iunit,'(I4.4,3(I2.2))') year, month, day, hour
 
   if ( advance_time_present ) then
      call get_date(advance_time, year, month, day, hour, minute, second)
      write(iunit, '(I4.4,3(I2.2))') year, month, day, hour
   endif

   if (advance_time_present) then

      ! calculate number of hours in forecast
   
      forecast_length = advance_time - model_time
      call get_time(forecast_length, second, day)
      hour = day*24 + second/3600
      write(iunit, '(I3.3)') hour

      minute = mod(second,3600)
      if (minute .ne. 0) then
         write(msgstring,*) 'forecast length ',hour,':',minute,'not integer number of hours.'
         call error_handler(E_ERR,'main',msgstring,source,revision,revdate)
      endif

   endif

   call close_file(iunit)

   !-------------------------------------------------------------------
   ! Print the times we found in the advance file.
   !-------------------------------------------------------------------

   call print_date(  model_time, str='dart_to_nogaps:DART   model date')
   call print_time(  model_time, str='dart_to_nogaps:DART   model time')
   if (advance_time_present) then
      call print_date(advance_time, str='dart_to_nogaps:DART advance date')
      call print_time(advance_time, str='dart_to_nogaps:DART advance time')
      call print_time(forecast_length, str='dart_to_nogaps:forecast length')
   endif

endif

! FIXME: we might need a broadcast here to spread the state vector
! to all the other tasks.  i'm not sure what dart_to_shist is expecting
! to find.   if it needs to be broadcast, this will do it:
 
!   call array_broadcast(statevector, 0)
 

!----------------------------------------------------------------------
! NOGAPS setup code, and then call dart_to_shist which does the work.
!----------------------------------------------------------------------

!  ********************************************************************
!

      call MPI_Comm_rank(mpi_comm_world,irank,istat)
      call MPI_Comm_size(mpi_comm_world,nproc,istat)

      al2al = .false.
      lmpi2 = .true.
      reprod = .true.
      jsplit = 1
      ir = irank + 1
      iproc = 0
      jproc = 0
!
#ifdef DP
      mpireal  = mpi_double_precision
      mpireal4 = mpi_real
#elif C90
      mpireal  = mpi_real
      mpireal4 = mpi_real
#else
      mpireal  = mpi_real
      mpireal4 = mpi_real
#endif
      mpicomm    = mpi_comm_world
      mpiint     = mpi_integer
      mpilogical = mpi_logical
!
!     call wnlline(1)
!
      if(irank.eq.0)then
!
        open(unit=4,file='rdifil',form='formatted',status='old', iostat=ioc)
        if(ioc.eq.0)read(4,rdilst)
!
        write(6,rdilst)
      endif
!
! allocate memeory and set up pointers to timing 'bins'
!
      allocate(rtime(2,30))
!
      call zilch(rtime,60)
!
      rjobtm  => rtime(:, 1)
      rtfcst  => rtime(:, 2)
      rtnm    => rtime(:, 3)
      rts2p   => rtime(:, 4)
      rtp2s   => rtime(:, 5)
      rtget   => rtime(:, 6)
      rtrs    => rtime(:, 7)
      rtsr    => rtime(:, 8)
      rtmisc  => rtime(:, 9)
      rtcubiy => rtime(:,10)
      rtcubij => rtime(:,11)
      rtcmp3  => rtime(:,12)
      rtdiag  => rtime(:,13)
      rtimpl  => rtime(:,14)
      rtcubv  => rtime(:,15)
      rtlsp   => rtime(:,16)
      rtcuml  => rtime(:,17)
      rtlwr   => rtime(:,18)
      rtswr   => rtime(:,19)
      rtpbl   => rtime(:,20)
      rtshcu  => rtime(:,21)
      rtfnio  => rtime(:,22)
      rtalloc => rtime(:,23)
      rtread  => rtime(:,24)
      rtwrit  => rtime(:,25)
      rthsort => rtime(:,26)
      rthgath => rtime(:,27)
      rthscat => rtime(:,28)
!
      rtgwd   => rtime(:,29)
      rtscmix => rtime(:,30)
!
!
! ********************************************************************
!
      call mpi_bcast(lgtrdy,1,mpilogical,0,mpicomm,istat)
      call mpi_bcast(lnmode,1,mpilogical,0,mpicomm,istat)
      call mpi_bcast( lfcst,1,mpilogical,0,mpicomm,istat)
      call mpi_bcast( loutp,1,mpilogical,0,mpicomm,istat)
      call mpi_bcast( ptime,1,mpilogical,0,mpicomm,istat)
      call mpi_bcast(reprod,1,mpilogical,0,mpicomm,istat)
      call mpi_bcast( al2al,1,mpilogical,0,mpicomm,istat)
      call mpi_bcast( lmpi2,1,mpilogical,0,mpicomm,istat)
      call mpi_bcast(jsplit,1,mpiint,    0,mpicomm,istat)
      call mpi_bcast( iproc,1,mpiint,    0,mpicomm,istat)
      call mpi_bcast( jproc,1,mpiint,    0,mpicomm,istat)
!
! check to see if geometry is OK
!
      if(iproc.gt.(jtrun/2+1))then
        if(irank.eq.0)then
          print *, 'IPROC should be .LE. (JTRUN/2+1)'
          print *, ' TERMINATING '
        endif
        go to 100
      endif
!
      if(iproc.gt.lm)then
        if(irank.eq.0)then
          print *, 'IPROC should be .LE. LM'
          print *, ' TERMINATING '
        endif
        go to 100
      endif
!
      if(jproc.gt.jm/2)then
        if(irank.eq.0)then
          print *, 'JPROC should be .LE. JM/2'
          print *, ' TERMINATING '
        endif
        go to 100
      endif
!
      if(jproc.gt.jtrun)then
        if(irank.eq.0)then
          print *, 'JPROC should be .LE. JTRUN'
          print *, ' TERMINATING '
        endif
        go to 100
      endif

!----------------------------------------------------------------------
!  This is where the magic happens.  This calls NOGAPS code which
!  uses mpi, and converts the state vector data back into NOGAPS
!  shist format.
!----------------------------------------------------------------------

      call dart_to_shist(lgtrdy, lnmode, lfcst, loutp, ptime, &
         size(statevector), statevector)


 100 continue

deallocate(statevector)

!----------------------------------------------------------------------
! This closes the files, shuts down MPI.  We shouldn't print or try to
! write anything after this call.
!----------------------------------------------------------------------

call finalize_mpi_utilities()

end program dart_to_nogaps


! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

