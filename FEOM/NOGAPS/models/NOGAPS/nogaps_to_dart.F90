! This code may (or may not) be part of the NOGAPS distribution,
! So it is not protected by the DART copyright agreement.
!
! DART $Id$

program nogaps_to_dart

!----------------------------------------------------------------------
! purpose: interface between nogaps and DART
!
! method: Read nogaps "restart" files of model state
!         Reform fields into a DART state vector (control vector).
!         Write out state vector in "proprietary" format for DART.
!         The output is a "DART restart file" format.
! 
! USAGE:  The nogaps filename is read from the nogaps_in namelist
!         <edit nogaps_to_dart_output_file in input.nml:nogaps_to_dart_nml>
!         nogaps_to_dart
!
!----------------------------------------------------------------------

use        types_mod, only  : r4, r8
use    utilities_mod, only  : E_ERR, E_WARN, E_MSG, error_handler, logfileunit, &
                              initialize_utilities, timestamp, &
                              find_namelist_in_file, check_namelist_read, &
                              open_file, close_file
use        model_mod, only  : static_init_model, get_model_size, get_gridsize
use  assim_model_mod, only  : awrite_state_restart, open_restart_write, close_restart
use time_manager_mod, only  : time_type, print_time, print_date, set_date, &
                              set_calendar_type, GREGORIAN
use mpi_utilities_mod, only : initialize_mpi_utilities, finalize_mpi_utilities, &
                              iam_task0

use gcvars
use times

use netcdf
implicit none


! This is from the original assembla server we used during collaboration.
! character(len=128), parameter :: &
!    source   = "$orgURL:https://svn2.assembla.com/svn/ngdart/nogaps_to_dart.F90 $", &
!    revision = "$orgrev: 112 $", &
!    revdate  = "$orgDate: 2010-06-11 13:02:17 -0600 (Fri, 11 Jun 2010) $"

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"


!-----------------------------------------------------------------------
! namelist parameters with default values.
!-----------------------------------------------------------------------

character (len = 128) :: nogaps_restart_filename     = 'specfiles/shist000000'
character (len = 128) :: nogaps_data_time_filename   = 'CRDATE.dat'
character (len = 128) :: nogaps_to_dart_output_file  = 'dart_new_vector'

namelist /nogaps_to_dart_nml/ nogaps_restart_filename, &
                              nogaps_data_time_filename, &
                              nogaps_to_dart_output_file

!----------------------------------------------------------------------
! global storage
!----------------------------------------------------------------------

integer               :: io, iunit, x_size
type(time_type)       :: model_time
real(r8), allocatable :: statevector(:), phis(:)
logical               :: verbose = .FALSE.
integer               :: numlons, numlats, numlevels
character (len = 10 ) :: vector_time
integer               :: year, month, day, hour

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

!----------------------------------------------------------------------
! This needs to be the first thing that executes, because it starts
! MPI and does all the setup for DART.  Don't try to print or write
! before executing this.
!----------------------------------------------------------------------

call initialize_mpi_utilities(progname='nogaps_to_dart')

!----------------------------------------------------------------------
! Read the namelist to get the input and output filenames.
!----------------------------------------------------------------------

call find_namelist_in_file("input.nml", "nogaps_to_dart_nml", iunit)
read(iunit, nml = nogaps_to_dart_nml, iostat = io)
call check_namelist_read(iunit, io, "nogaps_to_dart_nml") ! closes, too.

if (iam_task0()) then
   write(*,*)
   write(*,'(''nogaps_to_dart:converting nogaps restart file '',A, &
         &'' to DART file '',A)') &
          trim(nogaps_restart_filename), trim(nogaps_to_dart_output_file)
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

allocate(phis(numlons*numlats))

!----------------------------------------------------------------------
! NOGAPS setup code, and then call shist_to_dart which does the work.
!----------------------------------------------------------------------

!  ********************************************************************
 

      call MPI_Comm_rank(mpi_comm_world,irank,istat)
      call MPI_Comm_size(mpi_comm_world,nproc,istat)

      al2al = .false.
      lmpi2 = .true.
      reprod = .true.
      jsplit = 1
      ir = irank + 1
      iproc = 0
      jproc = 0
 
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
 
        open(unit=4,file='rdifil',form='formatted',status='old', iostat=ioc)
        if(ioc.eq.0)read(4,rdilst)
 
        write(6,rdilst)
      endif
!
! allocate memeory and set up pointers to timing 'bins'
!
      allocate(rtime(2,30))
 
      call zilch(rtime,60)
 
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
 
      rtgwd   => rtime(:,29)
      rtscmix => rtime(:,30)
 
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
 
      if(iproc.gt.lm)then
        if(irank.eq.0)then
          print *, 'IPROC should be .LE. LM'
          print *, ' TERMINATING '
        endif
        go to 100
      endif
 
      if(jproc.gt.jm/2)then
        if(irank.eq.0)then
          print *, 'JPROC should be .LE. JM/2'
          print *, ' TERMINATING '
        endif
        go to 100
      endif
 
      if(jproc.gt.jtrun)then
        if(irank.eq.0)then
          print *, 'JPROC should be .LE. JTRUN'
          print *, ' TERMINATING '
        endif
        go to 100
      endif

!----------------------------------------------------------------------
!  This is where the magic happens.  This calls NOGAPS code which
!  uses mpi, and computes both the contents of the state vector, and
!  fills in the phis array.
!----------------------------------------------------------------------


      call shist_to_dart(lgtrdy, lnmode, lfcst, loutp, ptime, &
         size(statevector), statevector, size(phis), phis)


! FIXME: do we need to gather the data back to task 0 if running
! with more than 1 task?   hopefully not - task 0 should have a
! full statevector with all the data.  at least, that's what this
! code is assuming.

!----------------------------------------------------------------------
! Only have task 0 write out the restart file and the timestamp file.
! FIXME: Do we have to write the phis array somewhere?  do it inside
! the task0 loop, so you only get 1 copy if running with tasks > 1.
!----------------------------------------------------------------------

if (iam_task0()) then

   iunit = open_file(nogaps_data_time_filename, 'formatted', 'read')
   read(iunit, "(A10)") vector_time
   call close_file(iunit)
   
   !print *, 'got this from time file: ', vector_time
   
   read(vector_time(1:4 ), "(I4)") year
   read(vector_time(5:6 ), "(I2)") month
   read(vector_time(7:8 ), "(I2)") day
   read(vector_time(9:10), "(I2)") hour
   model_time = set_date(year, month, day, hour, 0, 0)
   
   iunit = open_restart_write(nogaps_to_dart_output_file)
   
   call awrite_state_restart(model_time, statevector, iunit)
   call close_restart(iunit)

   call print_date(model_time, str='nogaps_to_dart:nogaps model date')
   call print_time(model_time, str='nogaps_to_dart:DART   model time')

endif

 100 continue

deallocate(statevector, phis)

!----------------------------------------------------------------------
! This closes the files, shuts down MPI.  We shouldn't print or try to
! write anything after this call.
!----------------------------------------------------------------------

call finalize_mpi_utilities()

end program nogaps_to_dart

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

