! Data Assimilation Research Testbed -- DART
! Copyright 2006, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module mpi_utilities_mod

! <next five lines automatically updated by CVS, do not edit>
! $Source$ 
! $Revision$ 
! $Date$ 
! $Author$ 
! $Name$ 

!-----------------------------------------------------------------------------
!
!   A collection of interfaces to the MPI (Message Passing Interface)
!   multi-processor communication library routines.
!
!  *** initialize_mpi_utilities()  Subroutine that initializes MPI and sets
!                                  local values needed later.  Must be called
!                                  before any other routine here.
!
!  *** finalize_mpi_utilities()  Subroutine that shuts down MPI cleanly.
!                                Must be called before program exits, and no
!                                other routines here can be used afterwards.
!
!  *** task_count()       Function that returns the total number of MPI tasks.
!
!  *** my_task_id()       Function that returns my task number.  Note that
!                         in the MPI world task numbers run from 0 to N-1.
!
!  *** send_to()          Subroutine which sends a 1D data array
!                         synchronously to another task (point-to-point).
!
!  *** receive_from()     Subroutine which receives a 1D data array
!                         synchronously from another task (point-to-point).
!
!  *** task_sync()        Subroutine that only returns after every task has
!                         reached this same location in the code.
!        
!    * transpose_array()  Subroutine that transposes a 2D array
!                         from column-major to row-major or back.
!
!  *** array_broadcast()  Subroutine that sends a copy of the entire data 
!                         array to all other tasks. 
!                
!   ** array_distribute() Subroutine that distributes a data array across the
!                         other tasks, so each task gets a non-overlapping 
!                         subset of the data.
!                
!   MPI cover routines more specific for DART and hopefully more useful.
!
!  *** iam_task0()        Function which returns .TRUE. if task id is 0,
!                         .FALSE. for anything else.
!
!  *** broadcast_send()   Subroutine which takes two r8 arrays and broadcasts
!                         them to all other tasks.  If sending ID is not the
!                         same as the local task ID, an error is returned.
!                         Does not return until all other tasks have called
!                         recv to pick up the data.
!
!  *** broadcast_recv()   Subroutine which receives two r8 arrays from the 
!                         sending task ID.  If the sending ID is the same as
!                         the local task ID, an error is returned.  All other
!                         tasks must call recv before this routine returns.
!
!   Lower level utility routines which interact with the utilities_mod.f90 
!   code to open a named pipe per MPI task, read and write from them, and
!   close and/or remove them.
!
!    * make_pipe()        Function that creates a named pipe (fifo), opens it,
!                         and returns the unit number.  Ok to call if the pipe
!                         already exists or is already open; it will skip
!                         those steps and just return the unit number.  The 
!                         name argument is used as a base and a filename
!                         in the form 'base.NNNN' is generated, where the N's
!                         are the MPI rank number, 0 padded.
!
!    * destroy_pipe()     The unit number is closed and the pipe file is 
!                         removed.
!
!    * read_pipe()        The character string is read from the pipe.
!                         (Can be overloaded to read ints if time or status
!                         info is useful to exchange between processes.) 
!                         This routine blocks until data is available.
!
!    * write_pipe()       The character string is written to the pipe.
!                         (Can be overloaded to write ints if time or status
!                         info is useful to exchange between processes.) 
!                         This routine writes and returns immediately.
!
!
! *** both code and interface are done (but untested so far)
!  ** interface with proposed arguments exists but code not complete
!   * interface name only; no arg list devised yet 
!
!-----------------------------------------------------------------------------
! 
! these do not exist - i believe a single transpose will work.  but if not,
! they can be separated into these two, which can either work on a real
! 2D array or a single linearized array which is logically 2D but in reality
! stored in a 1D fortran array:
!
!      transpose_row_major()  Subroutine which transposes a logical 2D array
!                             from column-major to row-major.  The source and
!                             destination arrays must be stored in 1D arrays 
!                             of length (nrows * ncols).
!
!      transpose_col_major()  Subroutine which transposes a logical 2D array
!                             from row-major to column-major.  The source and
!                             destination arrays must be stored in 1D arrays 
!                             of length (nrows * ncols).
!
!-----------------------------------------------------------------------------

use types_mod, only : r8
use utilities_mod, only : register_module, error_handler, & 
                          E_ERR, E_WARN, E_MSG, E_DBG
use time_manager_mod, only : time_type, get_time, set_time

!
! Some MPI installations have an MPI module; if one is present, use that.
! If not, there will be an MPI include file which defines the parameters.
! Use one but not both.  For help on compiling a module which uses MPI
! see the $DART/doc/mpi directory.

!use mpi

implicit none
private

include "mpif.h"


!   ---- private data for mpi_utilities ----

integer :: myrank          ! my mpi number
integer :: total_tasks     ! total mpi tasks/procs
integer :: my_local_comm   ! duplicate communicator private to this file
integer :: comm_size       ! if ens count < tasks, only the first N participate

!!! probably not needed; most of these are in the ensemble_handle.  but i just
!!! copied them all over from my test program "in case".  delete as it is clear
!!! that they aren't needed.
!!integer :: state_size    ! number of state vars in a single model
!!integer :: ens_size      ! number of models/ensembles running
!!integer :: obs_size      ! number of observations available
!!integer :: ens_per_task  ! if number of ensembles > tasks, how many per task
!!integer :: states_per_task ! state vars generally > tasks; how many per task
!!integer :: obs_per_task  ! not sure if this is needed; unused so far.
!!integer :: ec_total_size ! total array size for ensemble-complete data
!!integer :: sc_total_size ! total array size for state-complete data
!!integer :: max_print     ! limit for value dumps


public :: task_count, my_task_id, transpose_array, &
          initialize_mpi_utilities, finalize_mpi_utilities, &
          task_sync, array_broadcast, array_distribute, &
          send_to, receive_from, iam_task0, broadcast_send, broadcast_recv

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

logical, save :: module_initialized = .false.

character(len = 129) :: errstring


! Namelist input - placeholder for now.

!namelist /mpi_utilities_nml/ x

contains

!-----------------------------------------------------------------------------
! mpi cover routines
!-----------------------------------------------------------------------------

subroutine initialize_mpi_utilities()

! Initialize MPI and query it for global information.  Make a duplicate
! communicator so that any user code which wants to call MPI will not 
! interfere with any outstanding asynchronous requests, accidental tag
! matches, etc.  This routine must be called before any other routine in
! this file, and it should not be called more than once (but it does have
! defensive code in case that happens.)

integer :: errcode
logical :: already

if ( module_initialized ) then
   ! return without calling the code below multiple times
   write(errstring, *) 'initialize_mpi_utilities has already been called'
   call error_handler(E_WARN,'initialize_mpi_utilities', errstring, source, revision, revdate)
   return
endif

if ( .not. module_initialized ) then
   ! Initialize the module with utilities
   call register_module(source, revision, revdate)
   module_initialized = .true.
endif

! allow for the possibility that the user code has already initialized
! the MPI libs and do not try to call initialize again.

call MPI_Initialized(already, errcode)
if (errcode /= MPI_SUCCESS) then
   write(errstring, '(a,i8)') 'MPI_Initialized returned error code ', errcode
   call error_handler(E_ERR,'initialize_mpi_utilities', errstring, source, revision, revdate)
endif

! generally this will be the first use of mpi and this code will be executed.
if (.not.already) then
   errcode = -999
   call MPI_Init(errcode)
   if (errcode /= MPI_SUCCESS) then
      write(errstring, '(a,i8)') 'MPI_Init returned error code ', errcode
      call error_handler(E_ERR,'initialize_mpi_utilities', errstring, source, revision, revdate)
   endif
endif

! duplicate the world communicator to isolate us from any other user
! calls to MPI.  All subsequent mpi calls here will use the local communicator
! and not the global world comm.
call MPI_Comm_dup(MPI_COMM_WORLD, my_local_comm, errcode)
if (errcode /= MPI_SUCCESS) then
   write(errstring, '(a,i8)') 'MPI_Comm_dup returned error code ', errcode
   call error_handler(E_ERR,'initialize_mpi_utilities', errstring, source, revision, revdate)
endif

! find out who we are (0 to N-1).
call MPI_Comm_rank(my_local_comm, myrank, errcode)
if (errcode /= MPI_SUCCESS) then
   write(errstring, '(a,i8)') 'MPI_Comm_rank returned error code ', errcode
   call error_handler(E_ERR,'initialize_mpi_utilities', errstring, source, revision, revdate)
endif

! number of tasks (if 10, returns 10.  task id numbers go from 0-9)
call MPI_Comm_size(my_local_comm, total_tasks, errcode)
if (errcode /= MPI_SUCCESS) then
   write(errstring, '(a,i8)') 'MPI_Comm_size returned error code ', errcode
   call error_handler(E_ERR,'initialize_mpi_utilities', errstring, source, revision, revdate)
endif

! TODO: if there are fewer ensembles than tasks, all the collective routines
! need to take that into account and not participate if they are > comm_size.
comm_size = total_tasks

! MPI successfully initialized.

end subroutine initialize_mpi_utilities

!-----------------------------------------------------------------------------

subroutine finalize_mpi_utilities(callfinalize)
 logical, intent(in), optional :: callfinalize

! Shut down MPI cleanly.  This must be done before the program exits; on
! some implementations of MPI the final I/O flushes are not done until this
! is called.  The optional argument can prevent us from calling MPI_Finalize,
! so that user code can continue to use MPI after this returns.  For good
! coding practice you should not call any other routines in this file
! after calling this routine.

integer :: errcode
logical :: dofinalize

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'finalize_mpi_utilities', errstring, source, revision, revdate)
endif

! Release the private communicator we created at init time.
call MPI_Comm_free(my_local_comm, errcode)
if (errcode /= MPI_SUCCESS) then
   write(errstring, '(a,i8)') 'MPI_Comm_free returned error code ', errcode
   call error_handler(E_ERR,'finalize_mpi_utilities', errstring, source, revision, revdate)
endif

! If the optional argument is not given, or is given and is true, 
! shut down mpi.  Only if the argument is specified and is false do we
! skip the finalization.
if (.not.present(callfinalize)) then
   dofinalize = .TRUE.
else if (callfinalize) then
   dofinalize = .TRUE.
else
   dofinalize = .FALSE.
endif

! Normally we shut down MPI here.  If the user tells us not to shut down MPI
! they must call this routine from their own code before exiting.
if (dofinalize) then
   call MPI_Finalize(errcode)
   if (errcode /= MPI_SUCCESS) then
      write(errstring, '(a,i8)') 'MPI_Finalize returned error code ', errcode
      call error_handler(E_ERR,'finalize_mpi_utilities', errstring, source, revision, revdate)
   endif
endif


end subroutine finalize_mpi_utilities


!-----------------------------------------------------------------------------

function task_count()

! Return the total number of MPI tasks.  e.g. if the number of tasks is 4,
! it returns 4.  (The actual task numbers are 0-3.)

integer :: task_count

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'task_count', errstring, source, revision, revdate)
endif

task_count = total_tasks

end function task_count


!-----------------------------------------------------------------------------

function my_task_id()

! Return my unique task id.  Values run from 0 to N-1 (where N is the
! total number of MPI tasks.

integer :: my_task_id

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'my_task_id', errstring, source, revision, revdate)
endif

my_task_id = myrank

end function my_task_id


!-----------------------------------------------------------------------------

subroutine task_sync()

! Synchronize all tasks.  This subroutine does not return until all tasks
! execute this line of code.

integer :: errcode

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'task_sync', errstring, source, revision, revdate)
endif

call MPI_Barrier(my_local_comm, errcode)
if (errcode /= MPI_SUCCESS) then
   write(errstring, '(a,i8)') 'MPI_Barrier returned error code ', errcode
   call error_handler(E_ERR,'task_sync', errstring, source, revision, revdate)
endif

end subroutine task_sync


!-----------------------------------------------------------------------------

subroutine send_to(dest_id, srcarray, time)
 integer, intent(in) :: dest_id
 real(r8), intent(in) :: srcarray(:)
 type(time_type), intent(in), optional :: time

! Send the srcarray to the destination id.
! If time is specified, it is also sent in a separate communications call.  
! This is a synchronous call; it will not return until the destination has 
! called receive to accept the data.  If the send_to/receive_from calls are 
! not paired correctly the code will hang.

integer :: i, tag, errcode
integer :: datasize
integer :: itime(2)

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'send_to', errstring, source, revision, revdate)
endif

! simple idiotproofing
if ((dest_id < 0) .or. (dest_id >= total_tasks)) then
   write(errstring, '(a,i8,a,i8)') "destination task id ", dest_id, &
                                   "must be >= 0 and < ", total_tasks
   call error_handler(E_ERR,'send_to', errstring, source, revision, revdate)
endif

!print *, "kind = ", kind(srcarray(1))
!print *, "digits = ", digits(srcarray(1))
! TODO: FIX THIS based on what kind and digits say
datasize = MPI_DOUBLE_PRECISION

! use my task id as the tag; unused at this point.
tag = myrank

! call MPI to send the data to the remote task
call MPI_Ssend(srcarray, size(srcarray), datasize, dest_id, tag, &
              my_local_comm, errcode)
if (errcode /= MPI_SUCCESS) then
   write(errstring, '(a,i8)') 'MPI_Ssend returned error code ', errcode
   call error_handler(E_ERR,'send_to', errstring, source, revision, revdate)
endif

! if time specified, call MPI again to send the 2 time ints.
if (present(time)) then
   call get_time(time, itime(1), itime(2))
   call MPI_Ssend(itime, 2, MPI_INTEGER, dest_id, tag*2, my_local_comm, errcode)
   if (errcode /= MPI_SUCCESS) then
      write(errstring, '(a,i8)') 'MPI_Ssend returned error code ', errcode
      call error_handler(E_ERR,'send_to', errstring, source, revision, revdate)
   endif
endif


end subroutine send_to


!-----------------------------------------------------------------------------

subroutine receive_from(src_id, destarray, time)
 integer, intent(in) :: src_id
 real(r8), intent(out) :: destarray(:)
 type(time_type), intent(out), optional :: time

! Receive data into the destination array from the src task.
! If time is specified, it is received in a separate communications call.  
! This is a synchronous call; it will not return until the source has 
! sent the data.  If the send_to/receive_from calls are not paired correctly 
! the code will hang.

integer :: i, tag, errcode
integer :: datasize
integer :: itime(2)
real(r8), pointer :: temparray(:)
integer :: status(MPI_STATUS_SIZE)

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'receive_from', errstring, source, revision, revdate)
endif

! simple idiotproofing
if ((src_id < 0) .or. (src_id >= total_tasks)) then
   write(errstring, '(a,i8,a,i8)') "source task id ", src_id, &
                                   "must be >= 0 and < ", total_tasks
   call error_handler(E_ERR,'receive_from', errstring, source, revision, revdate)
endif

!print *, "kind = ", kind(destarray(1))
!print *, "digits = ", digits(destarray(1))
! TODO: FIX THIS based on what kind and digits say
datasize = MPI_DOUBLE_PRECISION

! send_to uses its own id as the tag.
tag = src_id

! call MPI to receive the data from the remote task
call MPI_Recv(destarray, size(destarray), datasize, src_id, MPI_ANY_TAG, &
              my_local_comm, status, errcode)
if (errcode /= MPI_SUCCESS) then
   write(errstring, '(a,i8)') 'MPI_Recv returned error code ', errcode
   call error_handler(E_ERR,'receive_from', errstring, source, revision, revdate)
endif

! if time specified, call MPI again to send the 2 time ints.
if (present(time)) then
   call MPI_Recv(itime, 2, MPI_INTEGER, src_id, tag*2, &
                 my_local_comm, status, errcode)
   if (errcode /= MPI_SUCCESS) then
      write(errstring, '(a,i8)') 'MPI_Recv returned error code ', errcode
      call error_handler(E_ERR,'receive_from', errstring, source, revision, revdate)
   endif
   time = set_time(itime(1), itime(2))
endif


end subroutine receive_from


!-----------------------------------------------------------------------------
! TODO: do i need to overload this for both integer and real?
!       do i need to handle 1D, 2D, 3D inputs?


subroutine transpose_array

! not implemented here yet.  will have arguments -- several of them.

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'transpose_array', errstring, source, revision, revdate)
endif

write(errstring, *) 'not implemented yet'
call error_handler(E_ERR,'transpose_array', errstring, source, revision, revdate)

end subroutine transpose_array


!-----------------------------------------------------------------------------
! TODO: do i need to overload this for both integer and real?
!       do i need to handle 2D inputs?

subroutine array_broadcast(array, root)
 real(r8), intent(inout) :: array(:)
 integer, intent(in) :: root

! The data array values on the root task will be broadcast to every other
! task.  When this routine returns, all tasks will have the contents of the
! root array in their own arrays.  Thus 'array' is intent(in) on root, and
! intent(out) on all other tasks.

integer :: itemcount, datasize, errcode

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'array_broadcast', errstring, source, revision, revdate)
endif

! simple idiotproofing
if ((root < 0) .or. (root >= total_tasks)) then
   write(errstring, '(a,i8,a,i8)') "root task id ", root, &
                                   "must be >= 0 and < ", total_tasks
   call error_handler(E_ERR,'array_broadcast', errstring, source, revision, revdate)
endif

itemcount = size(array)
!print *, "kind = ", kind(array(1))
!print *, "digits = ", digits(array(1))
! TODO: FIX THIS based on what kind and digits say
datasize = MPI_DOUBLE_PRECISION

call MPI_Bcast(array, itemcount, datasize, root, my_local_comm, errcode)
if (errcode /= MPI_SUCCESS) then
   write(errstring, '(a,i8)') 'MPI_Bcast returned error code ', errcode
   call error_handler(E_ERR,'array_broadcast', errstring, source, revision, revdate)
endif

end subroutine array_broadcast


!-----------------------------------------------------------------------------
! TODO: do i need to overload this for both integer and real?
!       do i need to handle 2D inputs?

subroutine array_distribute(srcarray, root, dstarray, dstcount, how, which)
 real(r8), intent(in) :: srcarray(:)
 integer, intent(in) :: root
 real(r8), intent(out) :: dstarray(:)
 integer, intent(out) :: dstcount
 integer, intent(in) :: how
 integer, intent(out) :: which(:)

! 'srcarray' on the root task will be distributed across all the tasks
! into 'dstarray'.  dstarray must be large enough to hold each task's share
! of the data.  The actual number of values returned on each task will be
! passed back in the 'count' argument.  'how' is a flag to select how to
! distribute the data (round-robin, contiguous chunks, etc).  'which' is an
! integer index array which lists which of the original values were selected
! and put into 'dstarray'.

real(r8), allocatable :: localchunk(:)
integer :: srccount, datasize, leftover
integer :: i, tag, errcode
logical :: iamroot
integer :: status(MPI_STATUS_SIZE)

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'array_distribute', errstring, source, revision, revdate)
endif

! simple idiotproofing
if ((root < 0) .or. (root >= total_tasks)) then
   write(errstring, '(a,i8,a,i8)') "root task id ", root, &
                                   "must be >= 0 and < ", total_tasks
   call error_handler(E_ERR,'array_broadcast', errstring, source, revision, revdate)
endif

iamroot = (root == myrank)
tag = 1

srccount = size(srcarray)
print *, "kind = ", kind(srcarray(1))
print *, "digits = ", digits(srcarray(1))
! TODO: FIX THIS based on what kind and digits say
datasize = MPI_DOUBLE_PRECISION

! TODO: right now this code does contig chunks only
! TODO: it should select on the 'how' argument
dstcount = srccount / total_tasks
leftover = srccount - (dstcount * total_tasks)
if (myrank == total_tasks-1) dstcount = dstcount + leftover


! idiotproofing, continued...
if (size(dstarray) < dstcount) then
   write(errstring, '(a,i8,a,i8)') "size of dstarray is", size(dstarray), & 
                      " but must be >= ", dstcount
   call error_handler(E_ERR,'array_broadcast', errstring, source, revision, revdate)
endif
if (size(which) < dstcount) then
   write(errstring, '(a,i8,a,i8)') "size of which is", size(which), & 
                      " but must be >= ", dstcount
   call error_handler(E_ERR,'array_broadcast', errstring, source, revision, revdate)
endif

! TODO: this code is separate from the 'dstcount' computation because we
! need to test to be sure the user has passed us in arrays large enough to
! hold the data, but then this section needs to have a select (how) and set
! the corresponding index numbers accordingly.
which(1:dstcount) = (/ (i, i= myrank *dstcount, (myrank+1)*dstcount - 1) /)
if (size(which) > dstcount) which(dstcount+1:) = -1
   

if (.not.iamroot) then

   ! my task is receiving data.
   call MPI_Recv(dstarray, dstcount, datasize, root, MPI_ANY_TAG, &
	          my_local_comm, status, errcode)
   if (errcode /= MPI_SUCCESS) then
      write(errstring, '(a,i8)') 'MPI_Recv returned error code ', errcode
      call error_handler(E_ERR,'array_broadcast', errstring, source, revision, revdate)
   endif

else
   ! my task must send to everyone else and copy to myself.
   allocate(localchunk(dstcount), stat=errcode)  
   if (errcode /= 0) then
      write(errstring, *) 'allocation error of allocatable array'
      call error_handler(E_ERR,'array_broadcast', errstring, source, revision, revdate)
   endif

   do i=0, total_tasks-1
      ! copy correct data from srcarray to localchunk for each destination
      if (i == myrank) then
         ! this is my task, so do a copy from localchunk to dstarray
         dstarray(1:dstcount) = localchunk(1:dstcount)
      else
         ! call MPI to send the data to the remote task
         call MPI_Ssend(localchunk, dstcount, datasize, i, tag, &
   	               my_local_comm, errcode)
         if (errcode /= MPI_SUCCESS) then
            write(errstring, '(a,i8)') 'MPI_Ssend returned error code ', errcode
            call error_handler(E_ERR,'array_broadcast', errstring, source, revision, revdate)
         endif
      endif
      tag = tag + 1
   enddo

   deallocate(localchunk, stat=errcode)  
   if (errcode /= 0) then
      write(errstring, *) 'deallocation error of allocatable array'
      call error_handler(E_ERR,'array_broadcast', errstring, source, revision, revdate)
   endif
endif
   
! set any additional space which wasn't filled with zeros.
if (size(dstarray) > dstcount) dstarray(dstcount+1:) = 0.0

end subroutine array_distribute

!-----------------------------------------------------------------------------
! DART-specific cover utilities
!-----------------------------------------------------------------------------

function iam_task0()

! Return .TRUE. if my local task id is 0, .FALSE. otherwise.
! (Task numbers in MPI start at 0, contrary to the rules of polite fortran.)

logical :: iam_task0

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'iam_task0', errstring, source, revision, revdate)
endif

iam_task0 = (myrank == 0)

end function iam_task0

!-----------------------------------------------------------------------------
subroutine broadcast_send(from, array1, array2)
 integer, intent(in) :: from
 ! really only intent(in) here, but must match array_broadcast() call.
 real(r8), intent(inout) :: array1(:), array2(:)

! cover routine for array broadcast.  one additional sanity check -- make 
! sure the 'from' matches my local task id.  also, these arrays are
! intent(in) here, but they call a routine which is intent(inout) so they
! must be the same here.

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'broadcast_send', errstring, source, revision, revdate)
endif

! simple idiotproofing
if (from /= myrank) then
   write(errstring, '(a,i8,a,i8)') "'from' task id ", from, &
                                   "must be same as current task id ", myrank
   call error_handler(E_ERR,'broadcast_send', errstring, source, revision, revdate)
endif

! this must be paired with broadcast_recv() on all other tasks. 
! it will not return until all tasks in the communications group have
! made the call.
call array_broadcast(array1, from)
call array_broadcast(array2, from)

end subroutine broadcast_send

!-----------------------------------------------------------------------------
subroutine broadcast_recv(from, array1, array2)
 integer, intent(in) :: from
 ! really only intent(out) here, but must match array_broadcast() call.
 real(r8), intent(inout) :: array1(:), array2(:)

! cover routine for array broadcast.  one additional sanity check -- make 
! sure the 'from' is not the same as my local task id.  these arrays are
! intent(out) here, but they call a routine which is intent(inout) so they
! must be the same here.

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'broadcast_recv', errstring, source, revision, revdate)
endif

! simple idiotproofing
if (from == myrank) then
   write(errstring, '(a,i8,a,i8)') "'from' task id ", from, &
                                   "cannot be same as current task id ", myrank
   call error_handler(E_ERR,'broadcast_recv', errstring, source, revision, revdate)
endif

! this must be paired with a single broadcast_send() on the 'from' task.
! it will not return until all tasks in the communications group have
! made the call.
call array_broadcast(array1, from)
call array_broadcast(array2, from)

end subroutine broadcast_recv


!-----------------------------------------------------------------------------
! pipe utilities
!-----------------------------------------------------------------------------
!    * make_pipe()        Function that creates a named pipe (fifo), opens it,
!                         and returns the unit number.  Ok to call if the pipe
!                         already exists or is already open; it will skip
!                         those steps and just return the unit number.  The 
!                         name argument is used as a base and a filename
!                         in the form 'base.NNNN' is generated, where the N's
!                         are the MPI rank number, 0 padded.
!
!-----------------------------------------------------------------------------
!    * destroy_pipe()     The unit number is closed and the pipe file is 
!                         removed.
!
!-----------------------------------------------------------------------------
!    * read_pipe()        The character string is read from the pipe.
!                         (Can be overloaded to read ints if time or status
!                         info is useful to exchange between processes.) 
!                         This routine blocks until data is available.
!
!-----------------------------------------------------------------------------
!    * write_pipe()       The character string is written to the pipe.
!                         (Can be overloaded to write ints if time or status
!                         info is useful to exchange between processes.) 
!                         This routine writes and returns immediately.
!
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

end module mpi_utilities_mod

