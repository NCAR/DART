! Data Assimilation Research Testbed -- DART
! Copyright 2004-2006, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module commtest_mod

! <next few lines automatically updated by version control software, do not edit>
! $Revision$
! $Date$
! $Id$
! 

implicit none
private

!! MPI -- message passing interface --
!!  if the module is available, use it.  otherwise, use the include file
!!  to define the MPI constants below.
!use mpi
include "mpif.h"

! subroutines exported to the driver program
public  setup, takedown, transpose_state_compl, transpose_ensemb_compl, &
        broadcast_increments, select_statevars, advance_model, &
        compute_priors, compute_qc, broadcast_obs
 

! module globals shared by most of the routines below.

! various sizes
integer :: myrank          ! my mpi number
integer :: total_tasks     ! total mpi tasks/procs
integer :: my_local_comm   ! copy of mpi world for now; dup eventually
integer :: comm_size       ! if ens count < tasks, only the first N participate
integer :: state_size      ! number of state vars in a single model
integer :: ens_size        ! number of models/ensembles running
integer :: obs_size        ! number of observations available
integer :: ens_per_task    ! if number of ensembles > tasks, how many per task
integer :: states_per_task ! state vars generally > tasks; how many per task
integer :: obs_per_task    ! not sure if this is needed; unused so far.
integer :: ec_total_size   ! total array size for ensemble-complete data
integer :: sc_total_size   ! total array size for state-complete data
integer :: max_print       ! limit for value dumps

! for simplicity, make these data areas global; they could end up being
! passed as arguments, or allocated, or ...
double precision, allocatable :: sc_vals(:)     ! (state_size * ens_per_task)
double precision, allocatable :: ec_vals(:)     ! (ens_size * states_per_task)
double precision, allocatable :: obs_vals(:)                   ! (obs_size)
double precision, allocatable :: stddev_vals(:)                ! (obs_size)
double precision, allocatable :: qc_vals(:)                    ! (obs_size)
double precision, allocatable :: prior_vals(:)                 ! (obs_size)
double precision, allocatable :: state_vals(:)                 ! (state_size)
double precision, allocatable :: single_increment_vals(:)      ! (ens_size)
double precision, allocatable :: single_prior_distrib_vals(:)  ! (ens_size)

contains

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
! called once at start to get MPI initialized and to set the sizes and
! module global numbers for later.
!-----------------------------------------------------------------------------
subroutine setup(state_vector_size, ensemble_size, observation_size)
   integer, intent(in) :: state_vector_size, ensemble_size, observation_size
   integer :: errcode

   print *, " "
   print *, "in setup"
 
   errcode = -999
   print *, "calling mpiinit"
   call MPI_Init(errcode)
   if (errcode /= MPI_SUCCESS) then
      print *, "init returned failure, errcode = ", errcode
      print *, "stopping on error"
      stop
   else
      print *, "init returned success"
   endif

   ! duplicate the world communicator to isolate us from any other user
   ! calls to MPI.
   call MPI_Comm_dup(MPI_COMM_WORLD, my_local_comm, errcode)
   if (errcode /= MPI_SUCCESS) then
      print *, "comm dup returned failure, errcode = ", errcode
      print *, "stopping on error"
      stop
   else
      print *, "comm dup returned success"
   endif
 
   errcode = -999
   print *, "calling commrank"
   call MPI_Comm_rank(my_local_comm, myrank, errcode)
   if (errcode /= MPI_SUCCESS) then
      print *, "commrank returned failure, errcode = ", errcode
      print *, "stopping on error"
      stop
   else
      print *, "commrank returned success, i am task ", myrank
   endif

   errcode = -999
   print *, "calling commsize"
   call MPI_Comm_size(my_local_comm, total_tasks, errcode)
   if (errcode /= MPI_SUCCESS) then
      print *, "commsize returned failure, errcode = ", errcode
      print *, "stopping on error"
      stop
   else
      print *, "commsize returned success, there are ", total_tasks, "tasks."
   endif

   print *, "MPI initialization completed successfully."
   print *, " "
   
   ! end of mpi setup.  now copy values into local variables for later.
   state_size = state_vector_size
   ens_size = ensemble_size
   obs_size = observation_size
      
   ! for now, refuse to run with more processors than there are ensembles.
   ! (this complicates the send/recv/broadcast loops; it can be handled
   ! but is not an interesting case for now).
   if (ens_size < total_tasks) then
      print *, "This version of the program will not run with a smaller number"
      print *, "of ensemble members than number of processors/tasks."
      print *, "Ensemble count, ", ens_size, &
               "is less than processor count, ", total_tasks
      print *, "stopping on error"
      stop
   endif

   ! this first test will not ever be true because above we have refused
   ! to continue in this case about 20 lines up.  but eventually we could
   ! support this by passing in a count of 0 for any communications involving
   ! tasks > ens size.
   if (total_tasks > ens_size) then
      comm_size = ens_size
      ens_per_task = 1
   else
      comm_size = total_tasks
      ens_per_task = (ens_size + total_tasks-1) / total_tasks
      ! also for now, refuse to run with an uneven number of ensembles per task.
      ! again, it should be handled eventually, but for now it complicates the
      ! loops unnecessarily for timing/testing.
      if ((ens_per_task * total_tasks) /= ens_size) then
         print *, "This version of the program will not run with an uneven"
         print *, "count of ensemble members per processors/tasks."
         print *, "Ensemble count, ", ens_size, &
                  "does not divide evenly across processor count, ", total_tasks
         print *, "stopping on error"
         stop
      endif
   endif
   print *, "each task computes ", ens_per_task, "ensembles"

   ! if ens_per_task is > 1, we have the option of either making the
   ! state complete arrays 2d:  data(state_size, ens_per_task) or leaving
   ! them 1d: data(state_size * ens_per_task).  precompute the total size and
   ! have the code handle it as 1d for now.
   sc_total_size = state_size * ens_per_task

   ! for now, refuse to run with more processors than there are state vars.
   ! (this complicates the send/recv/broadcast loops; it can be handled
   ! but is not an interesting case for now).
   if (state_size < total_tasks) then
      print *, "This version of the program will not run with a smaller number"
      print *, "of state variables than number of processors/tasks."
      print *, "State variables, ", state_size, &
               "is less than processor count, ", total_tasks
      print *, "stopping on error"
      stop
   endif

   ! this first test will not ever be true because above we have refused
   ! to continue in this case about 20 lines up.  but eventually we could
   ! support this by passing in a count of 0 for any communications involving
   ! tasks > state var size.
   if (total_tasks > state_size) then
      comm_size = state_size
      states_per_task = 1
   else
      comm_size = total_tasks
      states_per_task = (state_size + total_tasks-1) / total_tasks
      ! also for now, refuse to run with an uneven number of ensembles per task.
      ! again, it should be handled eventually, but for now it complicates the
      ! loops unnecessarily for timing/testing.
      if ((states_per_task * total_tasks) /= state_size) then
         print *, "This version of the program will not run with an uneven"
         print *, "count of state variables per processors/tasks."
         print *, "State Variable count, ", state_size, &
                  "does not divide evenly across processor count, ", total_tasks
         print *, "stopping on error"
         stop
      endif
   endif
   print *, "each task computes ", states_per_task, "state variables"

   ! if states_per_task is > 1, we have the option of either making the
   ! ensemble complete arrays 2d:  data(ens_size, states_per_task) or leaving
   ! them 1d: data(ens_size * states_per_task).  precompute the total size and
   ! have the code handle it as 1d for now.
   ec_total_size = ens_size * states_per_task

   ! limit printing of values from large arrays
   max_print = 100

   ! for the log, print out the values used.
   print *, " "
   print *, "debug dump: "
   print *, "myrank = ", myrank
   print *, "total_tasks = ", total_tasks
   print *, "my_local_comm = ", my_local_comm
   print *, "comm_size = ", comm_size
   print *, "state_size = ", state_size
   print *, "ens_size = ", ens_size
   print *, "obs_size = ", obs_size
   print *, "ens_per_task = ", ens_per_task
   print *, "states_per_task = ", states_per_task
   print *, "obs_per_task = ", obs_per_task
   print *, "ec_total_size = ", ec_total_size
   print *, "sc_total_size = ", sc_total_size
   print *, "max_print = ", max_print
   print *, " "

   ! now allocate the proper amount of space and use the same buffers
   ! for all the subsequent routines.
   allocate(sc_vals(state_size * ens_per_task), stat=errcode)
   if (errcode /= 0) stop
   allocate(ec_vals(ens_size * states_per_task), stat=errcode)
   if (errcode /= 0) stop
   allocate(obs_vals(obs_size), stat=errcode)
   if (errcode /= 0) stop
   allocate(stddev_vals(obs_size), stat=errcode)
   if (errcode /= 0) stop
   allocate(prior_vals(obs_size), stat=errcode)
   if (errcode /= 0) stop
   allocate(qc_vals(obs_size), stat=errcode)
   if (errcode /= 0) stop
   allocate(state_vals(state_size), stat=errcode)
   if (errcode /= 0) stop
   allocate(single_increment_vals(ens_size), stat=errcode)
   if (errcode /= 0) stop
   allocate(single_prior_distrib_vals(ens_size), stat=errcode)
   if (errcode /= 0) stop

   ! all the interesting values are now precomputed.
   print *, "setup finished successfully: i am", myrank, "of", total_tasks
   print *, " "

end subroutine setup

!-----------------------------------------------------------------------------
! shut down MPI cleanly.
!-----------------------------------------------------------------------------
subroutine takedown
   integer :: errcode

   print *, " "
   print *, "in takedown, iam ", myrank
  
   errcode = -999
   print *, "calling comm free"
   call MPI_Comm_free(my_local_comm, errcode)
   if (errcode /= MPI_SUCCESS) then
      print *, "comm free returned error, continuing.  code = ", errcode
   endif

   print *, "calling finalize"
   call MPI_Finalize(errcode)
   if (errcode /= MPI_SUCCESS) then
      print *, "finalize returned error, continuing.  code = ", errcode
   endif

   ! clean up some.
   deallocate(sc_vals, stat=errcode)
   if (errcode /= 0) stop
   deallocate(ec_vals, stat=errcode)
   if (errcode /= 0) stop
   deallocate(obs_vals, stat=errcode)
   if (errcode /= 0) stop
   deallocate(stddev_vals, stat=errcode)
   if (errcode /= 0) stop
   deallocate(prior_vals, stat=errcode)
   if (errcode /= 0) stop
   deallocate(qc_vals, stat=errcode)
   if (errcode /= 0) stop
   deallocate(state_vals, stat=errcode)
   if (errcode /= 0) stop
   deallocate(single_increment_vals, stat=errcode)
   if (errcode /= 0) stop
   deallocate(single_prior_distrib_vals, stat=errcode)
   if (errcode /= 0) stop

   print *, " "

end subroutine takedown

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
! "read in" obs on them master and broadcast them.
!-----------------------------------------------------------------------------
subroutine broadcast_obs(observation_size)
   integer, intent(in) :: observation_size

   integer :: errcode, root

   print *, " "
   print *, "in broadcast_obs, iam ", myrank
  
   root = 0
   if (myrank == root) then
      obs_vals(:) = 3.14159
   else
      obs_vals(:) = -666.0
   endif

   errcode = -999
   print *, "broadcasting observations"
   call MPI_Bcast(obs_vals, obs_size, MPI_DOUBLE_PRECISION, &
                  root, my_local_comm, errcode)
   if (errcode /= MPI_SUCCESS) then
      print *, "broadcast returned error, continuing.  code = ", errcode
      stop
   endif

   print *, " "

end subroutine broadcast_obs

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
! decide which state vars get collected together across the ensembles.
! the variations are on which state vars and observations are chosen.
!-----------------------------------------------------------------------------
subroutine select_statevars(state_vector_size, obs_size, index_list, algorithm)
   integer, intent(in) :: state_vector_size, obs_size
   ! intent(inout)
   integer, pointer :: index_list(:)
   integer, intent(in) :: algorithm

   if (algorithm == 1) then
      call select1(state_vector_size, index_list, obs_size)
   else if (algorithm == 2) then
      call select2(state_vector_size, index_list, obs_size)
   else
     print *, "bad value for select_statevar algorithm selection, ", algorithm
     print *, "stopping on error"
     stop
   endif

   call array_print("index selected for each observation/prior is: ", &
                iarray=index_list)

end subroutine select_statevars

!-----------------------------------------------------------------------------
! pick random values, which for efficiency means selecting contiguous blocks
! from the state vector and observations.
!-----------------------------------------------------------------------------
subroutine select1(state_vector_size, index_list, obs_size)
   integer, intent(in) :: state_vector_size, obs_size
   ! intent(inout)
   integer, pointer :: index_list(:)
 
   integer :: i, errcode, region

   print *, " "
   print *, "in select_1, iam ", myrank

   allocate(index_list(obs_size), stat=errcode)
   ! TODO: check error code

   region = obs_size / total_tasks

   ! agree that index_list() is 0 based, just like MPI task numbers
   do i=0, obs_size-1
      ! contig blocks
      index_list(i+1) = i / region
   enddo
   
   ! L contains the index values 
   print *, "Filled L with blocks of index values"
   print *, " "

end subroutine select1

!-----------------------------------------------------------------------------
! pick particular values, which means constructing an index list to use to
! select values from the state vector and observations.
!-----------------------------------------------------------------------------
subroutine select2(state_vector_size, index_list, obs_size)
   integer, intent(in) :: state_vector_size, obs_size
   ! intent(inout)
   integer, pointer :: index_list(:)

   integer :: i, errcode

   print *, " "
   print *, "in select_2, iam ", myrank

   allocate(index_list(obs_size), stat=errcode)
   ! TODO: check error code

   ! index_list numbers are 0 based, just like MPI task numbers
   do i=1, obs_size
      ! ok, not very random.  mix this up at some point.
      index_list(i) = mod(i, total_tasks)
   enddo
   
   ! L contains the index values 
   print *, "Filled L with cycling index values"
   print *, " "

end subroutine select2

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
! transpose from state complete to ensemble complete.
! the variations are in which MPI calls are used, and whether they are
! synchronous or asynchronous; in the end all the bits should be in the
! same places no matter which algorithm is specified.  it is a question of
! which method is fastest.
!-----------------------------------------------------------------------------
subroutine transpose_ensemb_compl(state_vector_size, ensemble_size, &
                                  index_list,  priors, qc, algorithm)
   integer, intent(in) :: state_vector_size, ensemble_size
   ! intent(in)
   integer, pointer :: index_list(:)
   double precision, pointer :: priors(:), qc(:)
   integer, intent(in) :: algorithm

   if (algorithm == 1) then
      call transpose1(state_vector_size, ensemble_size, index_list, priors, qc)
   else if (algorithm == 2) then
      call transpose2(state_vector_size, ensemble_size, index_list, priors, qc)
   else if (algorithm == 3) then
      call transpose3(state_vector_size, ensemble_size, index_list, priors, qc)
   else
     print *, "bad value for ensemb_compl algorithm selection, ", algorithm
     print *, "stopping on error"
     stop
   endif

   call array_print("ensemble_complete values are: ", darray=ec_vals)

end subroutine transpose_ensemb_compl

!-----------------------------------------------------------------------------
! simplest send/irecv version which will always work no matter what the
! size of the mpi buffers:  first call all non-blocking receives, and
! then call blocking sends.   then wait for all the receives to complete,
! one by one.   this is sending the right *number* of items, but not the
! right items (it is not doing a column-to-row transpose.)
!-----------------------------------------------------------------------------
subroutine transpose1(state_vector_size, ensemble_size, index_list, priors, qc)
   integer, intent(in) :: state_vector_size, ensemble_size
   ! intent(in)
   integer, pointer :: index_list(:)
   double precision, pointer :: priors(:), qc(:)

   ! TODO: consider whether making these be real 2d arrays would make
   ! the code simpler or more complicated.
   integer :: errcode, send_start, recv_start
   integer :: i, req(total_tasks)
   integer :: status(MPI_STATUS_SIZE)

   print *, " "
   print *, "in transpose_1, iam ", myrank

   ! the 9.0 intel compiler does not recognize real() as an intrinsic?
   !sc_vals = real(myrank)
   sc_vals = myrank * 1.0
   ec_vals = -1.0

   ! mpi rank numbers are 0 based, not 1 based.
   do i=0, total_tasks-1
      !print *, "ireceive from ", i
      recv_start = (i * states_per_task) + 1
      call MPI_Irecv(ec_vals(recv_start), ens_size, MPI_DOUBLE_PRECISION, &
                     i, MPI_ANY_TAG, my_local_comm, req(i+1), errcode)
      print *, "done with ireceive from ", i, " errcode = ", errcode
      if (errcode /= MPI_SUCCESS) stop
   enddo

   do i=0, total_tasks-1
      !print *, "sending to ", i
      send_start = (i * states_per_task) + 1
      call MPI_Send(sc_vals(send_start), ens_size, MPI_DOUBLE_PRECISION, &
                    i, 1, my_local_comm, errcode)
      print *, "done with sending to ", i, "errcode = ", errcode
      if (errcode /= MPI_SUCCESS) stop
   enddo
 
   do i=0, total_tasks-1
      !print *, "waiting for request num ", i
      call MPI_Wait(req(i+1), status, errcode)
      print *, "done with wait for ", i, "status = ", status
      if (errcode /= MPI_SUCCESS) stop
   enddo

end subroutine transpose1

!-----------------------------------------------------------------------------
! call all the non-blocking receives, then the non-blocking sends, then wait
! for all to complete at the end in a single combined wait call.
! this is sending the right *number* of items, but not the right items.
! (it is not doing a column-to-row transpose.)
!-----------------------------------------------------------------------------
subroutine transpose2(state_vector_size, ensemble_size, index_list, priors, qc)
   integer, intent(in) :: state_vector_size, ensemble_size
   ! intent(in)
   integer, pointer :: index_list(:)
   double precision, pointer :: priors(:), qc(:)

   integer :: errcode, send_start, recv_start
   integer :: i, req(total_tasks*2)
   integer :: status_array(MPI_STATUS_SIZE, total_tasks*2)

   print *, " "
   print *, "in transpose_2, iam ", myrank

   sc_vals = myrank * 1.0
   ec_vals = -1.0

   do i=0, total_tasks-1
      !print *, "ireceive from ", i
      recv_start = (i * states_per_task) + 1
      call MPI_Irecv(ec_vals(recv_start), ens_size, MPI_DOUBLE_PRECISION, &
                     i, MPI_ANY_TAG, my_local_comm, req(i+1), errcode)
      print *, "done with ireceive from ", i, " errcode = ", errcode
      if (errcode /= MPI_SUCCESS) stop
   enddo

   do i=0, total_tasks-1
      !print *, "sending to ", i
      send_start = (i * states_per_task) + 1
      call MPI_Isend(sc_vals(send_start), ens_size, MPI_DOUBLE_PRECISION, &
                    i, i, my_local_comm, req(total_tasks+i+1), errcode)
      print *, "done with sending to ", i, "errcode = ", errcode
      if (errcode /= MPI_SUCCESS) stop
   enddo
 
   call MPI_Waitall(total_tasks*2, req, status_array, errcode)
   print *, "back from waitall, iam ", myrank
   if (errcode /= MPI_SUCCESS) stop

end subroutine transpose2

!-----------------------------------------------------------------------------
! create an MPI vector type which has the strides set to copy a contiguous
! row from the source and place it in a strided column in the destination.
! this version is finally moving both the right number of items and putting
! them into the correct locations.  this routine handles multiple ensembles
! per task correctly (it does require the ensemble count divides evenly
! across the number of tasks.)
!-----------------------------------------------------------------------------
subroutine transpose3(state_vector_size, ensemble_size, index_list, priors, qc)
   integer, intent(in) :: state_vector_size, ensemble_size
   ! intent(in)
   integer, pointer :: index_list(:)
   double precision, pointer :: priors(:), qc(:)

   integer :: errcode
   integer :: i, j, req(ens_size*2)
   integer :: send_start, recv_start
   integer :: status_array(MPI_STATUS_SIZE, ens_size*2)
   integer :: stridetype, stride, region_size, base, offset

   print *, " "
   print *, "in transpose_3, iam ", myrank

   stride = ens_size
   region_size = state_size/total_tasks
  print *, "t3: region_size, stride = ", region_size, stride

   ! assign cell values so the 6th digit is the ens number; the ones
   ! are the state var number inside an ens.
   do j=1, ens_per_task
      do i=1, state_size
        sc_vals((j-1)*state_size + i) = (((myrank*ens_per_task)+j) &
                                         * 100000.0) + i
      enddo
   enddo
   ec_vals = -1.0

   ! create the mpi vector type with the striding.
   call MPI_Type_vector(region_size, 1, stride, MPI_DOUBLE_PRECISION, &
                        stridetype, errcode)
   call MPI_Type_commit(stridetype, errcode)

   base = 0
   do j=0, ens_per_task-1
      do i=0, total_tasks-1
         !print *, "ireceive from ", i
         offset = i + (j*total_tasks) + 1
         recv_start = offset
      print *, "t3:", i, j, "posting receive with starting offset ", &
                recv_start, " req index ", offset
         call MPI_Irecv(ec_vals(recv_start), 1, stridetype, i, MPI_ANY_TAG, &
                        my_local_comm, req(offset), errcode)
         print *, myrank, "done with ireceive from ", i, " errcode = ", errcode
         if (errcode /= MPI_SUCCESS) stop
      enddo
   enddo

   base = ens_size
   do j=0, ens_per_task-1
      do i=0, total_tasks-1
         !print *, "isending to ", i
         send_start = (j*state_size)+(i*region_size) + 1
         offset = i + (j*total_tasks) + 1
      print *, "t3:", i, j, "send with starting offset ", &
                send_start, " req index ", offset
         call MPI_Isend(sc_vals(send_start), region_size, MPI_DOUBLE_PRECISION, &
                        i, offset, my_local_comm, req(base+offset), errcode)
         print *, myrank, "done with isend to ", i, " errcode = ", errcode
         if (errcode /= MPI_SUCCESS) stop
      enddo
   enddo
 
   print *, "all sends and receives have been posted, starting to Wait"
   call MPI_Waitall(ens_size*2, req, status_array, errcode)
   print *, "back from waitall, iam ", myrank
   if (errcode /= MPI_SUCCESS) stop

   ! TODO: also send priors, qc, and selected subset of obs.

   call MPI_Type_free(stridetype, errcode)
   if (errcode /= MPI_SUCCESS) stop

   print *, " "

end subroutine transpose3

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
! transpose from ensemble complete to state complete.  as before, the
! variations are speed only; the same bits should arrive in the same
! places no matter which one is used.
!-----------------------------------------------------------------------------
subroutine transpose_state_compl(state_vector_size, ensemble_size, &
                                 index_list, priors, qc, algorithm)
   integer, intent(in) :: state_vector_size, ensemble_size
   ! intent(in)
   integer, pointer :: index_list(:)
   double precision, pointer :: priors(:), qc(:)
   integer, intent(in) :: algorithm

   !if (algorithm == 1) then
   !  call transpose1a(state_vector_size, ensemble_size, index_list, priors, qc)
   !else if (algorithm == 2) then
   !  call transpose2a(state_vector_size, ensemble_size, index_list, priors, qc)

   if (algorithm == 3) then
      call transpose3a(state_vector_size, ensemble_size, index_list, priors, qc)
   else
     print *, "bad value for state_compl algorithm selection, ", algorithm
     print *, "stopping on error"
     stop
   endif

   call array_print("state_complete values are: ", darray=sc_vals)

end subroutine transpose_state_compl

!-----------------------------------------------------------------------------
! create an MPI vector type which has the strides set to copy a strided column
! from the source and place it in a contiguous row in the destination.
! this version is finally moving both the right number of items and putting
! them into the correct locations.  this routine handles multiple ensembles
! per task correctly (it does require the ensemble count divides evenly
! across the number of tasks.)
!-----------------------------------------------------------------------------
subroutine transpose3a(state_vector_size, ensemble_size, index_list, priors, qc)
   integer, intent(in) :: state_vector_size, ensemble_size
   ! intent(in)
   integer, pointer :: index_list(:)
   double precision, pointer :: priors(:), qc(:)

   integer :: errcode
   integer :: i, j, req(ens_size*2)
   integer :: send_start, recv_start
   integer :: status_array(MPI_STATUS_SIZE, ens_size*2)
   integer :: stridetype, stride, region_size, base, offset

   print *, " "
   print *, "in transpose_3a, iam ", myrank

   stride = ens_size
   region_size = state_size/total_tasks
  print *, "t3a: region_size, stride = ", region_size, stride

   ! assign cell values so the 6th digit is the state number; the ones
   ! are the ens number.
   do j=1, states_per_task
      do i=1, ens_size
         ec_vals((j-1)*ens_size + i) = (((myrank*states_per_task)+j) &
                                       * 100000.0) + i
      enddo
   enddo
   sc_vals = -1.0

   ! create the mpi vector type with the striding.
   call MPI_Type_vector(region_size, 1, stride, MPI_DOUBLE_PRECISION, &
                        stridetype, errcode)
   call MPI_Type_commit(stridetype, errcode)

   base = 0
   do j=0, ens_per_task-1
      do i=0, total_tasks-1
         !print *, "ireceive from ", i
         recv_start = (j*state_size)+(i*region_size) + 1
         offset = i + (j*total_tasks) + 1
      print *, "t3a:", i, j, "posting ireceive with starting offset ", &
                recv_start, " req index ", offset
         call MPI_Irecv(sc_vals(recv_start), region_size, MPI_DOUBLE_PRECISION, &
                        i, MPI_ANY_TAG, my_local_comm, req(offset), errcode)
         print *, myrank, "done with irecv to ", i, " errcode = ", errcode
         if (errcode /= MPI_SUCCESS) stop
      enddo
   enddo

   base = ens_size
   do j=0, ens_per_task-1
      do i=0, total_tasks-1
         !print *, "isend to ", i
         offset = i + (j*total_tasks) + 1
         send_start = offset
      print *, "t3a:", i, j, "isend with starting offset ", &
                send_start, " req index ", offset
         call MPI_Isend(ec_vals(send_start), 1, stridetype, i, offset, &
                        my_local_comm, req(base+offset), errcode)
         print *, myrank, "done with ireceive from ", i, " errcode = ", errcode
         if (errcode /= MPI_SUCCESS) stop
      enddo
   enddo
 
   print *, "all sends and receives have been posted, starting to Wait"
   call MPI_Waitall(ens_size*2, req, status_array, errcode)
   print *, "back from waitall, iam ", myrank
   if (errcode /= MPI_SUCCESS) stop

   ! TODO: also send priors, qc, and selected subset of obs.

   call MPI_Type_free(stridetype, errcode)
   if (errcode /= MPI_SUCCESS) stop

   print *, " "

end subroutine transpose3a


!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
! broadcast the increments for the priors, plus the observation distribution.
! the variations in algorithm are which MPI calls to use; the results should
! be the same.
!-----------------------------------------------------------------------------
subroutine broadcast_increments(index, obs_count, index_list, ensemble_size, &
                                algorithm)
   integer, intent(in) :: index, obs_count
   ! intent(in)
   integer, pointer :: index_list(:)
   integer, intent(in) :: ensemble_size, algorithm

   if (algorithm == 1) then
      call broadcast1a(index, obs_count, index_list, ensemble_size)
   else if (algorithm == 2) then
      call broadcast2a(index, obs_count, index_list, ensemble_size)
   else if (algorithm == 3) then
      call broadcast3a(index, obs_count, index_list, ensemble_size)
   else
     print *, "bad value for broadcast algorithm selection, ", algorithm
     print *, "stopping on error"
     stop
   endif

   call array_print("broadcast increments are:", darray=single_increment_vals)
 
end subroutine broadcast_increments

!-----------------------------------------------------------------------------
! in this version, everyone has the index list and so everyone knows who
! the root is and can use the MPI Broadcast function.
!-----------------------------------------------------------------------------
subroutine broadcast1a(index, obs_count, index_list, ensemble_size)
   integer, intent(in) :: index, obs_count, ensemble_size
   ! intent(in)
   integer, pointer :: index_list(:)

   integer :: i, errcode, root

   print *, " "
   print *, "in broadcast1a, iam ", myrank


   ! can put timing loop here if needed

   ! if my observation, i am root
   root = index_list(index)
   print *, "b1a: index = ", index, "list = ", index_list(index), "root = ", root

   if (root == myrank) then
      single_increment_vals(:) = (myrank+1) * 1.1
      single_prior_distrib_vals(:) = 0.5
   else
      single_increment_vals(:) = -1
      single_prior_distrib_vals(:) = -1
   endif

   ! TODO: this moves the right number of bytes, but does not have an
   ! offset or multiple arrays for each set of increments.
   do i=1,ens_per_task
      print *, "b1a: ",  myrank, &
               "calling broadcast for increments and distribution, root = ", &
               root
      call MPI_Bcast(single_increment_vals, ens_size, &
                     MPI_DOUBLE_PRECISION, root, my_local_comm, errcode)
      if (errcode /= MPI_SUCCESS) stop
   
      call MPI_Bcast(single_prior_distrib_vals, ens_size, &
                     MPI_DOUBLE_PRECISION, root, my_local_comm, errcode)
      print *, "b1a: ", myrank, "back from broadcast"
      if (errcode /= MPI_SUCCESS) stop

   enddo

end subroutine broadcast1a

!-----------------------------------------------------------------------------
! in this version, everyone has the index list and so everyone knows who
! the next root is, but everyone is also computing at the same time.  so to
! maximize the compute/communication overlap time, each "not next root" will
! post a receive and then continue on computing until the "next root" has the
! value ready and sends it to everyone else.   everyone should know if they
! are the "next root" and expedite sending the values.
!-----------------------------------------------------------------------------
subroutine broadcast2a(index, obs_count, index_list, ensemble_size)
   integer, intent(in) :: index, obs_count, ensemble_size
   ! intent(in)
   integer, pointer :: index_list(:)

   integer :: i, j, k, errcode, root
   logical :: iam_root, done
   integer :: req(ens_size*2), tag
   integer :: status_array(MPI_STATUS_SIZE, ens_size*2)

   print *, " "
   print *, "in broadcast2a, iam ", myrank


   ! can put timing loop here if needed
   k = 0

   ! if my observation, i am root
   root = index_list(index)
   print *, "b2a: index = ", index, "list = ", index_list(index), "root = ", root

   if (root == myrank) then
      single_increment_vals(:) = (myrank+1) * 1.1
      single_prior_distrib_vals(:) = 0.5
      iam_root = .TRUE.
   else
      single_increment_vals(:) = -1
      single_prior_distrib_vals(:) = -1
      iam_root = .FALSE.
   endif

   ! we aren't using the tags right now but generate something which would
   ! be unique if it turns out to be necessary.  the 2 below is number of
   ! messages each send generates; if that changes it needs to be bumped up.
   tag = index * 2 * total_tasks

   ! if i am not root, post a receive and go about computing things.
   if (.not.iam_root) then

      ! these are non-blocking so they return immediately; they do not
      ! wait until the data arrives.
      do i=1,ens_per_task
         print *, "b2a: ", myrank, &
                  "calling irecv for increments and distribution, root = ", &
                  root
         call MPI_Irecv(single_increment_vals, ens_size, &
                        MPI_DOUBLE_PRECISION, root, MPI_ANY_TAG, &
                        my_local_comm, req((2*i)-1), errcode)
         if (errcode /= MPI_SUCCESS) stop
      
         call MPI_Irecv(single_prior_distrib_vals, ens_size, &
                        MPI_DOUBLE_PRECISION, root, MPI_ANY_TAG, &
                        my_local_comm, req(2*i), errcode)
         print *, "b2a: ", myrank, "back from irecvs"
         if (errcode /= MPI_SUCCESS) stop

      enddo
    
      done = .FALSE.
      do while (.not.done)
         ! alternate computing and checking to see if the next value is here

         ! compute here
         k = k + 2

         ! test to see if the data is back but if not, do not wait.  keep
         ! computing until you have to wait.  this version of the routine
         ! waits on the next value (which requires 2 receives).   if we 
         ! want to wait on multiple values, we can use Testsome() instead
         ! of Testall() and process them as they arrive.
         call MPI_Testall(ens_per_task, req, done, status_array, errcode)

         ! if all done, testall sets 'done' flag to true
         if (done) print *, "b2a: ", myrank, "all increments received"

         ! but prevent an infinite loop here on error
         if (errcode /= MPI_SUCCESS) done = .TRUE.

      enddo

   else  
      ! i am root this time -- compute my values and send them to the others.
      ! this is using a sync send -- it will not return until everyone
      ! else has received the data.  but since they should have already
      ! posted async receives, this should return quickly.

      ! TODO: this moves the right number of bytes, but does not have an
      ! offset or multiple arrays for each set of increments.
      do i=1,ens_per_task
         print *, "b2a: ", myrank, &
                  "calling send for increments and distribution, root = ", &
                  root
         do j=0,total_tasks-1
            ! skip ourselves
            if (j == myrank) cycle

            call MPI_Send(single_increment_vals, ens_size, &
                           MPI_DOUBLE_PRECISION, j, tag, my_local_comm, errcode)
            if (errcode /= MPI_SUCCESS) stop
         
            call MPI_Send(single_prior_distrib_vals, ens_size, &
                         MPI_DOUBLE_PRECISION, j, tag+1, my_local_comm, errcode)
            print *, "b2a: ", myrank, "back from sends"
            if (errcode /= MPI_SUCCESS) stop

            tag = tag + 2
         enddo
   
      enddo
   endif


end subroutine broadcast2a

!-----------------------------------------------------------------------------
! in this version, everyone has the index list and so everyone knows who
! the next root is, but everyone is also computing at the same time. 
! to maximize the compute/communication overlap time, each "not next root" will
! post a receive and then continue on computing until the "next root" has the
! value ready and sends it to everyone else.   everyone should know if they
! are the "next root" and expedite sending the values.  this differs from 2a
! above in that it tries to allow multiple outstanding receives per ensemble,
! so the send is async as well and this routine can return without the comms
! being done.  (this is a bit awkward in that the index number comes in as an
! argument instead of looping inside this routine -- that may have to change.)
!-----------------------------------------------------------------------------
subroutine broadcast3a(index, obs_count, index_list, ensemble_size)
   integer, intent(in) :: index, obs_count, ensemble_size
   ! intent(in)
   integer, pointer :: index_list(:)

   integer :: i, j, k, errcode, root
   logical :: iam_root, done
   integer :: max_reqs, inuse, tag, nextindex
   integer, allocatable :: req(:)
   integer, allocatable :: status_array(:, :)

   print *, " "
   print *, "in broadcast3a, iam ", myrank

   ! get the right amount of space for pending requests
   ! max_reqs = ens_size * 2 * obs_count
   ! will be that, but for now use:
   max_reqs = ens_size * 2 * total_tasks
   allocate(req(max_reqs), status_array(MPI_STATUS_SIZE, max_reqs), &
            stat=errcode)
   if (errcode /= 0) stop

   ! can put timing loop here if needed
   k = 0

   ! if my observation, i am root
   root = index_list(index)
   print *, "b3a: index = ", index, "list = ", index_list(index), "root = ", root

   if (root == myrank) then
      single_increment_vals(:) = (myrank+1) * 1.1
      single_prior_distrib_vals(:) = 0.5
      iam_root = .TRUE.
   else
      single_increment_vals(:) = -1
      single_prior_distrib_vals(:) = -1
      iam_root = .FALSE.
   endif

   ! we aren't using the tags right now but generate something which would
   ! be unique if it turns out to be necessary.  the 2 below is number of
   ! messages each send generates; if that changes it needs to be bumped up.
   tag = index * 2 * total_tasks

   ! if i am not root, post a receive and go about computing things.
   if (.not.iam_root) then

      nextindex = 0

      ! these are non-blocking so they return immediately; they do not
      ! wait until the data arrives.
      do i=1,ens_per_task
         print *, "b3a: ", myrank, &
                  "calling irecv for increments and distribution, root = ", &
                  root
         nextindex = nextindex + 1
         call MPI_Irecv(single_increment_vals, ens_size, &
                        MPI_DOUBLE_PRECISION, root, MPI_ANY_TAG, &
                        my_local_comm, req(nextindex), errcode)
         if (errcode /= MPI_SUCCESS) stop
      
         nextindex = nextindex + 1
         call MPI_Irecv(single_prior_distrib_vals, ens_size, &
                        MPI_DOUBLE_PRECISION, root, MPI_ANY_TAG, &
                        my_local_comm, req(nextindex), errcode)
         print *, "b3a: ", myrank, "back from irecvs"
         if (errcode /= MPI_SUCCESS) stop

      enddo
    
      inuse = nextindex

      done = .FALSE.
      do while (.not.done)
         ! alternate computing and checking to see if the next value is here

         ! compute here
         k = k + 2

         ! test to see if the data is back but if not, do not wait.  keep
         ! computing until you have to wait.  this version of the routine
         ! waits on the next value (which requires 2 receives).   if we 
         ! want to wait on multiple values, we can use Testsome() instead
         ! of Testall() and process them as they arrive.
         call MPI_Testall(inuse, req, done, status_array, errcode)

         ! if all done, testall sets 'done' flag to true
         if (done) print *, "b3a: ", myrank, "testall says all receives done"

         ! but prevent an infinite loop here on error
         if (errcode /= MPI_SUCCESS) done = .TRUE.

      enddo

   else  
      ! i am root this time -- compute my values and send them to the others.
      ! this is using a sync send -- it will not return until everyone
      ! else has received the data.  but since they should have already
      ! posted async receives, this should return quickly.

      nextindex = 0

      ! TODO: this moves the right number of bytes, but does not have an
      ! offset or multiple arrays for each set of increments.
      do i=1,ens_per_task
         print *, "b3a: ", myrank, &
                  "calling send for increments and distribution, root = ", &
                  root
         do j=0,total_tasks-1
            ! don't send to ourselves
            if (j == myrank) cycle

            nextindex = nextindex + 1
            call MPI_Isend(single_increment_vals, ens_size, &
                           MPI_DOUBLE_PRECISION, j, tag, my_local_comm, &
                           req(nextindex), errcode)
            if (errcode /= MPI_SUCCESS) stop
         
            nextindex = nextindex + 1
            call MPI_Isend(single_prior_distrib_vals, ens_size, &
                          MPI_DOUBLE_PRECISION, j, tag+1, my_local_comm, &
                          req(nextindex), errcode)
            print *, "b3a: ", myrank, "back from isends"
            if (errcode /= MPI_SUCCESS) stop

            tag = tag + 2
         enddo
   
      enddo

      inuse = nextindex

      ! now wait here -- if we move the index loop into this routine, then
      ! we can keep looping and overlap more.
      call MPI_Waitall(inuse, req, status_array, errcode)
      if (errcode /= MPI_SUCCESS) stop

      print *, "b3a: ", myrank, "waitall says all sends done"

   endif

   ! get rid of temp space
   deallocate(req, status_array, stat=errcode)
   if (errcode /= 0) stop

end subroutine broadcast3a

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
! for timing only.   do something to simulate the amount of time a model
! might take doing a real computation.
!-----------------------------------------------------------------------------
subroutine advance_model(state_vector_size, ensemble_size)
   integer, intent(in) :: state_vector_size, ensemble_size

   integer :: i, j
   real :: val

   val = 0.0
   do j=1, ensemble_size
     do i=1, state_vector_size
       val = val + i + j 
     enddo
   enddo

end subroutine advance_model

!-----------------------------------------------------------------------------
! this can be used for timing, but also does something real - allocates
! an array where the "priors" can be stored.
!-----------------------------------------------------------------------------
subroutine compute_priors(state_vector_size, obs_size, ensemble_size, priors)
   integer, intent(in) :: state_vector_size, obs_size, ensemble_size
   ! intent(inout)
   double precision, pointer :: priors(:)

   integer :: i, j
   real :: val

   val = 0.0
   do j=1, ensemble_size
     do i=1, state_vector_size
       val = val + i + j + obs_size
     enddo
   enddo

   prior_vals(:) = 100.0

end subroutine compute_priors

!-----------------------------------------------------------------------------
! same as above; can be used for timing, but also allocates a QC array.
!-----------------------------------------------------------------------------
subroutine compute_qc(obs_size, ensemble_size, quality_control)
   integer, intent(in) :: obs_size, ensemble_size
   ! intent(inout)
   double precision, pointer :: quality_control(:)

   integer :: i, j
   real :: val

   val = 0.0
   do j=1, ensemble_size
     do i=1, obs_size
       val = val + i + j
     enddo
   enddo

   qc_vals(:) = 1.0

end subroutine compute_qc

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
! util routines - print from the given count up to max_print values.
! one for integers, one for doubles
!-----------------------------------------------------------------------------
subroutine array_print(label, iarray, darray)
   character(len=*), intent(in) :: label
   integer, intent(in), optional :: iarray(:)
   double precision, intent(in), optional :: darray(:)

   integer :: i, to_print, count

   if (present(iarray)) then
      count = size(iarray)
   else if (present(darray)) then
      count = size(darray)
   else
      count = 0
   endif
  
   if (count < max_print) then
      to_print = count
   else
      to_print = max_print
   endif

   print *, " "
   print *, "iam ", myrank, label
   if (present(iarray)) then
      print "(10i7 /)", (iarray(i),i=1,to_print)
   else if (present(darray)) then
      print "(10f9.0 /)", (darray(i),i=1,to_print)
   else
      print *, "<no data>"
   endif
   print *, " "

end subroutine array_print

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

end module commtest_mod

