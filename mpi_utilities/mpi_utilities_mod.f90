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
!      initialize_mpi_utilities()  Subroutine that initializes MPI and sets
!                                  local values needed later.  Must be called
!                                  before any other routine here.
!
!      finalize_mpi_utilities()  Subroutine that shuts down MPI cleanly.
!                                Must be called before program exits, and no
!                                other routines here can be used afterwards.
!
!      task_count()    Function that returns the total number of MPI tasks.
!
!      my_task_id()    Function that returns my task number.  Note that
!                      in the MPI world task numbers run from 0 to N-1.
!
!      transpose_array()  Subroutine which transposes a 2D array
!                         from column-major to row-major or back.
!        
!-----------------------------------------------------------------------------
! 
! these do not exist - i believe a single transpose will work.  but if not,
! they can be separated into these two:
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

!
! Some MPI installations have an MPI module; if one is present, use that.
! If not, there will be an MPI include file which defines the parameters.
! Use one but not both.  For more help on compiling a module which uses MPI
! see the $DART/doc/mpi directory.

!use mpi

implicit none
private

include "mpif.h"


!   ---- private data for mpi_utilities ----

integer :: myrank          ! my mpi number
integer :: total_tasks     ! total mpi tasks/procs
integer :: my_local_comm   ! duplicate communicator private to this file
!!integer :: comm_size     ! if ens count < tasks, only the first N participate

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


public :: total_tasks, my_task_id, transpose_array, &
          initialize_mpi_utilities, finalize_mpi_utilities

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

subroutine initialize_mpi_utilities()

integer :: errcode

if ( .not. module_initialized ) then
   ! Initialize the module with utilities
   call register_module(source, revision, revdate)
   module_initialized = .true.
endif

errcode = -999
call MPI_Init(errcode)
if (errcode /= MPI_SUCCESS) then
   write(errstring, '(a, i8)') 'MPI_Init returned error, code = ', errcode
   call error_handler(E_ERR,'initialize_mpi_utilities', errstring, source, revision, revdate)
endif

! duplicate the world communicator to isolate us from any other user
! calls to MPI.  All subsequent mpi calls will use the local communicator
! and not world.
call MPI_Comm_dup(MPI_COMM_WORLD, my_local_comm, errcode)
if (errcode /= MPI_SUCCESS) then
   write(errstring, '(a, i8)') 'MPI_Comm_dup returned error code ', errcode
   call error_handler(E_ERR,'initialize_mpi_utilities', errstring, source, revision, revdate)
endif

call MPI_Comm_rank(my_local_comm, myrank, errcode)
if (errcode /= MPI_SUCCESS) then
   write(errstring, '(a, i8)') 'MPI_Comm_rank returned error code ', errcode
   call error_handler(E_ERR,'initialize_mpi_utilities', errstring, source, revision, revdate)
endif

call MPI_Comm_size(my_local_comm, total_tasks, errcode)
if (errcode /= MPI_SUCCESS) then
   write(errstring, '(a, i8)') 'MPI_Comm_size returned error code ', errcode
   call error_handler(E_ERR,'initialize_mpi_utilities', errstring, source, revision, revdate)
endif

! MPI successfully initialized.

end subroutine initialize_mpi_utilities

!-----------------------------------------------------------------------------

subroutine finalize_mpi_utilities

integer :: errcode

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'finalize_mpi_utilities', errstring, source, revision, revdate)
endif

call MPI_Comm_free(my_local_comm, errcode)
if (errcode /= MPI_SUCCESS) then
   write(errstring, '(a, i8)') 'MPI_Comm_free returned error code ', errcode
   call error_handler(E_ERR,'finalize_mpi_utilities', errstring, source, revision, revdate)
endif

call MPI_Finalize(errcode)
if (errcode /= MPI_SUCCESS) then
   write(errstring, '(a, i8)') 'MPI_Finalize returned error code ', errcode
   call error_handler(E_ERR,'finalize_mpi_utilities', errstring, source, revision, revdate)
endif


end subroutine finalize_mpi_utilities


!-----------------------------------------------------------------------------

function task_count

! returns total number of MPI tasks.  e.g. if the number of tasks is 4,
! it returns 4.  the actual task numbers are 0-3.

integer :: task_count

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'task_count', errstring, source, revision, revdate)
endif

task_count = total_tasks

end function task_count


!-----------------------------------------------------------------------------

function my_task_id

integer :: my_task_id

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'my_task_id', errstring, source, revision, revdate)
endif

my_task_id = myrank

end function my_task_id


!-----------------------------------------------------------------------------

subroutine transpose_array


if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'transpose_array', errstring, source, revision, revdate)
endif

write(errstring, *) 'not implemented yet'
call error_handler(E_ERR,'transpose_array', errstring, source, revision, revdate)

end subroutine transpose_array


!-----------------------------------------------------------------------------

end module mpi_utilities_mod

