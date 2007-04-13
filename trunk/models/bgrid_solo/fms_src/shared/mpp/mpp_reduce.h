!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                   !!
!!                   GNU General Public License                      !!
!!                                                                   !!
!! This file is part of the Flexible Modeling System (FMS).          !!
!!                                                                   !!
!! FMS is free software; you can redistribute it and/or modify       !!
!! it and are expected to follow the terms of the GNU General Public !!
!! License as published by the Free Software Foundation.             !!
!!                                                                   !!
!! FMS is distributed in the hope that it will be useful,            !!
!! but WITHOUT ANY WARRANTY; without even the implied warranty of    !!
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     !!
!! GNU General Public License for more details.                      !!
!!                                                                   !!
!! You should have received a copy of the GNU General Public License !!
!! along with FMS; if not, write to:                                 !!
!!          Free Software Foundation, Inc.                           !!
!!          59 Temple Place, Suite 330                               !!
!!          Boston, MA  02111-1307  USA                              !!
!! or see:                                                           !!
!!          http://www.gnu.org/licenses/gpl.txt                      !!
!!                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine MPP_REDUCE_( a, pelist )
!find the max of scalar a the PEs in pelist (all PEs if this argument is omitted)
!result is also automatically broadcast to all PEs
      MPP_TYPE_, intent(inout) :: a
      integer, intent(in), optional :: pelist(0:)
      integer :: n
#ifdef use_libSMA
!work holds pWrk array + 1 word for symmetric copy of a
      MPP_TYPE_ :: work(SHMEM_REDUCE_MIN_WRKDATA_SIZE+1)
      pointer( ptr, work )
      integer :: words
      character(len=8) :: text
#endif use_libSMA
#ifdef use_libMPI
      MPP_TYPE_ :: work
#endif
      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_REDUCE: You must first call mpp_init.' )
      n = get_peset(pelist); if( peset(n)%count.EQ.1 )return

      if( current_clock.NE.0 )call SYSTEM_CLOCK(start_tick)
#ifdef use_libSMA
!allocate space from the stack for pwrk and b
      ptr = LOC(mpp_stack)
      words = size(work)*size(transfer(work(1),word))
      if( words.GT.mpp_stack_size )then
          write( text, '(i8)' )words
          call mpp_error( FATAL, 'MPP_REDUCE user stack overflow: call mpp_set_stack_size('//text//') from all PEs.' )
      end if
      mpp_stack_hwm = max( words, mpp_stack_hwm )
      
      work(1) = a
      call SHMEM_REDUCE_( work, work, 1, peset(n)%start, peset(n)%log2stride, peset(n)%count, work(2), sync )
      call mpp_sync(pelist)
      a = work(1)
#endif use_libSMA
#ifdef use_libMPI
      if( verbose )call mpp_error( NOTE, 'MPP_REDUCE_: using MPI_ALLREDUCE...' )
      call MPI_ALLREDUCE( a, work, 1, MPI_TYPE_, MPI_REDUCE_, peset(n)%id, error )
      a = work
#endif
      if( current_clock.NE.0 )call increment_current_clock( EVENT_ALLREDUCE, MPP_TYPE_BYTELEN_ )
      return
    end subroutine MPP_REDUCE_
