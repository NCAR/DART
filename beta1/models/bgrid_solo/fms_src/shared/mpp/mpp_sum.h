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
    subroutine MPP_SUM_( a, length, pelist )
!sums array a over the PEs in pelist (all PEs if this argument is omitted)
!result is also automatically broadcast: all PEs have the sum in a at the end
!we are using f77-style call: array passed by address and not descriptor; further, the f90 conformance check is avoided.
      integer, intent(in) :: length
      integer, intent(in), optional :: pelist(:)
      MPP_TYPE_, intent(inout) :: a(*)
      integer :: n
#ifdef use_libSMA
!first <length> words are array, rest are pWrk
      MPP_TYPE_ :: work(length+length/2+1+SHMEM_REDUCE_MIN_WRKDATA_SIZE)
      pointer( ptr, work )
      integer :: words
      character(len=8) :: text
#else
      MPP_TYPE_ :: work(length)
#endif

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_SUM: You must first call mpp_init.' )
      n = get_peset(pelist); if( peset(n)%count.EQ.1 )return

      if( current_clock.NE.0 )call SYSTEM_CLOCK(start_tick)
#ifdef use_libSMA
!allocate space from the stack for pwrk and b
      ptr = LOC(mpp_stack)
      words = size(work)*size(transfer(work(1),word))
      if( words.GT.mpp_stack_size )then
          write( text, '(i8)' )words
          call mpp_error( FATAL, 'MPP_SUM user stack overflow: call mpp_set_stack_size('//text//') from all PEs.' )
      end if
      mpp_stack_hwm = max( words, mpp_stack_hwm )
      work(1:length) = a(1:length)
      call mpp_sync(pelist)
      call SHMEM_SUM_( work, work, length, peset(n)%start, peset(n)%log2stride, peset(n)%count, work(length+1), sync )
#endif use_libSMA
#ifdef use_libMPI
      if( verbose )call mpp_error( NOTE, 'MPP_SUM: using MPI_ALLREDUCE...' )
      if( debug )write( stderr(),* )'pe, n, peset(n)%id=', pe, n, peset(n)%id
      call MPI_ALLREDUCE( a, work, length, MPI_TYPE_, MPI_SUM, peset(n)%id, error )
#endif
      a(1:length) = work(1:length)
      if( current_clock.NE.0 )call increment_current_clock( EVENT_ALLREDUCE, length*MPP_TYPE_BYTELEN_ )
      return
    end subroutine MPP_SUM_

    subroutine MPP_SUM_SCALAR_( a, pelist )
!sums array a when only first element is passed: this routine just converts to a call to MPP_SUM_
      MPP_TYPE_, intent(inout) :: a
      integer, intent(in), optional :: pelist(:)
      MPP_TYPE_ :: b(1)

      b(1) = a
      if( debug )call mpp_error( NOTE, 'MPP_SUM_SCALAR_: calling MPP_SUM_ ...' )
      call MPP_SUM_( b, 1, pelist )
      a = b(1)
      return
    end subroutine MPP_SUM_SCALAR_

    subroutine MPP_SUM_2D_( a, length, pelist )
      MPP_TYPE_, intent(inout) :: a(:,:)
      integer, intent(in) :: length
      integer, intent(in), optional :: pelist(:)
      MPP_TYPE_ :: a1D(length)
#ifdef use_CRI_pointers
      pointer( ptr, a1D )
      ptr = LOC(a)
      call mpp_sum( a1D, length, pelist )
#else
      a1D = TRANSFER( a, a1D, length ) !faster than RESHAPE? length is probably redundant
!      a1D = RESHAPE( a, SHAPE(a1D) )
      call mpp_sum( a1D, length, pelist )
      a = RESHAPE( a1D, SHAPE(a) )
#endif
      return
    end subroutine MPP_SUM_2D_

    subroutine MPP_SUM_3D_( a, length, pelist )
      MPP_TYPE_, intent(inout) :: a(:,:,:)
      integer, intent(in) :: length
      integer, intent(in), optional :: pelist(:)
      MPP_TYPE_ :: a1D(length)
#ifdef use_CRI_pointers
      pointer( ptr, a1D )
      ptr = LOC(a)
      call mpp_sum( a1D, length, pelist )
#else
      a1D = TRANSFER( a, a1D, length ) !faster than RESHAPE? length is probably redundant
!      a1D = RESHAPE( a, SHAPE(a1D) )
      call mpp_sum( a1D, length, pelist )
      a = RESHAPE( a1D, SHAPE(a) )
#endif
      return
    end subroutine MPP_SUM_3D_

    subroutine MPP_SUM_4D_( a, length, pelist )
      MPP_TYPE_, intent(inout) :: a(:,:,:,:)
      integer, intent(in) :: length
      integer, intent(in), optional :: pelist(:)
      MPP_TYPE_ :: a1D(length)
#ifdef use_CRI_pointers
      pointer( ptr, a1D )
      ptr = LOC(a)
      call mpp_sum( a1D, length, pelist )
#else
      a1D = TRANSFER( a, a1D, length ) !faster than RESHAPE? length is probably redundant
!      a1D = RESHAPE( a, SHAPE(a1D) )
      call mpp_sum( a1D, length, pelist )
      a = RESHAPE( a1D, SHAPE(a) )
#endif
      return
    end subroutine MPP_SUM_4D_

    subroutine MPP_SUM_5D_( a, length, pelist )
      MPP_TYPE_, intent(inout) :: a(:,:,:,:,:)
      integer, intent(in) :: length
      integer, intent(in), optional :: pelist(:)
      MPP_TYPE_ :: a1D(length)
#ifdef use_CRI_pointers
      pointer( ptr, a1D )
      ptr = LOC(a)
      call mpp_sum( a1D, length, pelist )
#else
      a1D = TRANSFER( a, a1D, length ) !faster than RESHAPE? length is probably redundant
!      a1D = RESHAPE( a, SHAPE(a1D) )
      call mpp_sum( a1D, length, pelist )
      a = RESHAPE( a1D, SHAPE(a) )
#endif
      return
    end subroutine MPP_SUM_5D_
