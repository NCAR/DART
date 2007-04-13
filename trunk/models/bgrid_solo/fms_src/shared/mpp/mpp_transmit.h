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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                                  MPP_TRANSMIT                               !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine MPP_TRANSMIT_( put_data, put_len, to_pe, get_data, get_len, from_pe )
!a message-passing routine intended to be reminiscent equally of both MPI and SHMEM

!put_data and get_data are contiguous MPP_TYPE_ arrays

!at each call, your put_data array is put to   to_pe's get_data
!              your get_data array is got from from_pe's put_data
!i.e we assume that typically (e.g updating halo regions) each PE performs a put _and_ a get

!special PE designations:
!      NULL_PE: to disable a put or a get (e.g at boundaries)
!      ANY_PE:  if remote PE for the put or get is to be unspecific
!      ALL_PES: broadcast and collect operations (collect not yet implemented)

!ideally we would not pass length, but this f77-style call performs better (arrays passed by address, not descriptor)
!further, this permits <length> contiguous words from an array of any rank to be passed (avoiding f90 rank conformance check)

!caller is responsible for completion checks (mpp_sync_self) before and after

      integer, intent(in) :: put_len, to_pe, get_len, from_pe
      MPP_TYPE_, intent(in)  :: put_data(*)
      MPP_TYPE_, intent(out) :: get_data(*)
      integer :: i
#ifdef use_libSMA
      integer :: np
      integer(LONG_KIND) :: data_loc
!pointer to remote data
      MPP_TYPE_ :: remote_data(get_len)
      pointer( ptr_remote_data, remote_data )
      MPP_TYPE_ :: broadcast_data(get_len)
      pointer( ptr, broadcast_data )
      integer :: words
      character(len=8) :: text
#endif
      MPP_TYPE_, allocatable, save :: local_data(:) !local copy used by non-parallel code (no SHMEM or MPI)

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_TRANSMIT: You must first call mpp_init.' )
      if( to_pe.EQ.NULL_PE .AND. from_pe.EQ.NULL_PE )return
      
      if( debug )then
          call SYSTEM_CLOCK(tick)
          write( stdout(),'(a,i18,a,i5,a,2i5,2i8)' )&
               'T=',tick, ' PE=',pe, ' MPP_TRANSMIT begin: to_pe, from_pe, put_len, get_len=', to_pe, from_pe, put_len, get_len
      end if

!do put first and then get
      if( to_pe.GE.0 .AND. to_pe.LT.npes )then
#ifdef use_libSMA
!send data pointer to to_pe
#ifdef _CRAYT90
          call SHMEM_UDCFLUSH !invalidate data cache
#endif
          if( current_clock.NE.0 )call SYSTEM_CLOCK(start_tick)
          call SHMEM_INT8_WAIT( status(to_pe), MPP_WAIT )
          status(to_pe) = MPP_WAIT !prohibit puts to to_pe until it has retrieved this message
          if( current_clock.NE.0 )call increment_current_clock(EVENT_WAIT)
          data_loc = LOC(put_data)
          if( current_clock.NE.0 )call SYSTEM_CLOCK(start_tick)
          call SHMEM_INTEGER_PUT( mpp_from_pe, pe, 1, to_pe )
          call SHMEM_PUT8( remote_data_loc(pe), data_loc, 1, to_pe )
          if( current_clock.NE.0 )call increment_current_clock( EVENT_SEND, put_len*MPP_TYPE_BYTELEN_ )
#elif use_libMPI
!use non-blocking sends
          if( current_clock.NE.0 )call SYSTEM_CLOCK(start_tick)
          if( request(to_pe).NE.MPI_REQUEST_NULL )then !only one message from pe->to_pe in queue 
!              if( debug )write( stderr(),* )'PE waiting ', pe, to_pe
              call MPI_WAIT( request(to_pe), stat, error )
          end if
          call MPI_ISEND( put_data, put_len, MPI_TYPE_, to_pe, tag, MPI_COMM_WORLD, request(to_pe), error )
          if( current_clock.NE.0 )call increment_current_clock( EVENT_SEND, put_len*MPP_TYPE_BYTELEN_ )
#else                           !neither SHMEM nor MPI
          if( allocated(local_data) ) &
               call mpp_error( FATAL, 'MPP_TRANSMIT: local_data should have been deallocated by prior receive.' )
          allocate( local_data(put_len) )
          do i = 1,put_len
             local_data(i) = put_data(i)
          end do
#endif

      else if( to_pe.EQ.ALL_PES )then !this is a broadcast from from_pe
          if( from_pe.LT.0 .OR. from_pe.GE.npes )call mpp_error( FATAL, 'MPP_TRANSMIT: broadcasting from invalid PE.' )
          if( put_len.GT.get_len )call mpp_error( FATAL, 'MPP_TRANSMIT: size mismatch between put_data and get_data.' )
          if( pe.EQ.from_pe )then
#ifdef use_CRI_pointers
              if( LOC(get_data).NE.LOC(put_data) )then
!dir$ IVDEP
#endif
                  do i = 1,get_len
                     get_data(i) = put_data(i)
                  end do
#ifdef use_CRI_pointers
              end if
#endif
          end if
          call mpp_broadcast( get_data, get_len, from_pe )
          return

      else if( to_pe.EQ.ANY_PE )then !we don't have a destination to do puts to, so only do gets
#ifdef use_libSMA
          if( from_pe.LT.0 .OR. from_pe.GE.npes )call mpp_error( FATAL, 'MPP_TRANSMIT: invalid from_pe along with to_pe=ANY_PE.' )
          if( current_clock.NE.0 )call SYSTEM_CLOCK(start_tick)
          call SHMEM_GET_( get_data, put_data, get_len, from_pe )
          call SHMEM_PUT8( status(pe), MPP_READY, 1, from_pe ) !tell from_pe that you have retrieved this message
          if( current_clock.NE.0 )call increment_current_clock( EVENT_RECV, get_len*MPP_TYPE_BYTELEN_ )
          return
#endif
#ifdef use_libMPI
!...but you cannot have a pure get with MPI
          call mpp_error( FATAL, 'MPP_TRANSMIT: you cannot transmit to ANY_PE using MPI.' )
#endif

      else if( to_pe.NE.NULL_PE )then   !no other valid cases except NULL_PE
          call mpp_error( FATAL, 'MPP_TRANSMIT: invalid to_pe.' )
      end if

!do the get: for libSMA, a get means do a wait to ensure put on remote PE is complete
      if( from_pe.GE.0 .AND. from_pe.LT.npes )then
#ifdef use_libSMA
#ifdef _CRAYT90
          call SHMEM_UDCFLUSH !invalidate data cache
#endif
          if( current_clock.NE.0 )call SYSTEM_CLOCK(start_tick)
          if( debug )write( stderr(),* )'pe, from_pe, remote_data_loc(from_pe)=', pe, from_pe, remote_data_loc(from_pe)
          call SHMEM_INT8_WAIT( remote_data_loc(from_pe), MPP_WAIT )
          if( current_clock.NE.0 )call increment_current_clock(EVENT_WAIT)
          ptr_remote_data = remote_data_loc(from_pe)
          remote_data_loc(from_pe) = MPP_WAIT !reset
!          call SHMEM_PUT8( status(pe), MPP_READY, 1, from_pe ) !tell from_pe we have retrieved the location
          if( current_clock.NE.0 )call SYSTEM_CLOCK(start_tick)
#if defined(CRAYPVP) || defined(sgi_mipspro)
!since we have the pointer to remote data, just retrieve it with a simple copy
          if( LOC(get_data).NE.LOC(remote_data) )then
!dir$ IVDEP
              do i = 1,get_len
                 get_data(i) = remote_data(i)
              end do
          else
              call mpp_error(FATAL)
          end if
#else
          call SHMEM_GET_( get_data, remote_data, get_len, from_pe )
#endif
          call SHMEM_PUT8( status(pe), MPP_READY, 1, from_pe ) !tell from_pe we have retrieved the location
          if( current_clock.NE.0 )call increment_current_clock( EVENT_RECV, get_len*MPP_TYPE_BYTELEN_ )
#elif  use_libMPI
!receive from from_pe
          if( current_clock.NE.0 )call SYSTEM_CLOCK(start_tick)
          call MPI_RECV( get_data, get_len, MPI_TYPE_, from_pe, MPI_ANY_TAG, MPI_COMM_WORLD, stat, error )
          if( current_clock.NE.0 )call increment_current_clock( EVENT_RECV, get_len*MPP_TYPE_BYTELEN_ )
#else                           !neither SHMEM nor MPI
          if( .NOT.allocated(local_data) ) &
               call mpp_error( FATAL, 'MPP_TRANSMIT: local_data should have been allocated by prior send.' )
          do i = 1,get_len
             get_data(i) = local_data(i)
          end do
          deallocate(local_data)
!#else !neither use_libSMA nor use_libMPI
!          if( pe.EQ.from_pe )then
!#ifdef use_CRI_pointers
!!dir$ IVDEP
!              if( LOC(get_data).NE.LOC(put_data) ) &
!#endif
!                   get_data(1:put_len) = put_data(1:put_len)
!          end if
#endif

      else if( from_pe.EQ.ANY_PE )then
#ifdef use_libSMA
#ifdef _CRAYT90
          call SHMEM_UDCFLUSH !invalidate data cache
#endif
!since we don't know which PE is sending us data, we wait for remote PE to send us its ID
!this is only required for !CRAYPVP  && !sgi_mipspro, but is done there too, so that we can send put_is_done back.
          call shmem_integer_wait( mpp_from_pe, ANY_PE )
          if( current_clock.NE.0 )call SYSTEM_CLOCK(start_tick)
          call SHMEM_INT8_WAIT( remote_data_loc(mpp_from_pe), MPP_WAIT )
          if( current_clock.NE.0 )call increment_current_clock(EVENT_WAIT)
          ptr_remote_data = remote_data_loc(mpp_from_pe)
          remote_data_loc(mpp_from_pe) = MPP_WAIT !reset
          call SHMEM_PUT8( status(pe), MPP_READY, 1, mpp_from_pe ) !tell mpp_from_pe we have retrieved the location
#if defined(CRAYPVP) || defined(sgi_mipspro)
!since we have the pointer to remote data, just retrieve it with a simple copy
          if( current_clock.NE.0 )call SYSTEM_CLOCK(start_tick)
          if( LOC(get_data).NE.LOC(remote_data) )then
!dir$ IVDEP
              do i = 1,get_len
                 get_data(i) = remote_data(i)
              end do
          end if
#else
          call SHMEM_GET_( get_data, remote_data, get_len, mpp_from_pe )
#endif
          if( current_clock.NE.0 )call increment_current_clock( EVENT_RECV, get_len*MPP_TYPE_BYTELEN_ )
          mpp_from_pe = ANY_PE   !reset
#endif use_libSMA
#ifdef use_libMPI
!receive from MPI_ANY_SOURCE
          if( current_clock.NE.0 )call SYSTEM_CLOCK(start_tick)
          call MPI_RECV( get_data, get_len, MPI_TYPE_, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, stat, error )
          if( current_clock.NE.0 )call increment_current_clock( EVENT_RECV, get_len*MPP_TYPE_BYTELEN_ )
#endif

      else if( from_pe.EQ.ALL_PES )then
          call mpp_error( FATAL, 'MPP_TRANSMIT: from_pe=ALL_PES has ambiguous meaning, and hence is not implemented.' )

      else if( from_pe.NE.NULL_PE )then !only remaining valid choice is NULL_PE
          call mpp_error( FATAL, 'MPP_TRANSMIT: invalid from_pe.' )
      end if

      if( debug )then
          call SYSTEM_CLOCK(tick)
          write( stdout(),'(a,i18,a,i5,a,2i5,2i8)' )&
               'T=',tick, ' PE=',pe, ' MPP_TRANSMIT end: to_pe, from_pe, put_len, get_len=', to_pe, from_pe, put_len, get_len
      end if
      return
    end subroutine MPP_TRANSMIT_

    subroutine MPP_TRANSMIT_SCALAR_( put_data, to_pe, get_data, from_pe )
      integer, intent(in) :: to_pe, from_pe
      MPP_TYPE_, intent(in)  :: put_data
      MPP_TYPE_, intent(out) :: get_data
      MPP_TYPE_ :: put_data1D(1), get_data1D(1)
#ifdef use_CRI_pointers
      pointer( ptrp, put_data1D )
      pointer( ptrg, get_data1D )

      ptrp = LOC(put_data)
      ptrg = LOC(get_data)
      call MPP_TRANSMIT_ ( put_data1D, 1, to_pe, get_data1D, 1, from_pe )
#else
      put_data1D(1) = put_data
      call MPP_TRANSMIT_ ( put_data1D, 1, to_pe, get_data1D, 1, from_pe )
      get_data = get_data1D(1)
#endif
      return
    end subroutine MPP_TRANSMIT_SCALAR_

    subroutine MPP_TRANSMIT_2D_( put_data, put_len, to_pe, get_data, get_len, from_pe )
      integer, intent(in) :: put_len, to_pe, get_len, from_pe
      MPP_TYPE_, intent(in)  :: put_data(:,:)
      MPP_TYPE_, intent(out) :: get_data(:,:)
      MPP_TYPE_ :: put_data1D(put_len), get_data1D(get_len)
#ifdef use_CRI_pointers
      pointer( ptrp, put_data1D )
      pointer( ptrg, get_data1D )
      ptrp = LOC(put_data)
      ptrg = LOC(get_data)
      call mpp_transmit( put_data1D, put_len, to_pe, get_data1D, get_len, from_pe )
#else
      if( to_pe.NE.NULL_PE )put_data1D = TRANSFER( put_data, put_data1D, get_len ) !faster than RESHAPE? get_len is probably redundant
!      if( to_pe.NE.NULL_PE )put_data1D = RESHAPE( put_data, SHAPE(put_data1D) )
      call mpp_transmit( put_data1D, put_len, to_pe, get_data1D, get_len, from_pe )
      if( from_pe.NE.NULL_PE )get_data = RESHAPE( get_data1D, SHAPE(get_data) )
#endif
      return
    end subroutine MPP_TRANSMIT_2D_

    subroutine MPP_TRANSMIT_3D_( put_data, put_len, to_pe, get_data, get_len, from_pe )
      integer, intent(in) :: put_len, to_pe, get_len, from_pe
      MPP_TYPE_, intent(in)  :: put_data(:,:,:)
      MPP_TYPE_, intent(out) :: get_data(:,:,:)
      MPP_TYPE_ :: put_data1D(put_len), get_data1D(get_len)
#ifdef use_CRI_pointers
      pointer( ptrp, put_data1D )
      pointer( ptrg, get_data1D )
      ptrp = LOC(put_data)
      ptrg = LOC(get_data)
      call mpp_transmit( put_data1D, put_len, to_pe, get_data1D, get_len, from_pe )
#else
      if( to_pe.NE.NULL_PE )put_data1D = TRANSFER( put_data, put_data1D, get_len ) !faster than RESHAPE? get_len is probably redundant
!      if( to_pe.NE.NULL_PE )put_data1D = RESHAPE( put_data, SHAPE(put_data1D) )
      call mpp_transmit( put_data1D, put_len, to_pe, get_data1D, get_len, from_pe )
      if( from_pe.NE.NULL_PE )get_data = RESHAPE( get_data1D, SHAPE(get_data) )
#endif
      return
    end subroutine MPP_TRANSMIT_3D_

    subroutine MPP_TRANSMIT_4D_( put_data, put_len, to_pe, get_data, get_len, from_pe )
      integer, intent(in) :: put_len, to_pe, get_len, from_pe
      MPP_TYPE_, intent(in)  :: put_data(:,:,:,:)
      MPP_TYPE_, intent(out) :: get_data(:,:,:,:)
      MPP_TYPE_ :: put_data1D(put_len), get_data1D(get_len)
#ifdef use_CRI_pointers
      pointer( ptrp, put_data1D )
      pointer( ptrg, get_data1D )
      ptrp = LOC(put_data)
      ptrg = LOC(get_data)
      call mpp_transmit( put_data1D, put_len, to_pe, get_data1D, get_len, from_pe )
#else
      if( to_pe.NE.NULL_PE )put_data1D = TRANSFER( put_data, put_data1D, get_len ) !faster than RESHAPE? get_len is probably redundant
!      if( to_pe.NE.NULL_PE )put_data1D = RESHAPE( put_data, SHAPE(put_data1D) )
      call mpp_transmit( put_data1D, put_len, to_pe, get_data1D, get_len, from_pe )
      if( from_pe.NE.NULL_PE )get_data = RESHAPE( get_data1D, SHAPE(get_data) )
#endif
      return
    end subroutine MPP_TRANSMIT_4D_

    subroutine MPP_TRANSMIT_5D_( put_data, put_len, to_pe, get_data, get_len, from_pe )
      integer, intent(in) :: put_len, to_pe, get_len, from_pe
      MPP_TYPE_, intent(in)  :: put_data(:,:,:,:,:)
      MPP_TYPE_, intent(out) :: get_data(:,:,:,:,:)
      MPP_TYPE_ :: put_data1D(put_len), get_data1D(get_len)
#ifdef use_CRI_pointers
      pointer( ptrp, put_data1D )
      pointer( ptrg, get_data1D )
      ptrp = LOC(put_data)
      ptrg = LOC(get_data)
      call mpp_transmit( put_data1D, put_len, to_pe, get_data1D, get_len, from_pe )
#else
      if( to_pe.NE.NULL_PE )put_data1D = TRANSFER( put_data, put_data1D, get_len ) !faster than RESHAPE? get_len is probably redundant
!      if( to_pe.NE.NULL_PE )put_data1D = RESHAPE( put_data, SHAPE(put_data1D) )
      call mpp_transmit( put_data1D, put_len, to_pe, get_data1D, get_len, from_pe )
      if( from_pe.NE.NULL_PE )get_data = RESHAPE( get_data1D, SHAPE(get_data) )
#endif
      return
    end subroutine MPP_TRANSMIT_5D_
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                              MPP_SEND and RECV                              !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine MPP_RECV_( get_data, get_len, from_pe )
!a mpp_transmit with null arguments on the put side
      integer, intent(in) :: get_len, from_pe
      MPP_TYPE_, intent(out) :: get_data(*)
      MPP_TYPE_ :: dummy(1)
      call mpp_transmit( dummy, 1, NULL_PE, get_data, get_len, from_pe )
    end subroutine MPP_RECV_

    subroutine MPP_SEND_( put_data, put_len, to_pe )
!a mpp_transmit with null arguments on the get side
      integer, intent(in) :: put_len, to_pe
      MPP_TYPE_, intent(in) :: put_data(*)
      MPP_TYPE_ :: dummy(1)
      call mpp_transmit( put_data, put_len, to_pe, dummy, 1, NULL_PE )
    end subroutine MPP_SEND_

    subroutine MPP_RECV_SCALAR_( get_data, from_pe )
!a mpp_transmit with null arguments on the put side
      integer, intent(in) :: from_pe
      MPP_TYPE_, intent(out) :: get_data
      MPP_TYPE_ :: get_data1D(1)
      MPP_TYPE_ :: dummy(1)
#ifdef use_CRI_pointers
      pointer( ptr, get_data1D )
      ptr = LOC(get_data)
      call mpp_transmit( dummy, 1, NULL_PE, get_data1D, 1, from_pe )
#else
      call mpp_transmit( dummy, 1, NULL_PE, get_data1D, 1, from_pe )
      get_data = get_data1D(1)
#endif
    end subroutine MPP_RECV_SCALAR_

    subroutine MPP_SEND_SCALAR_( put_data, to_pe )
!a mpp_transmit with null arguments on the get side
      integer, intent(in) :: to_pe
      MPP_TYPE_, intent(in) :: put_data
      MPP_TYPE_ :: put_data1D(1)
      MPP_TYPE_ :: dummy(1)
#ifdef use_CRI_pointers
      pointer( ptr, put_data1D )
      ptr = LOC(put_data)
      call mpp_transmit( put_data1D, 1, to_pe, dummy, 1, NULL_PE )
#else
      put_data1D(1) = put_data
      call mpp_transmit( put_data1D, 1, to_pe, dummy, 1, NULL_PE )
#endif
    end subroutine MPP_SEND_SCALAR_

    subroutine MPP_RECV_2D_( get_data, get_len, from_pe )
!a mpp_transmit with null arguments on the put side
      integer, intent(in) :: get_len, from_pe
      MPP_TYPE_, intent(out) :: get_data(:,:)
      MPP_TYPE_ :: dummy(1,1)
      call mpp_transmit( dummy, 1, NULL_PE, get_data, get_len, from_pe )
    end subroutine MPP_RECV_2D_

    subroutine MPP_SEND_2D_( put_data, put_len, to_pe )
!a mpp_transmit with null arguments on the get side
      integer, intent(in) :: put_len, to_pe
      MPP_TYPE_, intent(in) :: put_data(:,:)
      MPP_TYPE_ :: dummy(1,1)
      call mpp_transmit( put_data, put_len, to_pe, dummy, 1, NULL_PE )
    end subroutine MPP_SEND_2D_

    subroutine MPP_RECV_3D_( get_data, get_len, from_pe )
!a mpp_transmit with null arguments on the put side
      integer, intent(in) :: get_len, from_pe
      MPP_TYPE_, intent(out) :: get_data(:,:,:)
      MPP_TYPE_ :: dummy(1,1,1)
      call mpp_transmit( dummy, 1, NULL_PE, get_data, get_len, from_pe )
    end subroutine MPP_RECV_3D_

    subroutine MPP_SEND_3D_( put_data, put_len, to_pe )
!a mpp_transmit with null arguments on the get side
      integer, intent(in) :: put_len, to_pe
      MPP_TYPE_, intent(in) :: put_data(:,:,:)
      MPP_TYPE_ :: dummy(1,1,1)
      call mpp_transmit( put_data, put_len, to_pe, dummy, 1, NULL_PE )
    end subroutine MPP_SEND_3D_

    subroutine MPP_RECV_4D_( get_data, get_len, from_pe )
!a mpp_transmit with null arguments on the put side
      integer, intent(in) :: get_len, from_pe
      MPP_TYPE_, intent(out) :: get_data(:,:,:,:)
      MPP_TYPE_ :: dummy(1,1,1,1)
      call mpp_transmit( dummy, 1, NULL_PE, get_data, get_len, from_pe )
    end subroutine MPP_RECV_4D_

    subroutine MPP_SEND_4D_( put_data, put_len, to_pe )
!a mpp_transmit with null arguments on the get side
      integer, intent(in) :: put_len, to_pe
      MPP_TYPE_, intent(in) :: put_data(:,:,:,:)
      MPP_TYPE_ :: dummy(1,1,1,1)
      call mpp_transmit( put_data, put_len, to_pe, dummy, 1, NULL_PE )
    end subroutine MPP_SEND_4D_

    subroutine MPP_RECV_5D_( get_data, get_len, from_pe )
!a mpp_transmit with null arguments on the put side
      integer, intent(in) :: get_len, from_pe
      MPP_TYPE_, intent(out) :: get_data(:,:,:,:,:)
      MPP_TYPE_ :: dummy(1,1,1,1,1)
      call mpp_transmit( dummy, 1, NULL_PE, get_data, get_len, from_pe )
    end subroutine MPP_RECV_5D_

    subroutine MPP_SEND_5D_( put_data, put_len, to_pe )
!a mpp_transmit with null arguments on the get side
      integer, intent(in) :: put_len, to_pe
      MPP_TYPE_, intent(in) :: put_data(:,:,:,:,:)
      MPP_TYPE_ :: dummy(1,1,1,1,1)
      call mpp_transmit( put_data, put_len, to_pe, dummy, 1, NULL_PE )
    end subroutine MPP_SEND_5D_
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                                MPP_BROADCAST                                !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine MPP_BROADCAST_( data, length, from_pe, pelist )
!this call was originally bundled in with mpp_transmit, but that doesn't allow
!broadcast to a subset of PEs. This version will, and mpp_transmit will remain
!backward compatible.
      MPP_TYPE_, intent(inout) :: data(*)
      integer, intent(in) :: length, from_pe
      integer, intent(in), optional :: pelist(:)
      integer :: n
#ifdef use_libSMA
      integer :: np, i
      integer(LONG_KIND) :: data_loc
!pointer to remote data
      MPP_TYPE_ :: bdata(length)
      pointer( ptr, bdata )
      integer :: words
      character(len=8) :: text
#endif

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_BROADCAST: You must first call mpp_init.' )
      n = get_peset(pelist); if( peset(n)%count.EQ.1 )return

      if( debug )then
          call SYSTEM_CLOCK(tick)
          write( stdout(),'(a,i18,a,i5,a,2i5,2i8)' )&
               'T=',tick, ' PE=',pe, ' MPP_BROADCAST begin: from_pe, length=', from_pe, length
      end if

      if( .NOT.ANY(from_pe.EQ.peset(current_peset_num)%list) ) &
           call mpp_error( FATAL, 'MPP_BROADCAST: broadcasting from invalid PE.' )

      if( current_clock.NE.0 )call SYSTEM_CLOCK(start_tick)
#ifdef use_libSMA
      ptr = LOC(mpp_stack)
      words = size(bdata)*size(transfer(bdata(1),word))
      if( words.GT.mpp_stack_size )then
          write( text, '(i8)' )words
          call mpp_error( FATAL, 'MPP_BROADCAST user stack overflow: call mpp_set_stack_size('//text//') from all PEs.' )
      end if
      mpp_stack_hwm = max( words, mpp_stack_hwm )
      if( mpp_npes().GT.1 )then
!dir$ IVDEP
          do i = 1,length
             bdata(i) = data(i)
          end do
          call mpp_sync(pelist) !eliminate?
#ifdef _CRAYT90
          call SHMEM_UDCFLUSH !invalidate data cache
#endif
          call SHMEM_BROADCAST_( bdata, bdata, length, from_pe, peset(n)%start, peset(n)%log2stride, peset(n)%count, sync )
          call mpp_sync(pelist) !eliminate?
!dir$ IVDEP
          do i = 1,length
             data(i) = bdata(i)
          end do
      end if
#endif
#ifdef use_libMPI
      if( mpp_npes().GT.1 )call MPI_BCAST( data, length, MPI_TYPE_, from_pe, peset(n)%id, error )
#endif
      if( current_clock.NE.0 )call increment_current_clock( EVENT_BROADCAST, length*MPP_TYPE_BYTELEN_ )
      return
    end subroutine MPP_BROADCAST_

    subroutine MPP_BROADCAST_SCALAR_( data, from_pe, pelist )
      MPP_TYPE_, intent(inout) :: data
      integer, intent(in) :: from_pe
      integer, intent(in), optional :: pelist(:)
      MPP_TYPE_ :: data1D(1)
#ifdef use_CRI_pointers
      pointer( ptr, data1D )

      ptr = LOC(data)
      call MPP_BROADCAST_( data1D, 1, from_pe, pelist )
#else
      data1D(1) = data
      call MPP_BROADCAST_( data1D, 1, from_pe, pelist )
      data = data1D(1)
#endif
      return
    end subroutine MPP_BROADCAST_SCALAR_
    
    subroutine MPP_BROADCAST_2D_( data, length, from_pe, pelist )
!this call was originally bundled in with mpp_transmit, but that doesn't allow
!broadcast to a subset of PEs. This version will, and mpp_transmit will remain
!backward compatible.
      MPP_TYPE_, intent(inout) :: data(:,:)
      integer, intent(in) :: length, from_pe
      integer, intent(in), optional :: pelist(:)
      MPP_TYPE_ :: data1D(length)
#ifdef use_CRI_pointers
      pointer( ptr, data1D )
      ptr = LOC(data)
      call mpp_broadcast( data1D, length, from_pe, pelist )
#else
      data1D = TRANSFER( data, data1D, length ) !faster than RESHAPE? length is probably redundant
!      data1D = RESHAPE( data, SHAPE(data1D) )
      call mpp_broadcast( data1D, length, from_pe, pelist )
      data = RESHAPE( data1D, SHAPE(data) )
#endif
      return
    end subroutine MPP_BROADCAST_2D_
    
    subroutine MPP_BROADCAST_3D_( data, length, from_pe, pelist )
!this call was originally bundled in with mpp_transmit, but that doesn't allow
!broadcast to a subset of PEs. This version will, and mpp_transmit will remain
!backward compatible.
      MPP_TYPE_, intent(inout) :: data(:,:,:)
      integer, intent(in) :: length, from_pe
      integer, intent(in), optional :: pelist(:)
      MPP_TYPE_ :: data1D(length)
#ifdef use_CRI_pointers
      pointer( ptr, data1D )
      ptr = LOC(data)
      call mpp_broadcast( data1D, length, from_pe, pelist )
#else
      data1D = TRANSFER( data, data1D, length ) !faster than RESHAPE? length is probably redundant
!      data1D = RESHAPE( data, SHAPE(data1D) )
      call mpp_broadcast( data1D, length, from_pe, pelist )
      data = RESHAPE( data1D, SHAPE(data) )
#endif
      return
    end subroutine MPP_BROADCAST_3D_
    
    subroutine MPP_BROADCAST_4D_( data, length, from_pe, pelist )
!this call was originally bundled in with mpp_transmit, but that doesn't allow
!broadcast to a subset of PEs. This version will, and mpp_transmit will remain
!backward compatible.
      MPP_TYPE_, intent(inout) :: data(:,:,:,:)
      integer, intent(in) :: length, from_pe
      integer, intent(in), optional :: pelist(:)
      MPP_TYPE_ :: data1D(length)
#ifdef use_CRI_pointers
      pointer( ptr, data1D )
      ptr = LOC(data)
      call mpp_broadcast( data1D, length, from_pe, pelist )
#else
      data1D = TRANSFER( data, data1D, length ) !faster than RESHAPE? length is probably redundant
!      data1D = RESHAPE( data, SHAPE(data1D) )
      call mpp_broadcast( data1D, length, from_pe, pelist )
      data = RESHAPE( data1D, SHAPE(data) )
#endif
      return
    end subroutine MPP_BROADCAST_4D_
    
    subroutine MPP_BROADCAST_5D_( data, length, from_pe, pelist )
!this call was originally bundled in with mpp_transmit, but that doesn't allow
!broadcast to a subset of PEs. This version will, and mpp_transmit will remain
!backward compatible.
      MPP_TYPE_, intent(inout) :: data(:,:,:,:,:)
      integer, intent(in) :: length, from_pe
      integer, intent(in), optional :: pelist(:)
      MPP_TYPE_ :: data1D(length)
#ifdef use_CRI_pointers
      pointer( ptr, data1D )
      ptr = LOC(data)
      call mpp_broadcast( data1D, length, from_pe, pelist )
#else
      data1D = TRANSFER( data, data1D, length ) !faster than RESHAPE? length is probably redundant
!      data1D = RESHAPE( data, SHAPE(data1D) )
      call mpp_broadcast( data1D, length, from_pe, pelist )
      data = RESHAPE( data1D, SHAPE(data) )
#endif
      return
    end subroutine MPP_BROADCAST_5D_
