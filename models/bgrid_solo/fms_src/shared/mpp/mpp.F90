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
!-----------------------------------------------------------------------
#ifdef __aix
#define FLUSH FLUSH_
#define flush flush_
#endif

!                 Communication for message-passing codes
!
! AUTHOR: V. Balaji (vb@gfdl.gov)
!         SGI/GFDL Princeton University
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! For the full text of the GNU General Public License,
! write to: Free Software Foundation, Inc.,
!           675 Mass Ave, Cambridge, MA 02139, USA.  
!-----------------------------------------------------------------------

!these are used to determine hardware/OS/compiler
#include <os.h>
#ifdef __aix
#define FLUSH FLUSH_
#define flush flush_
#endif


!only one of SMA or MPI can be used
!(though mixing calls is allowed, this module will not)
#ifdef use_libSMA
#undef use_libMPI
#endif

!shmalloc is used on MPP SGI/Cray systems for shmem
#if defined(use_libSMA) && defined(SGICRAY_MPP)
#define use_shmalloc
#endif

module mpp_mod
!string BWA is used to tag lines that are bug workarounds and will disappear
!when offending compiler bug is fixed
!a generalized communication package for use with shmem and MPI
!will add: co_array_fortran, MPI2
!Balaji (vb@gfdl.gov) 11 May 1998

! <CONTACT EMAIL="vb@gfdl.noaa.gov">
!   V. Balaji
! </CONTACT>

! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>
! <RCSLOG SRC="http://www.gfdl.noaa.gov/~vb/changes_mpp.html"/>

! <OVERVIEW>
!   <TT>mpp_mod</TT>, is a set of simple calls to provide a uniform interface
!   to different message-passing libraries. It currently can be
!   implemented either in the SGI/Cray native SHMEM library or in the MPI
!   standard. Other libraries (e.g MPI-2, Co-Array Fortran) can be
!   incorporated as the need arises.
! </OVERVIEW>

! <DESCRIPTION>
!   The data transfer between a processor and its own memory is based
!   on <TT>load</TT> and <TT>store</TT> operations upon
!   memory. Shared-memory systems (including distributed shared memory
!   systems) have a single address space and any processor can acquire any
!   data within the memory by <TT>load</TT> and
!   <TT>store</TT>. The situation is different for distributed
!   parallel systems. Specialized MPP systems such as the T3E can simulate
!   shared-memory by direct data acquisition from remote memory. But if
!   the parallel code is distributed across a cluster, or across the Net,
!   messages must be sent and received using the protocols for
!   long-distance communication, such as TCP/IP. This requires a
!   ``handshaking'' between nodes of the distributed system. One can think
!   of the two different methods as involving <TT>put</TT>s or
!   <TT>get</TT>s (e.g the SHMEM library), or in the case of
!   negotiated communication (e.g MPI), <TT>send</TT>s and
!   <TT>recv</TT>s.
!   
!   The difference between SHMEM and MPI is that SHMEM uses one-sided
!   communication, which can have very low-latency high-bandwidth
!   implementations on tightly coupled systems. MPI is a standard
!   developed for distributed computing across loosely-coupled systems,
!   and therefore incurs a software penalty for negotiating the
!   communication. It is however an open industry standard whereas SHMEM
!   is a proprietary interface. Besides, the <TT>put</TT>s or
!   <TT>get</TT>s on which it is based cannot currently be implemented in
!   a cluster environment (there are recent announcements from Compaq that
!   occasion hope).
!   
!   The message-passing requirements of climate and weather codes can be
!   reduced to a fairly simple minimal set, which is easily implemented in
!   any message-passing API. <TT>mpp_mod</TT> provides this API.
!
!    Features of <TT>mpp_mod</TT> include:
!   
!    1) Simple, minimal API, with free access to underlying API for
!       more complicated stuff.<BR/>
!    2) Design toward typical use in climate/weather CFD codes.<BR/>
!    3) Performance to be not significantly lower than any native API.
!   
!   This module is used to develop higher-level calls for <LINK 
!   SRC="mpp_domains.html">domain decomposition</LINK> and <LINK
!   SRC="mpp_io.html">parallel I/O</LINK>.
!   
!   Parallel computing is initially daunting, but it soon becomes
!   second nature, much the way many of us can now write vector code
!   without much effort. The key insight required while reading and
!   writing parallel code is in arriving at a mental grasp of several
!   independent parallel execution streams through the same code (the SPMD
!   model). Each variable you examine may have different values for each
!   stream, the processor ID being an obvious example. Subroutines and
!   function calls are particularly subtle, since it is not always obvious
!   from looking at a call what synchronization between execution streams
!   it implies. An example of erroneous code would be a global barrier
!   call (see <LINK SRC="#mpp_sync"><TT>mpp_sync</TT></LINK> below) placed
!   within a code block that not all PEs will execute, e.g:
!   
!   <PRE>
!   if( pe.EQ.0 )call mpp_sync()
!   </PRE>
!   
!   Here only PE 0 reaches the barrier, where it will wait
!   indefinitely. While this is a particularly egregious example to
!   illustrate the coding flaw, more subtle versions of the same are
!   among the most common errors in parallel code.
!   
!   It is therefore important to be conscious of the context of a
!   subroutine or function call, and the implied synchronization. There
!   are certain calls here (e.g <TT>mpp_declare_pelist, mpp_init,
!   mpp_malloc, mpp_set_stack_size</TT>) which must be called by all
!   PEs. There are others which must be called by a subset of PEs (here
!   called a <TT>pelist</TT>) which must be called by all the PEs in the
!   <TT>pelist</TT> (e.g <TT>mpp_max, mpp_sum, mpp_sync</TT>). Still
!   others imply no synchronization at all. I will make every effort to
!   highlight the context of each call in the MPP modules, so that the
!   implicit synchronization is spelt out.  
!   
!   For performance it is necessary to keep synchronization as limited
!   as the algorithm being implemented will allow. For instance, a single
!   message between two PEs should only imply synchronization across the
!   PEs in question. A <I>global</I> synchronization (or <I>barrier</I>)
!   is likely to be slow, and is best avoided. But codes first
!   parallelized on a Cray T3E tend to have many global syncs, as very
!   fast barriers were implemented there in hardware.
!   
!   Another reason to use pelists is to run a single program in MPMD
!   mode, where different PE subsets work on different portions of the
!   code. A typical example is to assign an ocean model and atmosphere
!   model to different PE subsets, and couple them concurrently instead of
!   running them serially. The MPP module provides the notion of a
!   <I>current pelist</I>, which is set when a group of PEs branch off
!   into a subset. Subsequent calls that omit the <TT>pelist</TT> optional
!   argument (seen below in many of the individual calls) assume that the
!   implied synchronization is across the current pelist. The calls
!   <TT>mpp_root_pe</TT> and <TT>mpp_npes</TT> also return the values
!   appropriate to the current pelist. The <TT>mpp_set_current_pelist</TT>
!   call is provided to set the current pelist.

! </DESCRIPTION>
! <PUBLIC>
!  F90 is a strictly-typed language, and the syntax pass of the
!  compiler requires matching of type, kind and rank (TKR). Most calls
!  listed here use a generic type, shown here as <TT>MPP_TYPE_</TT>. This
!  is resolved in the pre-processor stage to any of a variety of
!  types. In general the MPP operations work on 4-byte and 8-byte
!  variants of <TT>integer, real, complex, logical</TT> variables, of
!  rank 0 to 5, leading to 48 specific module procedures under the same
!  generic interface. Any of the variables below shown as
!  <TT>MPP_TYPE_</TT> is treated in this way.
! </PUBLIC>
#ifdef sgi_mipspro
#ifdef use_libSMA
  use shmem_interface
#endif
#ifdef use_libMPI
  use mpi
#endif
#endif
  implicit none
  private
  character(len=128), private :: version= &
       '$Revision$'
  character(len=128), private :: tagname= &
       '$Id$'

!various lengths (see shpalloc) are estimated in "words" which are 32bit on SGI, 64bit on Cray
!these are also the expected sizeof of args to MPI/shmem libraries
#ifdef _CRAY
  integer(LONG_KIND), private :: word(1)
#endif
#ifdef sgi_mipspro
  integer(INT_KIND), private :: word(1)
#endif

#ifdef SGICRAY
!see intro_io(3F): to see why these values are used rather than 5,6,0
  integer, private :: in_unit=100, out_unit=101, err_unit=102
#else
  integer, private :: in_unit=5, out_unit=6, err_unit=0
#endif
  integer :: log_unit, etc_unit
  logical, private :: module_is_initialized=.FALSE.
  integer, private :: pe=0, node=0, npes=1, root_pe=0
  integer, private :: error
  integer, parameter, private :: MAXPES=2048 !used for dimensioning stuff that might be indexed by pe
  character(len=32) :: configfile='logfile.out'
  character(len=32) :: etcfile='._mpp.nonrootpe.stdout'

!initialization flags
  integer, parameter, public :: MPP_VERBOSE=1, MPP_DEBUG=2
  logical, private :: verbose=.FALSE., debug=.FALSE.

!flags to transmit routines
  integer, parameter, public :: ALL_PES=-1, ANY_PE=-2, NULL_PE=-3

!errortype flags
  integer, parameter, public :: NOTE=0, WARNING=1, FATAL=2
  logical, private :: warnings_are_fatal = .FALSE.
  integer, private :: error_state=0

  integer(LONG_KIND), parameter, private :: MPP_WAIT=-1, MPP_READY=-2
#ifdef use_libSMA
#include <mpp/shmem.fh>
  integer :: sync(SHMEM_REDUCE_SYNC_SIZE+SHMEM_BCAST_SYNC_SIZE+SHMEM_BARRIER_SYNC_SIZE)
!status and remote_data_loc are used to synchronize communication is MPP_TRANSMIT
#ifdef use_shmalloc
  integer(LONG_KIND), private, dimension(0:MAXPES) :: status, remote_data_loc
#else
  integer(LONG_KIND), private, allocatable, dimension(:) :: status, remote_data_loc
#endif
  integer, private :: mpp_from_pe !used to announce from where data is coming from
#ifdef use_shmalloc
!we call shpalloc in mpp_init() to ensure all these are remotely accessible
!on PVP where shpalloc doesn't exist, module variables are automatically
!guaranteed to be remotely accessible
  pointer( ptr_sync, sync )
  pointer( ptr_status, status )
  pointer( ptr_from, mpp_from_pe )
  pointer( ptr_remote, remote_data_loc )
#endif
#endif use_libSMA
#ifdef use_libMPI
#ifndef sgi_mipspro
!sgi_mipspro gets this from 'use mpi'
#include <mpif.h>
#endif
!tag is never used, but must be initialized to non-negative integer
  integer, private :: tag=1, stat(MPI_STATUS_SIZE)
!  integer, private, allocatable :: request(:)
  integer, public, allocatable :: request(:)
#ifdef _CRAYT3E
!BWA: mpif.h on t3e currently does not contain MPI_INTEGER8 datatype
!(O2k and t90 do)
!(t3e: fixed on 3.3 I believe)
  integer, parameter :: MPI_INTEGER8=MPI_INTEGER
#endif
#endif use_libMPI

!mpp_stack is used by SHMEM collective ops
!must be SHPALLOC'd on SGICRAY_MPP, but is allocatable on PVP
#ifdef use_shmalloc
  real(DOUBLE_KIND), private :: mpp_stack(1)
  pointer( ptr_stack, mpp_stack )
#else
  real(DOUBLE_KIND), private, allocatable :: mpp_stack(:)
#endif
  integer, private :: mpp_stack_size=0, mpp_stack_hwm=0

!peset hold communicators as SHMEM-compatible triads (start, log2(stride), num)
  type, private :: communicator
     character(len=32) :: name
     integer, pointer :: list(:)
     integer :: count
#ifdef use_libSMA
     integer :: start, log2stride
#elif use_libMPI
     integer :: id, group    !MPI communicator and group id for this PE set
#endif
  end type
  integer, parameter :: PESET_MAX=32 !should be .LE. max num of MPI communicators
  type(communicator) :: peset(0:PESET_MAX) !0 is a dummy used to hold single-PE "self" communicator
  integer :: peset_num=0, current_peset_num=0
  integer :: world_peset_num !the world communicator

!performance profiling
!  This profiles every type of MPI/SHMEM call within
!    a specified region of high-level code
!  Initialize or retrieve a clock with
!  id = mpp_clock_id( 'Region identifier name' )
!  Then set caliper points around the region using:
!  call mpp_clock_begin(id)
!  ...
!  call mpp_clock_end(id)
!  mpp_exit will print out the results.
#ifdef __sgi
#define SYSTEM_CLOCK system_clock_sgi
#endif

#ifdef use_libMPI
#define SYSTEM_CLOCK system_clock_mpi
#endif

#if defined(__sgi) || defined(use_libMPI)
  integer(LONG_KIND), private :: tick, ticks_per_sec, max_ticks, start_tick, end_tick, tick0=0
#else
  integer, private :: tick, ticks_per_sec, max_ticks, start_tick, end_tick, tick0=0
#endif
  real, private :: tick_rate
  integer, private, parameter :: MAX_CLOCKS=100, MAX_EVENT_TYPES=5, MAX_EVENTS=40000
!event types
  integer, private, parameter :: EVENT_ALLREDUCE=1, EVENT_BROADCAST=2, EVENT_RECV=3, EVENT_SEND=4, EVENT_WAIT=5
  integer, private :: clock_num=0, current_clock=0
  integer, private :: clock0    !measures total runtime from mpp_init to mpp_exit
  integer, private :: clock_grain=HUGE(1)
!the event contains information for each type of event (e.g SHMEM_PUT)
  type, private :: event
     character(len=16)  :: name
     integer(LONG_KIND) :: ticks(MAX_EVENTS), bytes(MAX_EVENTS)
     integer            :: calls
  end type event
!a clock contains an array of event profiles for a region
  integer, parameter, public :: MPP_CLOCK_SYNC=1, MPP_CLOCK_DETAILED=2
  type, private :: clock
     character(len=32) :: name
#if defined(__sgi) || defined(use_libMPI)
     integer(LONG_KIND) :: tick
#else
     integer :: tick
#endif
     integer(LONG_KIND) :: total_ticks
     integer :: peset_num
     logical :: sync_on_begin, detailed
     type(event), pointer :: events(:) !if needed, allocate to MAX_EVENT_TYPES
  end type
  type(clock) :: clocks(MAX_CLOCKS)

  integer,parameter :: MAX_BINS=20
  TYPE :: Clock_Data_Summary
    character(len=16) :: name
    real(DOUBLE_KIND) :: msg_size_sums(MAX_BINS)
    real(DOUBLE_KIND) :: msg_time_sums(MAX_BINS)
    real(DOUBLE_KIND) :: total_data
    real(DOUBLE_KIND) :: total_time
    integer(LONG_KIND) :: msg_size_cnts(MAX_BINS)
    integer(LONG_KIND) :: total_cnts
  END TYPE Clock_Data_Summary

  TYPE :: Summary_Struct
    character(len=16)         :: name
    type (Clock_Data_Summary) :: event(MAX_EVENT_TYPES)
  END TYPE Summary_Struct
  type(Summary_Struct) :: clock_summary(MAX_CLOCKS)
  
!public interfaces
! <INTERFACE NAME="mpp_max">
!  <OVERVIEW>
!    Reduction operations.
!  </OVERVIEW>
!  <DESCRIPTION>
!    Find the max of scalar a the PEs in pelist
!    result is also automatically broadcast to all PEs
!  </DESCRIPTION>
!  <TEMPLATE>
!    call  mpp_max( a, pelist )
!  </TEMPLATE>
!  <IN NAME="a">
!    <TT>real</TT> or <TT>integer</TT>, of 4-byte of 8-byte kind.
!  </IN>
!  <IN NAME="pelist">
!    If <TT>pelist</TT> is omitted, the context is assumed to be the
!    current pelist. This call implies synchronization across the PEs in
!    <TT>pelist</TT>, or the current pelist if <TT>pelist</TT> is absent.
!  </IN>
! </INTERFACE>
  interface mpp_max
     module procedure mpp_max_real8
#ifndef no_8byte_integers
     module procedure mpp_max_int8
#endif
#ifndef no_4byte_reals
     module procedure mpp_max_real4
#endif
     module procedure mpp_max_int4
  end interface
  interface mpp_min
     module procedure mpp_min_real8
#ifndef no_8byte_integers
     module procedure mpp_min_int8
#endif
#ifndef no_4byte_reals
     module procedure mpp_min_real4
#endif
     module procedure mpp_min_int4
  end interface

! <INTERFACE NAME="mpp_sum">
!  <OVERVIEW>
!    Reduction operation.
!  </OVERVIEW>
!  <DESCRIPTION>
!    <TT>MPP_TYPE_</TT> corresponds to any 4-byte and 8-byte variant of
!    <TT>integer, real, complex</TT> variables, of rank 0 or 1. A
!    contiguous block from a multi-dimensional array may be passed by its
!    starting address and its length, as in <TT>f77</TT>.
!
!    Library reduction operators are not required or guaranteed to be
!    bit-reproducible. In any case, changing the processor count changes
!    the data layout, and thus very likely the order of operations. For
!    bit-reproducible sums of distributed arrays, consider using the
!    <TT>mpp_global_sum</TT> routine provided by the <LINK
!    SRC="mpp_domains.html"><TT>mpp_domains</TT></LINK> module.
!
!    The <TT>bit_reproducible</TT> flag provided in earlier versions of
!    this routine has been removed.
!
!
!    If <TT>pelist</TT> is omitted, the context is assumed to be the
!    current pelist. This call implies synchronization across the PEs in
!    <TT>pelist</TT>, or the current pelist if <TT>pelist</TT> is absent.
!  </DESCRIPTION>
!  <TEMPLATE>
!    call mpp_sum( a, length, pelist )
!  </TEMPLATE>
!  <IN NAME="length"></IN>
!  <IN NAME="pelist"></IN>
!  <INOUT NAME="a"></INOUT>
! </INTERFACE>
  interface mpp_sum
#ifndef no_8byte_integers
     module procedure mpp_sum_int8
     module procedure mpp_sum_int8_scalar
     module procedure mpp_sum_int8_2d
     module procedure mpp_sum_int8_3d
     module procedure mpp_sum_int8_4d
     module procedure mpp_sum_int8_5d
#endif
     module procedure mpp_sum_real8
     module procedure mpp_sum_real8_scalar
     module procedure mpp_sum_real8_2d
     module procedure mpp_sum_real8_3d
     module procedure mpp_sum_real8_4d
     module procedure mpp_sum_real8_5d
     module procedure mpp_sum_cmplx8
     module procedure mpp_sum_cmplx8_scalar
     module procedure mpp_sum_cmplx8_2d
     module procedure mpp_sum_cmplx8_3d
     module procedure mpp_sum_cmplx8_4d
     module procedure mpp_sum_cmplx8_5d
     module procedure mpp_sum_int4
     module procedure mpp_sum_int4_scalar
     module procedure mpp_sum_int4_2d
     module procedure mpp_sum_int4_3d
     module procedure mpp_sum_int4_4d
     module procedure mpp_sum_int4_5d
#ifndef no_4byte_reals
     module procedure mpp_sum_real4
     module procedure mpp_sum_real4_scalar
     module procedure mpp_sum_real4_2d
     module procedure mpp_sum_real4_3d
     module procedure mpp_sum_real4_4d
     module procedure mpp_sum_real4_5d
     module procedure mpp_sum_cmplx4
     module procedure mpp_sum_cmplx4_scalar
     module procedure mpp_sum_cmplx4_2d
     module procedure mpp_sum_cmplx4_3d
     module procedure mpp_sum_cmplx4_4d
     module procedure mpp_sum_cmplx4_5d
#endif
  end interface

! <INTERFACE NAME="mpp_transmit">
!  <OVERVIEW>
!    Basic message-passing call.
!  </OVERVIEW>
!  <DESCRIPTION>
!    <TT>MPP_TYPE_</TT> corresponds to any 4-byte and 8-byte variant of
!    <TT>integer, real, complex, logical</TT> variables, of rank 0 or 1. A
!    contiguous block from a multi-dimensional array may be passed by its
!    starting address and its length, as in <TT>f77</TT>.
!    
!    <TT>mpp_transmit</TT> is currently implemented as asynchronous
!    outward transmission and synchronous inward transmission. This follows
!    the behaviour of <TT>shmem_put</TT> and <TT>shmem_get</TT>. In MPI, it
!    is implemented as <TT>mpi_isend</TT> and <TT>mpi_recv</TT>. For most
!    applications, transmissions occur in pairs, and are here accomplished
!    in a single call.
!    
!    The special PE designations <TT>NULL_PE</TT>,
!    <TT>ANY_PE</TT> and <TT>ALL_PES</TT> are provided by use
!    association.
!    
!    <TT>NULL_PE</TT>: is used to disable one of the pair of
!    transmissions.<BR/>
!    <TT>ANY_PE</TT>: is used for unspecific remote
!    destination. (Please note that <TT>put_pe=ANY_PE</TT> has no meaning
!    in the MPI context, though it is available in the SHMEM invocation. If
!    portability is a concern, it is best avoided).<BR/>
!    <TT>ALL_PES</TT>: is used for broadcast operations.
!    
!    It is recommended that <LINK
!    SRC="#mpp_broadcast"><TT>mpp_broadcast</TT></LINK> be used for
!    broadcasts.
!    
!    The following example illustrates the use of
!    <TT>NULL_PE</TT> and <TT>ALL_PES</TT>:
!    
!    <PRE>
!    real, dimension(n) :: a
!    if( pe.EQ.0 )then
!        do p = 1,npes-1
!           call mpp_transmit( a, n, p, a, n, NULL_PE )
!        end do
!    else
!        call mpp_transmit( a, n, NULL_PE, a, n, 0 )
!    end if
!    
!    call mpp_transmit( a, n, ALL_PES, a, n, 0 )
!    </PRE>
!    
!    The do loop and the broadcast operation above are equivalent.
!    
!    Two overloaded calls <TT>mpp_send</TT> and
!     <TT>mpp_recv</TT> have also been
!    provided. <TT>mpp_send</TT> calls <TT>mpp_transmit</TT>
!    with <TT>get_pe=NULL_PE</TT>. <TT>mpp_recv</TT> calls
!    <TT>mpp_transmit</TT> with <TT>put_pe=NULL_PE</TT>. Thus
!    the do loop above could be written more succinctly:
!    
!    <PRE>
!    if( pe.EQ.0 )then
!        do p = 1,npes-1
!           call mpp_send( a, n, p )
!        end do
!    else
!        call mpp_recv( a, n, 0 )
!    end if
!    </PRE>
!  </DESCRIPTION>
!  <TEMPLATE>
!    call mpp_transmit( put_data, put_len, put_pe, get_data, get_len, get_pe )
!  </TEMPLATE>
! </INTERFACE>
  interface mpp_transmit
     module procedure mpp_transmit_real8
     module procedure mpp_transmit_real8_scalar
     module procedure mpp_transmit_real8_2d
     module procedure mpp_transmit_real8_3d
     module procedure mpp_transmit_real8_4d
     module procedure mpp_transmit_real8_5d
     module procedure mpp_transmit_cmplx8
     module procedure mpp_transmit_cmplx8_scalar
     module procedure mpp_transmit_cmplx8_2d
     module procedure mpp_transmit_cmplx8_3d
     module procedure mpp_transmit_cmplx8_4d
     module procedure mpp_transmit_cmplx8_5d
#ifndef no_8byte_integers
     module procedure mpp_transmit_int8
     module procedure mpp_transmit_int8_scalar
     module procedure mpp_transmit_int8_2d
     module procedure mpp_transmit_int8_3d
     module procedure mpp_transmit_int8_4d
     module procedure mpp_transmit_int8_5d
     module procedure mpp_transmit_logical8
     module procedure mpp_transmit_logical8_scalar
     module procedure mpp_transmit_logical8_2d
     module procedure mpp_transmit_logical8_3d
     module procedure mpp_transmit_logical8_4d
     module procedure mpp_transmit_logical8_5d
#endif
#ifndef no_4byte_reals
     module procedure mpp_transmit_real4
     module procedure mpp_transmit_real4_scalar
     module procedure mpp_transmit_real4_2d
     module procedure mpp_transmit_real4_3d
     module procedure mpp_transmit_real4_4d
     module procedure mpp_transmit_real4_5d
     module procedure mpp_transmit_cmplx4
     module procedure mpp_transmit_cmplx4_scalar
     module procedure mpp_transmit_cmplx4_2d
     module procedure mpp_transmit_cmplx4_3d
     module procedure mpp_transmit_cmplx4_4d
     module procedure mpp_transmit_cmplx4_5d
#endif
     module procedure mpp_transmit_int4
     module procedure mpp_transmit_int4_scalar
     module procedure mpp_transmit_int4_2d
     module procedure mpp_transmit_int4_3d
     module procedure mpp_transmit_int4_4d
     module procedure mpp_transmit_int4_5d
     module procedure mpp_transmit_logical4
     module procedure mpp_transmit_logical4_scalar
     module procedure mpp_transmit_logical4_2d
     module procedure mpp_transmit_logical4_3d
     module procedure mpp_transmit_logical4_4d
     module procedure mpp_transmit_logical4_5d
  end interface
  interface mpp_recv
     module procedure mpp_recv_real8
     module procedure mpp_recv_real8_scalar
     module procedure mpp_recv_real8_2d
     module procedure mpp_recv_real8_3d
     module procedure mpp_recv_real8_4d
     module procedure mpp_recv_real8_5d
     module procedure mpp_recv_cmplx8
     module procedure mpp_recv_cmplx8_scalar
     module procedure mpp_recv_cmplx8_2d
     module procedure mpp_recv_cmplx8_3d
     module procedure mpp_recv_cmplx8_4d
     module procedure mpp_recv_cmplx8_5d
#ifndef no_8byte_integers
     module procedure mpp_recv_int8
     module procedure mpp_recv_int8_scalar
     module procedure mpp_recv_int8_2d
     module procedure mpp_recv_int8_3d
     module procedure mpp_recv_int8_4d
     module procedure mpp_recv_int8_5d
     module procedure mpp_recv_logical8
     module procedure mpp_recv_logical8_scalar
     module procedure mpp_recv_logical8_2d
     module procedure mpp_recv_logical8_3d
     module procedure mpp_recv_logical8_4d
     module procedure mpp_recv_logical8_5d
#endif
#ifndef no_4byte_reals
     module procedure mpp_recv_real4
     module procedure mpp_recv_real4_scalar
     module procedure mpp_recv_real4_2d
     module procedure mpp_recv_real4_3d
     module procedure mpp_recv_real4_4d
     module procedure mpp_recv_real4_5d
     module procedure mpp_recv_cmplx4
     module procedure mpp_recv_cmplx4_scalar
     module procedure mpp_recv_cmplx4_2d
     module procedure mpp_recv_cmplx4_3d
     module procedure mpp_recv_cmplx4_4d
     module procedure mpp_recv_cmplx4_5d
#endif
     module procedure mpp_recv_int4
     module procedure mpp_recv_int4_scalar
     module procedure mpp_recv_int4_2d
     module procedure mpp_recv_int4_3d
     module procedure mpp_recv_int4_4d
     module procedure mpp_recv_int4_5d
     module procedure mpp_recv_logical4
     module procedure mpp_recv_logical4_scalar
     module procedure mpp_recv_logical4_2d
     module procedure mpp_recv_logical4_3d
     module procedure mpp_recv_logical4_4d
     module procedure mpp_recv_logical4_5d
  end interface
  interface mpp_send
     module procedure mpp_send_real8
     module procedure mpp_send_real8_scalar
     module procedure mpp_send_real8_2d
     module procedure mpp_send_real8_3d
     module procedure mpp_send_real8_4d
     module procedure mpp_send_real8_5d
     module procedure mpp_send_cmplx8
     module procedure mpp_send_cmplx8_scalar
     module procedure mpp_send_cmplx8_2d
     module procedure mpp_send_cmplx8_3d
     module procedure mpp_send_cmplx8_4d
     module procedure mpp_send_cmplx8_5d
#ifndef no_8byte_integers
     module procedure mpp_send_int8
     module procedure mpp_send_int8_scalar
     module procedure mpp_send_int8_2d
     module procedure mpp_send_int8_3d
     module procedure mpp_send_int8_4d
     module procedure mpp_send_int8_5d
     module procedure mpp_send_logical8
     module procedure mpp_send_logical8_scalar
     module procedure mpp_send_logical8_2d
     module procedure mpp_send_logical8_3d
     module procedure mpp_send_logical8_4d
     module procedure mpp_send_logical8_5d
#endif
#ifndef no_4byte_reals
     module procedure mpp_send_real4
     module procedure mpp_send_real4_scalar
     module procedure mpp_send_real4_2d
     module procedure mpp_send_real4_3d
     module procedure mpp_send_real4_4d
     module procedure mpp_send_real4_5d
     module procedure mpp_send_cmplx4
     module procedure mpp_send_cmplx4_scalar
     module procedure mpp_send_cmplx4_2d
     module procedure mpp_send_cmplx4_3d
     module procedure mpp_send_cmplx4_4d
     module procedure mpp_send_cmplx4_5d
#endif
     module procedure mpp_send_int4
     module procedure mpp_send_int4_scalar
     module procedure mpp_send_int4_2d
     module procedure mpp_send_int4_3d
     module procedure mpp_send_int4_4d
     module procedure mpp_send_int4_5d
     module procedure mpp_send_logical4
     module procedure mpp_send_logical4_scalar
     module procedure mpp_send_logical4_2d
     module procedure mpp_send_logical4_3d
     module procedure mpp_send_logical4_4d
     module procedure mpp_send_logical4_5d
  end interface


! <INTERFACE NAME="mpp_broadcast">

!   <OVERVIEW>
!     Parallel broadcasts.
!   </OVERVIEW>
!   <DESCRIPTION>
!     The <TT>mpp_broadcast</TT> call has been added because the original
!     syntax (using <TT>ALL_PES</TT> in <TT>mpp_transmit</TT>) did not
!     support a broadcast across a pelist.
!
!     <TT>MPP_TYPE_</TT> corresponds to any 4-byte and 8-byte variant of
!     <TT>integer, real, complex, logical</TT> variables, of rank 0 or 1. A
!     contiguous block from a multi-dimensional array may be passed by its
!     starting address and its length, as in <TT>f77</TT>.
!
!     Global broadcasts through the <TT>ALL_PES</TT> argument to <LINK
!     SRC="#mpp_transmit"><TT>mpp_transmit</TT></LINK> are still provided for
!     backward-compatibility.
!
!     If <TT>pelist</TT> is omitted, the context is assumed to be the
!     current pelist. <TT>from_pe</TT> must belong to the current
!     pelist. This call implies synchronization across the PEs in
!     <TT>pelist</TT>, or the current pelist if <TT>pelist</TT> is absent.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call mpp_broadcast( data, length, from_pe, pelist )
!   </TEMPLATE>
!   <IN NAME="length"> </IN>
!   <IN NAME="from_pe"> </IN>
!   <IN NAME="pelist"> </IN>
!   <INOUT NAME="data(*)"> </INOUT>
! </INTERFACE>
  interface mpp_broadcast
     module procedure mpp_broadcast_real8
     module procedure mpp_broadcast_real8_scalar
     module procedure mpp_broadcast_real8_2d
     module procedure mpp_broadcast_real8_3d
     module procedure mpp_broadcast_real8_4d
     module procedure mpp_broadcast_real8_5d
     module procedure mpp_broadcast_cmplx8
     module procedure mpp_broadcast_cmplx8_scalar
     module procedure mpp_broadcast_cmplx8_2d
     module procedure mpp_broadcast_cmplx8_3d
     module procedure mpp_broadcast_cmplx8_4d
     module procedure mpp_broadcast_cmplx8_5d
#ifndef no_8byte_integers
     module procedure mpp_broadcast_int8
     module procedure mpp_broadcast_int8_scalar
     module procedure mpp_broadcast_int8_2d
     module procedure mpp_broadcast_int8_3d
     module procedure mpp_broadcast_int8_4d
     module procedure mpp_broadcast_int8_5d
     module procedure mpp_broadcast_logical8
     module procedure mpp_broadcast_logical8_scalar
     module procedure mpp_broadcast_logical8_2d
     module procedure mpp_broadcast_logical8_3d
     module procedure mpp_broadcast_logical8_4d
     module procedure mpp_broadcast_logical8_5d
#endif
#ifndef no_4byte_reals
     module procedure mpp_broadcast_real4
     module procedure mpp_broadcast_real4_scalar
     module procedure mpp_broadcast_real4_2d
     module procedure mpp_broadcast_real4_3d
     module procedure mpp_broadcast_real4_4d
     module procedure mpp_broadcast_real4_5d
     module procedure mpp_broadcast_cmplx4
     module procedure mpp_broadcast_cmplx4_scalar
     module procedure mpp_broadcast_cmplx4_2d
     module procedure mpp_broadcast_cmplx4_3d
     module procedure mpp_broadcast_cmplx4_4d
     module procedure mpp_broadcast_cmplx4_5d
#endif
     module procedure mpp_broadcast_int4
     module procedure mpp_broadcast_int4_scalar
     module procedure mpp_broadcast_int4_2d
     module procedure mpp_broadcast_int4_3d
     module procedure mpp_broadcast_int4_4d
     module procedure mpp_broadcast_int4_5d
     module procedure mpp_broadcast_logical4
     module procedure mpp_broadcast_logical4_scalar
     module procedure mpp_broadcast_logical4_2d
     module procedure mpp_broadcast_logical4_3d
     module procedure mpp_broadcast_logical4_4d
     module procedure mpp_broadcast_logical4_5d
  end interface

! <INTERFACE NAME="mpp_chksum">

!   <OVERVIEW>
!     Parallel checksums.
!   </OVERVIEW>
!   <DESCRIPTION>
!     <TT>mpp_chksum</TT> is a parallel checksum routine that returns an
!     identical answer for the same array irrespective of how it has been
!     partitioned across processors. <TT>LONG_KIND</TT>is the <TT>KIND</TT>
!     parameter corresponding to long integers (see discussion on
!     OS-dependent preprocessor directives) defined in
!     the header file <TT>os.h</TT>. <TT>MPP_TYPE_</TT> corresponds to any
!     4-byte and 8-byte variant of <TT>integer, real, complex, logical</TT>
!     variables, of rank 0 to 5.
!
!     Integer checksums on FP data use the F90 <TT>TRANSFER()</TT>
!     intrinsic.
!
!     The <LINK SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/shared/chksum/chksum.html">serial checksum module</LINK> is superseded
!     by this function, and is no longer being actively maintained. This
!     provides identical results on a single-processor job, and to perform
!     serial checksums on a single processor of a parallel job, you only
!     need to use the optional <TT>pelist</TT> argument.
!     <PRE>
!     use mpp_mod
!     integer :: pe, chksum
!     real :: a(:)
!     pe = mpp_pe()
!     chksum = mpp_chksum( a, (/pe/) )
!     </PRE>
!
!     The additional functionality of <TT>mpp_chksum</TT> over
!     serial checksums is to compute the checksum across the PEs in
!     <TT>pelist</TT>. The answer is guaranteed to be the same for
!     the same distributed array irrespective of how it has been
!     partitioned.
!
!     If <TT>pelist</TT> is omitted, the context is assumed to be the
!     current pelist. This call implies synchronization across the PEs in
!     <TT>pelist</TT>, or the current pelist if <TT>pelist</TT> is absent.
!   </DESCRIPTION>
!   <TEMPLATE>
!     mpp_chksum( var, pelist )
!   </TEMPLATE>
!   <IN NAME="pelist" TYPE="integer" DIM="(:)"> </IN>
!   <IN NAME="var" TYPE="MPP_TYPE_"> </IN>
! </INTERFACE>
  interface mpp_chksum
#ifndef no_8byte_integers
     module procedure mpp_chksum_i8_1d
     module procedure mpp_chksum_i8_2d
     module procedure mpp_chksum_i8_3d
     module procedure mpp_chksum_i8_4d
#endif
     module procedure mpp_chksum_i4_1d
     module procedure mpp_chksum_i4_2d
     module procedure mpp_chksum_i4_3d
     module procedure mpp_chksum_i4_4d
     module procedure mpp_chksum_r8_0d
     module procedure mpp_chksum_r8_1d
     module procedure mpp_chksum_r8_2d
     module procedure mpp_chksum_r8_3d
     module procedure mpp_chksum_r8_4d
     module procedure mpp_chksum_r8_5d
     module procedure mpp_chksum_c8_0d
     module procedure mpp_chksum_c8_1d
     module procedure mpp_chksum_c8_2d
     module procedure mpp_chksum_c8_3d
     module procedure mpp_chksum_c8_4d
     module procedure mpp_chksum_c8_5d
#ifndef no_4byte_reals
     module procedure mpp_chksum_r4_0d
     module procedure mpp_chksum_r4_1d
     module procedure mpp_chksum_r4_2d
     module procedure mpp_chksum_r4_3d
     module procedure mpp_chksum_r4_4d
     module procedure mpp_chksum_r4_5d
     module procedure mpp_chksum_c4_0d
     module procedure mpp_chksum_c4_1d
     module procedure mpp_chksum_c4_2d
     module procedure mpp_chksum_c4_3d
     module procedure mpp_chksum_c4_4d
     module procedure mpp_chksum_c4_5d
#endif
  end interface

! <INTERFACE NAME="mpp_error">
!  <OVERVIEW>
!    Error handler.
!  </OVERVIEW>
!  <DESCRIPTION>
!    It is strongly recommended that all error exits pass through
!    <TT>mpp_error</TT> to assure the program fails cleanly. An individual
!    PE encountering a <TT>STOP</TT> statement, for instance, can cause the
!    program to hang. The use of the <TT>STOP</TT> statement is strongly
!    discouraged.
!    
!    Calling mpp_error with no arguments produces an immediate error
!    exit, i.e:
!    <PRE>
!    call mpp_error
!    call mpp_error(FATAL)
!    </PRE>
!    are equivalent.
!    
!    The argument order
!    <PRE>
!    call mpp_error( routine, errormsg, errortype )
!    </PRE>
!    is also provided to support legacy code. In this version of the
!    call, none of the arguments may be omitted.
!    
!    The behaviour of <TT>mpp_error</TT> for a <TT>WARNING</TT> can be
!    controlled with an additional call <TT>mpp_set_warn_level</TT>.
!    <PRE>
!    call mpp_set_warn_level(ERROR)
!    </PRE>
!    causes <TT>mpp_error</TT> to treat <TT>WARNING</TT>
!    exactly like <TT>FATAL</TT>.
!    <PRE>
!    call mpp_set_warn_level(WARNING)
!    </PRE>
!    resets to the default behaviour described above.
!    
!    <TT>mpp_error</TT> also has an internal error state which
!    maintains knowledge of whether a warning has been issued. This can be
!    used at startup in a subroutine that checks if the model has been
!    properly configured. You can generate a series of warnings using
!    <TT>mpp_error</TT>, and then check at the end if any warnings has been
!    issued using the function <TT>mpp_error_state()</TT>. If the value of
!    this is <TT>WARNING</TT>, at least one warning has been issued, and
!    the user can take appropriate action:
!    
!    <PRE>
!    if( ... )call mpp_error( WARNING, '...' )
!    if( ... )call mpp_error( WARNING, '...' )
!    if( ... )call mpp_error( WARNING, '...' )
!    ...
!    if( mpp_error_state().EQ.WARNING )call mpp_error( FATAL, '...' )
!    </PRE>
!  </DESCRIPTION>
!  <TEMPLATE>
!    call mpp_error( errortype, routine, errormsg )
!  </TEMPLATE>
!  <IN NAME="errortype">
!    One of <TT>NOTE</TT>, <TT>WARNING</TT> or <TT>FATAL</TT> 
!    (these definitions are acquired by use association).
!    <TT>NOTE</TT> writes <TT>errormsg</TT> to <TT>STDOUT</TT>. 
!    <TT>WARNING</TT> writes <TT>errormsg</TT> to <TT>STDERR</TT>.
!    <TT>FATAL</TT> writes <TT>errormsg</TT> to <TT>STDERR</TT>,
!    and induces a clean error exit with a call stack traceback.
!  </IN>
! </INTERFACE>
  interface mpp_error
     module procedure mpp_error_basic
     module procedure mpp_error_mesg
     module procedure mpp_error_noargs
  end interface

#ifdef use_libSMA
!currently SMA contains no generic shmem_wait for different integer kinds:
!I have inserted one here
  interface shmem_integer_wait
     module procedure shmem_int4_wait_local
     module procedure shmem_int8_wait_local
  end interface
#endif
  public :: mpp_chksum, mpp_max, mpp_min, mpp_sum
  public :: mpp_exit, mpp_init
  public :: mpp_pe, mpp_node, mpp_npes, mpp_root_pe, mpp_set_root_pe, mpp_set_stack_size
  public :: mpp_clock_begin, mpp_clock_end, mpp_clock_id, mpp_clock_set_grain
  public :: mpp_error, mpp_error_state, mpp_set_warn_level
  public :: mpp_sync, mpp_sync_self
  public :: mpp_transmit, mpp_send, mpp_recv, mpp_broadcast
  public :: stdin, stdout, stderr, stdlog
  public :: mpp_declare_pelist, mpp_get_current_pelist, mpp_set_current_pelist
#ifdef use_shmalloc
  public :: mpp_malloc
#endif

  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!       ROUTINES TO INITIALIZE/FINALIZE MPP MODULE: mpp_init, mpp_exit        !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! <SUBROUTINE NAME="mpp_init">
!  <OVERVIEW>
!   Initialize <TT>mpp_mod</TT>.
!  </OVERVIEW>
!  <DESCRIPTION>
!   Called to initialize the <TT>mpp_mod</TT> package. It is recommended
!   that this call be the first executed line in your program. It sets the
!   number of PEs assigned to this run (acquired from the command line, or
!   through the environment variable <TT>NPES</TT>), and associates an ID
!   number to each PE. These can be accessed by calling <LINK
!   SRC="#mpp_npes"><TT>mpp_npes</TT></LINK> and <LINK
!   SRC="#mpp_pe"><TT>mpp_pe</TT></LINK>.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call mpp_init( flags )
!  </TEMPLATE>
!  <IN NAME="flags" TYPE="integer">
!   <TT>flags</TT> can be set to <TT>MPP_VERBOSE</TT> to
!   have <TT>mpp_mod</TT> keep you informed of what it's up to.
!  </IN>
! </SUBROUTINE>
      subroutine mpp_init( flags )
      integer, optional, intent(in) :: flags
!    subroutine mpp_init( flags, in, out, err, log )
!      integer, optional, intent(in) :: flags, in, out, err, log
      integer :: my_pe, num_pes, lngth
      integer :: i
      logical :: opened
#ifdef _CRAYT3E
      intrinsic my_pe
#endif

      if( module_is_initialized )return

#ifdef use_libSMA
      call START_PES(0)         !the argument 0 means extract from environment variable NPES on PVP/SGI, from mpprun -n on t3e
      pe = my_pe()
      node = pe                 !on an SMP this should return node ID rather than PE number.
      npes = num_pes()
#elif use_libMPI
      call MPI_INITIALIZED( opened, error ) !in case called from another MPI package
      if( .NOT.opened )call MPI_INIT(error)
      call MPI_COMM_RANK( MPI_COMM_WORLD, pe,   error )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, npes, error )
      allocate( request(0:npes-1) )
      request(:) = MPI_REQUEST_NULL
#endif
      module_is_initialized = .TRUE.

!PEsets: make defaults illegal
      peset(:)%count = -1
#ifdef use_libSMA
      peset(:)%start = -1
      peset(:)%log2stride = -1
#elif use_libMPI
      peset(:)%id = -1
      peset(:)%group = -1
#endif
!0=single-PE, initialized so that count returns 1
      peset(0)%count = 1
      allocate( peset(0)%list(1) )
      peset(0)%list = pe
#ifdef use_libMPI
      current_peset_num = 0
      peset(0)%id = MPI_COMM_WORLD
      call MPI_COMM_GROUP( MPI_COMM_WORLD, peset(0)%group, error )
#endif
      world_peset_num = get_peset( (/(i,i=0,npes-1)/) )
      current_peset_num = world_peset_num !initialize current PEset to world

!initialize clocks
      call SYSTEM_CLOCK( count=tick0, count_rate=ticks_per_sec, count_max=max_ticks )
      tick_rate = 1./ticks_per_sec
      clock0 = mpp_clock_id( 'Total runtime', flags=MPP_CLOCK_SYNC )

      if( PRESENT(flags) )then
          debug   = flags.EQ.MPP_DEBUG
          verbose = flags.EQ.MPP_VERBOSE .OR. debug
      end if

#ifdef use_libSMA
#ifdef use_shmalloc
!we use shpalloc to ensure all these are remotely accessible
      lngth=0; ptr_sync = LOC(pe)   !null initialization
      call mpp_malloc( ptr_sync,        size(TRANSFER(sync,word)),            lngth )
      lngth=0; ptr_status = LOC(pe)  !null initialization
      call mpp_malloc( ptr_status, npes*size(TRANSFER(status(0),word)),   lngth )
      lngth=0; ptr_remote = LOC(pe) !null initialization
      call mpp_malloc( ptr_remote, npes*size(TRANSFER(remote_data_loc(0),word)), lngth )
      lngth=0; ptr_from = LOC(pe)   !null initialization
      call mpp_malloc( ptr_from,        size(TRANSFER(mpp_from_pe,word)),     lngth )
#else
      allocate( status(0:npes-1) )
      allocate( remote_data_loc(0:npes-1) )
#endif
      sync(:) = SHMEM_SYNC_VALUE
      status(0:npes-1) = MPP_READY
      remote_data_loc(0:npes-1) = MPP_WAIT
      call mpp_set_stack_size(32768) !default initial value
#endif
!logunit: log messages are written to configfile.out by default
      etc_unit=get_unit()
!      write( etcfile,'(a,i4.4)' )trim(etcfile)//'.', pe
      if( pe.EQ.root_pe )open( unit=etc_unit, file=trim(etcfile), status='REPLACE' )
      call mpp_sync()                         
      if( pe.NE.root_pe )open( unit=etc_unit, file=trim(etcfile), status='OLD' )
!if optional argument logunit=stdout, write messages to stdout instead.
!if specifying non-defaults, you must specify units not yet in use.
!      if( PRESENT(in) )then
!          inquire( unit=in, opened=opened )
!          if( opened )call mpp_error( FATAL, 'MPP_INIT: unable to open stdin.' )
!          in_unit=in
!      end if
!      if( PRESENT(out) )then
!          inquire( unit=out, opened=opened )
!          if( opened )call mpp_error( FATAL, 'MPP_INIT: unable to open stdout.' )
!          out_unit=out
!      end if
!      if( PRESENT(err) )then
!          inquire( unit=err, opened=opened )
!          if( opened )call mpp_error( FATAL, 'MPP_INIT: unable to open stderr.' )
!          err_unit=err
!      end if
!      log_unit=get_unit()
!      if( PRESENT(log) )then
!          inquire( unit=log, opened=opened )
!          if( opened .AND. log.NE.out_unit )call mpp_error( FATAL, 'MPP_INIT: unable to open stdlog.' )
!          log_unit=log
!      end if
!!log_unit can be written to only from root_pe, all others write to stdout
!      if( log_unit.NE.out_unit )then
!          inquire( unit=log_unit, opened=opened )
!          if( opened )call mpp_error( FATAL, 'MPP_INIT: specified unit for stdlog already in use.' )
!          if( pe.EQ.root_pe )open( unit=log_unit, file=trim(configfile), status='REPLACE' )
!          call mpp_sync()
!          if( pe.NE.root_pe )open( unit=log_unit, file=trim(configfile), status='OLD' )
!      end if
      if( pe.EQ.root_pe )then
          log_unit = get_unit()
          open( unit=log_unit, file=trim(configfile), status='REPLACE' )
          close(log_unit)
      end if
!messages
      if( verbose )call mpp_error( NOTE, 'MPP_INIT: initializing MPP module...' )
      if( pe.EQ.root_pe )then
          write( stdlog(),'(/a)' )'MPP module '//trim(version)//trim(tagname)
          write( stdlog(),'(a,i4)' )'MPP started with NPES=', npes
#ifdef use_libSMA
          write( stdlog(),'(a)' )'Using SMA (shmem) library for message passing...'
#endif
#ifdef use_libMPI
          write( stdlog(),'(a)' )'Using MPI library for message passing...'
#endif
          write( stdlog(), '(a,es12.4,a,i10,a)' ) &
               'Realtime clock resolution=', tick_rate, ' sec (', ticks_per_sec, ' ticks/sec)'
          write( stdlog(), '(a,es12.4,a,i20,a)' ) &
               'Clock rolls over after ', max_ticks*tick_rate, ' sec (', max_ticks, ' ticks)'
      end if
      call mpp_clock_begin(clock0)

      return
    end subroutine mpp_init

! <FUNCTION NAME="stdin">
!  <OVERVIEW>
!    Standard fortran unit numbers.
!  </OVERVIEW>
!  <DESCRIPTION>
!    This function, as well as stdout(), stderr(), stdlog(), returns the current 
!    standard fortran unit numbers for
!    input, output, error messages and log messages. Log messages, by
!    convention, are written to the file <TT>logfile.out</TT>.
!  </DESCRIPTION>
!  <TEMPLATE>
!   stdin()
!  </TEMPLATE>
! </FUNCTION>
    function stdin()
      integer :: stdin
      stdin = in_unit
      return
    end function stdin

    function stdout()
      integer :: stdout
      stdout = out_unit
      if( pe.NE.root_pe )stdout = etc_unit
      return
    end function stdout

    function stderr()
      integer :: stderr
      stderr = err_unit
      return
    end function stderr

    function stdlog()
      integer :: stdlog
      logical :: opened
      if( pe.EQ.root_pe )then
          inquire( file=trim(configfile), opened=opened )
          if( opened )then
              call FLUSH(log_unit)
          else
              log_unit=get_unit()
              open( unit=log_unit, status='OLD', file=trim(configfile), position='APPEND', err=10 )
          end if
          stdlog = log_unit
      else
          stdlog = etc_unit
      end if
      return
   10 call mpp_error( FATAL, 'STDLOG: unable to open '//trim(configfile)//'.' )
    end function stdlog

! <SUBROUTINE NAME="mpp_exit">
!  <OVERVIEW>
!   Exit <TT>mpp_mod</TT>.
!  </OVERVIEW>
!  <DESCRIPTION>
!   Called at the end of the run, or to re-initialize <TT>mpp_mod</TT>,
!   should you require that for some odd reason.
!
!   This call implies synchronization across all PEs.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call mpp_exit()
!  </TEMPLATE>
! </SUBROUTINE>
    subroutine mpp_exit()
!to be called at the end of a run
      integer :: i, j, k, n, nmax
      real :: t, tmin, tmax, tavg, tstd
      real :: m, mmin, mmax, mavg, mstd
      real :: t_total

      if( .NOT.module_is_initialized )return
      call mpp_clock_end(clock0)
      t_total = clocks(clock0)%total_ticks*tick_rate
      if( clock_num.GT.0 )then
          if( ANY(clocks(1:clock_num)%detailed) )then
              call sum_clock_data; call dump_clock_summary
          end if
          if( pe.EQ.root_pe )then
              write( stdout(),'(/a,i4,a)' ) 'Tabulating mpp_clock statistics across ', npes, ' PEs...'
              if( ANY(clocks(1:clock_num)%detailed) ) &
                   write( stdout(),'(a)' )'   ... see mpp_clock.out.#### for details on individual PEs.'
              write( stdout(),'(/32x,a)' ) '          tmin          tmax          tavg          tstd  tfrac'
          end if
          do i = 1,clock_num
             call mpp_set_current_pelist() !implied global barrier
             current_peset_num = clocks(i)%peset_num
             if( .NOT.ANY(peset(current_peset_num)%list(:).EQ.pe) )cycle
!times between mpp_clock ticks
             t = clocks(i)%total_ticks*tick_rate
             tmin = t; call mpp_min(tmin)
             tmax = t; call mpp_max(tmax)
             tavg = t; call mpp_sum(tavg); tavg = tavg/mpp_npes()
             tstd = (t-tavg)**2; call mpp_sum(tstd); tstd = sqrt( tstd/mpp_npes() )
             if( pe.EQ.root_pe )write( stdout(),'(a32,4f14.6,f7.3)' ) &
                  clocks(i)%name, tmin, tmax, tavg, tstd, tavg/t_total
          end do
          if( ANY(clocks(1:clock_num)%detailed) .AND. pe.EQ.root_pe )write( stdout(),'(/32x,a)' ) &
               '       tmin       tmax       tavg       tstd       mmin       mmax       mavg       mstd  mavg/tavg'
          do i = 1,clock_num
!messages: bytelengths and times
             if( .NOT.clocks(i)%detailed )cycle
             do j = 1,MAX_EVENT_TYPES
                n = clocks(i)%events(j)%calls; nmax = n
                call mpp_max(nmax)
                if( nmax.NE.0 )then
!don't divide by n because n might be 0
                    m = 0
                    if( n.GT.0 )m = sum(clocks(i)%events(j)%bytes(1:n))
                    mmin = m; call mpp_min(mmin)
                    mmax = m; call mpp_max(mmax)
                    mavg = m; call mpp_sum(mavg); mavg = mavg/mpp_npes()
                    mstd = (m-mavg)**2; call mpp_sum(mstd); mstd = sqrt( mstd/mpp_npes() )
                    t = 0
                    if( n.GT.0 )t = sum(clocks(i)%events(j)%ticks(1:n))*tick_rate
                    tmin = t; call mpp_min(tmin)
                    tmax = t; call mpp_max(tmax)
                    tavg = t; call mpp_sum(tavg); tavg = tavg/mpp_npes()
                    tstd = (t-tavg)**2; call mpp_sum(tstd); tstd = sqrt( tstd/mpp_npes() )
                    if( pe.EQ.root_pe )write( stdout(),'(a32,4f11.3,5es11.3)' ) &
                         trim(clocks(i)%name)//' '//trim(clocks(i)%events(j)%name), &
                         tmin, tmax, tavg, tstd, mmin, mmax, mavg, mstd, mavg/tavg
                end if
             end do
          end do
      end if
      call mpp_set_current_pelist()
      call mpp_sync()
      call mpp_max(mpp_stack_hwm)
      if( pe.EQ.root_pe )write( stdout(),* )'MPP_STACK high water mark=', mpp_stack_hwm
#ifdef use_libMPI
      call MPI_FINALIZE(error)
#endif

      return
    end subroutine mpp_exit

! <FUNCTION NAME="mpp_pe">
!  <OVERVIEW>
!    Returns processor ID.
!  </OVERVIEW>
!  <DESCRIPTION>
!    This returns the unique ID associated with a PE. This number runs
!    between 0 and <TT>npes-1</TT>, where <TT>npes</TT> is the total
!    processor count, returned by <TT>mpp_npes</TT>. For a uniprocessor
!    application this will always return 0.
!  </DESCRIPTION>
!  <TEMPLATE>
!    mpp_pe()
!  </TEMPLATE>
! </FUNCTION>
    function mpp_pe()
      integer :: mpp_pe

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_PE: You must first call mpp_init.' )
      mpp_pe = pe
      return
    end function mpp_pe

    function mpp_node()
      integer :: mpp_node

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_NODE: You must first call mpp_init.' )
      mpp_node = node
      return
    end function mpp_node

! <FUNCTION NAME="mpp_npes">
!  <OVERVIEW>
!    Returns processor count for current pelist.
!  </OVERVIEW>
!  <DESCRIPTION>
!    This returns the number of PEs in the current pelist. For a
!    uniprocessor application, this will always return 1.
!  </DESCRIPTION>
!  <TEMPLATE>
!    mpp_npes()
!  </TEMPLATE>
! </FUNCTION>
    function mpp_npes()
      integer :: mpp_npes

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_NPES: You must first call mpp_init.' )
!      mpp_npes = npes
      mpp_npes = size(peset(current_peset_num)%list)
      return
    end function mpp_npes

    function mpp_root_pe()
      integer :: mpp_root_pe

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_ROOT_PE: You must first call mpp_init.' )
      mpp_root_pe = root_pe
      return
    end function mpp_root_pe

    subroutine mpp_set_root_pe(num)
      integer, intent(in) :: num
      logical :: opened

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_SET_ROOT_PE: You must first call mpp_init.' )
      if( .NOT.(ANY(num.EQ.peset(current_peset_num)%list)) ) &
           call mpp_error( FATAL, 'MPP_SET_ROOT_PE: you cannot set a root PE outside the current pelist.' )
!actions to take if root_pe has changed:
! open log_unit on new root_pe, close it on old root_pe and point its log_unit to stdout.
!      if( num.NE.root_pe )then  !root_pe has changed
!          if( pe.EQ.num )then
!!on the new root_pe
!              if( log_unit.NE.out_unit )then
!                  inquire( unit=log_unit, opened=opened )
!                  if( .NOT.opened )open( unit=log_unit, status='OLD', file=trim(configfile), position='APPEND' )
!              end if
!          else if( pe.EQ.root_pe )then
!!on the old root_pe
!              if( log_unit.NE.out_unit )then
!                  inquire( unit=log_unit, opened=opened )
!                  if( opened )close(log_unit)
!                  log_unit = out_unit
!              end if
!          end if
!      end if
      root_pe = num
      return
    end subroutine mpp_set_root_pe

! <SUBROUTINE NAME="mpp_declare_pelist">
!  <OVERVIEW>
!    Declare a pelist.
!  </OVERVIEW>
!  <DESCRIPTION>
!    This call is written specifically to accommodate a MPI restriction
!    that requires a parent communicator to create a child communicator, In
!    other words: a pelist cannot go off and declare a communicator, but
!    every PE in the parent, including those not in pelist(:), must get
!    together for the <TT>MPI_COMM_CREATE</TT> call. The parent is
!    typically <TT>MPI_COMM_WORLD</TT>, though it could also be a subset
!    that includes all PEs in <TT>pelist</TT>.
!  
!    The restriction does not apply to SMA but to have uniform code, you
!    may as well call it.
!  
!    This call implies synchronization across the PEs in the current
!    pelist, of which <TT>pelist</TT> is a subset.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call mpp_declare_pelist( pelist,name )
!  </TEMPLATE>
!  <IN NAME="pelist" DIM="(:)" TYPE="integer"></IN>
! </SUBROUTINE>

    subroutine mpp_declare_pelist( pelist, name )
!this call is written specifically to accommodate a brain-dead MPI restriction
!that requires a parent communicator to create a child communicator:
!in other words: a pelist cannot go off and declare a communicator, but every PE
!in the parent, including those not in pelist(:), must get together for the
!MPI_COMM_CREATE call. The parent is typically MPI_COMM_WORLD, though it could also
!be a subset that includes all PEs in pelist.
!This restriction does not apply to SMA but to have uniform code,
!you may as well call it. It must be placed in a context where all PEs call it.
!Subsequent calls that use the pelist should be called PEs in the pelist only.
      integer, intent(in) :: pelist(:)
      character(len=*), optional :: name
      integer :: i

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_DECLARE_PELIST: You must first call mpp_init.' )
      i = get_peset(pelist)
      write( peset(i)%name,'(a,i2.2)' ) 'PElist', i !default name
      if( PRESENT(name) )peset(i)%name = name
      return
    end subroutine mpp_declare_pelist

! <SUBROUTINE NAME="mpp_set_current_pelist">
!  <OVERVIEW>
!    Set context pelist.
!  </OVERVIEW>
!  <DESCRIPTION>
!    This call sets the value of the current pelist, which is the
!    context for all subsequent "global" calls where the optional
!    <TT>pelist</TT> argument is omitted. All the PEs that are to be in the
!    current pelist must call it.
!  
!    In MPI, this call may hang unless <TT>pelist</TT> has been previous
!    declared using <LINK
!    SRC="#mpp_declare_pelist"><TT>mpp_declare_pelist</TT></LINK>.
!  
!    If the argument <TT>pelist</TT> is absent, the current pelist is
!    set to the "world" pelist, of all PEs in the job.
!  </DESCRIPTION>
!  <TEMPLATE>
!    call mpp_set_current_pelist( pelist )
!  </TEMPLATE>
!  <IN NAME="pliest" TYPE="integer"></IN>
! </SUBROUTINE>

    subroutine mpp_set_current_pelist( pelist )
!Once we branch off into a PE subset, we want subsequent "global" calls to
!sync only across this subset. This is declared as the current pelist (peset(current_peset_num)%list)
!when current_peset all pelist ops with no pelist should apply the current pelist.
!also, we set the start PE in this pelist to be the root_pe.
!unlike mpp_declare_pelist, this is called by the PEs in the pelist only
!so if the PEset has not been previously declared, this will hang in MPI.
!if pelist is omitted, we reset pelist to the world pelist.
      integer, intent(in), optional :: pelist(:)
      integer :: i

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_SET_CURRENT_PELIST: You must first call mpp_init.' )
      if( PRESENT(pelist) )then
          if( .NOT.ANY(pe.EQ.pelist) )call mpp_error( FATAL, 'MPP_SET_CURRENT_PELIST: pe must be in pelist.' )
          current_peset_num = get_peset(pelist)
      else
          current_peset_num = world_peset_num
      end if
      call mpp_set_root_pe( MINVAL(peset(current_peset_num)%list) )
      call mpp_sync()           !this is called to make sure everyone in the current pelist is here.
!      npes = mpp_npes()
      return
    end subroutine mpp_set_current_pelist

    subroutine mpp_get_current_pelist( pelist, name )
!this is created for use by mpp_define_domains within a pelist
!will be published but not publicized
      integer, intent(out) :: pelist(:)
      character(len=*), intent(out), optional :: name

      if( size(pelist).NE.size(peset(current_peset_num)%list) ) &
           call mpp_error( FATAL, 'MPP_GET_CURRENT_PELIST: size(pelist) is wrong.' )
      pelist(:) = peset(current_peset_num)%list(:)
      if( PRESENT(name) )name = peset(current_peset_num)%name

      return
    end subroutine mpp_get_current_pelist

    function get_peset(pelist)
      integer :: get_peset
!makes a PE set out of a PE list
!a PE list is an ordered list of PEs
!a PE set is a triad (start,log2stride,size) for SHMEM, an a communicator for MPI
!if stride is non-uniform or not a power of 2, will return error (not required for MPI but enforced for uniformity)
      integer, intent(in), optional :: pelist(:)
      integer :: group
      integer :: i, n, stride
      integer, allocatable :: sorted(:)

      if( .NOT.PRESENT(pelist) )then !set it to current_peset_num
          get_peset = current_peset_num; return
      end if
      if( size(pelist).EQ.1 .AND. npes.GT.1 )then    !collective ops on single PEs should return
          get_peset = 0; return
      end if
!make a sorted list
      n = 1
      if( ascend_sort(pelist).NE.1 )call mpp_error( FATAL, 'GET_PESET: sort error.' )   !result is the array sorted(:)
      if( debug )write( stderr(),* )'pelist=', pelist, ' sorted=', sorted
!find if this array matches any existing peset
      do i = 1,peset_num
         if( debug )write( stderr(),'(a,3i4)' )'pe, i, peset_num=', pe, i, peset_num
         if( size(sorted).EQ.size(peset(i)%list) )then
             if( ALL(sorted.EQ.peset(i)%list) )then
                 deallocate(sorted)
                 get_peset = i; return
             end if
         end if
      end do
!not found, so create new peset
      peset_num = peset_num + 1
      if( peset_num.GE.PESET_MAX )call mpp_error( FATAL, 'GET_PESET: number of PE sets exceeds PESET_MAX.' )
      i = peset_num             !shorthand
!create list
      allocate( peset(i)%list(size(sorted)) )
      peset(i)%list(:) = sorted(:)
      peset(i)%count = size(sorted)
#ifdef use_libSMA
      peset(i)%start = sorted(1)
      if( size(sorted).GT.1 )then
          stride = sorted(2)-sorted(1)
          if( ANY(sorted(2:n)-sorted(1:n-1).NE.stride) ) &
               call mpp_error( WARNING, 'GET_PESET: pelist must have constant stride.' )
          peset(i)%log2stride = nint( log(real(stride))/log(2.) )
          if( 2**peset(i)%log2stride.NE.stride )call mpp_error( WARNING, 'GET_PESET: pelist must have power-of-2 stride.' )
      else
          peset(i)%log2stride = 0
      end if
#elif use_libMPI
      call MPI_GROUP_INCL( peset(current_peset_num)%group, size(sorted), sorted, peset(i)%group, error )
      call MPI_COMM_CREATE( peset(current_peset_num)%id, peset(i)%group, peset(i)%id, error )
#endif
      deallocate(sorted)
      get_peset = i

      return

      contains
        
        recursive function ascend_sort(a) result(a_sort)
          integer :: a_sort
          integer, intent(in) :: a(:)
          integer :: b, i
          if( size(a).EQ.1 .OR. ALL(a.EQ.a(1)) )then
              allocate( sorted(n) )
              sorted(n) = a(1)
              a_sort = n
              return
          end if
          b = minval(a)
          n = n + 1
          i = ascend_sort( pack(a,mask=a.NE.b) )
          a_sort = i - 1
          sorted(i-1) = b
          return
        end function ascend_sort

    end function get_peset

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                        PERFORMANCE PROFILING CALLS                          !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! <SUBROUTINE NAME="mpp_clock_set_grain">
!  <OVERVIEW>
!    Set the level of granularity of timing measurements.   
!  </OVERVIEW>
!  <DESCRIPTION>
!    This routine and three other routines, mpp_clock_id, mpp_clock_begin(id),
!    and mpp_clock_end(id) may be used to time parallel code sections, and
!    extract parallel statistics. Clocks are identified by names, which
!    should be unique in the first 32 characters. The <TT>mpp_clock_id</TT>
!    call initializes a clock of a given name and returns an integer
!    <TT>id</TT>. This <TT>id</TT> can be used by subsequent
!    <TT>mpp_clock_begin</TT> and <TT>mpp_clock_end</TT> calls set around a
!    code section to be timed. Example:
!    <PRE>
!    integer :: id
!    id = mpp_clock_id( 'Atmosphere' )
!    call mpp_clock_begin(id)
!    call atmos_model()
!    call mpp_clock_end()
!    </PRE>  
!     Two flags may be used to alter the behaviour of
!     <TT>mpp_clock</TT>. If the flag <TT>MPP_CLOCK_SYNC</TT> is turned on
!     by <TT>mpp_clock_id</TT>, the clock calls <TT>mpp_sync</TT> across all
!     the PEs in the current pelist at the top of the timed code section,
!     but allows each PE to complete the code section (and reach
!     <TT>mpp_clock_end</TT>) at different times. This allows us to measure
!     load imbalance for a given code section. Statistics are written to
!     <TT>stdout</TT> by <TT>mpp_exit</TT>.
!     
!     The flag <TT>MPP_CLOCK_DETAILED</TT> may be turned on by
!     <TT>mpp_clock_id</TT> to get detailed communication
!     profiles. Communication events of the types <TT>SEND, RECV, BROADCAST,
!     REDUCE</TT> and <TT>WAIT</TT> are separately measured for data volume
!     and time. Statistics are written to <TT>stdout</TT> by
!     <TT>mpp_exit</TT>, and individual PE info is also written to the file
!     <TT>mpp_clock.out.####</TT> where <TT>####</TT> is the PE id given by
!     <TT>mpp_pe</TT>.
!     
!     The flags <TT>MPP_CLOCK_SYNC</TT> and <TT>MPP_CLOCK_DETAILED</TT> are
!     integer parameters available by use association, and may be summed to
!     turn them both on.
!     
!     While the nesting of clocks is allowed, please note that turning on
!     the non-optional flags on inner clocks has certain subtle issues.
!     Turning on <TT>MPP_CLOCK_SYNC</TT> on an inner
!     clock may distort outer clock measurements of load imbalance. Turning
!     on <TT>MPP_CLOCK_DETAILED</TT> will stop detailed measurements on its
!     outer clock, since only one detailed clock may be active at one time.
!     Also, detailed clocks only time a certain number of events per clock
!     (currently 40000) to conserve memory. If this array overflows, a
!     warning message is printed, and subsequent events for this clock are
!     not timed.
!     
!     Timings are done using the <TT>f90</TT> standard
!     <TT>SYSTEM_CLOCK</TT> intrinsic.
!     
!     The resolution of SYSTEM_CLOCK is often too coarse for use except
!     across large swaths of code. On SGI systems this is transparently
!     overloaded with a higher resolution clock made available in a
!     non-portable fortran interface made available by
!     <TT>nsclock.c</TT>. This approach will eventually be extended to other
!     platforms.
!     
!     New behaviour added at the Havana release allows the user to embed
!     profiling calls at varying levels of granularity all over the code,
!     and for any particular run, set a threshold of granularity so that
!     finer-grained clocks become dormant.
!     
!     The threshold granularity is held in the private module variable
!     <TT>clock_grain</TT>. This value may be modified by the call
!     <TT>mpp_clock_set_grain</TT>, and affect clocks initiated by
!     subsequent calls to <TT>mpp_clock_id</TT>. The value of
!     <TT>clock_grain</TT> is set to an arbitrarily large number initially.
!     
!     Clocks initialized by <TT>mpp_clock_id</TT> can set a new optional
!     argument <TT>grain</TT> setting their granularity level. Clocks check
!     this level against the current value of <TT>clock_grain</TT>, and are
!     only triggered if they are <I>at or below ("coarser than")</I> the
!     threshold. Finer-grained clocks are dormant for that run.
!     
!     Note that subsequent changes to <TT>clock_grain</TT> do not
!     change the status of already initiated clocks, and that if the
!     optional <TT>grain</TT> argument is absent, the clock is always
!     triggered. This guarantees backward compatibility.
!  </DESCRIPTION>
!  <TEMPLATE>
!     call mpp_clock_set_grain( grain )
!  </TEMPLATE>
!  <IN NAME="grain" TYPE="integer"></IN>
! </SUBROUTINE>

    subroutine mpp_clock_set_grain( grain )
      integer, intent(in) :: grain
!set the granularity of times: only clocks whose grain is lower than
!clock_grain are triggered, finer-grained clocks are dormant.
!clock_grain is initialized to HUGE(1), so all clocks are triggered if
!this is never called.   
      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_CLOCK_SET_GRAIN: You must first call mpp_init.' )

      clock_grain = grain
      return
    end subroutine mpp_clock_set_grain

    function mpp_clock_id( name, flags, grain )
!return an ID for a new or existing clock
      integer :: mpp_clock_id
      character(len=*), intent(in) :: name
      integer, intent(in), optional :: flags, grain
      integer :: i

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_CLOCK_ID: You must first call mpp_init.' )
!if grain is present, the clock is only triggered if it
!is low ("coarse") enough: compared to clock_grain
!finer-grained clocks are dormant.
!if grain is absent, clock is triggered.
      if( PRESENT(grain) )then
          if( grain.GT.clock_grain )then
              mpp_clock_id = 0
              return
          end if
      end if
      mpp_clock_id = 1
      if( clock_num.EQ.0 )then  !first
!         allocate( clocks(MAX_CLOCKS) )
          clock_num = mpp_clock_id
          call clock_init(mpp_clock_id,name,flags)
      else
          FIND_CLOCK: do while( trim(name).NE.trim(clocks(mpp_clock_id)%name) )
             mpp_clock_id = mpp_clock_id + 1
             if( mpp_clock_id.GT.clock_num )then
                 if( mpp_clock_id.GT.MAX_CLOCKS )then
                     call mpp_error( WARNING, 'MPP_CLOCK_ID: too many clock requests, this one is ignored.' )
                 else               !new clock: initialize
                     clock_num = mpp_clock_id
                     call clock_init(mpp_clock_id,name,flags)
                     exit FIND_CLOCK
                 end if
             end if
          end do FIND_CLOCK
      endif
      return
    end function mpp_clock_id

    subroutine mpp_clock_begin(id)
      integer, intent(in) :: id

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_CLOCK_BEGIN: You must first call mpp_init.' )
      if( id.EQ.0 )return
      if( id.LT.0 .OR. id.GT.clock_num )call mpp_error( FATAL, 'MPP_CLOCK_BEGIN: invalid id.' )

      if( clocks(id)%peset_num.EQ.0 )clocks(id)%peset_num = current_peset_num
      if( clocks(id)%peset_num.NE.current_peset_num ) &
           call mpp_error( FATAL, 'MPP_CLOCK_BEGIN: cannot change pelist context of a clock.' )
      if( clocks(id)%sync_on_begin )then
!do an untimed sync at the beginning of the clock
!this puts all PEs in the current pelist on par, so that measurements begin together
!ending time will be different, thus measuring load imbalance for this clock.
          current_clock = 0; call mpp_sync()
      end if
      current_clock = id
      call SYSTEM_CLOCK( clocks(id)%tick )
      return
    end subroutine mpp_clock_begin

    subroutine mpp_clock_end(id)
!the id argument is not used for anything at present
      integer, intent(in), optional :: id
      integer(LONG_KIND) :: delta

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_CLOCK_END: You must first call mpp_init.' )
      if( id.EQ.0 )return
      if( id.LT.0 .OR. id.GT.clock_num )call mpp_error( FATAL, 'MPP_CLOCK_BEGIN: invalid id.' )
      call SYSTEM_CLOCK(end_tick)
      if( clocks(id)%peset_num.NE.current_peset_num ) &
           call mpp_error( FATAL, 'MPP_CLOCK_END: cannot change pelist context of a clock.' )
      delta = end_tick - clocks(id)%tick
      if( delta.LT.0 )then
          write( stderr(),* )'pe, id, start_tick, end_tick, delta, max_ticks=', pe, id, clocks(id)%tick, end_tick, delta, max_ticks
          delta = delta + max_ticks + 1
          call mpp_error( WARNING, 'MPP_CLOCK_END: Clock rollover, assumed single roll.' )
      end if
      clocks(id)%total_ticks = clocks(id)%total_ticks + delta
      current_clock = 0
      return
    end subroutine mpp_clock_end

    subroutine increment_current_clock( event_id, bytes )
      integer, intent(in) :: event_id
      integer, intent(in), optional :: bytes
      integer :: n
      integer(LONG_KIND) :: delta

      if( current_clock.EQ.0 )return
      if( current_clock.LT.0 .OR. current_clock.GT.clock_num )call mpp_error( FATAL, 'MPP_CLOCK_BEGIN: invalid current_clock.' )
      if( .NOT.clocks(current_clock)%detailed )return
      call SYSTEM_CLOCK(end_tick)
      n = clocks(current_clock)%events(event_id)%calls + 1

      if( n.EQ.MAX_EVENTS )call mpp_error( WARNING, &
           'MPP_CLOCK: events exceed MAX_EVENTS, ignore detailed profiling data for clock '//trim(clocks(current_clock)%name) )
      if( n.GT.MAX_EVENTS )return

      clocks(current_clock)%events(event_id)%calls = n
      delta = end_tick - start_tick
      if( delta.LT.0 )then
          delta = delta + max_ticks + 1
          call mpp_error( WARNING, 'MPP_CLOCK_END: Clock rollover, assumed single roll.' )
      end if
      clocks(current_clock)%events(event_id)%ticks(n) = delta
      if( PRESENT(bytes) )clocks(current_clock)%events(event_id)%bytes(n) = bytes
      return
    end subroutine increment_current_clock

  subroutine dump_clock_summary()
    implicit none

    real :: total_time,total_time_all,total_data
    real :: msg_size,eff_BW,s
    integer :: SD_UNIT
    integer :: total_calls
    integer :: i,j,k,ct
    integer :: msg_cnt
    character(len=2)  :: u
    character(len=18) :: filename
    character(len=20),dimension(MAX_BINS),save :: bin

    data bin( 1)  /'  0   -    8    B:  '/
    data bin( 2)  /'  8   -   16    B:  '/
    data bin( 3)  /' 16   -   32    B:  '/
    data bin( 4)  /' 32   -   64    B:  '/
    data bin( 5)  /' 64   -  128    B:  '/
    data bin( 6)  /'128   -  256    B:  '/
    data bin( 7)  /'256   -  512    B:  '/
    data bin( 8)  /'512   - 1024    B:  '/
    data bin( 9)  /'  1.0 -    2.1 KB:  '/
    data bin(10)  /'  2.1 -    4.1 KB:  '/
    data bin(11)  /'  4.1 -    8.2 KB:  '/
    data bin(12)  /'  8.2 -   16.4 KB:  '/
    data bin(13)  /' 16.4 -   32.8 KB:  '/
    data bin(14)  /' 32.8 -   65.5 KB:  '/
    data bin(15)  /' 65.5 -  131.1 KB:  '/
    data bin(16)  /'131.1 -  262.1 KB:  '/
    data bin(17)  /'262.1 -  524.3 KB:  '/
    data bin(18)  /'524.3 - 1048.6 KB:  '/
    data bin(19)  /'  1.0 -    2.1 MB:  '/
    data bin(20)  /' >2.1          MB:  '/

    if( .NOT.ANY(clocks(1:clock_num)%detailed) )return
    write( filename,'(a,i4.4)' )'mpp_clock.out.', pe

    SD_UNIT = get_unit()
    open(SD_UNIT,file=trim(filename),form='formatted')

    COMM_TYPE: do ct = 1,clock_num

      if( .NOT.clocks(ct)%detailed )cycle
      write(SD_UNIT,*) &
          clock_summary(ct)%name(1:15),' Communication Data for PE ',pe

      write(SD_UNIT,*) ' '
      write(SD_UNIT,*) ' '

      total_time_all = 0.0
      EVENT_TYPE: do k = 1,MAX_EVENT_TYPES-1

        if(clock_summary(ct)%event(k)%total_time == 0.0)cycle

        total_time = clock_summary(ct)%event(k)%total_time
        total_time_all = total_time_all + total_time
        total_data = clock_summary(ct)%event(k)%total_data
        total_calls = clock_summary(ct)%event(k)%total_cnts

        write(SD_UNIT,1000) clock_summary(ct)%event(k)%name(1:9) // ':'

        write(SD_UNIT,1001) 'Total Data: ',total_data*1.0e-6, &
                            'MB; Total Time: ', total_time, &
                            'secs; Total Calls: ',total_calls

        write(SD_UNIT,*) ' '
        write(SD_UNIT,1002) '     Bin            Counts      Avg Size        Eff B/W'
        write(SD_UNIT,*) ' '

        BIN_LOOP: do j=1,MAX_BINS

          if(clock_summary(ct)%event(k)%msg_size_cnts(j)==0)cycle

          if(j<=8)then
            s = 1.0
            u = ' B'
          elseif(j<=18)then
            s = 1.0e-3
            u = 'KB'
          else
            s = 1.0e-6
            u = 'MB'
          endif

          msg_cnt = clock_summary(ct)%event(k)%msg_size_cnts(j)
          msg_size = &
            s*(clock_summary(ct)%event(k)%msg_size_sums(j)/real(msg_cnt))
          eff_BW = (1.0e-6)*( clock_summary(ct)%event(k)%msg_size_sums(j) / &
                                  clock_summary(ct)%event(k)%msg_time_sums(j) )

          write(SD_UNIT,1003) bin(j),msg_cnt,msg_size,u,eff_BW

        end do BIN_LOOP

        write(SD_UNIT,*) ' '
        write(SD_UNIT,*) ' '
      end do EVENT_TYPE

   ! "Data-less" WAIT

      if(clock_summary(ct)%event(MAX_EVENT_TYPES)%total_time>0.0)then

        total_time = clock_summary(ct)%event(MAX_EVENT_TYPES)%total_time
        total_time_all = total_time_all + total_time
        total_calls = clock_summary(ct)%event(MAX_EVENT_TYPES)%total_cnts

        write(SD_UNIT,1000) clock_summary(ct)%event(MAX_EVENT_TYPES)%name(1:9) // ':'

        write(SD_UNIT,1004) 'Total Calls: ',total_calls,'; Total Time: ', &
                             total_time,'secs'

      endif

      write(SD_UNIT,*) ' '
      write(SD_UNIT,1005) 'Total communication time spent for ' // &
                      clock_summary(ct)%name(1:9) // ': ',total_time_all,'secs'
      write(SD_UNIT,*) ' '
      write(SD_UNIT,*) ' '
      write(SD_UNIT,*) ' '

    end do COMM_TYPE

    close(SD_UNIT)

1000  format(a)
1001  format(a,f8.2,a,f8.2,a,i6)
1002  format(a)
1003  format(a,i6,'    ','  ',f6.1,a,'    ',f7.3,'MB/sec')
1004  format(a,i8,a,f9.2,a)
1005  format(a,f9.2,a)
    return
  end subroutine dump_clock_summary

      integer function get_unit()
        implicit none

        integer,save :: i
        logical :: l_open

        i = 10
        do i=10,99
           inquire(unit=i,opened=l_open)
           if(.not.l_open)exit
        end do

        if(i==100)then
            call mpp_error(FATAL,'Unable to get I/O unit')
        else
            get_unit = i
        endif

        return
      end function get_unit

  subroutine sum_clock_data()
    implicit none

    integer :: i,j,k,ct,event_size,event_cnt
    real    :: msg_time

    CLOCK_TYPE: do ct=1,clock_num
      if( .NOT.clocks(ct)%detailed )cycle
      EVENT_TYPE: do j=1,MAX_EVENT_TYPES-1
        event_cnt = clocks(ct)%events(j)%calls
        EVENT_SUMMARY: do i=1,event_cnt

        clock_summary(ct)%event(j)%total_cnts = &
              clock_summary(ct)%event(j)%total_cnts + 1

        event_size = clocks(ct)%events(j)%bytes(i)

        k = find_bin(event_size)

        clock_summary(ct)%event(j)%msg_size_cnts(k) = &
              clock_summary(ct)%event(j)%msg_size_cnts(k) + 1

        clock_summary(ct)%event(j)%msg_size_sums(k) = &
              clock_summary(ct)%event(j)%msg_size_sums(k) &
            + clocks(ct)%events(j)%bytes(i)

        clock_summary(ct)%event(j)%total_data = &
              clock_summary(ct)%event(j)%total_data &
            + clocks(ct)%events(j)%bytes(i)

        msg_time = clocks(ct)%events(j)%ticks(i)
        msg_time = tick_rate * real( clocks(ct)%events(j)%ticks(i) )

        clock_summary(ct)%event(j)%msg_time_sums(k) = &
              clock_summary(ct)%event(j)%msg_time_sums(k) + msg_time

        clock_summary(ct)%event(j)%total_time = &
              clock_summary(ct)%event(j)%total_time + msg_time

        end do EVENT_SUMMARY
      end do EVENT_TYPE

      j = MAX_EVENT_TYPES ! WAITs
           ! "msg_size_cnts" doesn't really mean anything for WAIT
           ! but position will be used to store number of counts for now.

      event_cnt = clocks(ct)%events(j)%calls
      clock_summary(ct)%event(j)%msg_size_cnts(1) = event_cnt
      clock_summary(ct)%event(j)%total_cnts       = event_cnt

      msg_time = tick_rate * real( sum ( clocks(ct)%events(j)%ticks(1:event_cnt) ) )
      clock_summary(ct)%event(j)%msg_time_sums(1) = &
              clock_summary(ct)%event(j)%msg_time_sums(1) + msg_time

      clock_summary(ct)%event(j)%total_time = clock_summary(ct)%event(j)%msg_time_sums(1)

    end do CLOCK_TYPE

    return
    contains
      integer function find_bin(event_size)
        implicit none

        integer,intent(in) :: event_size
        integer :: k,msg_size

        msg_size = 8
        k = 1
        do while(event_size>msg_size .and. k<MAX_BINS)
           k = k+1
           msg_size = msg_size*2
        end do
        find_bin = k
        return
      end function find_bin

  end subroutine sum_clock_data

  subroutine clock_init(id,name,flags)
    integer, intent(in) :: id
    character(len=*), intent(in) :: name
    integer, intent(in), optional :: flags
    integer :: i

    clocks(id)%name = name
    clocks(id)%tick = 0
    clocks(id)%total_ticks = 0
    clocks(id)%sync_on_begin = .FALSE.
    clocks(id)%detailed      = .FALSE.
    clocks(id)%peset_num = 0
    if( PRESENT(flags) )then
        if( BTEST(flags,0) )clocks(id)%sync_on_begin = .TRUE.
        if( BTEST(flags,1) )clocks(id)%detailed      = .TRUE.
    end if
    if( clocks(id)%detailed )then
        allocate( clocks(id)%events(MAX_EVENT_TYPES) )
        clocks(id)%events(EVENT_ALLREDUCE)%name = 'ALLREDUCE'
        clocks(id)%events(EVENT_BROADCAST)%name = 'BROADCAST'
        clocks(id)%events(EVENT_RECV)%name = 'RECV'
        clocks(id)%events(EVENT_SEND)%name = 'SEND'
        clocks(id)%events(EVENT_WAIT)%name = 'WAIT'
        do i=1,MAX_EVENT_TYPES
           clocks(id)%events(i)%ticks(:) = 0
           clocks(id)%events(i)%bytes(:) = 0
           clocks(id)%events(i)%calls = 0
        end do
        clock_summary(id)%name = name
        clock_summary(id)%event(EVENT_ALLREDUCE)%name = 'ALLREDUCE'
        clock_summary(id)%event(EVENT_BROADCAST)%name = 'BROADCAST'
        clock_summary(id)%event(EVENT_RECV)%name = 'RECV'
        clock_summary(id)%event(EVENT_SEND)%name = 'SEND'
        clock_summary(id)%event(EVENT_WAIT)%name = 'WAIT'
        do i=1,MAX_EVENT_TYPES
           clock_summary(id)%event(i)%msg_size_sums(:) = 0.0
           clock_summary(id)%event(i)%msg_time_sums(:) = 0.0
           clock_summary(id)%event(i)%total_data = 0.0
           clock_summary(id)%event(i)%total_time = 0.0
           clock_summary(id)%event(i)%msg_size_cnts(:) = 0
           clock_summary(id)%event(i)%total_cnts = 0
        end do
    end if
    return
  end subroutine clock_init

#ifdef __sgi
    subroutine system_clock_sgi( count, count_rate, count_max )
!mimics F90 SYSTEM_CLOCK intrinsic
      integer(LONG_KIND), intent(out), optional :: count, count_rate, count_max
      integer(LONG_KIND) :: sgi_tick, sgi_ticks_per_sec, sgi_max_tick
!sgi_max_tick currently returns 64
!count must return a number between 0 and count_max
      integer(LONG_KIND), save :: maxtick=0
      if( maxtick.EQ.0 )then
          maxtick = sgi_max_tick() !actually reports #bits in maxtick
          if( maxtick.LT.BIT_SIZE(maxtick) )then
              maxtick = 2**maxtick
          else
              maxtick = huge(maxtick)
          end if
      end if
      if( PRESENT(count) )then
          count = modulo( sgi_tick()-tick0, maxtick )
!          count = sgi_tick()
      end if
      if( PRESENT(count_rate) )then
          count_rate = sgi_ticks_per_sec()
      end if
      if( PRESENT(count_max) )then
          count_max = maxtick-1
      end if
      return
    end subroutine system_clock_sgi
#endif

#ifdef use_libMPI
    subroutine system_clock_mpi( count, count_rate, count_max )
!mimics F90 SYSTEM_CLOCK intrinsic
      integer(LONG_KIND), intent(out), optional :: count, count_rate, count_max
!count must return a number between 0 and count_max
      integer(LONG_KIND), parameter :: maxtick=HUGE(count_max)
      if( PRESENT(count) )then
          count = MPI_WTime()/MPI_WTick()
      end if
      if( PRESENT(count_rate) )then
          count_rate = MPI_Wtick()**(-1)
      end if
      if( PRESENT(count_max) )then
          count_max = maxtick-1
      end if
      return
    end subroutine system_clock_mpi
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                BASIC MESSAGE PASSING ROUTINE: mpp_transmit                  !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define MPP_TRANSMIT_ mpp_transmit_real8
#define MPP_TRANSMIT_SCALAR_ mpp_transmit_real8_scalar
#define MPP_TRANSMIT_2D_ mpp_transmit_real8_2d
#define MPP_TRANSMIT_3D_ mpp_transmit_real8_3d
#define MPP_TRANSMIT_4D_ mpp_transmit_real8_4d
#define MPP_TRANSMIT_5D_ mpp_transmit_real8_5d
#define MPP_RECV_ mpp_recv_real8
#define MPP_RECV_SCALAR_ mpp_recv_real8_scalar
#define MPP_RECV_2D_ mpp_recv_real8_2d
#define MPP_RECV_3D_ mpp_recv_real8_3d
#define MPP_RECV_4D_ mpp_recv_real8_4d
#define MPP_RECV_5D_ mpp_recv_real8_5d
#define MPP_SEND_ mpp_send_real8
#define MPP_SEND_SCALAR_ mpp_send_real8_scalar
#define MPP_SEND_2D_ mpp_send_real8_2d
#define MPP_SEND_3D_ mpp_send_real8_3d
#define MPP_SEND_4D_ mpp_send_real8_4d
#define MPP_SEND_5D_ mpp_send_real8_5d
#define MPP_BROADCAST_ mpp_broadcast_real8
#define MPP_BROADCAST_SCALAR_ mpp_broadcast_real8_scalar
#define MPP_BROADCAST_2D_ mpp_broadcast_real8_2d
#define MPP_BROADCAST_3D_ mpp_broadcast_real8_3d
#define MPP_BROADCAST_4D_ mpp_broadcast_real8_4d
#define MPP_BROADCAST_5D_ mpp_broadcast_real8_5d
#define MPP_TYPE_ real(DOUBLE_KIND)
#define MPP_TYPE_BYTELEN_ 8
#define MPI_TYPE_ MPI_REAL8
#define SHMEM_BROADCAST_ SHMEM_BROADCAST8
#define SHMEM_GET_ SHMEM_GET8
#include <mpp_transmit.h>

#define MPP_TRANSMIT_ mpp_transmit_cmplx8
#define MPP_TRANSMIT_SCALAR_ mpp_transmit_cmplx8_scalar
#define MPP_TRANSMIT_2D_ mpp_transmit_cmplx8_2d
#define MPP_TRANSMIT_3D_ mpp_transmit_cmplx8_3d
#define MPP_TRANSMIT_4D_ mpp_transmit_cmplx8_4d
#define MPP_TRANSMIT_5D_ mpp_transmit_cmplx8_5d
#define MPP_RECV_ mpp_recv_cmplx8
#define MPP_RECV_SCALAR_ mpp_recv_cmplx8_scalar
#define MPP_RECV_2D_ mpp_recv_cmplx8_2d
#define MPP_RECV_3D_ mpp_recv_cmplx8_3d
#define MPP_RECV_4D_ mpp_recv_cmplx8_4d
#define MPP_RECV_5D_ mpp_recv_cmplx8_5d
#define MPP_SEND_ mpp_send_cmplx8
#define MPP_SEND_SCALAR_ mpp_send_cmplx8_scalar
#define MPP_SEND_2D_ mpp_send_cmplx8_2d
#define MPP_SEND_3D_ mpp_send_cmplx8_3d
#define MPP_SEND_4D_ mpp_send_cmplx8_4d
#define MPP_SEND_5D_ mpp_send_cmplx8_5d
#define MPP_BROADCAST_ mpp_broadcast_cmplx8
#define MPP_BROADCAST_SCALAR_ mpp_broadcast_cmplx8_scalar
#define MPP_BROADCAST_2D_ mpp_broadcast_cmplx8_2d
#define MPP_BROADCAST_3D_ mpp_broadcast_cmplx8_3d
#define MPP_BROADCAST_4D_ mpp_broadcast_cmplx8_4d
#define MPP_BROADCAST_5D_ mpp_broadcast_cmplx8_5d
#define MPP_TYPE_ complex(DOUBLE_KIND)
#define MPP_TYPE_BYTELEN_ 16
#define MPI_TYPE_ MPI_DOUBLE_COMPLEX
#define SHMEM_BROADCAST_ SHMEM_BROADCAST8
#define SHMEM_GET_ SHMEM_GET128
#include <mpp_transmit.h>

#ifndef no_4byte_reals
#define MPP_TRANSMIT_ mpp_transmit_real4
#define MPP_TRANSMIT_SCALAR_ mpp_transmit_real4_scalar
#define MPP_TRANSMIT_2D_ mpp_transmit_real4_2d
#define MPP_TRANSMIT_3D_ mpp_transmit_real4_3d
#define MPP_TRANSMIT_4D_ mpp_transmit_real4_4d
#define MPP_TRANSMIT_5D_ mpp_transmit_real4_5d
#define MPP_RECV_ mpp_recv_real4
#define MPP_RECV_SCALAR_ mpp_recv_real4_scalar
#define MPP_RECV_2D_ mpp_recv_real4_2d
#define MPP_RECV_3D_ mpp_recv_real4_3d
#define MPP_RECV_4D_ mpp_recv_real4_4d
#define MPP_RECV_5D_ mpp_recv_real4_5d
#define MPP_SEND_ mpp_send_real4
#define MPP_SEND_SCALAR_ mpp_send_real4_scalar
#define MPP_SEND_2D_ mpp_send_real4_2d
#define MPP_SEND_3D_ mpp_send_real4_3d
#define MPP_SEND_4D_ mpp_send_real4_4d
#define MPP_SEND_5D_ mpp_send_real4_5d
#define MPP_BROADCAST_ mpp_broadcast_real4
#define MPP_BROADCAST_SCALAR_ mpp_broadcast_real4_scalar
#define MPP_BROADCAST_2D_ mpp_broadcast_real4_2d
#define MPP_BROADCAST_3D_ mpp_broadcast_real4_3d
#define MPP_BROADCAST_4D_ mpp_broadcast_real4_4d
#define MPP_BROADCAST_5D_ mpp_broadcast_real4_5d
#define MPP_TYPE_ real(FLOAT_KIND)
#define MPP_TYPE_BYTELEN_ 4
#define MPI_TYPE_ MPI_REAL4
#define SHMEM_BROADCAST_ SHMEM_BROADCAST4
#define SHMEM_GET_ SHMEM_GET4
#include <mpp_transmit.h>

#define MPP_TRANSMIT_ mpp_transmit_cmplx4
#define MPP_TRANSMIT_SCALAR_ mpp_transmit_cmplx4_scalar
#define MPP_TRANSMIT_2D_ mpp_transmit_cmplx4_2d
#define MPP_TRANSMIT_3D_ mpp_transmit_cmplx4_3d
#define MPP_TRANSMIT_4D_ mpp_transmit_cmplx4_4d
#define MPP_TRANSMIT_5D_ mpp_transmit_cmplx4_5d
#define MPP_RECV_ mpp_recv_cmplx4
#define MPP_RECV_SCALAR_ mpp_recv_cmplx4_scalar
#define MPP_RECV_2D_ mpp_recv_cmplx4_2d
#define MPP_RECV_3D_ mpp_recv_cmplx4_3d
#define MPP_RECV_4D_ mpp_recv_cmplx4_4d
#define MPP_RECV_5D_ mpp_recv_cmplx4_5d
#define MPP_SEND_ mpp_send_cmplx4
#define MPP_SEND_SCALAR_ mpp_send_cmplx4_scalar
#define MPP_SEND_2D_ mpp_send_cmplx4_2d
#define MPP_SEND_3D_ mpp_send_cmplx4_3d
#define MPP_SEND_4D_ mpp_send_cmplx4_4d
#define MPP_SEND_5D_ mpp_send_cmplx4_5d
#define MPP_BROADCAST_ mpp_broadcast_cmplx4
#define MPP_BROADCAST_SCALAR_ mpp_broadcast_cmplx4_scalar
#define MPP_BROADCAST_2D_ mpp_broadcast_cmplx4_2d
#define MPP_BROADCAST_3D_ mpp_broadcast_cmplx4_3d
#define MPP_BROADCAST_4D_ mpp_broadcast_cmplx4_4d
#define MPP_BROADCAST_5D_ mpp_broadcast_cmplx4_5d
#define MPP_TYPE_ complex(FLOAT_KIND)
#define MPP_TYPE_BYTELEN_ 8
#define MPI_TYPE_ MPI_COMPLEX
#define SHMEM_BROADCAST_ SHMEM_BROADCAST4
#define SHMEM_GET_ SHMEM_GET64
#include <mpp_transmit.h>
#endif

#ifndef no_8byte_integers
#define MPP_TRANSMIT_ mpp_transmit_int8
#define MPP_TRANSMIT_SCALAR_ mpp_transmit_int8_scalar
#define MPP_TRANSMIT_2D_ mpp_transmit_int8_2d
#define MPP_TRANSMIT_3D_ mpp_transmit_int8_3d
#define MPP_TRANSMIT_4D_ mpp_transmit_int8_4d
#define MPP_TRANSMIT_5D_ mpp_transmit_int8_5d
#define MPP_RECV_ mpp_recv_int8
#define MPP_RECV_SCALAR_ mpp_recv_int8_scalar
#define MPP_RECV_2D_ mpp_recv_int8_2d
#define MPP_RECV_3D_ mpp_recv_int8_3d
#define MPP_RECV_4D_ mpp_recv_int8_4d
#define MPP_RECV_5D_ mpp_recv_int8_5d
#define MPP_SEND_ mpp_send_int8
#define MPP_SEND_SCALAR_ mpp_send_int8_scalar
#define MPP_SEND_2D_ mpp_send_int8_2d
#define MPP_SEND_3D_ mpp_send_int8_3d
#define MPP_SEND_4D_ mpp_send_int8_4d
#define MPP_SEND_5D_ mpp_send_int8_5d
#define MPP_BROADCAST_ mpp_broadcast_int8
#define MPP_BROADCAST_SCALAR_ mpp_broadcast_int8_scalar
#define MPP_BROADCAST_2D_ mpp_broadcast_int8_2d
#define MPP_BROADCAST_3D_ mpp_broadcast_int8_3d
#define MPP_BROADCAST_4D_ mpp_broadcast_int8_4d
#define MPP_BROADCAST_5D_ mpp_broadcast_int8_5d
#define MPP_TYPE_ integer(LONG_KIND)
#define MPP_TYPE_BYTELEN_ 8
#define MPI_TYPE_ MPI_INTEGER8
#define SHMEM_BROADCAST_ SHMEM_BROADCAST8
#define SHMEM_GET_ SHMEM_GET8
#include <mpp_transmit.h>
#endif

#define MPP_TRANSMIT_ mpp_transmit_int4
#define MPP_TRANSMIT_SCALAR_ mpp_transmit_int4_scalar
#define MPP_TRANSMIT_2D_ mpp_transmit_int4_2d
#define MPP_TRANSMIT_3D_ mpp_transmit_int4_3d
#define MPP_TRANSMIT_4D_ mpp_transmit_int4_4d
#define MPP_TRANSMIT_5D_ mpp_transmit_int4_5d
#define MPP_RECV_ mpp_recv_int4
#define MPP_RECV_SCALAR_ mpp_recv_int4_scalar
#define MPP_RECV_2D_ mpp_recv_int4_2d
#define MPP_RECV_3D_ mpp_recv_int4_3d
#define MPP_RECV_4D_ mpp_recv_int4_4d
#define MPP_RECV_5D_ mpp_recv_int4_5d
#define MPP_SEND_ mpp_send_int4
#define MPP_SEND_SCALAR_ mpp_send_int4_scalar
#define MPP_SEND_2D_ mpp_send_int4_2d
#define MPP_SEND_3D_ mpp_send_int4_3d
#define MPP_SEND_4D_ mpp_send_int4_4d
#define MPP_SEND_5D_ mpp_send_int4_5d
#define MPP_BROADCAST_ mpp_broadcast_int4
#define MPP_BROADCAST_SCALAR_ mpp_broadcast_int4_scalar
#define MPP_BROADCAST_2D_ mpp_broadcast_int4_2d
#define MPP_BROADCAST_3D_ mpp_broadcast_int4_3d
#define MPP_BROADCAST_4D_ mpp_broadcast_int4_4d
#define MPP_BROADCAST_5D_ mpp_broadcast_int4_5d
#define MPP_TYPE_ integer(INT_KIND)
#define MPP_TYPE_BYTELEN_ 4
#define MPI_TYPE_ MPI_INTEGER4
#define SHMEM_BROADCAST_ SHMEM_BROADCAST4
#define SHMEM_GET_ SHMEM_GET4
#include <mpp_transmit.h>

#ifndef no_8byte_integers
#define MPP_TRANSMIT_ mpp_transmit_logical8
#define MPP_TRANSMIT_SCALAR_ mpp_transmit_logical8_scalar
#define MPP_TRANSMIT_2D_ mpp_transmit_logical8_2d
#define MPP_TRANSMIT_3D_ mpp_transmit_logical8_3d
#define MPP_TRANSMIT_4D_ mpp_transmit_logical8_4d
#define MPP_TRANSMIT_5D_ mpp_transmit_logical8_5d
#define MPP_RECV_ mpp_recv_logical8
#define MPP_RECV_SCALAR_ mpp_recv_logical8_scalar
#define MPP_RECV_2D_ mpp_recv_logical8_2d
#define MPP_RECV_3D_ mpp_recv_logical8_3d
#define MPP_RECV_4D_ mpp_recv_logical8_4d
#define MPP_RECV_5D_ mpp_recv_logical8_5d
#define MPP_SEND_ mpp_send_logical8
#define MPP_SEND_SCALAR_ mpp_send_logical8_scalar
#define MPP_SEND_2D_ mpp_send_logical8_2d
#define MPP_SEND_3D_ mpp_send_logical8_3d
#define MPP_SEND_4D_ mpp_send_logical8_4d
#define MPP_SEND_5D_ mpp_send_logical8_5d
#define MPP_BROADCAST_ mpp_broadcast_logical8
#define MPP_BROADCAST_SCALAR_ mpp_broadcast_logical8_scalar
#define MPP_BROADCAST_2D_ mpp_broadcast_logical8_2d
#define MPP_BROADCAST_3D_ mpp_broadcast_logical8_3d
#define MPP_BROADCAST_4D_ mpp_broadcast_logical8_4d
#define MPP_BROADCAST_5D_ mpp_broadcast_logical8_5d
#define MPP_TYPE_ logical(LONG_KIND)
#define MPP_TYPE_BYTELEN_ 8
#define MPI_TYPE_ MPI_INTEGER8
#define SHMEM_BROADCAST_ SHMEM_BROADCAST8
#define SHMEM_GET_ SHMEM_GET8
#include <mpp_transmit.h>
#endif

#define MPP_TRANSMIT_ mpp_transmit_logical4
#define MPP_TRANSMIT_SCALAR_ mpp_transmit_logical4_scalar
#define MPP_TRANSMIT_2D_ mpp_transmit_logical4_2d
#define MPP_TRANSMIT_3D_ mpp_transmit_logical4_3d
#define MPP_TRANSMIT_4D_ mpp_transmit_logical4_4d
#define MPP_TRANSMIT_5D_ mpp_transmit_logical4_5d
#define MPP_RECV_ mpp_recv_logical4
#define MPP_RECV_SCALAR_ mpp_recv_logical4_scalar
#define MPP_RECV_2D_ mpp_recv_logical4_2d
#define MPP_RECV_3D_ mpp_recv_logical4_3d
#define MPP_RECV_4D_ mpp_recv_logical4_4d
#define MPP_RECV_5D_ mpp_recv_logical4_5d
#define MPP_SEND_ mpp_send_logical4
#define MPP_SEND_SCALAR_ mpp_send_logical4_scalar
#define MPP_SEND_2D_ mpp_send_logical4_2d
#define MPP_SEND_3D_ mpp_send_logical4_3d
#define MPP_SEND_4D_ mpp_send_logical4_4d
#define MPP_SEND_5D_ mpp_send_logical4_5d
#define MPP_BROADCAST_ mpp_broadcast_logical4
#define MPP_BROADCAST_SCALAR_ mpp_broadcast_logical4_scalar
#define MPP_BROADCAST_2D_ mpp_broadcast_logical4_2d
#define MPP_BROADCAST_3D_ mpp_broadcast_logical4_3d
#define MPP_BROADCAST_4D_ mpp_broadcast_logical4_4d
#define MPP_BROADCAST_5D_ mpp_broadcast_logical4_5d
#define MPP_TYPE_ logical(INT_KIND)
#define MPP_TYPE_BYTELEN_ 4
#define MPI_TYPE_ MPI_INTEGER4
#define SHMEM_BROADCAST_ SHMEM_BROADCAST4
#define SHMEM_GET_ SHMEM_GET4
#include <mpp_transmit.h>

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!            GLOBAL REDUCTION ROUTINES: mpp_max, mpp_sum, mpp_min             !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define MPP_REDUCE_ mpp_max_real8
#define MPP_TYPE_ real(DOUBLE_KIND)
#define SHMEM_REDUCE_ SHMEM_REAL8_MAX_TO_ALL
#define MPI_TYPE_ MPI_REAL8
#define MPI_REDUCE_ MPI_MAX
#include <mpp_reduce.h>

#ifndef no_4byte_reals
#define MPP_REDUCE_ mpp_max_real4
#define MPP_TYPE_ real(FLOAT_KIND)
#define SHMEM_REDUCE_ SHMEM_REAL4_MAX_TO_ALL
#define MPI_TYPE_ MPI_REAL4
#define MPI_REDUCE_ MPI_MAX
#include <mpp_reduce.h>
#endif

#ifndef no_8byte_integers   
#define MPP_REDUCE_ mpp_max_int8
#define MPP_TYPE_ integer(LONG_KIND)
#define SHMEM_REDUCE_ SHMEM_INT8_MAX_TO_ALL
#define MPI_TYPE_ MPI_INTEGER8
#define MPI_REDUCE_ MPI_MAX
#include <mpp_reduce.h>
#endif

#define MPP_REDUCE_ mpp_max_int4
#define MPP_TYPE_ integer(INT_KIND)
#define SHMEM_REDUCE_ SHMEM_INT4_MAX_TO_ALL
#define MPI_TYPE_ MPI_INTEGER4
#define MPI_REDUCE_ MPI_MAX
#include <mpp_reduce.h>

#define MPP_REDUCE_ mpp_min_real8
#define MPP_TYPE_ real(DOUBLE_KIND)
#define SHMEM_REDUCE_ SHMEM_REAL8_MIN_TO_ALL
#define MPI_TYPE_ MPI_REAL8
#define MPI_REDUCE_ MPI_MIN
#include <mpp_reduce.h>

#ifndef no_4byte_reals
#define MPP_REDUCE_ mpp_min_real4
#define MPP_TYPE_ real(FLOAT_KIND)
#define SHMEM_REDUCE_ SHMEM_REAL4_MIN_TO_ALL
#define MPI_TYPE_ MPI_REAL4
#define MPI_REDUCE_ MPI_MIN
#include <mpp_reduce.h>
#endif

#ifndef no_8byte_integers   
#define MPP_REDUCE_ mpp_min_int8
#define MPP_TYPE_ integer(LONG_KIND)
#define SHMEM_REDUCE_ SHMEM_INT8_MIN_TO_ALL
#define MPI_TYPE_ MPI_INTEGER8
#define MPI_REDUCE_ MPI_MIN
#include <mpp_reduce.h>
#endif

#define MPP_REDUCE_ mpp_min_int4
#define MPP_TYPE_ integer(INT_KIND)
#define SHMEM_REDUCE_ SHMEM_INT4_MIN_TO_ALL
#define MPI_TYPE_ MPI_INTEGER4
#define MPI_REDUCE_ MPI_MIN
#include <mpp_reduce.h>

#define MPP_SUM_ mpp_sum_real8
#define MPP_SUM_SCALAR_ mpp_sum_real8_scalar
#define MPP_SUM_2D_ mpp_sum_real8_2d
#define MPP_SUM_3D_ mpp_sum_real8_3d
#define MPP_SUM_4D_ mpp_sum_real8_4d
#define MPP_SUM_5D_ mpp_sum_real8_5d
#define MPP_TYPE_ real(DOUBLE_KIND)
#define SHMEM_SUM_ SHMEM_REAL8_SUM_TO_ALL
#define MPI_TYPE_ MPI_REAL8
#define MPP_TYPE_BYTELEN_ 8
#include <mpp_sum.h>

#define MPP_SUM_ mpp_sum_cmplx8
#define MPP_SUM_SCALAR_ mpp_sum_cmplx8_scalar
#define MPP_SUM_2D_ mpp_sum_cmplx8_2d
#define MPP_SUM_3D_ mpp_sum_cmplx8_3d
#define MPP_SUM_4D_ mpp_sum_cmplx8_4d
#define MPP_SUM_5D_ mpp_sum_cmplx8_5d
#define MPP_TYPE_ complex(DOUBLE_KIND)
#define SHMEM_SUM_ SHMEM_COMP8_SUM_TO_ALL
#define MPI_TYPE_ MPI_DOUBLE_COMPLEX
#define MPP_TYPE_BYTELEN_ 16
#include <mpp_sum.h>

#ifndef no_4byte_reals
#define MPP_SUM_ mpp_sum_real4
#define MPP_SUM_SCALAR_ mpp_sum_real4_scalar
#define MPP_SUM_2D_ mpp_sum_real4_2d
#define MPP_SUM_3D_ mpp_sum_real4_3d
#define MPP_SUM_4D_ mpp_sum_real4_4d
#define MPP_SUM_5D_ mpp_sum_real4_5d
#define MPP_TYPE_ real(FLOAT_KIND)
#define SHMEM_SUM_ SHMEM_REAL4_SUM_TO_ALL
#define MPI_TYPE_ MPI_REAL4
#define MPP_TYPE_BYTELEN_ 4
#include <mpp_sum.h>

#define MPP_SUM_ mpp_sum_cmplx4
#define MPP_SUM_SCALAR_ mpp_sum_cmplx4_scalar
#define MPP_SUM_2D_ mpp_sum_cmplx4_2d
#define MPP_SUM_3D_ mpp_sum_cmplx4_3d
#define MPP_SUM_4D_ mpp_sum_cmplx4_4d
#define MPP_SUM_5D_ mpp_sum_cmplx4_5d
#define MPP_TYPE_ complex(FLOAT_KIND)
#define SHMEM_SUM_ SHMEM_COMP4_SUM_TO_ALL
#define MPI_TYPE_ MPI_COMPLEX
#define MPP_TYPE_BYTELEN_ 8
#include <mpp_sum.h>
#endif

#ifndef no_8byte_integers
#define MPP_SUM_ mpp_sum_int8
#define MPP_SUM_SCALAR_ mpp_sum_int8_scalar
#define MPP_SUM_2D_ mpp_sum_int8_2d
#define MPP_SUM_3D_ mpp_sum_int8_3d
#define MPP_SUM_4D_ mpp_sum_int8_4d
#define MPP_SUM_5D_ mpp_sum_int8_5d
#define MPP_TYPE_ integer(LONG_KIND)
#define SHMEM_SUM_ SHMEM_INT8_SUM_TO_ALL
#define MPI_TYPE_ MPI_INTEGER8
#define MPP_TYPE_BYTELEN_ 8
#include <mpp_sum.h>
#endif

#define MPP_SUM_ mpp_sum_int4
#define MPP_SUM_SCALAR_ mpp_sum_int4_scalar
#define MPP_SUM_2D_ mpp_sum_int4_2d
#define MPP_SUM_3D_ mpp_sum_int4_3d
#define MPP_SUM_4D_ mpp_sum_int4_4d
#define MPP_SUM_5D_ mpp_sum_int4_5d
#define MPP_TYPE_ integer(INT_KIND)
#define SHMEM_SUM_ SHMEM_INT4_SUM_TO_ALL
#define MPI_TYPE_ MPI_INTEGER4
#define MPP_TYPE_BYTELEN_ 4
#include <mpp_sum.h>

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!           SYNCHRONIZATION ROUTINES: mpp_sync, mpp_sync_self                 !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! <SUBROUTINE NAME="mpp_sync">
!  <OVERVIEW>
!    Global synchronization.
!  </OVERVIEW>
!  <DESCRIPTION>
!    Synchronizes PEs at this point in the execution. If
!   <TT>pelist</TT> is omitted all PEs are synchronized. This can
!   be expensive on many systems, and should be avoided if possible. Under
!   MPI, we do not call <TT>MPI_BARRIER</TT>, as you might
!   expect. This is because this call can be prohibitively slow on many
!   systems. Instead, we perform the same operation as
!   <TT>mpp_sync_self</TT>, i.e all participating PEs wait for
!   completion of all their outstanding non-blocking operations.
! 
!   If <TT>pelist</TT> is omitted, the context is assumed to be the
!   current pelist. This call implies synchronization across the PEs in
!   <TT>pelist</TT>, or the current pelist if <TT>pelist</TT> is absent.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call mpp_sync( pelist )
!  </TEMPLATE>
!  <IN NAME="pelist" TYPE="integer" DIM="(:)"></IN>
! </SUBROUTINE>
    subroutine mpp_sync( pelist )
!synchronize PEs in list
      integer, intent(in), optional :: pelist(:)
      integer :: n

      call mpp_sync_self(pelist)

      n = get_peset(pelist); if( peset(n)%count.EQ.1 )return

      if( current_clock.NE.0 )call SYSTEM_CLOCK(start_tick)
#ifdef use_libSMA
      if( n.EQ.world_peset_num )then
          call SHMEM_BARRIER_ALL() !special call is faster
      else
          call SHMEM_BARRIER( peset(n)%start, peset(n)%log2stride, peset(n)%count, sync )
      end if
#endif
#ifdef use_libMPI
      call MPI_BARRIER( peset(n)%id, error )
#endif
      if( current_clock.NE.0 )call increment_current_clock(EVENT_WAIT)

      return
    end subroutine mpp_sync

! <SUBROUTINE NAME="mpp_sync_self">
!  <OVERVIEW>
!    Local synchronization.
!  </OVERVIEW>
!  <DESCRIPTION>
!   <TT>mpp_transmit</TT> is implemented as asynchronous
!   <TT>put/send</TT> and synchronous
!   <TT>get/recv</TT>. <TT>mpp_sync_self</TT> guarantees that outstanding
!   asynchronous operations from the calling PE are complete. If
!   <TT>pelist</TT> is supplied, <TT>mpp_sync_self</TT> checks only for
!   outstanding puts to the PEs in <TT>pelist</TT>.
!
!   If <TT>pelist</TT> is omitted, the context is assumed to be the
!   current pelist. This call implies synchronization across the PEs in
!   <TT>pelist</TT>, or the current pelist if <TT>pelist</TT> is absent.
!  </DESCRIPTION>
!  <IN NAME="pelist" TYPE="integer" DIM="(:)"></IN>
! </SUBROUTINE>
    subroutine mpp_sync_self( pelist )
!this is to check if current PE's outstanding puts are complete
!but we can't use shmem_fence because we are actually waiting for
!a remote PE to complete its get
      integer, intent(in), optional :: pelist(:)
      integer :: i, m, n, stride

      n = get_peset(pelist); if( peset(n)%count.EQ.1 )return

      if( current_clock.NE.0 )call SYSTEM_CLOCK(start_tick)
#ifdef use_libSMA
#ifdef _CRAYT90
      call SHMEM_UDCFLUSH !invalidate data cache
#endif
#endif
      do m = 1,peset(n)%count
         i = peset(n)%list(m)
#ifdef use_libSMA
         call SHMEM_INT8_WAIT( status(i), MPP_WAIT ) !wait for status.NE.MPP_WAIT
#endif
#ifdef use_libMPI
         if( request(i).NE.MPI_REQUEST_NULL )call MPI_WAIT( request(i), stat, error )
#endif
      end do
      if( current_clock.NE.0 )call increment_current_clock(EVENT_WAIT)
      return
    end subroutine mpp_sync_self

#ifdef use_libSMA
!these local versions are written for grouping into shmem_integer_wait
    subroutine shmem_int4_wait_local( ivar, cmp_value )
!dir$ INLINEALWAYS shmem_int4_wait_local
      integer(INT_KIND), intent(in) :: cmp_value
      integer(INT_KIND), intent(inout) :: ivar
      call SHMEM_INT4_WAIT( ivar, cmp_value )
      return
    end subroutine shmem_int4_wait_local
    subroutine shmem_int8_wait_local( ivar, cmp_value )
!dir$ INLINEALWAYS shmem_int8_wait_local
      integer(LONG_KIND), intent(in) :: cmp_value
      integer(LONG_KIND), intent(inout) :: ivar
      call SHMEM_INT8_WAIT( ivar, cmp_value )
      return
    end subroutine shmem_int8_wait_local
#endif
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!         MISCELLANEOUS UTILITIES: mpp_error, mpp_chksum, mpp_malloc          !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine mpp_error_basic( errortype, errormsg )
!a very basic error handler
!uses ABORT and FLUSH calls, may need to use cpp to rename
      integer, intent(in) :: errortype
      character(len=*), intent(in), optional :: errormsg
      character(len=128) :: text
      logical :: opened
      
      if( .NOT.module_is_initialized )call ABORT()

      select case( errortype )
      case(NOTE)
          text = 'NOTE'         !just FYI
      case(WARNING)
          text = 'WARNING'      !probable error
      case(FATAL)
          text = 'FATAL'        !fatal error
      case default
          text = 'WARNING: non-existent errortype (must be NOTE|WARNING|FATAL)'
      end select

      if( npes.GT.1 )write( text,'(a,i5)' )trim(text)//' from PE', pe  !this is the mpp part
      if( PRESENT(errormsg) )text = trim(text)//': '//trim(errormsg)

      select case( errortype )
      case(NOTE)
          write( stdout(),'(a)' )trim(text)
      case default
          write( stderr(),'(/a/)' )trim(text)
          if( errortype.EQ.FATAL .OR. warnings_are_fatal )then
              call FLUSH(stdout())
#ifdef sgi_mipspro
              call TRACE_BACK_STACK_AND_PRINT()
#endif
#ifdef use_libMPI
#ifndef sgi_mipspro
!the call to MPI_ABORT is not trapped by TotalView on sgi
              call MPI_ABORT( MPI_COMM_WORLD, 1, error )
#endif
#endif
              call ABORT()      !automatically calls traceback on Cray systems
          end if
      end select

      error_state = errortype
      return
    end subroutine mpp_error_basic
!overloads to mpp_error_basic
    subroutine mpp_error_mesg( routine, errormsg, errortype )
!support for error_mesg routine in FMS
      character(len=*), intent(in) :: routine, errormsg
      integer, intent(in) :: errortype
      call mpp_error( errortype, trim(routine)//': '//trim(errormsg) )
      return
    end subroutine mpp_error_mesg
    subroutine mpp_error_noargs()
      call mpp_error(FATAL)
    end subroutine mpp_error_noargs
      
    subroutine mpp_set_warn_level(flag)
      integer, intent(in) :: flag

      if( flag.EQ.WARNING )then
          warnings_are_fatal = .FALSE.
      else if( flag.EQ.FATAL )then
          warnings_are_fatal = .TRUE.
      else
          call mpp_error( FATAL, 'MPP_SET_WARN_LEVEL: warning flag must be set to WARNING or FATAL.' )
      end if
      return
    end subroutine mpp_set_warn_level

    function mpp_error_state()
      integer :: mpp_error_state
      mpp_error_state = error_state
      return
    end function mpp_error_state

#ifdef use_shmalloc
! <SUBROUTINE NAME="mpp_malloc">
!  <OVERVIEW>
!    Symmetric memory allocation.
!  </OVERVIEW>
!  <DESCRIPTION>
!    This routine is used on SGI systems when <TT>mpp_mod</TT> is
!    invoked in the SHMEM library. It ensures that dynamically allocated
!    memory can be used with <TT>shmem_get</TT> and
!    <TT>shmem_put</TT>. This is called <I>symmetric
!    allocation</I> and is described in the
!    <TT>intro_shmem</TT> man page. <TT>ptr</TT> is a <I>Cray
!    pointer</I> (see the section on <LINK
!    SRC="#PORTABILITY">portability</LINK>).  The operation can be expensive
!    (since it requires a global barrier). We therefore attempt to re-use
!    existing allocation whenever possible. Therefore <TT>len</TT>
!    and <TT>ptr</TT> must have the <TT>SAVE</TT> attribute
!    in the calling routine, and retain the information about the last call
!    to <TT>mpp_malloc</TT>. Additional memory is symmetrically
!    allocated if and only if <TT>newlen</TT> exceeds
!    <TT>len</TT>.
!
!    This is never required on Cray PVP or MPP systems. While the T3E
!    manpages do talk about symmetric allocation, <TT>mpp_mod</TT>
!    is coded to remove this restriction.
!
!    It is never required if <TT>mpp_mod</TT> is invoked in MPI.
!
!   This call implies synchronization across all PEs.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call mpp_malloc( ptr, newlen, len )
!  </TEMPLATE>
!  <IN NAME="ptr">
!     a cray pointer, points to a dummy argument in this routine.
!  </IN>
!  <IN NAME="newlen" TYPE="integer">
!     the required allocation length for the pointer ptr
!  </IN>
!  <IN NAME="len" TYPE="integer">
!     the current allocation (0 if unallocated).
!  </IN>
! </SUBROUTINE>

    subroutine mpp_malloc( ptr, newlen, len )
!routine to perform symmetric allocation:
!this is required on the t3e/O2k for variables that will be non-local arguments
!to a shmem call (see man intro_shmem(3F)).
!newlen is the required allocation length for the pointer ptr
!   len is the current allocation (0 if unallocated)
      integer, intent(in) :: newlen
      integer, intent(inout) :: len
      real :: dummy
      integer :: words_per_long
      integer(LONG_KIND) :: long
!argument ptr is a cray pointer, points to a dummy argument in this routine
      pointer( ptr, dummy )
!      integer(LONG_KIND) :: error_8

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_MALLOC: You must first call mpp_init.' )
!use existing allocation if it is enough
      if( newlen.LE.len )return

      call SHMEM_BARRIER_ALL()
!if the pointer is already allocated, deallocate
!      if( len.NE.0 )call SHPDEALLC( ptr, error_8, -1 ) !BWA: error_8 instead of error, see PV 682618 (fixed in mpt.1.3.0.1)
      if( len.NE.0 )call SHPDEALLC( ptr, error, -1 )
!allocate new length: assume that the array is KIND=8
      words_per_long = size(transfer(long,word))
      call SHPALLOC( ptr, newlen*words_per_long, error, -1 )
      len = newlen
      call SHMEM_BARRIER_ALL()

      if( debug )then
          call SYSTEM_CLOCK(tick)
          write( stdout(),'(a,i18,a,i5,a,2i8,i16)' )'T=', tick, ' PE=', pe, ' MPP_MALLOC: len, newlen, ptr=', len, newlen, ptr
      end if
      return
    end subroutine mpp_malloc
#endif use_shmalloc

! <SUBROUTINE NAME="mpp_set_stack_size">
!  <OVERVIEW>
!    Allocate module internal workspace.
!  </OVERVIEW>
!  <DESCRIPTION>
!    <TT>mpp_mod</TT> maintains a private internal array called
!    <TT>mpp_stack</TT> for private workspace. This call sets the length,
!    in words, of this array. 
!
!    The <TT>mpp_init</TT> call sets this
!    workspace length to a default of 32768, and this call may be used if a
!    longer workspace is needed.
!    
!    This call implies synchronization across all PEs.
!    
!    This workspace is symmetrically allocated, as required for
!    efficient communication on SGI and Cray MPP systems. Since symmetric
!    allocation must be performed by <I>all</I> PEs in a job, this call
!    must also be called by all PEs, using the same value of
!    <TT>n</TT>. Calling <TT>mpp_set_stack_size</TT> from a subset of PEs,
!    or with unequal argument <TT>n</TT>, may cause the program to hang.
!    
!    If any MPP call using <TT>mpp_stack</TT> overflows the declared
!    stack array, the program will abort with a message specifying the
!    stack length that is required. Many users wonder why, if the required
!    stack length can be computed, it cannot also be specified at that
!    point. This cannot be automated because there is no way for the
!    program to know if all PEs are present at that call, and with equal
!    values of <TT>n</TT>. The program must be rerun by the user with the
!    correct argument to <TT>mpp_set_stack_size</TT>, called at an
!    appropriate point in the code where all PEs are known to be present.
!  </DESCRIPTION>
!  <TEMPLATE>
!    call mpp_set_stack_size(n)
!  </TEMPLATE>
!  <IN NAME="n" TYPE="integer"></IN>
! </SUBROUTINE>
    subroutine mpp_set_stack_size(n)
!set the mpp_stack variable to be at least n LONG words long
      integer, intent(in) :: n
      character(len=8) :: text
#ifdef use_shmalloc
      call mpp_malloc( ptr_stack, n, mpp_stack_size )
#else
      if( n.GT.mpp_stack_size .AND. allocated(mpp_stack) )deallocate(mpp_stack)
      if( .NOT.allocated(mpp_stack) )then
          allocate( mpp_stack(n) )
          mpp_stack_size = n
      end if
#endif
      write( text,'(i8)' )n
      if( pe.EQ.root_pe )call mpp_error( NOTE, 'MPP_SET_STACK_SIZE: stack size set to '//text//'.' )

      return
    end subroutine mpp_set_stack_size

#ifndef no_8byte_integers
#define MPP_CHKSUM_INT_ mpp_chksum_i8_1d
#define MPP_TYPE_ integer(LONG_KIND)
#define MPP_RANK_  (:)
#include <mpp_chksum_int.h>

#define MPP_CHKSUM_INT_ mpp_chksum_i8_2d
#define MPP_TYPE_ integer(LONG_KIND)
#define MPP_RANK_  (:,:)
#include <mpp_chksum_int.h>

#define MPP_CHKSUM_INT_ mpp_chksum_i8_3d
#define MPP_TYPE_ integer(LONG_KIND)
#define MPP_RANK_  (:,:,:)
#include <mpp_chksum_int.h>

#define MPP_CHKSUM_INT_ mpp_chksum_i8_4d
#define MPP_TYPE_ integer(LONG_KIND)
#define MPP_RANK_  (:,:,:,:)
#include <mpp_chksum_int.h>

#define MPP_CHKSUM_INT_ mpp_chksum_i8_5d
#define MPP_TYPE_ integer(LONG_KIND)
#define MPP_RANK_  (:,:,:,:,:)
#include <mpp_chksum_int.h>
#endif

#define MPP_CHKSUM_INT_ mpp_chksum_i4_1d
#define MPP_TYPE_ integer(INT_KIND)
#define MPP_RANK_  (:)
#include <mpp_chksum_int.h>

#define MPP_CHKSUM_INT_ mpp_chksum_i4_2d
#define MPP_TYPE_ integer(INT_KIND)
#define MPP_RANK_  (:,:)
#include <mpp_chksum_int.h>

#define MPP_CHKSUM_INT_ mpp_chksum_i4_3d
#define MPP_TYPE_ integer(INT_KIND)
#define MPP_RANK_  (:,:,:)
#include <mpp_chksum_int.h>

#define MPP_CHKSUM_INT_ mpp_chksum_i4_4d
#define MPP_TYPE_ integer(INT_KIND)
#define MPP_RANK_  (:,:,:,:)
#include <mpp_chksum_int.h>

#define MPP_CHKSUM_INT_ mpp_chksum_i4_5d
#define MPP_TYPE_ integer(INT_KIND)
#define MPP_RANK_  (:,:,:,:,:)
#include <mpp_chksum_int.h>

#define MPP_CHKSUM_ mpp_chksum_r8_0d
#define MPP_TYPE_ real(DOUBLE_KIND)
#define MPP_RANK_ !
#include <mpp_chksum.h>

#define MPP_CHKSUM_ mpp_chksum_r8_1d
#define MPP_TYPE_ real(DOUBLE_KIND)
#define MPP_RANK_ (:)
#include <mpp_chksum.h>

#define MPP_CHKSUM_ mpp_chksum_r8_2d
#define MPP_TYPE_ real(DOUBLE_KIND)
#define MPP_RANK_ (:,:)
#include <mpp_chksum.h>

#define MPP_CHKSUM_ mpp_chksum_r8_3d
#define MPP_TYPE_ real(DOUBLE_KIND)
#define MPP_RANK_ (:,:,:)
#include <mpp_chksum.h>

#define MPP_CHKSUM_ mpp_chksum_r8_4d
#define MPP_TYPE_ real(DOUBLE_KIND)
#define MPP_RANK_ (:,:,:,:)
#include <mpp_chksum.h>

#define MPP_CHKSUM_ mpp_chksum_r8_5d
#define MPP_TYPE_ real(DOUBLE_KIND)
#define MPP_RANK_ (:,:,:,:,:)
#include <mpp_chksum.h>

#define MPP_CHKSUM_ mpp_chksum_c8_0d
#define MPP_TYPE_ complex(DOUBLE_KIND)
#define MPP_RANK_ !
#include <mpp_chksum.h>

#define MPP_CHKSUM_ mpp_chksum_c8_1d
#define MPP_TYPE_ complex(DOUBLE_KIND)
#define MPP_RANK_ (:)
#include <mpp_chksum.h>

#define MPP_CHKSUM_ mpp_chksum_c8_2d
#define MPP_TYPE_ complex(DOUBLE_KIND)
#define MPP_RANK_ (:,:)
#include <mpp_chksum.h>

#define MPP_CHKSUM_ mpp_chksum_c8_3d
#define MPP_TYPE_ complex(DOUBLE_KIND)
#define MPP_RANK_ (:,:,:)
#include <mpp_chksum.h>

#define MPP_CHKSUM_ mpp_chksum_c8_4d
#define MPP_TYPE_ complex(DOUBLE_KIND)
#define MPP_RANK_ (:,:,:,:)
#include <mpp_chksum.h>

#define MPP_CHKSUM_ mpp_chksum_c8_5d
#define MPP_TYPE_ complex(DOUBLE_KIND)
#define MPP_RANK_ (:,:,:,:,:)
#include <mpp_chksum.h>

#ifndef no_4byte_reals
!CAUTION: the r4 versions of these may produce
!unpredictable results: I'm not sure what the result
!of the TRANSFER() to integer(8) is from an odd number of real(4)s?
!However the complex(4) will work, since it is guaranteed even.
#define MPP_CHKSUM_ mpp_chksum_r4_0d
#define MPP_TYPE_ real(FLOAT_KIND)
#define MPP_RANK_ !
#include <mpp_chksum.h>

#define MPP_CHKSUM_ mpp_chksum_r4_1d
#define MPP_TYPE_ real(FLOAT_KIND)
#define MPP_RANK_ (:)
#include <mpp_chksum.h>

#define MPP_CHKSUM_ mpp_chksum_r4_2d
#define MPP_TYPE_ real(FLOAT_KIND)
#define MPP_RANK_ (:,:)
#include <mpp_chksum.h>

#define MPP_CHKSUM_ mpp_chksum_r4_3d
#define MPP_TYPE_ real(FLOAT_KIND)
#define MPP_RANK_ (:,:,:)
#include <mpp_chksum.h>

#define MPP_CHKSUM_ mpp_chksum_r4_4d
#define MPP_TYPE_ real(FLOAT_KIND)
#define MPP_RANK_ (:,:,:,:)
#include <mpp_chksum.h>

#define MPP_CHKSUM_ mpp_chksum_r4_5d
#define MPP_TYPE_ real(FLOAT_KIND)
#define MPP_RANK_ (:,:,:,:,:)
#include <mpp_chksum.h>

#define MPP_CHKSUM_ mpp_chksum_c4_0d
#define MPP_TYPE_ complex(FLOAT_KIND)
#define MPP_RANK_ !
#include <mpp_chksum.h>

#define MPP_CHKSUM_ mpp_chksum_c4_1d
#define MPP_TYPE_ complex(FLOAT_KIND)
#define MPP_RANK_ (:)
#include <mpp_chksum.h>

#define MPP_CHKSUM_ mpp_chksum_c4_2d
#define MPP_TYPE_ complex(FLOAT_KIND)
#define MPP_RANK_ (:,:)
#include <mpp_chksum.h>

#define MPP_CHKSUM_ mpp_chksum_c4_3d
#define MPP_TYPE_ complex(FLOAT_KIND)
#define MPP_RANK_ (:,:,:)
#include <mpp_chksum.h>

#define MPP_CHKSUM_ mpp_chksum_c4_4d
#define MPP_TYPE_ complex(FLOAT_KIND)
#define MPP_RANK_ (:,:,:,:)
#include <mpp_chksum.h>

#define MPP_CHKSUM_ mpp_chksum_c4_5d
#define MPP_TYPE_ complex(FLOAT_KIND)
#define MPP_RANK_ (:,:,:,:,:)
#include <mpp_chksum.h>
#endif

  end module mpp_mod

#ifdef test_mpp
#ifdef SYSTEM_CLOCK
#undef SYSTEM_CLOCK
#endif
  program test
!test various aspects of mpp_mod
#ifdef sgi_mipspro
    use shmem_interface
#endif
    use mpp_mod
    implicit none
    integer :: pe, npes, root
    integer, parameter :: n=1048576
    real, allocatable, dimension(:) :: a, b, c
    integer :: tick, tick0, ticks_per_sec, id
    integer :: i, j, k, l, m
    real :: dt

    call mpp_init()
    call mpp_set_stack_size(3145746)
    pe = mpp_pe()
    npes = mpp_npes()
    root = mpp_root_pe()

    call SYSTEM_CLOCK( count_rate=ticks_per_sec )
    allocate( a(n), b(n) )
    id = mpp_clock_id( 'Random number' )
    call mpp_clock_begin(id)
    call random_number(a)
    call mpp_clock_end  (id)
!time transmit, compare against shmem_put and get
    if( pe.EQ.root )then
        print *, 'Time mpp_transmit for various lengths...'
#ifdef SGICRAY
        print *, 'For comparison, times for shmem_get and shmem_put are also provided.'
#endif
        print *
    end if
    id = mpp_clock_id( 'mpp_transmit' )
    call mpp_clock_begin(id)
!timing is done for cyclical pass (more useful than ping-pong etc)
    l = n
    do while( l.GT.0 )
!mpp_transmit
       call mpp_sync()
       call SYSTEM_CLOCK(tick0)
       do i = 1,npes
          call mpp_transmit( a, l, modulo(pe+npes-i,npes), b, l, modulo(pe+i,npes) )
!          call mpp_sync_self( (/modulo(pe+npes-i,npes)/) )
       end do
       call mpp_sync()
       call SYSTEM_CLOCK(tick)
       dt = real(tick-tick0)/(npes*ticks_per_sec)
       dt = max( dt, epsilon(dt) )
       if( pe.EQ.root )write( stdout(),'(/a,i8,f13.6,f8.2)' )'MPP_TRANSMIT length, time, bw(Mb/s)=', l, dt, l*8e-6/dt
#ifdef SGICRAY
!shmem_put
       call mpp_sync()
       call SYSTEM_CLOCK(tick0)
       do i = 1,npes
          call shmem_real_put( b, a, l, modulo(pe+1,npes) )
       end do
       call mpp_sync()
       call SYSTEM_CLOCK(tick)
       dt = real(tick-tick0)/(npes*ticks_per_sec)
       dt = max( dt, epsilon(dt) )
       if( pe.EQ.root )write( stdout(),'( a,i8,f13.6,f8.2)' )'SHMEM_PUT    length, time, bw(Mb/s)=', l, dt, l*8e-6/dt
!shmem_get
       call mpp_sync()
       call SYSTEM_CLOCK(tick0)
       do i = 1,npes
          call shmem_real_get( b, a, l, modulo(pe+1,npes) )
       end do
       call SYSTEM_CLOCK(tick)
       dt = real(tick-tick0)/(npes*ticks_per_sec)
       dt = max( dt, epsilon(dt) )
       if( pe.EQ.root )write( stdout(),'( a,i8,f13.6,f8.2)' )'SHMEM_GET    length, time, bw(Mb/s)=', l, dt, l*8e-6/dt
#endif
       l = l/2
    end do

!test mpp_sum
    if( pe.EQ.root )then
        print '(/a)', 'Time mpp_sum...'
    end if
    a = real(pe+1)
    call mpp_sync()
    call SYSTEM_CLOCK(tick0)
    call mpp_sum(a,n)
    call SYSTEM_CLOCK(tick)
    dt = real(tick-tick0)/ticks_per_sec
    dt = max( dt, epsilon(dt) )
    if( pe.EQ.root )write( stdout(),'(a,2i4,f9.1,i8,f13.6,f8.2/)' ) &
         'mpp_sum: pe, npes, sum(pe+1), length, time, bw(Mb/s)=', pe, npes, a(1), n, dt, n*8e-6/dt
    call mpp_clock_end(id)

!test mpp_max
    if( pe.EQ.root )then
        print *
        print *, 'Test mpp_max...'
    end if
    a = real(pe+1)
    print *, 'pe,     pe+1 =', pe, a(1)
    call mpp_max( a(1) )
    print *, 'pe, max(pe+1)=', pe, a(1)

!pelist check
    call mpp_sync()
    call flush(stdout())
    if( npes.GE.2 )then
        if( pe.EQ.root )print *, 'Test of pelists: bcast, sum and max using PEs 0...npes-2 (excluding last PE)'
        call mpp_declare_pelist( (/(i,i=0,npes-2)/) )
            
        a = real(pe+1)
        if( pe.NE.npes-1 )call mpp_broadcast( a, n, npes-2, (/(i,i=0,npes-2)/) )
        print *, 'bcast(npes-1) from 0 to npes-2=', pe, a(1)
        a = real(pe+1)
        if( pe.NE.npes-1 )call mpp_sum( a, n, (/(i,i=0,npes-2)/) )
        if( pe.EQ.root )print *, 'sum(pe+1) from 0 to npes-2=', a(1)
        a = real(pe+1)
        if( pe.NE.npes-1 )call mpp_max( a(1), (/(i,i=0,npes-2)/) )
        if( pe.EQ.root )print *, 'max(pe+1) from 0 to npes-2=', a(1)
    end if
#ifdef use_CRI_pointers
!test mpp_chksum
    if( modulo(n,npes).EQ.0 )then  !only set up for even division
        if( pe.EQ.root )call random_number(a)
        call mpp_sync()
        call mpp_transmit( a, n, ALL_PES, a, n, root )
        m= n/npes
        allocate( c(m) )
        c = a(pe*m+1:pe*m+m)

        if( pe.EQ.root )then
            print *
            print *, 'Test mpp_chksum...'
            print *, 'This test shows that a whole array and a distributed array give identical checksums.'
        end if
        print *, 'chksum(a)=', mpp_chksum(a,(/pe/))
        print *, 'chksum(c)=', mpp_chksum(c)
    end if
#endif
    call mpp_exit()
  end program test
#endif test_mpp

! <INFO>

!  <COMPILER NAME="">     
!    Any module or program unit using <TT>mpp_mod</TT> must contain the line
    
!    <PRE>
!    use mpp_mod
!    </PRE>
    
!    The source file for <TT>mpp_mod</TT> is <LINK
!    SRC="ftp://ftp.gfdl.gov/pub/vb/mpp/mpp.F90"><TT>mpp.F90</TT></LINK>.
!    Activate the preprocessor flag <TT>-Duse_libSMA</TT> to invoke
!    the SHMEM library, or <TT>-Duse_libMPI</TT> to invoke the MPI
!    library. Global translation of preprocessor macros is required. This
!    required the activation of the <TT>-F</TT> flag on Cray systems
!    and the <TT>-ftpp -macro_expand</TT> flags on SGI systems. On
!    non-SGI/Cray systems, please consult the f90 manpage for the
!    equivalent flag.
!    
!    On Cray PVP systems, <I>all</I> routines in a message-passing
!    program must be compiled with <TT>-a taskcommon</TT>.
!    
!    On SGI systems, it is required to use 4-byte integers and 8-byte
!    reals, and the 64-bit ABI (<TT>-i4 -r8 -64 -mips4</TT>). It is also
!    required on SGI systems to link the following libraries explicitly:
!    one of <TT>-lmpi</TT> and <TT>-lsma</TT>, depending on
!    whether you wish to use the SHMEM or MPI implementations; and
!    <TT>-lexc</TT>). On Cray systems, all the required flags are
!    default.
!    
!    On SGI, use MIPSPro f90 7.3.1.2 or higher.
!    
!    On Cray, use cf90 3.0.0.0 or higher.
!    
!    On either, use the message-passing toolkit MPT 1.2 or higher.
!    
!    The declaration <TT>MPI_INTEGER8</TT> for 8-byte integers
!    was provided by <TT>mpp_mod</TT> because it was absent in early
!    releases of the Message Passing Toolkit. It has since been included
!    there, and the declaration in <TT>mpp_mod</TT> commented
!    out. This declaration may need to be reinstated if you get a compiler
!    error from this (i.e you are using a superseded version of the MPT).
!    
!    By turning on the cpp flag <TT>-Dtest_mpp</TT> and compiling
!    <TT>mpp_mod</TT> by itself, you may create a test program to
!    exercise certain aspects of <TT>mpp_mod</TT>, e.g
!    
!    <PRE>
!    f90 -F -Duse_libSMA -Dtest_mpp mpp.F90
!    mpprun -n4 a.out
!    </PRE>
!    
!    runs a 4-PE test on a t3e.
!  </COMPILER>
!  <PRECOMP FLAG="">      
!    While the SHMEM library is currently available only on SGI/Cray
!    systems, <TT>mpp_mod</TT> can be used on any other system with
!    a standard-compliant f90 compiler and MPI library. SHMEM is now
!    becoming available on other systems as well.
!    
!    There are some <LINK SRC="os.html">OS-dependent
!    pre-processor directives</LINK> that you might need to modify on
!    non-SGI/Cray systems and compilers.
!    
!    On SGI systems, the <TT>f90</TT> standard <TT>SYSTEM_CLOCK</TT>
!    intrinsic is overloaded with a non-portable fortran interface to a
!    higher-precision clock. This is distributed with the MPP package as
!    <TT>nsclock.c</TT>. This approach will eventually be extended to other
!    platforms, since the resolution of the default clock is often too
!    coarse for our needs.
!  </PRECOMP> 
!  <LOADER FLAG="">       
!    The <TT>mpp</TT> source consists of the main source file
!   <TT>mpp.F90</TT> and also requires the following include files:
!
!   <TT>mpp/shmem.fh</TT> (when compiled with <TT>-Duse_libSMA</TT>)<BR/>
!   <TT>mpif.h</TT> (when compiled with <TT>-Duse_libMPI</TT>)<BR/>
!   <TT>os.h</TT><BR/>
!   <TT>mpp_transmit.h</TT><BR/>
!   <TT>mpp_reduce.h</TT><BR/>
!   <TT>mpp_sum.h</TT><BR/>
!   <TT>mpp_chksum.h</TT><BR/>
!   <TT>mpp_chksum_int.h</TT>
!
!   GFDL users can check it out of the main CVS repository as part of
!   the <TT>mpp</TT> CVS module. The current public tag is <TT>fez</TT>.
!   External users can download the latest <TT>mpp</TT> package <LINK
!   SRC="ftp://ftp.gfdl.gov/pub/vb/mpp/mpp.tar.Z">here</LINK>. Public access
!   to the GFDL CVS repository will soon be made available.
!  </LOADER>
!  <BUG>
!   The <TT>SYSTEM_CLOCK</TT> intrinsic has a limited range before the
!   clock rolls over. The maximum time interval that may be measured
!   before rollover depends on the default integer precision, and is
!   <TT>COUNT_MAX/COUNT_RATE</TT> seconds. Timing a code section longer
!   than this interval will give incorrect results. The <TT>mpp</TT>
!   entry in the logfile reports the rollover time interval. Note that
!   this is a limitation, or "feature" of the <TT>f90 SYSTEM_CLOCK</TT>
!   intrinsic.
!  </BUG>
! </INFO>
