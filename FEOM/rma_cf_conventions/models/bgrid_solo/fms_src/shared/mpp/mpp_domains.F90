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
!   Domain decomposition and domain update for message-passing codes
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
#include <os.h>

!shmalloc is used on MPP SGI/Cray systems for shmem
#if defined(use_libSMA) && defined(SGICRAY_MPP)
#define use_shmalloc
#endif

! <CONTACT EMAIL="vb@gfdl.noaa.gov">
!   V. Balaji
! </CONTACT>

! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>
! <RCSLOG SRC="http://www.gfdl.noaa.gov/~vb/changes_mpp_domains.html"/>

! <OVERVIEW>
!   <TT>mpp_domains_mod</TT> is a set of simple calls for domain
!   decomposition and domain updates on rectilinear grids. It requires the
!   module <LINK SRC="mpp.html"><TT>mpp_mod</TT></LINK>, upon which it is built.
! </OVERVIEW>

! <DESCRIPTION>
!   Scalable implementations of finite-difference codes are generally
!   based on decomposing the model domain into subdomains that are
!   distributed among processors. These domains will then be obliged to
!   exchange data at their boundaries if data dependencies are merely
!   nearest-neighbour, or may need to acquire information from the global
!   domain if there are extended data dependencies, as in the spectral
!   transform. The domain decomposition is a key operation in the
!   development of parallel codes.
!   
!   <TT>mpp_domains_mod</TT> provides a domain decomposition and domain
!   update API for <I>rectilinear</I> grids, built on top of the <LINK
!   SRC="mpp.html"><TT>mpp_mod</TT></LINK> API for message passing. Features
!   of <TT>mpp_domains_mod</TT> include:
! 
!   Simple, minimal API, with free access to underlying API for more complicated stuff.
!
!   Design toward typical use in climate/weather CFD codes.
!  
!   <H4>Domains</H4>
! 
!   I have assumed that domain decomposition will mainly be in 2
!   horizontal dimensions, which will in general be the two
!   fastest-varying indices. There is a separate implementation of 1D
!   decomposition on the fastest-varying index, and 1D decomposition on
!   the second index, treated as a special case of 2D decomposition, is
!   also possible. We define <I>domain</I> as the grid associated with a <I>task</I>.
!   We define the <I>compute domain</I> as the set of gridpoints that are
!   computed by a task, and the <I>data domain</I> as the set of points
!   that are required by the task for the calculation. There can in
!   general be more than 1 task per PE, though often
!   the number of domains is the same as the processor count. We define
!   the <I>global domain</I> as the global computational domain of the
!   entire model (i.e, the same as the computational domain if run on a
!   single processor). 2D domains are defined using a derived type <TT>domain2D</TT>,
!   constructed as follows (see comments in code for more details):
!   
!   <PRE>
!     type, public :: domain_axis_spec
!        private
!        integer :: begin, end, size, max_size
!        logical :: is_global
!     end type domain_axis_spec
!     type, public :: domain1D
!        private
!        type(domain_axis_spec) :: compute, data, global, active
!        logical :: mustputb, mustgetb, mustputf, mustgetf, folded
!        type(domain1D), pointer, dimension(:) :: list
!        integer :: pe              !PE to which this domain is assigned
!        integer :: pos
!     end type domain1D
!domaintypes of higher rank can be constructed from type domain1D
!typically we only need 1 and 2D, but could need higher (e.g 3D LES)
!some elements are repeated below if they are needed once per domain
!     type, public :: domain2D
!        private
!        type(domain1D) :: x
!        type(domain1D) :: y
!        type(domain2D), pointer, dimension(:) :: list
!        integer :: pe              !PE to which this domain is assigned
!        integer :: pos
!     end type domain2D
!     type(domain1D), public :: NULL_DOMAIN1D
!     type(domain2D), public :: NULL_DOMAIN2D
!   </PRE>

!   The <TT>domain2D</TT> type contains all the necessary information to
!   define the global, compute and data domains of each task, as well as the PE
!   associated with the task. The PEs from which remote data may be
!   acquired to update the data domain are also contained in a linked list
!   of neighbours.
! </DESCRIPTION>

module mpp_domains_mod
!a generalized domain decomposition package for use with mpp_mod
!Balaji (vb@gfdl.gov) 15 March 1999
  use mpp_mod
  implicit none
  private
  character(len=128), private :: version= &
       '$Revision$'
  character(len=128), private :: tagname= &
       '$Id$'
  character(len=128), private :: version_update_domains2D, version_global_reduce, version_global_sum, version_global_field

!parameters used to define domains: these are passed to the flags argument of mpp_define_domains
!  if data domain is to have global extent, set GLOBAL_DATA_DOMAIN
!  if global domain has periodic boundaries, set CYCLIC_GLOBAL_DOMAIN
!  sum flags together if more than one of the above conditions is to be met.
  integer, parameter, private :: GLOBAL=0, CYCLIC=1
  integer, parameter, private :: WEST=2, EAST=3, SOUTH=4, NORTH=5
  integer, parameter, private :: SEND=1, RECV=2
  integer, parameter, public :: GLOBAL_DATA_DOMAIN=2**GLOBAL, CYCLIC_GLOBAL_DOMAIN=2**CYCLIC
!gridtypes
  integer, parameter, private :: AGRID=0, BGRID=1, CGRID=2
  integer, parameter, public :: BGRID_NE=BGRID+2**NORTH+2**EAST
  integer, parameter, public :: BGRID_SW=BGRID+2**SOUTH+2**WEST
  integer, parameter, public :: CGRID_NE=CGRID+2**NORTH+2**EAST
  integer, parameter, public :: CGRID_SW=CGRID+2**SOUTH+2**WEST
  integer, private :: grid_offset_type=AGRID
!folds
  integer, parameter, public :: FOLD_WEST_EDGE = 2**WEST, FOLD_EAST_EDGE = 2**EAST
  integer, parameter, public :: FOLD_SOUTH_EDGE=2**SOUTH, FOLD_NORTH_EDGE=2**NORTH
!update
  integer, parameter, public :: WUPDATE=2**WEST, EUPDATE=2**EAST, SUPDATE=2**SOUTH, NUPDATE=2**NORTH
  integer, parameter, public :: XUPDATE=WUPDATE+EUPDATE, YUPDATE=SUPDATE+NUPDATE
!used by mpp_global_sum
  integer, parameter, public :: BITWISE_EXACT_SUM=1

  type, public :: domain_axis_spec        !type used to specify index limits along an axis of a domain
     private
     integer :: begin, end, size, max_size      !start, end of domain axis, size, max size in set
     logical :: is_global       !TRUE if domain axis extent covers global domain
  end type domain_axis_spec
  type, public :: domain1D
     private
     type(domain_axis_spec) :: compute, data, global
     logical :: cyclic
     type(domain1D), pointer :: list(:)
     integer :: pe              !PE to which this domain is assigned
     integer :: pos             !position of this PE within link list, i.e domain%list(pos)%pe = pe
  end type domain1D
!domaintypes of higher rank can be constructed from type domain1D
!typically we only need 1 and 2D, but could need higher (e.g 3D LES)
!some elements are repeated below if they are needed once per domain, not once per axis
  type, private :: rectangle
     integer :: is, ie, js, je
     logical :: overlap, folded
  end type rectangle
  type, public :: domain2D
     private
     type(domain1D) :: x
     type(domain1D) :: y
     type(domain2D), pointer :: list(:)
     integer :: pe              !PE to which this domain is assigned
     integer :: pos             !position of this PE within link list, i.e domain%list(pos)%pe = pe
     integer :: fold, gridtype
     logical :: overlap
     type(rectangle) :: recv_e, recv_se, recv_s, recv_sw, &
                        recv_w, recv_nw, recv_n, recv_ne
     type(rectangle) :: send_e, send_se, send_s, send_sw, &
                        send_w, send_nw, send_n, send_ne
     logical :: remote_domains_initialized
     type(rectangle) :: recv_e_off, recv_se_off, recv_s_off, recv_sw_off, &
                        recv_w_off, recv_nw_off, recv_n_off, recv_ne_off
     type(rectangle) :: send_e_off, send_se_off, send_s_off, send_sw_off, &
                        send_w_off, send_nw_off, send_n_off, send_ne_off
     logical :: remote_off_domains_initialized
  end type domain2D
  type(domain1D), public :: NULL_DOMAIN1D
  type(domain2D), public :: NULL_DOMAIN2D

  integer, private :: pe

  integer, private :: tk
  logical, private :: verbose=.FALSE., debug=.FALSE., domain_clocks_on=.FALSE.
  logical, private :: module_is_initialized=.FALSE.
  integer, parameter, public :: MPP_DOMAIN_TIME=MPP_DEBUG+1
  integer :: send_clock=0, recv_clock=0, unpk_clock=0, wait_clock=0, pack_clock=0, pack_loop_clock=0

!stack for internal buffers
!allocated differently if use_shmalloc
#ifdef use_shmalloc
  real(DOUBLE_KIND), private :: mpp_domains_stack(1)
  pointer( ptr_stack, mpp_domains_stack )
#else
  real(DOUBLE_KIND), private, allocatable :: mpp_domains_stack(:)
#endif
  integer, private :: mpp_domains_stack_size=0, mpp_domains_stack_hwm=0

!used by mpp_define_domains2D_new to transmit data
  integer :: domain_info_buf(16)
#ifdef use_shmalloc
  pointer( ptr_info, domain_info_buf )
#endif

!public interfaces

! this interface can add halo to a no-halo domain and copy everything else.

  interface mpp_copy_domains
     module procedure mpp_copy_domains1D
     module procedure mpp_copy_domains2D
  end interface

! <INTERFACE NAME="mpp_define_domains">

!   <OVERVIEW>
!     Set up a domain decomposition.
!   </OVERVIEW>
!   <DESCRIPTION>
!     There are two forms for the <TT>mpp_define_domains</TT> call. The 2D
!     version is generally to be used but is built by repeated calls to the
!     1D version, also provided.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call mpp_define_domains( global_indices, ndivs, domain, &
!                                   pelist, flags, halo, extent, maskmap )
!   </TEMPLATE>
!  <TEMPLATE>
!    call mpp_define_domains( global_indices, layout, domain, pelist, &
!                                   xflags, yflags, xhalo, yhalo,           &
!                                   xextent, yextent, maskmap, name )
!  </TEMPLATE>
!   <IN NAME="global_indices" >
!     Defines the global domain.
!   </IN>
!   <IN NAME="ndivs">
!     Is the number of domain divisions required.
!   </IN>
!   <INOUT NAME="domain">
!     Holds the resulting domain decomposition.
!   </INOUT>
!   <IN NAME="pelist">
!     List of PEs to which the domains are to be assigned.
!   </IN>
!   <IN NAME="flags">
!      An optional flag to pass additional information
!      about the desired domain topology. Useful flags in a 1D decomposition
!      include <TT>GLOBAL_DATA_DOMAIN</TT> and
!      <TT>CYCLIC_GLOBAL_DOMAIN</TT>. Flags are integers: multiple flags may
!      be added together. The flag values are public parameters available by
!      use association.
!   </IN>
!   <IN NAME="halo">
!     Width of the halo.
!   </IN>
!   <IN NAME="extent">
!      Normally <TT>mpp_define_domains</TT> attempts
!      an even division of the global domain across <TT>ndivs</TT>
!      domains. The <TT>extent</TT> array can be used by the user to pass a
!      custom domain division. The <TT>extent</TT> array has <TT>ndivs</TT>
!      elements and holds the compute domain widths, which should add up to
!      cover the global domain exactly.
!   </IN>
!   <IN NAME="maskmap">
!     Some divisions may be masked
!     (<TT>maskmap=.FALSE.</TT>) to exclude them from the computation (e.g
!     for ocean model domains that are all land). The <TT>maskmap</TT> array
!     is dimensioned <TT>ndivs</TT> and contains <TT>.TRUE.</TT> values for
!     any domain that must be <I>included</I> in the computation (default
!     all). The <TT>pelist</TT> array length should match the number of
!     domains included in the computation.
!    </IN>   

!  <IN NAME="layout"></IN>
!  <IN NAME="xflags, yflags"></IN>
!  <IN NAME="xhalo, yhalo"></IN>
!  <IN NAME="xextent, yextent"></IN>
!  <IN NAME="name" ></IN>

!  <NOTE>    
!    For example:
!    
!    <PRE>
!    call mpp_define_domains( (/1,100/), 10, domain, &
!         flags=GLOBAL_DATA_DOMAIN+CYCLIC_GLOBAL_DOMAIN, halo=2 )
!    </PRE>
!    
!    defines 10 compute domains spanning the range [1,100] of the global
!    domain. The compute domains are non-overlapping blocks of 10. All the data
!    domains are global, and with a halo of 2 span the range [-1:102]. And
!    since the global domain has been declared to be cyclic,
!    <TT>domain(9)%next => domain(0)</TT> and <TT>domain(0)%prev =>
!    domain(9)</TT>. A field is allocated on the data domain, and computations proceed on
!    the compute domain. A call to <LINK
!    SRC="#mpp_update_domains"><TT>mpp_update_domains</TT></LINK> would fill in
!    the values in the halo region:
    
!    <PRE>
!    call mpp_get_data_domain( domain, isd, ied ) !returns -1 and 102
!    call mpp_get_compute_domain( domain, is, ie ) !returns (1,10) on PE 0 ...
!    allocate( a(isd:ied) )
!    do i = is,ie
!       a(i) = &lt;perform computations&gt;
!    end do
!    call mpp_update_domains( a, domain )
!    </PRE>

!    The call to <TT>mpp_update_domains</TT> fills in the regions outside
!    the compute domain. Since the global domain is cyclic, the values at
!    <TT>i=(-1,0)</TT> are the same as at <TT>i=(99,100)</TT>; and
!    <TT>i=(101,102)</TT> are the same as <TT>i=(1,2)</TT>.
!    
!    The 2D version is just an extension of this syntax to two
!    dimensions.
!
!    The 2D version of the above should generally be used in
!    codes, including 1D-decomposed ones, if there is a possibility of
!    future evolution toward 2D decomposition. The arguments are similar to
!    the 1D case, except that now we have optional arguments
!    <TT>flags</TT>, <TT>halo</TT>, <TT>extent</TT> and <TT>maskmap</TT>
!    along two axes.
!    
!    <TT>flags</TT> can now take an additional possible value to fold
!    one or more edges. This is done by using flags
!    <TT>FOLD_WEST_EDGE</TT>, <TT>FOLD_EAST_EDGE</TT>,
!    <TT>FOLD_SOUTH_EDGE</TT> or <TT>FOLD_NORTH_EDGE</TT>. When a fold
!    exists (e.g cylindrical domain), vector fields reverse sign upon
!    crossing the fold. This parity reversal is performed only in the
!    vector version of <LINK
!    SRC="#mpp_update_domains"><TT>mpp_update_domains</TT></LINK>. In
!    addition, shift operations may need to be applied to vector fields on
!    staggered grids, also described in the vector interface to
!    <TT>mpp_update_domains</TT>.
!    
!    <TT>name</TT> is the name associated with the decomposition,
!    e.g <TT>'Ocean model'</TT>. If this argument is present,
!    <TT>mpp_define_domains</TT> will print the domain decomposition
!    generated to <TT>stdlog</TT>.
!    
!    Examples:
!    
!    <PRE>
!    call mpp_define_domains( (/1,100,1,100/), (/2,2/), domain, xhalo=1 )
!    </PRE>
!    
!    will create the following domain layout:
!    <PRE>
!                   |---------|-----------|-----------|-------------|
!                   |domain(1)|domain(2)  |domain(3)  |domain(4)    |
!    |--------------|---------|-----------|-----------|-------------|
!    |Compute domain|1,50,1,50|51,100,1,50|1,50,51,100|51,100,51,100|
!    |--------------|---------|-----------|-----------|-------------|
!    |Data domain   |0,51,1,50|50,101,1,50|0,51,51,100|50,101,51,100|
!    |--------------|---------|-----------|-----------|-------------|
!    </PRE>
!    
!    Again, we allocate arrays on the data domain, perform computations
!    on the compute domain, and call <TT>mpp_update_domains</TT> to update
!    the halo region.
!    
!    If we wished to perfom a 1D decomposition along <TT>Y</TT>
!    on the same global domain, we could use:
    
!    <PRE>
!    call mpp_define_domains( (/1,100,1,100/), layout=(/4,1/), domain, xhalo=1 )
!    </PRE>
    
!    This will create the following domain layout:
!    <PRE>
!                   |----------|-----------|-----------|------------|
!                   |domain(1) |domain(2)  |domain(3)  |domain(4)   |
!    |--------------|----------|-----------|-----------|------------|
!    |Compute domain|1,100,1,25|1,100,26,50|1,100,51,75|1,100,76,100|
!    |--------------|----------|-----------|-----------|------------|
!    |Data domain   |0,101,1,25|0,101,26,50|0,101,51,75|1,101,76,100|
!    |--------------|----------|-----------|-----------|------------|
!    </PRE>
!   </NOTE>
! </INTERFACE>
  interface mpp_define_domains
     module procedure mpp_define_domains1D
     module procedure mpp_define_domains2D
  end interface

! <INTERFACE NAME="mpp_update_domains">
!  <OVERVIEW>
!     Halo updates.
!  </OVERVIEW>
!  <DESCRIPTION>
!    <TT>mpp_update_domains</TT> is used to perform a halo update of a
!    domain-decomposed array on each PE. <TT>MPP_TYPE_</TT> can be of type
!    <TT>complex</TT>, <TT>integer</TT>, <TT>logical</TT> or <TT>real</TT>;
!    of 4-byte or 8-byte kind; of rank up to 5. The vector version (with
!    two input data fields) is only present for <TT>real</TT> types.
!    
!    For 2D domain updates, if there are halos present along both
!    <TT>x</TT> and <TT>y</TT>, we can choose to update one only, by
!    specifying <TT>flags=XUPDATE</TT> or <TT>flags=YUPDATE</TT>. In
!    addition, one-sided updates can be performed by setting <TT>flags</TT>
!    to any combination of <TT>WUPDATE</TT>, <TT>EUPDATE</TT>,
!    <TT>SUPDATE</TT> and <TT>NUPDATE</TT>, to update the west, east, north
!    and south halos respectively. Any combination of halos may be used by
!    adding the requisite flags, e.g: <TT>flags=XUPDATE+SUPDATE</TT> or
!    <TT>flags=EUPDATE+WUPDATE+SUPDATE</TT> will update the east, west and
!    south halos.
!    
!    If a call to <TT>mpp_update_domains</TT> involves at least one E-W
!    halo and one N-S halo, the corners involved will also be updated, i.e,
!    in the example above, the SE and SW corners will be updated.
!    
!    If <TT>flags</TT> is not supplied, that is
!    equivalent to <TT>flags=XUPDATE+YUPDATE</TT>.
!    
!    The vector version is passed the <TT>x</TT> and <TT>y</TT>
!    components of a vector field in tandem, and both are updated upon
!    return. They are passed together to treat parity issues on various
!    grids. For example, on a cubic sphere projection, the <TT>x</TT> and
!    <TT>y</TT> components may be interchanged when passing from an
!    equatorial cube face to a polar face. For grids with folds, vector
!    components change sign on crossing the fold.
!    
!    Special treatment at boundaries such as folds is also required for
!    staggered grids. The following types of staggered grids are
!    recognized:
!    
!    1) <TT>AGRID</TT>: values are at grid centers.<BR/>
!    2) <TT>BGRID_NE</TT>: vector fields are at the NE vertex of a grid
!    cell, i.e: the array elements <TT>u(i,j)</TT> and <TT>v(i,j)</TT> are
!    actually at (i+&#189;,j+&#189;) with respect to the grid centers.<BR/>
!    3) <TT>BGRID_SW</TT>: vector fields are at the SW vertex of a grid
!    cell, i.e: the array elements <TT>u(i,j)</TT> and <TT>v(i,j)</TT> are
!    actually at (i-&#189;,j-&#189;) with respect to the grid centers.<BR/>
!    4) <TT>CGRID_NE</TT>: vector fields are at the N and E faces of a
!    grid cell, i.e: the array elements <TT>u(i,j)</TT> and <TT>v(i,j)</TT>
!    are actually at (i+&#189;,j) and (i,j+&#189;) with respect to the
!    grid centers.<BR/>
!    5) <TT>CGRID_SW</TT>: vector fields are at the S and W faces of a
!    grid cell, i.e: the array elements <TT>u(i,j)</TT> and <TT>v(i,j)</TT>
!    are actually at (i-&#189;,j) and (i,j-&#189;) with respect to the
!    grid centers.
!
!    The gridtypes listed above are all available by use association as
!    integer parameters. The scalar version of <TT>mpp_update_domains</TT>
!    assumes that the values of a scalar field are always at <TT>AGRID</TT>
!    locations, and no special boundary treatment is required. If vector
!    fields are at staggered locations, the optional argument
!    <TT>gridtype</TT> must be appropriately set for correct treatment at
!    boundaries.
!    
!    It is safe to apply vector field updates to the appropriate arrays
!    irrespective of the domain topology: if the topology requires no
!    special treatment of vector fields, specifying <TT>gridtype</TT> will
!    do no harm.
!
!    <TT>mpp_update_domains</TT> internally buffers the date being sent
!    and received into single messages for efficiency. A turnable internal
!    buffer area in memory is provided for this purpose by
!    <TT>mpp_domains_mod</TT>. The size of this buffer area can be set by
!    the user by calling <LINK SRC="mpp_domains.html#mpp_domains_set_stack_size">
!    <TT>mpp_domains_set_stack_size</TT></LINK>.
!  </DESCRIPTION>
!  <TEMPLATE>
!    call mpp_update_domains( field, domain, flags )
!  </TEMPLATE>
!  <TEMPLATE>
!    call mpp_update_domains( fieldx, fieldy, domain, flags, gridtype )
!  </TEMPLATE>
! </INTERFACE>
  interface mpp_update_domains
     module procedure mpp_update_domain2D_r8_2d
     module procedure mpp_update_domain2D_r8_3d
     module procedure mpp_update_domain2D_r8_4d
     module procedure mpp_update_domain2D_r8_5d
     module procedure mpp_update_domain2D_r8_2dv
     module procedure mpp_update_domain2D_r8_3dv
     module procedure mpp_update_domain2D_r8_4dv
     module procedure mpp_update_domain2D_r8_5dv
     module procedure mpp_update_domain2D_c8_2d
     module procedure mpp_update_domain2D_c8_3d
     module procedure mpp_update_domain2D_c8_4d
     module procedure mpp_update_domain2D_c8_5d
#ifndef no_8byte_integers
     module procedure mpp_update_domain2D_i8_2d
     module procedure mpp_update_domain2D_i8_3d
     module procedure mpp_update_domain2D_i8_4d
     module procedure mpp_update_domain2D_i8_5d
     module procedure mpp_update_domain2D_l8_2d
     module procedure mpp_update_domain2D_l8_3d
     module procedure mpp_update_domain2D_l8_4d
     module procedure mpp_update_domain2D_l8_5d
#endif
#ifndef no_4byte_reals
     module procedure mpp_update_domain2D_r4_2d
     module procedure mpp_update_domain2D_r4_3d
     module procedure mpp_update_domain2D_r4_4d
     module procedure mpp_update_domain2D_r4_5d
     module procedure mpp_update_domain2D_c4_2d
     module procedure mpp_update_domain2D_c4_3d
     module procedure mpp_update_domain2D_c4_4d
     module procedure mpp_update_domain2D_c4_5d
     module procedure mpp_update_domain2D_r4_2dv
     module procedure mpp_update_domain2D_r4_3dv
     module procedure mpp_update_domain2D_r4_4dv
     module procedure mpp_update_domain2D_r4_5dv
#endif
     module procedure mpp_update_domain2D_i4_2d
     module procedure mpp_update_domain2D_i4_3d
     module procedure mpp_update_domain2D_i4_4d
     module procedure mpp_update_domain2D_i4_5d
     module procedure mpp_update_domain2D_l4_2d
     module procedure mpp_update_domain2D_l4_3d
     module procedure mpp_update_domain2D_l4_4d
     module procedure mpp_update_domain2D_l4_5d

!     module procedure mpp_update_domain1D_r8_2d
!     module procedure mpp_update_domain1D_r8_3d
!     module procedure mpp_update_domain1D_r8_4d
!     module procedure mpp_update_domain1D_r8_5d
!     module procedure mpp_update_domain1D_c8_2d
!     module procedure mpp_update_domain1D_c8_3d
!     module procedure mpp_update_domain1D_c8_4d
!     module procedure mpp_update_domain1D_c8_5d
!#ifndef no_8byte_integers
!     module procedure mpp_update_domain1D_i8_2d
!     module procedure mpp_update_domain1D_i8_3d
!     module procedure mpp_update_domain1D_i8_4d
!     module procedure mpp_update_domain1D_i8_5d
!     module procedure mpp_update_domain1D_l8_2d
!     module procedure mpp_update_domain1D_l8_3d
!     module procedure mpp_update_domain1D_l8_4d
!     module procedure mpp_update_domain1D_l8_5d
!#endif
!     module procedure mpp_update_domain1D_r4_2d
!     module procedure mpp_update_domain1D_r4_3d
!     module procedure mpp_update_domain1D_r4_4d
!     module procedure mpp_update_domain1D_r4_5d
!     module procedure mpp_update_domain1D_c4_2d
!     module procedure mpp_update_domain1D_c4_3d
!     module procedure mpp_update_domain1D_c4_4d
!     module procedure mpp_update_domain1D_c4_5d
!     module procedure mpp_update_domain1D_i4_2d
!     module procedure mpp_update_domain1D_i4_3d
!     module procedure mpp_update_domain1D_i4_4d
!     module procedure mpp_update_domain1D_i4_5d
!     module procedure mpp_update_domain1D_l4_2d
!     module procedure mpp_update_domain1D_l4_3d
!     module procedure mpp_update_domain1D_l4_4d
!     module procedure mpp_update_domain1D_l4_5d
  end interface

! <INTERFACE NAME="mpp_redistribute">
!  <OVERVIEW>
!    Reorganization of distributed global arrays.
!  </OVERVIEW>
!  <DESCRIPTION>
!    <TT>mpp_redistribute</TT> is used to reorganize a distributed
!    array.  <TT>MPP_TYPE_</TT> can be of type <TT>integer</TT>,
!    <TT>complex</TT>, or <TT>real</TT>; of 4-byte or 8-byte kind; of rank
!    up to 5.
!  </DESCRIPTION>
!  <TEMPLATE>
!    call mpp_redistribute( domain_in, field_in, domain_out, field_out )
!  </TEMPLATE>
!  <IN NAME="field_in" TYPE="MPP_TYPE_">
!    <TT>field_in</TT> is dimensioned on the data domain of <TT>domain_in</TT>.
!  </IN>
!  <OUT NAME="field_out" TYPE="MPP_TYPE_">
!    <TT>field_out</TT> on the data domain of <TT>domain_out</TT>.
!  </OUT>
! </INTERFACE>
  interface mpp_redistribute
     module procedure mpp_redistribute_r8_2D
     module procedure mpp_redistribute_r8_3D
     module procedure mpp_redistribute_r8_4D
     module procedure mpp_redistribute_r8_5D
     module procedure mpp_redistribute_c8_2D
     module procedure mpp_redistribute_c8_3D
     module procedure mpp_redistribute_c8_4D
     module procedure mpp_redistribute_c8_5D
#ifndef no_8byte_integers
     module procedure mpp_redistribute_i8_2D
     module procedure mpp_redistribute_i8_3D
     module procedure mpp_redistribute_i8_4D
     module procedure mpp_redistribute_i8_5D
     module procedure mpp_redistribute_l8_2D
     module procedure mpp_redistribute_l8_3D
     module procedure mpp_redistribute_l8_4D
     module procedure mpp_redistribute_l8_5D
#endif
#ifndef no_4byte_reals
     module procedure mpp_redistribute_r4_2D
     module procedure mpp_redistribute_r4_3D
     module procedure mpp_redistribute_r4_4D
     module procedure mpp_redistribute_r4_5D
     module procedure mpp_redistribute_c4_2D
     module procedure mpp_redistribute_c4_3D
     module procedure mpp_redistribute_c4_4D
     module procedure mpp_redistribute_c4_5D
#endif
     module procedure mpp_redistribute_i4_2D
     module procedure mpp_redistribute_i4_3D
     module procedure mpp_redistribute_i4_4D
     module procedure mpp_redistribute_i4_5D
     module procedure mpp_redistribute_l4_2D
     module procedure mpp_redistribute_l4_3D
     module procedure mpp_redistribute_l4_4D
     module procedure mpp_redistribute_l4_5D
  end interface

! <INTERFACE NAME="mpp_global_field">
!  <OVERVIEW>
!    Fill in a global array from domain-decomposed arrays.
!  </OVERVIEW>
!  <DESCRIPTION>
!    <TT>mpp_global_field</TT> is used to get an entire
!    domain-decomposed array on each PE. <TT>MPP_TYPE_</TT> can be of type
!    <TT>complex</TT>, <TT>integer</TT>, <TT>logical</TT> or <TT>real</TT>;
!    of 4-byte or 8-byte kind; of rank up to 5.
!    
!    All PEs in a domain decomposition must call
!    <TT>mpp_global_field</TT>, and each will have a complete global field
!    at the end. Please note that a global array of rank 3 or higher could
!    occupy a lot of memory.
!  </DESCRIPTION>
!  <TEMPLATE>
!    call mpp_global_field( domain, local, global, flags )
!  </TEMPLATE>
!  <IN NAME="domain" TYPE="type(domain2D)"></IN>
!  <IN NAME="local" TYPE="MPP_TYPE_">
!    <TT>local</TT> is dimensioned on either the compute domain or the
!    data domain of <TT>domain</TT>.
!  </IN>
!  <OUT NAME="global" TYPE="MPP_TYPE_">
!    <TT>global</TT> is dimensioned on the corresponding global domain.
!  </OUT>
!  <IN NAME="flags" TYPE="integer">
!    <TT>flags</TT> can be given the value <TT>XONLY</TT> or
!    <TT>YONLY</TT>, to specify a globalization on one axis only.
!  </IN>
! </INTERFACE>
  interface mpp_global_field
     module procedure mpp_global_field2D_r8_2d
     module procedure mpp_global_field2D_r8_3d
     module procedure mpp_global_field2D_r8_4d
     module procedure mpp_global_field2D_r8_5d
     module procedure mpp_global_field2D_c8_2d
     module procedure mpp_global_field2D_c8_3d
     module procedure mpp_global_field2D_c8_4d
     module procedure mpp_global_field2D_c8_5d
#ifndef no_8byte_integers
     module procedure mpp_global_field2D_i8_2d
     module procedure mpp_global_field2D_i8_3d
     module procedure mpp_global_field2D_i8_4d
     module procedure mpp_global_field2D_i8_5d
     module procedure mpp_global_field2D_l8_2d
     module procedure mpp_global_field2D_l8_3d
     module procedure mpp_global_field2D_l8_4d
     module procedure mpp_global_field2D_l8_5d
#endif
#ifndef no_4byte_reals
     module procedure mpp_global_field2D_r4_2d
     module procedure mpp_global_field2D_r4_3d
     module procedure mpp_global_field2D_r4_4d
     module procedure mpp_global_field2D_r4_5d
     module procedure mpp_global_field2D_c4_2d
     module procedure mpp_global_field2D_c4_3d
     module procedure mpp_global_field2D_c4_4d
     module procedure mpp_global_field2D_c4_5d
#endif
     module procedure mpp_global_field2D_i4_2d
     module procedure mpp_global_field2D_i4_3d
     module procedure mpp_global_field2D_i4_4d
     module procedure mpp_global_field2D_i4_5d
     module procedure mpp_global_field2D_l4_2d
     module procedure mpp_global_field2D_l4_3d
     module procedure mpp_global_field2D_l4_4d
     module procedure mpp_global_field2D_l4_5d

     module procedure mpp_global_field1D_r8_2d
     module procedure mpp_global_field1D_c8_2d
#ifndef no_8byte_integers
     module procedure mpp_global_field1D_i8_2d
     module procedure mpp_global_field1D_l8_2d
#endif
#ifndef no_4byte_reals
     module procedure mpp_global_field1D_r4_2d
     module procedure mpp_global_field1D_c4_2d
#endif
     module procedure mpp_global_field1D_i4_2d
     module procedure mpp_global_field1D_l4_2d
  end interface

! <INTERFACE NAME="mpp_global_max">
!  <OVERVIEW>
!    Global max/min of domain-decomposed arrays.
!  </OVERVIEW>
!  <DESCRIPTION>
!    <TT>mpp_global_max</TT> is used to get the maximum value of a
!    domain-decomposed array on each PE. <TT>MPP_TYPE_</TT> can be of type
!    <TT>integer</TT> or <TT>real</TT>; of 4-byte or 8-byte kind; of rank
!    up to 5. The dimension of <TT>locus</TT> must equal the rank of
!    <TT>field</TT>.
!    
!    All PEs in a domain decomposition must call
!    <TT>mpp_global_max</TT>, and each will have the result upon exit.
!    
!    The function <TT>mpp_global_min</TT>, with an identical syntax. is
!    also available.
!  </DESCRIPTION>
!  <TEMPLATE>
!    mpp_global_max( domain, field, locus )
!  </TEMPLATE>
!  <IN NAME="domain" TYPE="type(domain2D)"></IN>
!  <IN NAME="field" TYPE="MPP_TYPE_">  
!    <TT>field</TT> is dimensioned on either the compute domain or the
!    data domain of <TT>domain</TT>.
!  </IN>
!  <OUT NAME="locus" TYPE="integer" DIM="(:)">
!    <TT>locus</TT>, if present, can be used to retrieve the location of
!    the maximum (as in the <TT>MAXLOC</TT> intrinsic of f90).
!  </OUT>
! </INTERFACE>

  interface mpp_global_max
     module procedure mpp_global_max_r8_2d
     module procedure mpp_global_max_r8_3d
     module procedure mpp_global_max_r8_4d
     module procedure mpp_global_max_r8_5d
#ifndef no_4byte_reals
     module procedure mpp_global_max_r4_2d
     module procedure mpp_global_max_r4_3d
     module procedure mpp_global_max_r4_4d
     module procedure mpp_global_max_r4_5d
#endif
#ifndef no_8byte_integers
     module procedure mpp_global_max_i8_2d
     module procedure mpp_global_max_i8_3d
     module procedure mpp_global_max_i8_4d
     module procedure mpp_global_max_i8_5d
#endif
     module procedure mpp_global_max_i4_2d
     module procedure mpp_global_max_i4_3d
     module procedure mpp_global_max_i4_4d
     module procedure mpp_global_max_i4_5d
  end interface

  interface mpp_global_min
     module procedure mpp_global_min_r8_2d
     module procedure mpp_global_min_r8_3d
     module procedure mpp_global_min_r8_4d
     module procedure mpp_global_min_r8_5d
#ifndef no_4byte_reals
     module procedure mpp_global_min_r4_2d
     module procedure mpp_global_min_r4_3d
     module procedure mpp_global_min_r4_4d
     module procedure mpp_global_min_r4_5d
#endif
#ifndef no_8byte_integers
     module procedure mpp_global_min_i8_2d
     module procedure mpp_global_min_i8_3d
     module procedure mpp_global_min_i8_4d
     module procedure mpp_global_min_i8_5d
#endif
     module procedure mpp_global_min_i4_2d
     module procedure mpp_global_min_i4_3d
     module procedure mpp_global_min_i4_4d
     module procedure mpp_global_min_i4_5d
  end interface

! <INTERFACE NAME="mpp_global_sum">
!  <OVERVIEW>
!    Global sum of domain-decomposed arrays.
!  </OVERVIEW>
!  <DESCRIPTION>
!    <TT>mpp_global_sum</TT> is used to get the sum of a
!    domain-decomposed array on each PE. <TT>MPP_TYPE_</TT> can be of type
!    <TT>integer</TT>, <TT>complex</TT>, or <TT>real</TT>; of 4-byte or
!    8-byte kind; of rank up to 5.
!  </DESCRIPTION>
!  <TEMPLATE>
!    call mpp_global_sum( domain, field, flags )
!  </TEMPLATE>
!  <IN NAME="domain" TYPE="type(domain2D)"></IN>
!  <IN NAME="field" TYPE="MPP_TYPE_">
!    <TT>field</TT> is dimensioned on either the compute domain or the
!    data domain of <TT>domain</TT>.
!  </IN>
!  <IN NAME="flags" TYPE="integer">
!    <TT>flags</TT>, if present, must have the value
!    <TT>BITWISE_EXACT_SUM</TT>. This produces a sum that is guaranteed to
!    produce the identical result irrespective of how the domain is
!    decomposed. This method does the sum first along the ranks beyond 2,
!    and then calls <LINK
!    SRC="#mpp_global_field"><TT>mpp_global_field</TT></LINK> to produce a
!    global 2D array which is then summed. The default method, which is
!    considerably faster, does a local sum followed by <LINK
!    SRC="mpp.html#mpp_sum"><TT>mpp_sum</TT></LINK> across the domain
!    decomposition.
!  </IN>
!  <NOTE>
!    All PEs in a domain decomposition must call
!    <TT>mpp_global_sum</TT>, and each will have the result upon exit.
!  </NOTE>
! </INTERFACE>

  interface mpp_global_sum
     module procedure mpp_global_sum_r8_2d
     module procedure mpp_global_sum_r8_3d
     module procedure mpp_global_sum_r8_4d
     module procedure mpp_global_sum_r8_5d
     module procedure mpp_global_sum_c8_2d
     module procedure mpp_global_sum_c8_3d
     module procedure mpp_global_sum_c8_4d
     module procedure mpp_global_sum_c8_5d
#ifndef no_4byte_reals
     module procedure mpp_global_sum_r4_2d
     module procedure mpp_global_sum_r4_3d
     module procedure mpp_global_sum_r4_4d
     module procedure mpp_global_sum_r4_5d
     module procedure mpp_global_sum_c4_2d
     module procedure mpp_global_sum_c4_3d
     module procedure mpp_global_sum_c4_4d
     module procedure mpp_global_sum_c4_5d
#endif
#ifndef no_8byte_integers
     module procedure mpp_global_sum_i8_2d
     module procedure mpp_global_sum_i8_3d
     module procedure mpp_global_sum_i8_4d
     module procedure mpp_global_sum_i8_5d
#endif
     module procedure mpp_global_sum_i4_2d
     module procedure mpp_global_sum_i4_3d
     module procedure mpp_global_sum_i4_4d
     module procedure mpp_global_sum_i4_5d
  end interface

! <INTERFACE NAME="operator">
!  <OVERVIEW>
!    Equality/inequality operators for domaintypes.
!  </OVERVIEW>
!  <DESCRIPTION>
!    The module provides public operators to check for
!    equality/inequality of domaintypes, e.g:
!    
!    <PRE>
!    type(domain1D) :: a, b
!    type(domain2D) :: c, d
!    ...
!    if( a.NE.b )then
!        ...
!    end if
!    if( c==d )then
!        ...
!    end if
!    </PRE>
!    
!    Domains are considered equal if and only if the start and end
!    indices of each of their component global, data and compute domains
!    are equal.
!  </DESCRIPTION>
! </INTERFACE>
  interface operator(.EQ.)
     module procedure mpp_domain1D_eq
     module procedure mpp_domain2D_eq
  end interface

  interface operator(.NE.)
     module procedure mpp_domain1D_ne
     module procedure mpp_domain2D_ne
  end interface

! <INTERFACE NAME="mpp_get_compute_domain">
!  <OVERVIEW>
!    These routines retrieve the axis specifications associated with the compute domains.
!  </OVERVIEW>
!  <DESCRIPTION>
!    The domain is a derived type with private elements. These routines 
!    retrieve the axis specifications associated with the compute domains
!    The 2D version of these is a simple extension of 1D.
!  </DESCRIPTION>
!  <TEMPLATE>
!    call mpp_get_compute_domain
!  </TEMPLATE>
! </INTERFACE>
  interface mpp_get_compute_domain
     module procedure mpp_get_compute_domain1D
     module procedure mpp_get_compute_domain2D
  end interface

! <INTERFACE NAME="mpp_get_compute_domains">
!  <OVERVIEW>
!    Retrieve the entire array of compute domain extents associated with a decomposition.
!  </OVERVIEW>
!  <DESCRIPTION>
!    Retrieve the entire array of compute domain extents associated with a decomposition.
!  </DESCRIPTION>
!  <TEMPLATE>
!    call mpp_get_compute_domains( domain, xbegin, xend, xsize, &
!                                                ybegin, yend, ysize )
!  </TEMPLATE>
!  <IN NAME="domain" TYPE="type(domain2D)"></IN>
!  <OUT NAME="xbegin,ybegin" TYPE="integer" DIM="(:)"></OUT>
!  <OUT NAME="xend,yend" TYPE="integer" DIM="(:)"></OUT>
!  <OUT NAME="xsize,ysize" TYPE="integer" DIM="(:)"></OUT>
! </INTERFACE>
  interface mpp_get_compute_domains
     module procedure mpp_get_compute_domains1D
     module procedure mpp_get_compute_domains2D
  end interface

! <INTERFACE NAME="mpp_get_data_domain">
!  <OVERVIEW>
!    These routines retrieve the axis specifications associated with the data domains.
!  </OVERVIEW>
!  <DESCRIPTION>
!    The domain is a derived type with private elements. These routines 
!    retrieve the axis specifications associated with the data domains.
!    The 2D version of these is a simple extension of 1D.
!  </DESCRIPTION>
!  <TEMPLATE>
!    call mpp_get_data_domain
!  </TEMPLATE>
! </INTERFACE>
  interface mpp_get_data_domain
     module procedure mpp_get_data_domain1D
     module procedure mpp_get_data_domain2D
  end interface

! <INTERFACE NAME="mpp_get_global_domain">
!  <OVERVIEW>
!    These routines retrieve the axis specifications associated with the global domains.
!  </OVERVIEW>
!  <DESCRIPTION>
!    The domain is a derived type with private elements. These routines 
!    retrieve the axis specifications associated with the global domains.
!    The 2D version of these is a simple extension of 1D.
!  </DESCRIPTION>
!  <TEMPLATE>
!    call mpp_get_global_domain
!  </TEMPLATE>
! </INTERFACE>
  interface mpp_get_global_domain
     module procedure mpp_get_global_domain1D
     module procedure mpp_get_global_domain2D
  end interface

! <INTERFACE NAME="mpp_define_layout">
!  <OVERVIEW>
!    Retrieve layout associated with a domain decomposition.
!  </OVERVIEW>
!  <DESCRIPTION>
!    Given a global 2D domain and the number of divisions in the
!    decomposition (<TT>ndivs</TT>: usually the PE count unless some
!    domains are masked) this calls returns a 2D domain layout.
!    
!    By default, <TT>mpp_define_layout</TT> will attempt to divide the
!    2D index space into domains that maintain the aspect ratio of the
!    global domain. If this cannot be done, the algorithm favours domains
!    that are longer in <TT>x</TT> than <TT>y</TT>, a preference that could
!    improve vector performance.
!  </DESCRIPTION>
!  <TEMPLATE>
!    call mpp_define_layout( global_indices, ndivs, layout )
!  </TEMPLATE>
!  <IN NAME="global_indices"></IN>
!  <IN NAME="ndivs"></IN>
!  <OUT NAME="layout"></OUT>
! </INTERFACE>

  interface mpp_define_layout
     module procedure mpp_define_layout2D
  end interface

! <INTERFACE NAME="mpp_get_pelist">
!  <OVERVIEW>
!    Retrieve list of PEs associated with a domain decomposition.
!  </OVERVIEW>
!  <DESCRIPTION>
!    The 1D version of this call returns an array of the PEs assigned to this 1D domain
!    decomposition. In addition the optional argument <TT>pos</TT> may be
!    used to retrieve the 0-based position of the domain local to the
!    calling PE, i.e <TT>domain%list(pos)%pe</TT> is the local PE,
!    as returned by <LINK SRC="mpp.html#mpp_pe"><TT>mpp_pe()</TT></LINK>.
!    The 2D version of this call is identical to 1D version.
!  </DESCRIPTION>
!  <IN NAME="domain"></IN>
!  <OUT NAME="pelist"></OUT>
!  <OUT NAME="pos"></OUT>
! </INTERFACE>
  interface mpp_get_pelist
     module procedure mpp_get_pelist1D
     module procedure mpp_get_pelist2D
  end interface

! <INTERFACE NAME="mpp_get_layout">
!  <OVERVIEW>
!    Retrieve layout associated with a domain decomposition.
!  </OVERVIEW>
!  <DESCRIPTION>
!    The 1D version of this call returns the number of divisions that was assigned to this
!    decomposition axis. The 2D version of this call returns an array of
!    dimension 2 holding the results on two axes.
!  </DESCRIPTION>
!  <TEMPLATE>
!    call mpp_get_layout( domain, layout )
!  </TEMPLATE>
!  <IN NAME="domain"></IN>
!  <OUT NAME="layout"></OUT>
! </INTERFACE>
  interface mpp_get_layout
     module procedure mpp_get_layout1D
     module procedure mpp_get_layout2D
  end interface

  public :: mpp_broadcast_domain, mpp_define_layout, mpp_define_domains, mpp_domains_init, mpp_domains_set_stack_size, &
            mpp_domains_exit, mpp_get_compute_domain, mpp_get_compute_domains, mpp_get_data_domain, mpp_get_global_domain, &
            mpp_get_domain_components, mpp_get_layout, mpp_get_pelist, mpp_redistribute, mpp_update_domains, &
            mpp_global_field, mpp_global_max, mpp_global_min, mpp_global_sum, operator(.EQ.), operator(.NE.), &
            mpp_copy_domains

  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                 MPP_DOMAINS: initialization and termination                 !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! <SUBROUTINE NAME="mpp_domains_init">
!  <OVERVIEW>
!    Initialize domain decomp package.
!  </OVERVIEW>
!  <DESCRIPTION>
!    Called to initialize the <TT>mpp_domains_mod</TT> package.
!    
!    <TT>flags</TT> can be set to <TT>MPP_VERBOSE</TT> to have
!    <TT>mpp_domains_mod</TT> keep you informed of what it's up
!    to. <TT>MPP_DEBUG</TT> returns even more information for debugging.
!    
!    <TT>mpp_domains_init</TT> will call <TT>mpp_init</TT>, to make sure
!    <LINK SRC="mpp.html"><TT>mpp_mod</TT></LINK> is initialized. (Repeated
!    calls to <TT>mpp_init</TT> do no harm, so don't worry if you already
!    called it).
!  </DESCRIPTION>
!  <TEMPLATE>
!    call mpp_domains_init(flags)
!  </TEMPLATE>
!  <IN NAME="flags" TYPE="integer"></IN>
! </SUBROUTINE>
    subroutine mpp_domains_init(flags)
!initialize domain decomp package
      integer, intent(in), optional :: flags
      integer :: l=0

      if( module_is_initialized )return
      call mpp_init(flags)           !this is a no-op if already initialized
      pe = mpp_pe()
      module_is_initialized = .TRUE.
      if( pe.EQ.mpp_root_pe() )then
          write( stdlog(),'(/a)' )'MPP_DOMAINS module '//trim(version) &
               //trim(tagname)
!          write( stdlog(),'(a)' )trim(version_update_domains2D)
      end if

      if( PRESENT(flags) )then
          debug   = flags.EQ.MPP_DEBUG
          verbose = flags.EQ.MPP_VERBOSE .OR. debug
          domain_clocks_on = flags.EQ.MPP_DOMAIN_TIME
      end if

      call mpp_domains_set_stack_size(32768) !default, pretty arbitrary
#ifdef use_shmalloc
      call mpp_malloc( ptr_info, 16, l )
#endif

!NULL_DOMAIN is a domaintype that can be used to initialize to undef
      NULL_DOMAIN1D%global%begin  = -1; NULL_DOMAIN1D%global%end  = -1; NULL_DOMAIN1D%global%size = 0
      NULL_DOMAIN1D%data%begin    = -1; NULL_DOMAIN1D%data%end    = -1; NULL_DOMAIN1D%data%size = 0
      NULL_DOMAIN1D%compute%begin = -1; NULL_DOMAIN1D%compute%end = -1; NULL_DOMAIN1D%compute%size = 0
      NULL_DOMAIN1D%pe = ANY_PE
      NULL_DOMAIN2D%x = NULL_DOMAIN1D
      NULL_DOMAIN2D%y = NULL_DOMAIN1D
      NULL_DOMAIN2D%pe = ANY_PE

      if( domain_clocks_on )then
          pack_clock = mpp_clock_id( 'Halo pack' )
          pack_loop_clock = mpp_clock_id( 'Halo pack loop' )
          send_clock = mpp_clock_id( 'Halo send' )
          recv_clock = mpp_clock_id( 'Halo recv' )
          unpk_clock = mpp_clock_id( 'Halo unpk' )
          wait_clock = mpp_clock_id( 'Halo wait' )
      end if
      return
    end subroutine mpp_domains_init

! <SUBROUTINE NAME="mpp_domains_set_stack_size">
!  <OVERVIEW>
!    Set user stack size.
! </OVERVIEW>
! <DESCRIPTION>
!    This sets the size of an array that is used for internal storage by
!    <TT>mpp_domains</TT>. This array is used, for instance, to buffer the
!    data sent and received in halo updates.
!    
!    This call has implied global synchronization. It should be
!    placed somewhere where all PEs can call it.
!  </DESCRIPTION>
!  <TEMPLATE>
!    call mpp_domains_set_stack_size(n)
!  </TEMPLATE>
!  <IN NAME="n" TYPE="integer"></IN>
! </SUBROUTINE>
    subroutine mpp_domains_set_stack_size(n)
!set the mpp_domains_stack variable to be at least n LONG words long
      integer, intent(in) :: n
      character(len=8) :: text

      if( n.LE.mpp_domains_stack_size )return
#ifdef use_shmalloc
      call mpp_malloc( ptr_stack, n, mpp_domains_stack_size )
#else
      if( allocated(mpp_domains_stack) )deallocate(mpp_domains_stack)
      allocate( mpp_domains_stack(n) )
      mpp_domains_stack_size = n
#endif
      write( text,'(i8)' )n
      if( pe.EQ.mpp_root_pe() )call mpp_error( NOTE, 'MPP_DOMAINS_SET_STACK_SIZE: stack size set to '//text//'.' )

      return
    end subroutine mpp_domains_set_stack_size

! <SUBROUTINE NAME="mpp_domains_exit">
!  <OVERVIEW>
!    Exit <TT>mpp_domains_mod</TT>.
!  </OVERVIEW>
!  <DESCRIPTION>
!    Serves no particular purpose, but is provided should you require to
!    re-initialize <TT>mpp_domains_mod</TT>, for some odd reason.
!  </DESCRIPTION>
!  <TEMPLATE>
!    call mpp_domains_exit()
!  </TEMPLATE>
! </SUBROUTINE>
    subroutine mpp_domains_exit()
!currently does not have much to do, but provides the possibility of re-initialization
      if( .NOT.module_is_initialized )return
      call mpp_max(mpp_domains_stack_hwm)
      if( pe.EQ.mpp_root_pe() )write( stdout(),* )'MPP_DOMAINS_STACK high water mark=', mpp_domains_stack_hwm
      module_is_initialized = .FALSE.
      return
    end subroutine mpp_domains_exit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                MPP_DOMAINS: overloaded operators (==, /=)                   !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function mpp_domain1D_eq( a, b )
      logical :: mpp_domain1D_eq
      type(domain1D), intent(in) :: a, b

      mpp_domain1D_eq = ( a%compute%begin.EQ.b%compute%begin .AND. &
                          a%compute%end  .EQ.b%compute%end   .AND. &
                          a%data%begin   .EQ.b%data%begin    .AND. &
                          a%data%end     .EQ.b%data%end      .AND. & 
                          a%global%begin .EQ.b%global%begin  .AND. &
                          a%global%end   .EQ.b%global%end    )
!compare pelists
!      if( mpp_domain1D_eq )mpp_domain1D_eq = ASSOCIATED(a%list) .AND. ASSOCIATED(b%list)
!      if( mpp_domain1D_eq )mpp_domain1D_eq = size(a%list).EQ.size(b%list)
!      if( mpp_domain1D_eq )mpp_domain1D_eq = ALL(a%list%pe.EQ.b%list%pe)
      
      return
    end function mpp_domain1D_eq

    function mpp_domain1D_ne( a, b )
      logical :: mpp_domain1D_ne
      type(domain1D), intent(in) :: a, b

      mpp_domain1D_ne = .NOT. ( a.EQ.b )
      return
    end function mpp_domain1D_ne

    function mpp_domain2D_eq( a, b )
      logical :: mpp_domain2D_eq
      type(domain2D), intent(in) :: a, b

      mpp_domain2D_eq = a%x.EQ.b%x .AND. a%y.EQ.b%y
      if( mpp_domain2D_eq .AND. ((a%pe.EQ.ANY_PE).OR.(b%pe.EQ.ANY_PE)) )return !NULL_DOMAIN2D
!compare pelists
      if( mpp_domain2D_eq )mpp_domain2D_eq = ASSOCIATED(a%list) .AND. ASSOCIATED(b%list)
      if( mpp_domain2D_eq )mpp_domain2D_eq = size(a%list).EQ.size(b%list)
      if( mpp_domain2D_eq )mpp_domain2D_eq = ALL(a%list%pe.EQ.b%list%pe)
      return
    end function mpp_domain2D_eq

    function mpp_domain2D_ne( a, b )
      logical :: mpp_domain2D_ne
      type(domain2D), intent(in) :: a, b

      mpp_domain2D_ne = .NOT. ( a.EQ.b )
      return
    end function mpp_domain2D_ne


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!              MPP_COPY_DOMAINS: Copy domain and add halo                     !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine mpp_copy_domains1D(domain_in, domain_out, halo)
! ARGUMENTS
! domain_in is the input domain
! domain_out is the returned domain
! halo (optional) defines halo width (currently the same on both sides)

      type(domain1D), intent(in)    :: domain_in
      type(domain1D), intent(out)   :: domain_out
      integer, intent(in), optional :: halo

      integer, dimension(:), allocatable :: pelist, extent
      integer :: ndivs, global_indices(2) !(/ isg, ieg /)
      integer :: isg, ieg, halosz

      
! get the global indices of the input domain
      call mpp_get_global_domain(domain_in, isg, ieg )
      global_indices(1) = isg;       global_indices(2) = ieg

! get the number of divisions of domain
      call mpp_get_layout(domain_in, ndivs)

      allocate(pelist(0:ndivs-1),extent(0:ndivs-1))

! get the pe list of the input domain
      call mpp_get_pelist( domain_in, pelist)

! get the halo
      halosz = 0
      if(present(halo)) halosz = halo

! get the extent
      call mpp_get_compute_domains(Domain_in, size = extent(0:ndivs-1))      

      call mpp_define_domains( global_indices, ndivs, domain_out, pelist = pelist, &
                                 halo = halosz, extent = extent )

    end subroutine mpp_copy_domains1D


!----------------------------------------------------------------------------------


    subroutine mpp_copy_domains2D(domain_in, domain_out, xhalo, yhalo)
! ARGUMENTS
! domain_in is the input domain
! domain_out is the returned domain
! xhalo, yhalo (optional) defines halo width. 

      type(domain2D), intent(in)    :: domain_in
      type(domain2D), intent(out)   :: domain_out
      integer, intent(in), optional :: xhalo, yhalo

      integer, dimension(:), allocatable :: pelist, xextent, yextent, xbegin, xend, &
                                            xsize, ybegin, yend, ysize
      integer :: ndivx, ndivy, ndivs, npes, global_indices(4), layout(2)
      integer :: isg, ieg, jsg, jeg, xhalosz, yhalosz, i, j, m, n, pe
      
! get the global indices of the input domain
      call mpp_get_global_domain(domain_in, isg, ieg, jsg, jeg )
      global_indices(1) = isg;       global_indices(2) = ieg
      global_indices(3) = jsg;       global_indices(4) = jeg

! get the number of divisions of domain
      call mpp_get_layout(domain_in, layout)
      ndivx = layout(1); ndivy = layout(2)
      ndivs = ndivx * ndivy
      npes = mpp_npes()

      if( ndivs.NE.npes )call mpp_error( 'mpp_domains_mod', &
                 'mpp_copy_domains: number of PEs is not consistent with the layout', FATAL )
      allocate(pelist(0:npes-1), xsize(0:npes-1), ysize(0:npes-1), &
               xbegin(0:npes-1), xend(0:npes-1),                   &
               ybegin(0:npes-1), yend(0:npes-1),                   &
               xextent(0:ndivx), yextent(0:ndivy))

! get the pe list of the input domain
      call mpp_get_pelist2D( domain_in, pelist)

! get the halo
      xhalosz = 0; yhalosz = 0
      if(present(xhalo)) xhalosz = xhalo
      if(present(yhalo)) yhalosz = yhalo

! get the extent
      call mpp_get_compute_domains(Domain_in, xbegin, xend, xsize, &
                                          ybegin, yend, ysize  )
            
!  compute the exact x-axis decomposition
      i = isg
      m = 0
      xextent = 0
      do pe = 0, npes-1
        if ( xbegin(pe) == i ) then
          m = m+1
          xextent (m) = xsize (pe)
          i = i + xsize(pe)
          if ( m == size(xextent) .or. i > ieg ) exit
        endif
      enddo

!  compute the exact y-axis decomposition
      j = jsg 
      n = 0   
      yextent = 0
      do pe = 0, npes-1
        if ( ybegin(pe) == j ) then
          n = n+1 
          yextent (n) = ysize (pe)
          j = j + ysize(pe)
          if ( n == size(yextent) .or. j > jeg ) exit
        endif   
      enddo   

      call mpp_define_domains( global_indices, layout, domain_out, pelist = pelist, &
                               xhalo = xhalosz, yhalo = yhalosz, xextent = xextent, & 
                               yextent = yextent)

    end subroutine mpp_copy_domains2D



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!              MPP_DEFINE_DOMAINS: define layout and decomposition            !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! <SUBROUTINE NAME="mpp_define_domains1D" INTERFACE="mpp_define_domains">>
!   <IN NAME="global_indices" TYPE="integer" DIM="(2)"> </IN>
!   <IN NAME="ndivs" TYPE="integer">  </IN>
!   <INOUT NAME="domain" TYPE="type(domain1D)"> </INOUT>
!   <IN NAME="pelist" TYPE="integer" DIM="(0:)">  </IN>
!   <IN NAME="flags" TYPE="integer">  </IN>
!   <IN NAME="halo" TYPE="integer">  </IN>
!   <IN NAME="extent" TYPE="integer" DIM="(0:)">  </IN>
!   <IN NAME="maskmap" TYPE="logical" DIM="(0:)"> </IN>
! </SUBROUTINE>
    subroutine mpp_define_domains1D( global_indices, ndivs, domain, pelist, flags, halo, extent, maskmap )
!routine to divide global array indices among domains, and assign domains to PEs
!domain is of type domain1D
!ARGUMENTS:
!      global_indices(2)=(isg,ieg) gives the extent of global domain
!      ndivs is number of divisions of domain: even divisions unless extent is present.
!      domain is the returned domain1D
!      pelist (optional) list of PEs to which domains are to be assigned (default 0...npes-1)
!                 size of pelist must correspond to number of mask=.TRUE. divisions
!      flags define whether compute and data domains are global (undecomposed) and whether global domain has periodic boundaries
!      halo (optional) defines halo width (currently the same on both sides)
!      extent (optional) array defines width of each division (used for non-uniform domain decomp, for e.g load-balancing)
!      maskmap (optional) a division whose maskmap=.FALSE. is not assigned to any domain
!  By default we assume decomposition of compute and data domains, non-periodic boundaries, no halo, as close to uniform extents
!  as the input parameters permit
      integer, intent(in) :: global_indices(2) !(/ isg, ieg /)
      integer, intent(in) :: ndivs
      type(domain1D), intent(inout) :: domain !declared inout so that existing links, if any, can be nullified
      integer, intent(in), optional :: pelist(0:)
      integer, intent(in), optional :: flags, halo
      integer, intent(in), optional :: extent(0:)
      logical, intent(in), optional :: maskmap(0:)

      logical :: compute_domain_is_global, data_domain_is_global
      integer :: ndiv, n, isg, ieg, is, ie, i, off, pos, hs, he
      integer, allocatable :: pes(:)
      logical, allocatable :: mask(:)
      integer :: halosz
!used by symmetry algorithm
      integer :: imax, ndmax, ndmirror
      logical :: symmetrize
!statement functions
      logical :: even, odd
      even(n) = (mod(n,2).EQ.0)
      odd (n) = (mod(n,2).EQ.1)
      
      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS1D: You must first call mpp_domains_init.' )
!get global indices
      isg = global_indices(1)
      ieg = global_indices(2)
      if( ndivs.GT.ieg-isg+1 )call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS1D: more divisions requested than rows available.' )
!get the list of PEs on which to assign domains; if pelist is absent use 0..npes-1
      if( PRESENT(pelist) )then
          if( .NOT.any(pelist.EQ.pe) )then
              write( stderr(),* )'pe=', pe, ' pelist=', pelist
              call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS1D: pe must be in pelist.' )
          end if
          allocate( pes(0:size(pelist)-1) )
          pes(:) = pelist(:)
      else
          allocate( pes(0:mpp_npes()-1) )
          pes(:) = (/ (i,i=0,mpp_npes()-1) /)
      end if

!get number of real domains: 1 mask domain per PE in pes
      allocate( mask(0:ndivs-1) )
      mask = .TRUE.                 !default mask
      if( PRESENT(maskmap) )then
          if( size(maskmap).NE.ndivs ) &
               call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS1D: maskmap array size must equal number of domain divisions.' )
          mask(:) = maskmap(:)
      end if
      if( count(mask).NE.size(pes) ) &
           call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS1D: number of TRUEs in maskmap array must match PE count.' )
      if( PRESENT(extent) )then
          if( size(extent).NE.ndivs ) &
               call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS1D: extent array size must equal number of domain divisions.' )
      end if
!get halosize
      halosz = 0
      if( PRESENT(halo) )halosz = halo

!get flags
      compute_domain_is_global = .FALSE.
      data_domain_is_global    = .FALSE.
      domain%cyclic = .FALSE.
      if( PRESENT(flags) )then
!NEW: obsolete flag global_compute_domain, since ndivs is non-optional and you cannot have global compute and ndivs.NE.1 
          compute_domain_is_global = ndivs.EQ.1
!if compute domain is global, data domain must also be
          data_domain_is_global    = BTEST(flags,GLOBAL) .OR. compute_domain_is_global
          domain%cyclic  = BTEST(flags,CYCLIC) .AND. halosz.NE.0
      end if

!set up links list
      allocate( domain%list(0:ndivs-1) )

!set global domain
      domain%list(:)%global%begin     = isg
      domain%list(:)%global%end       = ieg
      domain%list(:)%global%size      = ieg-isg+1
      domain%list(:)%global%max_size  = ieg-isg+1
      domain%list(:)%global%is_global = .TRUE. !always

!get compute domain
      if( compute_domain_is_global )then
          domain%list(:)%compute%begin = isg
          domain%list(:)%compute%end   = ieg
          domain%list(:)%compute%is_global = .TRUE.
          domain%list(:)%pe = pes(:)
          domain%pos = 0
      else
          domain%list(:)%compute%is_global = .FALSE.
          is = isg
          n = 0
          do ndiv=0,ndivs-1
             if( PRESENT(extent) )then
                 ie = is + extent(ndiv) - 1
                 if( ndiv.EQ.ndivs-1 .AND. ie.NE.ieg ) &
                      call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS: extent array limits do not match global domain.' )
             else
!modified for mirror-symmetry
!original line
!                 ie = is + CEILING( float(ieg-is+1)/(ndivs-ndiv) ) - 1

!problem of dividing nx points into n domains maintaining symmetry
!i.e nx=18 n=4 4554 and 5445 are solutions but 4455 is not.
!this will always work for nx even n even or odd
!this will always work for nx odd, n odd
!this will never  work for nx odd, n even: for this case we supersede the mirror calculation
!                 symmetrize = .NOT. ( mod(ndivs,2).EQ.0 .AND. mod(ieg-isg+1,2).EQ.1 )
!nx even n odd fails if n>nx/2
                 symmetrize = ( even(ndivs) .AND. even(ieg-isg+1) ) .OR. &
                              (  odd(ndivs) .AND.  odd(ieg-isg+1) ) .OR. &
                              (  odd(ndivs) .AND. even(ieg-isg+1) .AND. ndivs.LT.(ieg-isg+1)/2 )

!mirror domains are stored in the list and retrieved if required.
                 if( ndiv.EQ.0 )then
!initialize max points and max domains
                     imax = ieg
                     ndmax = ndivs
                 end if
!do bottom half of decomposition, going over the midpoint for odd ndivs
                 if( ndiv.LT.(ndivs-1)/2+1 )then
!domain is sized by dividing remaining points by remaining domains
                     ie = is + CEILING( REAL(imax-is+1)/(ndmax-ndiv) ) - 1
                     ndmirror = (ndivs-1) - ndiv !mirror domain
                     if( ndmirror.GT.ndiv .AND. symmetrize )then !only for domains over the midpoint
!mirror extents, the max(,) is to eliminate overlaps
                         domain%list(ndmirror)%compute%begin = max( isg+ieg-ie, ie+1 )
                         domain%list(ndmirror)%compute%end   = max( isg+ieg-is, ie+1 )
                         imax = domain%list(ndmirror)%compute%begin - 1
                         ndmax = ndmax - 1
                     end if
                 else
                     if( symmetrize )then
!do top half of decomposition by retrieving saved values
                         is = domain%list(ndiv)%compute%begin
                         ie = domain%list(ndiv)%compute%end
                     else
                         ie = is + CEILING( REAL(imax-is+1)/(ndmax-ndiv) ) - 1
                     end if
                 end if
             end if
             if( ie.LT.is )call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS: domain extents must be positive definite.' )
             domain%list(ndiv)%compute%begin = is
             domain%list(ndiv)%compute%end   = ie
             if( ndiv.GT.0 .AND. is.NE.domain%list(ndiv-1)%compute%end+1 ) &
                  call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS: domain extents do not span space completely.' )
             if( ndiv.EQ.ndivs-1 .AND. domain%list(ndiv)%compute%end.NE.ieg ) &
                  call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS: domain extents do not span space completely.' )
             if( mask(ndiv) )then
                 domain%list(ndiv)%pe = pes(n)
                 if( pe.EQ.pes(n) )domain%pos = ndiv
                 n = n + 1
             end if
             is = ie + 1
          end do
      end if
          
      domain%list(:)%compute%size  = domain%list(:)%compute%end - domain%list(:)%compute%begin + 1

!get data domain
!data domain is at least equal to compute domain
      domain%list(:)%data%begin = domain%list(:)%compute%begin
      domain%list(:)%data%end   = domain%list(:)%compute%end
      domain%list(:)%data%is_global = .FALSE.
!apply global flags
      if( data_domain_is_global )then
          domain%list(:)%data%begin  = isg
          domain%list(:)%data%end    = ieg
          domain%list(:)%data%is_global = .TRUE.
      end if
!apply margins
      domain%list(:)%data%begin = domain%list(:)%data%begin - halosz
      domain%list(:)%data%end   = domain%list(:)%data%end   + halosz  
      domain%list(:)%data%size  = domain%list(:)%data%end - domain%list(:)%data%begin + 1

!      domain = domain%list(pos) !load domain from domain%list(pos)
      domain%compute = domain%list(domain%pos)%compute
      domain%data = domain%list(domain%pos)%data
      domain%global = domain%list(domain%pos)%global
      domain%compute%max_size = MAXVAL( domain%list(:)%compute%size )
      domain%data%max_size    = MAXVAL( domain%list(:)%data%size )
      domain%global%max_size  = domain%global%size

!PV786667: the deallocate stmts can be removed when fixed (7.3.1.3m)
      deallocate( pes, mask )
      return

      contains

        function if_overlap( hs, he, cs, ce, os, oe )
!function to look if a portion of the halo [hs,he] lies with in the compute region [cs,ce]
!if yes, if_overlap returns true, and the overlap region is returned in [os,oe]
          logical :: if_overlap
          integer, intent(in) :: hs, he, cs, ce
          integer, intent(out) :: os, oe
          os = max(hs,cs)
          oe = min(he,ce)
          if( debug )write( stderr(),'(a,7i4)' ) &
               'MPP_DEFINE_DOMAINS1D: pe, hs, he, cs, ce, os, oe=', pe, hs, he, cs, ce, os, oe
          if_overlap = oe.GE.os
          return
        end function if_overlap
          
    end subroutine mpp_define_domains1D
! <SUBROUTINE NAME="mpp_define_domains2D" INTERFACE="mpp_define_domains">
!  <IN NAME="global_indices" TYPE="integer" DIM="(4)"> </IN>
!  <IN NAME="layout" TYPE="integer" DIM="(2)"></IN>
!  <INOUT NAME="domain" TYPE="type(domain2D)"></INOUT>
!  <IN NAME="pelist" TYPE="integer" DIM="(0:)"></IN>
!  <IN NAME="xflags, yflags" TYPE="integer"></IN>
!  <IN NAME="xhalo, yhalo" TYPE="integer"></IN>
!  <IN NAME="xextent, yextent" TYPE="integer" DIM="(0:)"></IN>
!  <IN NAME="maskmap" TYPE="logical" DIM="(:,:)"></IN>
!  <IN NAME="name" TYPE="character(len=*)"></IN>
! </SUBROUTINE>
    subroutine mpp_define_domains2D( global_indices, layout, domain, pelist, &
                                                     xflags, yflags, xhalo, yhalo, xextent, yextent, maskmap, name )
!define 2D data and computational domain on global rectilinear cartesian domain (isg:ieg,jsg:jeg) and assign them to PEs
      integer, intent(in) :: global_indices(4) !(/ isg, ieg, jsg, jeg /)
      integer, intent(in) :: layout(2)
      type(domain2D), intent(inout) :: domain
      integer, intent(in), optional :: pelist(0:)
      integer, intent(in), optional :: xflags, yflags, xhalo, yhalo
      integer, intent(in), optional :: xextent(0:), yextent(0:)
      logical, intent(in), optional :: maskmap(0:,0:)
      character(len=*), intent(in), optional :: name
      integer :: i, j, m, n
      integer :: ipos, jpos, pos
      integer :: ndivx, ndivy, isg, ieg, jsg, jeg, isd, ied, jsd, jed
      
      logical, allocatable :: mask(:,:)
      integer, allocatable :: pes(:), pearray(:,:)
      character(len=8) :: text

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS2D: You must first call mpp_domains_init.' )
      ndivx = layout(1); ndivy = layout(2)
      isg = global_indices(1); ieg = global_indices(2); jsg = global_indices(3); jeg = global_indices(4)

      if( PRESENT(pelist) )then
          allocate( pes(0:size(pelist)-1) )
          pes = pelist
      else
          allocate( pes(0:mpp_npes()-1) )
          call mpp_get_current_pelist(pes)
!          pes = (/ (i,i=0,mpp_npes()-1) /)
      end if

      allocate( mask(0:ndivx-1,0:ndivy-1) )
      mask = .TRUE.
      if( PRESENT(maskmap) )then
          if( size(maskmap,1).NE.ndivx .OR. size(maskmap,2).NE.ndivy ) &
               call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS2D: maskmap array does not match layout.' )
          mask(:,:) = maskmap(:,:)
      end if
!number of unmask domains in layout must equal number of PEs assigned
      n = count(mask)
      if( n.NE.size(pes) )then
          write( text,'(i8)' )n
          call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS2D: incorrect number of PEs assigned for this layout and maskmap. Use ' &
               //text//' PEs for this domain decomposition.' )
      end if

!place on PE array; need flag to assign them to j first and then i
      allocate( pearray(0:ndivx-1,0:ndivy-1) )
      pearray(:,:) = NULL_PE
      ipos = NULL_PE; jpos = NULL_PE; pos = NULL_PE
      n = 0
      do j = 0,ndivy-1
         do i = 0,ndivx-1
            if( mask(i,j) )then
                pearray(i,j) = pes(n)
                if( pes(n).EQ.pe )then
                    pos = n
                    ipos = i
                    jpos = j
                end if
                n = n + 1
            end if
         end do
      end do
      if( ipos.EQ.NULL_PE .OR. jpos.EQ.NULL_PE .or. pos.EQ.NULL_PE ) &
           call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS2D: pelist must include this PE.' )
      if( debug )write( stderr(), * )'pe, ipos, jpos=', pe, ipos, jpos, ' pearray(:,jpos)=', pearray(:,jpos), &
                                                                        ' pearray(ipos,:)=', pearray(ipos,:)
      
!do domain decomposition using 1D versions in X and Y
      call mpp_define_domains( global_indices(1:2), ndivx, domain%x, &
           pack(pearray(:,jpos),mask(:,jpos)), xflags, xhalo, xextent, mask(:,jpos) )
      call mpp_define_domains( global_indices(3:4), ndivy, domain%y, &
           pack(pearray(ipos,:),mask(ipos,:)), yflags, yhalo, yextent, mask(ipos,:) )
      if( domain%x%list(domain%x%pos)%pe.NE.domain%y%list(domain%y%pos)%pe ) &
           call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS2D: domain%x%list(ipos)%pe.NE.domain%y%list(jpos)%pe.' ) 
      domain%pos = pos
      domain%pe  = pe

!set up fold
      domain%fold = 0
      if( PRESENT(xflags) )then
          if( BTEST(xflags,WEST) )domain%fold = domain%fold + FOLD_WEST_EDGE
          if( BTEST(xflags,EAST) )domain%fold = domain%fold + FOLD_EAST_EDGE
      end if
      if( PRESENT(yflags) )then
          if( BTEST(yflags,SOUTH) )domain%fold = domain%fold + FOLD_SOUTH_EDGE
          if( BTEST(yflags,NORTH) )domain%fold = domain%fold + FOLD_NORTH_EDGE
      end if
          
      if( BTEST(domain%fold,SOUTH) .OR. BTEST(domain%fold,NORTH) )then
          if( domain%y%cyclic )call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS: an axis cannot be both folded and cyclic.' )
          if( modulo(domain%x%global%size,2).NE.0 ) &
               call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS: number of points in X must be even when there is a fold in Y.' )
!check if folded domain boundaries line up in X: compute domains lining up is a sufficient condition for symmetry
          n = ndivx - 1
          do i = 0,n/2
             if( domain%x%list(i)%compute%size.NE.domain%x%list(n-i)%compute%size ) &
                  call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS: Folded domain boundaries must line up (mirror-symmetric extents).' )
          end do
      end if
      if( BTEST(domain%fold,WEST) .OR. BTEST(domain%fold,EAST) )then
          if( domain%x%cyclic )call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS: an axis cannot be both folded and cyclic.' )
          if( modulo(domain%y%global%size,2).NE.0 ) &
               call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS: number of points in Y must be even when there is a fold in X.' )
!check if folded domain boundaries line up in Y: compute domains lining up is a sufficient condition for symmetry
          n = ndivy - 1
          do i = 0,n/2
             if( domain%y%list(i)%compute%size.NE.domain%y%list(n-i)%compute%size ) &
                  call mpp_error( FATAL, 'MPP_DEFINE_DOMAINS: Folded domain boundaries must line up (mirror-symmetric extents).' )
          end do
      end if

!set up domain%list
      if( debug )write( stderr(),'(a,9i4)' )'pe, domain=', pe, domain_info_buf(1:8)
      if( pe.EQ.mpp_root_pe() .AND. PRESENT(name) )then
          write( stdlog(), '(/a,i3,a,i3)' )trim(name)//' domain decomposition: ', ndivx, ' X', ndivy
          write( stdlog(), '(3x,a)' )'pe,   is,  ie,  js,  je,    isd, ied, jsd, jed'
      end if
      call mpp_sync()
      call mpp_get_compute_domain( domain, domain_info_buf(1), domain_info_buf(2), domain_info_buf(3), domain_info_buf(4) )
      call mpp_get_data_domain   ( domain, domain_info_buf(5), domain_info_buf(6), domain_info_buf(7), domain_info_buf(8) )
      n = size(pes)
      allocate( domain%list(0:n-1) ) !this is only used for storage of remote compute and data domain info
      do i = 0,n-1
         m = mod(pos+i,n)
         domain%list(m)%pe = pes(m)
         call mpp_transmit( domain_info_buf(1:8), 8, pes(mod(pos+n-i,n)), domain_info_buf(9:16), 8, pes(m) )
         domain%list(m)%x%compute%begin = domain_info_buf(9)
         domain%list(m)%x%compute%end   = domain_info_buf(10)
         domain%list(m)%y%compute%begin = domain_info_buf(11)
         domain%list(m)%y%compute%end   = domain_info_buf(12)
         domain%list(m)%x%data%begin = domain_info_buf(13)
         domain%list(m)%x%data%end   = domain_info_buf(14)
         domain%list(m)%y%data%begin = domain_info_buf(15)
         domain%list(m)%y%data%end   = domain_info_buf(16)
         if( pe.EQ.mpp_root_pe() .AND. PRESENT(name) )write( stdlog(), '(2x,i3,x,4i5,3x,4i5)' )pes(m), domain_info_buf(9:)
      end do
      call mpp_sync_self(pes)
      domain%list(:)%x%compute%size = domain%list(:)%x%compute%end - domain%list(:)%x%compute%begin + 1
      domain%list(:)%y%compute%size = domain%list(:)%y%compute%end - domain%list(:)%y%compute%begin + 1
      domain%list(:)%x%data%size = domain%list(:)%x%data%end - domain%list(:)%x%data%begin + 1
      domain%list(:)%y%data%size = domain%list(:)%y%data%end - domain%list(:)%y%data%begin + 1

      domain%remote_domains_initialized = .FALSE.
      domain%remote_off_domains_initialized = .FALSE.
      call compute_overlaps(domain)
!PV786667: the deallocate stmts can be removed when fixed (7.3.1.3m)
      deallocate( pes, mask, pearray )
          
      return
    end subroutine mpp_define_domains2D

    subroutine mpp_broadcast_domain( domain )
!broadcast domain (useful only outside the context of its own pelist)
      type(domain2D), intent(inout) :: domain
      integer, allocatable :: pes(:)
      logical :: native         !true if I'm on the pelist of this domain
      integer :: listsize, listpos
      integer :: n
      integer, dimension(5) :: msg, info         !pe and compute domain of each item in list

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_BROADCAST_DOMAIN: You must first call mpp_domains_init.' )

!get the current pelist
      allocate( pes(0:mpp_npes()-1) )
      call mpp_get_current_pelist(pes)

!am I part of this domain?
      native = ASSOCIATED(domain%list)

!set local list size
      if( native )then
          listsize = size(domain%list)
      else
          listsize = 0
      end if
      call mpp_max(listsize)

      if( .NOT.native )then
!initialize domain%list and set null values in message
          allocate( domain%list(0:listsize-1) )
          domain%pe = NULL_PE
          domain%pos = -1
          domain%x%compute%begin =  HUGE(1)
          domain%x%compute%end   = -HUGE(1)
          domain%y%compute%begin =  HUGE(1)
          domain%y%compute%end   = -HUGE(1)
      end if
!initialize values in info
      info(1) = domain%pe
      call mpp_get_compute_domain( domain, info(2), info(3), info(4), info(5) )

!broadcast your info across current pelist and unpack if needed
      listpos = 0
      do n = 0,mpp_npes()-1
         msg = info
         if( pe.EQ.pes(n) .AND. debug )write( stderr(),* )'PE ', pe, 'broadcasting msg ', msg
         call mpp_broadcast( msg, 5, pes(n) )
!no need to unpack message if native
!no need to unpack message from non-native PE
         if( .NOT.native .AND. msg(1).NE.NULL_PE )then
             domain%list(listpos)%pe = msg(1)
             domain%list(listpos)%x%compute%begin = msg(2)
             domain%list(listpos)%x%compute%end   = msg(3)
             domain%list(listpos)%y%compute%begin = msg(4)
             domain%list(listpos)%y%compute%end   = msg(5)
             listpos = listpos + 1
             if( debug )write( stderr(),* )'PE ', pe, 'received domain from PE ', msg(1), 'is,ie,js,je=', msg(2:5)
         end if
      end do
    end subroutine mpp_broadcast_domain

    subroutine compute_overlaps( domain )
!computes remote domain overlaps
!assumes only one in each direction
      type(domain2D), intent(inout) :: domain
      integer :: i, j, k, m, n
      integer :: is, ie, js, je, isc, iec, jsc, jec, isd, ied, jsd, jed, isg, ieg, jsg, jeg, ioff, joff
      integer :: list

      if( grid_offset_type.EQ.AGRID .AND. domain%remote_domains_initialized     )return
      if( grid_offset_type.NE.AGRID .AND. domain%remote_off_domains_initialized )return
      domain%gridtype = grid_offset_type
      n = size(domain%list)
!send
      call mpp_get_compute_domain( domain, isc, iec, jsc, jec )
      call mpp_get_global_domain ( domain, isg, ieg, jsg, jeg, xsize=ioff, ysize=joff ) !cyclic offsets
      domain%list(:)%overlap = .FALSE.
      do list = 0,n-1
         m = mod( domain%pos+list, n )
!to_pe's eastern halo
         is = domain%list(m)%x%compute%end+1; ie = domain%list(m)%x%data%end
         js = domain%list(m)%y%compute%begin; je = domain%list(m)%y%compute%end
         if( is.GT.ieg )then
             if( domain%x%cyclic )then !try cyclic offset
                 is = is-ioff; ie = ie-ioff
             else if( BTEST(domain%fold,EAST) )then
                 i=is; is = 2*ieg-ie+1; ie = 2*ieg-i+1
                 j=js; js = jsg+jeg-je; je = jsg+jeg-j
                 if( BTEST(grid_offset_type,EAST) )then
                     is = is - 1; ie = ie - 1
                 end if
             end if
         end if
         is = max(is,isc); ie = min(ie,iec)
         js = max(js,jsc); je = min(je,jec)
         if( ie.GE.is .AND. je.GE.js )then
             domain%list(m)%overlap = .TRUE.
             if( grid_offset_type.NE.AGRID )then
                 domain%list(m)%send_w_off%overlap = .TRUE.
                 domain%list(m)%send_w_off%is = is
                 domain%list(m)%send_w_off%ie = ie
                 domain%list(m)%send_w_off%js = js
                 domain%list(m)%send_w_off%je = je
             else
                 domain%list(m)%send_w%overlap = .TRUE.
                 domain%list(m)%send_w%is = is
                 domain%list(m)%send_w%ie = ie
                 domain%list(m)%send_w%js = js
                 domain%list(m)%send_w%je = je
             end if
         else
             domain%list(m)%send_w%overlap = .FALSE.
             domain%list(m)%send_w_off%overlap = .FALSE.
         end if
!to_pe's SE halo
         is = domain%list(m)%x%compute%end+1; ie = domain%list(m)%x%data%end
         js = domain%list(m)%y%data%begin; je = domain%list(m)%y%compute%begin-1
         if( is.GT.ieg )then
             if( domain%x%cyclic )then !try cyclic offset
                 is = is-ioff; ie = ie-ioff
             else if( BTEST(domain%fold,EAST) )then
                 i=is; is = 2*ieg-ie+1; ie = 2*ieg-i+1
                 j=js; js = jsg+jeg-je; je = jsg+jeg-j
                 if( BTEST(grid_offset_type,EAST) )then
                     is = is - 1; ie = ie - 1
                 end if
             end if
         end if
         if( jsg.GT.je )then
             if( domain%y%cyclic )then !try cyclic offset
                 js = js+joff; je = je+joff
             else if( BTEST(domain%fold,SOUTH) )then
                 i=is; is = isg+ieg-ie; ie = isg+ieg-i
                 j=js; js = 2*jsg-je-1; je = 2*jsg-j-1
                 if( BTEST(grid_offset_type,SOUTH) )then
                     js = js + 1; je = je + 1
                 end if
             end if
         end if
         is = max(is,isc); ie = min(ie,iec)
         js = max(js,jsc); je = min(je,jec)
         if( ie.GE.is .AND. je.GE.js )then
             domain%list(m)%overlap = .TRUE.
             if( grid_offset_type.NE.AGRID )then
                 domain%list(m)%send_nw_off%overlap = .TRUE.
                 domain%list(m)%send_nw_off%is = is
                 domain%list(m)%send_nw_off%ie = ie
                 domain%list(m)%send_nw_off%js = js
                 domain%list(m)%send_nw_off%je = je
             else
                 domain%list(m)%send_nw%overlap = .TRUE.
                 domain%list(m)%send_nw%is = is
                 domain%list(m)%send_nw%ie = ie
                 domain%list(m)%send_nw%js = js
                 domain%list(m)%send_nw%je = je
             end if
         else
             domain%list(m)%send_nw%overlap = .FALSE.
             domain%list(m)%send_nw_off%overlap = .FALSE.
         end if
!to_pe's southern halo
         is = domain%list(m)%x%compute%begin; ie = domain%list(m)%x%compute%end
         js = domain%list(m)%y%data%begin; je = domain%list(m)%y%compute%begin-1
         if( jsg.GT.je )then
             if( domain%y%cyclic )then !try cyclic offset
                 js = js+joff; je = je+joff
             else if( BTEST(domain%fold,SOUTH) )then
                 i=is; is = isg+ieg-ie; ie = isg+ieg-i
                 j=js; js = 2*jsg-je-1; je = 2*jsg-j-1
                 if( BTEST(grid_offset_type,SOUTH) )then
                     js = js + 1; je = je + 1
                 end if
             end if
         end if
         is = max(is,isc); ie = min(ie,iec)
         js = max(js,jsc); je = min(je,jec)
         if( ie.GE.is .AND. je.GE.js )then
             domain%list(m)%overlap = .TRUE.
             if( grid_offset_type.NE.AGRID )then
                 domain%list(m)%send_n_off%overlap = .TRUE.
                 domain%list(m)%send_n_off%is = is
                 domain%list(m)%send_n_off%ie = ie
                 domain%list(m)%send_n_off%js = js
                 domain%list(m)%send_n_off%je = je
             else
                 domain%list(m)%send_n%overlap = .TRUE.
                 domain%list(m)%send_n%is = is
                 domain%list(m)%send_n%ie = ie
                 domain%list(m)%send_n%js = js
                 domain%list(m)%send_n%je = je
             end if
         else
             domain%list(m)%send_n%overlap = .FALSE.
             domain%list(m)%send_n_off%overlap = .FALSE.
         end if
!to_pe's SW halo
         is = domain%list(m)%x%data%begin; ie = domain%list(m)%x%compute%begin-1
         js = domain%list(m)%y%data%begin; je = domain%list(m)%y%compute%begin-1
         if( isg.GT.ie )then
             if( domain%x%cyclic )then !try cyclic offset
                 is = is+ioff; ie = ie+ioff
             else if( BTEST(domain%fold,WEST) )then
                 i=is; is = 2*isg-ie-1; ie = 2*isg-i-1
                 j=js; js = jsg+jeg-je; je = jsg+jeg-j
                 if( BTEST(grid_offset_type,WEST) )then
                     is = is + 1; ie = ie + 1
                 end if
             end if
         end if
         if( jsg.GT.je )then
             if( domain%y%cyclic )then !try cyclic offset
                 js = js+joff; je = je+joff
             else if( BTEST(domain%fold,SOUTH) )then
                 i=is; is = isg+ieg-ie; ie = isg+ieg-i
                 j=js; js = 2*jsg-je-1; je = 2*jsg-j-1
                 if( BTEST(grid_offset_type,SOUTH) )then
                     js = js + 1; je = je + 1
                 end if
             end if
         end if
         is = max(is,isc); ie = min(ie,iec)
         js = max(js,jsc); je = min(je,jec)
         if( ie.GE.is .AND. je.GE.js )then
             domain%list(m)%overlap = .TRUE.
             if( grid_offset_type.NE.AGRID )then
                 domain%list(m)%send_ne_off%overlap = .TRUE.
                 domain%list(m)%send_ne_off%is = is
                 domain%list(m)%send_ne_off%ie = ie
                 domain%list(m)%send_ne_off%js = js
                 domain%list(m)%send_ne_off%je = je
             else
                 domain%list(m)%send_ne%overlap = .TRUE.
                 domain%list(m)%send_ne%is = is
                 domain%list(m)%send_ne%ie = ie
                 domain%list(m)%send_ne%js = js
                 domain%list(m)%send_ne%je = je
             end if
         else
             domain%list(m)%send_ne%overlap = .FALSE.
             domain%list(m)%send_ne_off%overlap = .FALSE.
         end if
!to_pe's western halo
         is = domain%list(m)%x%data%begin; ie = domain%list(m)%x%compute%begin-1
         js = domain%list(m)%y%compute%begin; je = domain%list(m)%y%compute%end
         if( isg.GT.ie )then
             if( domain%x%cyclic )then !try cyclic offset
                 is = is+ioff; ie = ie+ioff
             else if( BTEST(domain%fold,WEST) )then
                 i=is; is = 2*isg-ie-1; ie = 2*isg-i-1
                 j=js; js = jsg+jeg-je; je = jsg+jeg-j
                 if( BTEST(grid_offset_type,WEST) )then
                     is = is + 1; ie = ie + 1
                 end if
             end if
         end if
         is = max(is,isc); ie = min(ie,iec)
         js = max(js,jsc); je = min(je,jec)
         if( ie.GE.is .AND. je.GE.js )then
             domain%list(m)%overlap = .TRUE.
             if( grid_offset_type.NE.AGRID )then
                 domain%list(m)%send_e_off%overlap = .TRUE.
                 domain%list(m)%send_e_off%is = is
                 domain%list(m)%send_e_off%ie = ie
                 domain%list(m)%send_e_off%js = js
                 domain%list(m)%send_e_off%je = je
             else
                 domain%list(m)%send_e%overlap = .TRUE.
                 domain%list(m)%send_e%is = is
                 domain%list(m)%send_e%ie = ie
                 domain%list(m)%send_e%js = js
                 domain%list(m)%send_e%je = je
             end if
         else
             domain%list(m)%send_e%overlap = .FALSE.
             domain%list(m)%send_e_off%overlap = .FALSE.
         end if
!to_pe's NW halo
         is = domain%list(m)%x%data%begin; ie = domain%list(m)%x%compute%begin-1
         js = domain%list(m)%y%compute%end+1; je = domain%list(m)%y%data%end
         if( isg.GT.ie )then
             if( domain%x%cyclic )then !try cyclic offset
                 is = is+ioff; ie = ie+ioff
             else if( BTEST(domain%fold,WEST) )then
                 i=is; is = 2*isg-ie-1; ie = 2*isg-i-1
                 j=js; js = jsg+jeg-je; je = jsg+jeg-j
                 if( BTEST(grid_offset_type,WEST) )then
                     is = is + 1; ie = ie + 1
                 end if
             end if
         end if
         if( js.GT.jeg )then
             if( domain%y%cyclic )then !try cyclic offset
                 js = js-joff; je = je-joff
             else if( BTEST(domain%fold,NORTH) )then
                 i=is; is = isg+ieg-ie; ie = isg+ieg-i
                 j=js; js = 2*jeg-je+1; je = 2*jeg-j+1
                 if( BTEST(grid_offset_type,NORTH) )then
                     js = js - 1; je = je - 1
                 end if
             end if
         end if
         is = max(is,isc); ie = min(ie,iec)
         js = max(js,jsc); je = min(je,jec)
         if( ie.GE.is .AND. je.GE.js )then
             domain%list(m)%overlap = .TRUE.
             if( grid_offset_type.NE.AGRID )then
                 domain%list(m)%send_se_off%overlap = .TRUE.
                 domain%list(m)%send_se_off%is = is
                 domain%list(m)%send_se_off%ie = ie
                 domain%list(m)%send_se_off%js = js
                 domain%list(m)%send_se_off%je = je
             else
                 domain%list(m)%send_se%overlap = .TRUE.
                 domain%list(m)%send_se%is = is
                 domain%list(m)%send_se%ie = ie
                 domain%list(m)%send_se%js = js
                 domain%list(m)%send_se%je = je
             end if
         else
             domain%list(m)%send_se%overlap = .FALSE.
             domain%list(m)%send_se_off%overlap = .FALSE.
         end if
!to_pe's northern halo
         is = domain%list(m)%x%compute%begin; ie = domain%list(m)%x%compute%end
         js = domain%list(m)%y%compute%end+1; je = domain%list(m)%y%data%end
         if( js.GT.jeg )then
             if( domain%y%cyclic )then !try cyclic offset
                 js = js-joff; je = je-joff
             else if( BTEST(domain%fold,NORTH) )then
                 i=is; is = isg+ieg-ie; ie = isg+ieg-i
                 j=js; js = 2*jeg-je+1; je = 2*jeg-j+1
                 if( BTEST(grid_offset_type,NORTH) )then
                     js = js - 1; je = je - 1
                 end if
             end if
         end if
         is = max(is,isc); ie = min(ie,iec)
         js = max(js,jsc); je = min(je,jec)
         if( ie.GE.is .AND. je.GE.js )then
             domain%list(m)%overlap = .TRUE.
             if( grid_offset_type.NE.AGRID )then
                 domain%list(m)%send_s_off%overlap = .TRUE.
                 domain%list(m)%send_s_off%is = is
                 domain%list(m)%send_s_off%ie = ie
                 domain%list(m)%send_s_off%js = js
                 domain%list(m)%send_s_off%je = je
             else
                 domain%list(m)%send_s%overlap = .TRUE.
                 domain%list(m)%send_s%is = is
                 domain%list(m)%send_s%ie = ie
                 domain%list(m)%send_s%js = js
                 domain%list(m)%send_s%je = je
             end if
         else
             domain%list(m)%send_s%overlap = .FALSE.
             domain%list(m)%send_s_off%overlap = .FALSE.
         end if
!to_pe's NE halo
         is = domain%list(m)%x%compute%end+1; ie = domain%list(m)%x%data%end
         js = domain%list(m)%y%compute%end+1; je = domain%list(m)%y%data%end
         if( is.GT.ieg )then
             if( domain%x%cyclic )then !try cyclic offset
                 is = is-ioff; ie = ie-ioff
             else if( BTEST(domain%fold,EAST) )then
                 i=is; is = 2*ieg-ie+1; ie = 2*ieg-i+1
                 j=js; js = jsg+jeg-je; je = jsg+jeg-j
             end if
         end if
         if( js.GT.jeg )then
             if( domain%y%cyclic )then !try cyclic offset
                 js = js-joff; je = je-joff
             else if( BTEST(domain%fold,NORTH) )then
                 i=is; is = isg+ieg-ie; ie = isg+ieg-i
                 j=js; js = 2*jeg-je+1; je = 2*jeg-j+1
                 if( BTEST(grid_offset_type,NORTH) )then
                     js = js - 1; je = je - 1
                 end if
             end if
         end if
         is = max(is,isc); ie = min(ie,iec)
         js = max(js,jsc); je = min(je,jec)
         if( ie.GE.is .AND. je.GE.js )then
             domain%list(m)%overlap = .TRUE.
             if( grid_offset_type.NE.AGRID )then
                 domain%list(m)%send_sw_off%overlap = .TRUE.
                 domain%list(m)%send_sw_off%is = is
                 domain%list(m)%send_sw_off%ie = ie
                 domain%list(m)%send_sw_off%js = js
                 domain%list(m)%send_sw_off%je = je
             else
                 domain%list(m)%send_sw%overlap = .TRUE.
                 domain%list(m)%send_sw%is = is
                 domain%list(m)%send_sw%ie = ie
                 domain%list(m)%send_sw%js = js
                 domain%list(m)%send_sw%je = je
             end if
         else
             domain%list(m)%send_sw%overlap = .FALSE.
             domain%list(m)%send_sw_off%overlap = .FALSE.
         end if
      end do
             
!recv
      do list = 0,n-1
         m = mod( domain%pos+n-list, n )
         call mpp_get_compute_domain( domain%list(m), isc, iec, jsc, jec )
!recv_e
         isd = domain%x%compute%end+1; ied = domain%x%data%end
         jsd = domain%y%compute%begin; jed = domain%y%compute%end
         is=isc; ie=iec; js=jsc; je=jec
         domain%list(m)%recv_e%folded = .FALSE.
         if( isd.GT.ieg )then
             if( domain%x%cyclic )then !try cyclic offset
                 is = is+ioff; ie = ie+ioff
             else if( BTEST(domain%fold,EAST) )then
                 domain%list(m)%recv_e%folded = .TRUE.
                 i=is; is = 2*ieg-ie+1; ie = 2*ieg-i+1
                 j=js; js = jsg+jeg-je; je = jsg+jeg-j
                 if( BTEST(grid_offset_type,EAST) )then
                     is = is - 1; ie = ie - 1
                 end if
             end if
         end if
         is = max(isd,is); ie = min(ied,ie)
         js = max(jsd,js); je = min(jed,je)
         if( ie.GE.is .AND. je.GE.js )then
             domain%list(m)%overlap = .TRUE.
             if( grid_offset_type.NE.AGRID )then
                 domain%list(m)%recv_e_off%overlap = .TRUE.
                 domain%list(m)%recv_e_off%is = is
                 domain%list(m)%recv_e_off%ie = ie
                 domain%list(m)%recv_e_off%js = js
                 domain%list(m)%recv_e_off%je = je
             else
                 domain%list(m)%recv_e%overlap = .TRUE.
                 domain%list(m)%recv_e%is = is
                 domain%list(m)%recv_e%ie = ie
                 domain%list(m)%recv_e%js = js
                 domain%list(m)%recv_e%je = je
             endif
         else
             domain%list(m)%recv_e%overlap = .FALSE.
             domain%list(m)%recv_e_off%overlap = .FALSE.
         end if
!recv_se
         isd = domain%x%compute%end+1; ied = domain%x%data%end
         jsd = domain%y%data%begin; jed = domain%y%compute%begin-1
         is=isc; ie=iec; js=jsc; je=jec
         domain%list(m)%recv_se%folded = .FALSE.
         if( jed.LT.jsg )then
             if( domain%y%cyclic )then !try cyclic offset
                 js = js-joff; je = je-joff
             else if( BTEST(domain%fold,SOUTH) )then
                 i=is; is = isg+ieg-ie; ie = isg+ieg-i
                 j=js; js = 2*jsg-je-1; je = 2*jsg-j-1
                 domain%list(m)%recv_se%folded = .TRUE.
                 if( BTEST(grid_offset_type,SOUTH) )then
                     js = js + 1; je = je + 1
                 end if
             end if
         end if
         if( isd.GT.ieg )then
             if( domain%x%cyclic )then !try cyclic offset
                 is = is+ioff; ie = ie+ioff
             else if( BTEST(domain%fold,EAST) )then
                 i=is; is = 2*ieg-ie+1; ie = 2*ieg-i+1
                 j=js; js = jsg+jeg-je; je = jsg+jeg-j
                 domain%list(m)%recv_se%folded = .TRUE.
                 if( BTEST(grid_offset_type,EAST) )then
                     is = is - 1; ie = ie - 1
                 end if
             end if
         end if
         is = max(isd,is); ie = min(ied,ie)
         js = max(jsd,js); je = min(jed,je)
         if( ie.GE.is .AND. je.GE.js )then
             domain%list(m)%overlap = .TRUE.
             if( grid_offset_type.NE.AGRID )then
                 domain%list(m)%recv_se_off%overlap = .TRUE.
                 domain%list(m)%recv_se_off%is = is
                 domain%list(m)%recv_se_off%ie = ie
                 domain%list(m)%recv_se_off%js = js
                 domain%list(m)%recv_se_off%je = je
             else
                 domain%list(m)%recv_se%overlap = .TRUE.
                 domain%list(m)%recv_se%is = is
                 domain%list(m)%recv_se%ie = ie
                 domain%list(m)%recv_se%js = js
                 domain%list(m)%recv_se%je = je
             endif
         else
             domain%list(m)%recv_se%overlap = .FALSE.
             domain%list(m)%recv_se_off%overlap = .FALSE.
         end if
!recv_s
         isd = domain%x%compute%begin; ied = domain%x%compute%end
         jsd = domain%y%data%begin; jed = domain%y%compute%begin-1
         is=isc; ie=iec; js=jsc; je=jec
         domain%list(m)%recv_s%folded = .FALSE.
         if( jed.LT.jsg )then
             if( domain%y%cyclic )then !try cyclic offset
                 js = js-joff; je = je-joff
             else if( BTEST(domain%fold,SOUTH) )then
                 i=is; is = isg+ieg-ie; ie = isg+ieg-i
                 j=js; js = 2*jsg-je-1; je = 2*jsg-j-1
                 domain%list(m)%recv_s%folded = .TRUE.
                 if( BTEST(grid_offset_type,SOUTH) )then
                     js = js + 1; je = je + 1
                 end if
             end if
         end if
         is = max(isd,is); ie = min(ied,ie)
         js = max(jsd,js); je = min(jed,je)
         if( ie.GE.is .AND. je.GE.js )then
             domain%list(m)%overlap = .TRUE.
             if( grid_offset_type.NE.AGRID )then
                 domain%list(m)%recv_s_off%overlap = .TRUE.
                 domain%list(m)%recv_s_off%is = is
                 domain%list(m)%recv_s_off%ie = ie
                 domain%list(m)%recv_s_off%js = js
                 domain%list(m)%recv_s_off%je = je
             else
                 domain%list(m)%recv_s%overlap = .TRUE.
                 domain%list(m)%recv_s%is = is
                 domain%list(m)%recv_s%ie = ie
                 domain%list(m)%recv_s%js = js
                 domain%list(m)%recv_s%je = je
             endif
         else
             domain%list(m)%recv_s%overlap = .FALSE.
             domain%list(m)%recv_s_off%overlap = .FALSE.

         end if
!recv_sw
         isd = domain%x%data%begin; ied = domain%x%compute%begin-1
         jsd = domain%y%data%begin; jed = domain%y%compute%begin-1
         is=isc; ie=iec; js=jsc; je=jec
         domain%list(m)%recv_sw%folded = .FALSE.
         if( jed.LT.jsg )then
             if( domain%y%cyclic )then !try cyclic offset
                 js = js-joff; je = je-joff
             else if( BTEST(domain%fold,SOUTH) )then
                 i=is; is = isg+ieg-ie; ie = isg+ieg-i
                 j=js; js = 2*jsg-je-1; je = 2*jsg-j-1
                 domain%list(m)%recv_sw%folded = .TRUE.
                 if( BTEST(grid_offset_type,SOUTH) )then
                     js = js + 1; je = je + 1
                 end if
             end if
         end if
         if( ied.LT.isg )then
             if( domain%x%cyclic )then !try cyclic offset
                 is = is-ioff; ie = ie-ioff
             else if( BTEST(domain%fold,WEST) )then
                 i=is; is = 2*isg-ie-1; ie = 2*isg-i-1
                 j=js; js = jsg+jeg-je; je = jsg+jeg-j
                 domain%list(m)%recv_sw%folded = .TRUE.
                 if( BTEST(grid_offset_type,WEST) )then
                     is = is + 1; ie = ie + 1
                 end if
             end if
         end if
         is = max(isd,is); ie = min(ied,ie)
         js = max(jsd,js); je = min(jed,je)
         if( ie.GE.is .AND. je.GE.js )then
             domain%list(m)%overlap = .TRUE.
             if( grid_offset_type.NE.AGRID )then
                 domain%list(m)%recv_sw_off%overlap = .TRUE.
                 domain%list(m)%recv_sw_off%is = is
                 domain%list(m)%recv_sw_off%ie = ie
                 domain%list(m)%recv_sw_off%js = js
                 domain%list(m)%recv_sw_off%je = je
             else
                 domain%list(m)%recv_sw%overlap = .TRUE.
                 domain%list(m)%recv_sw%is = is
                 domain%list(m)%recv_sw%ie = ie
                 domain%list(m)%recv_sw%js = js
                 domain%list(m)%recv_sw%je = je
             endif
         else
             domain%list(m)%recv_sw%overlap = .FALSE.
             domain%list(m)%recv_sw_off%overlap = .FALSE.
         end if
!recv_w
         isd = domain%x%data%begin; ied = domain%x%compute%begin-1
         jsd = domain%y%compute%begin; jed = domain%y%compute%end
         is=isc; ie=iec; js=jsc; je=jec
         domain%list(m)%recv_w%folded = .FALSE.
         if( ied.LT.isg )then
             if( domain%x%cyclic )then !try cyclic offset
                 is = is-ioff; ie = ie-ioff
             else if( BTEST(domain%fold,WEST) )then
                 i=is; is = 2*isg-ie-1; ie = 2*isg-i-1
                 j=js; js = jsg+jeg-je; je = jsg+jeg-j
                 domain%list(m)%recv_w%folded = .TRUE.
                 if( BTEST(grid_offset_type,WEST) )then
                     is = is + 1; ie = ie + 1
                 end if
             end if
         end if
         is = max(isd,is); ie = min(ied,ie)
         js = max(jsd,js); je = min(jed,je)
         if( ie.GE.is .AND. je.GE.js )then
             domain%list(m)%overlap = .TRUE.
             if( grid_offset_type.NE.AGRID )then
                 domain%list(m)%recv_w_off%overlap = .TRUE.
                 domain%list(m)%recv_w_off%is = is
                 domain%list(m)%recv_w_off%ie = ie
                 domain%list(m)%recv_w_off%js = js
                 domain%list(m)%recv_w_off%je = je
             else
                 domain%list(m)%recv_w%overlap = .TRUE.
                 domain%list(m)%recv_w%is = is
                 domain%list(m)%recv_w%ie = ie
                 domain%list(m)%recv_w%js = js
                 domain%list(m)%recv_w%je = je
             endif
         else
             domain%list(m)%recv_w%overlap = .FALSE.
             domain%list(m)%recv_w_off%overlap = .FALSE.
         end if
!recv_nw
         isd = domain%x%data%begin; ied = domain%x%compute%begin-1
         jsd = domain%y%compute%end+1; jed = domain%y%data%end
         is=isc; ie=iec; js=jsc; je=jec
         domain%list(m)%recv_nw%folded = .FALSE.
         if( jsd.GT.jeg )then
             if( domain%y%cyclic )then !try cyclic offset
                 js = js+joff; je = je+joff
             else if( BTEST(domain%fold,NORTH) )then
                 i=is; is = isg+ieg-ie; ie = isg+ieg-i
                 j=js; js = 2*jeg-je+1; je = 2*jeg-j+1
                 domain%list(m)%recv_nw%folded = .TRUE.
                 if( BTEST(grid_offset_type,NORTH) )then
                     js = js - 1; je = je - 1
                 end if
             end if
         end if
         if( ied.LT.isg )then
             if( domain%x%cyclic )then !try cyclic offset
                 is = is-ioff; ie = ie-ioff
             else if( BTEST(domain%fold,WEST) )then
                 i=is; is = 2*isg-ie-1; ie = 2*isg-i-1
                 j=js; js = jsg+jeg-je; je = jsg+jeg-j
                 domain%list(m)%recv_nw%folded = .TRUE.
                 if( BTEST(grid_offset_type,WEST) )then
                     is = is + 1; ie = ie + 1
                 end if
             end if
         end if
         is = max(isd,is); ie = min(ied,ie)
         js = max(jsd,js); je = min(jed,je)
         if( ie.GE.is .AND. je.GE.js )then
             domain%list(m)%overlap = .TRUE.
             if( grid_offset_type.NE.AGRID )then
                 domain%list(m)%recv_nw_off%overlap = .TRUE.
                 domain%list(m)%recv_nw_off%is = is
                 domain%list(m)%recv_nw_off%ie = ie
                 domain%list(m)%recv_nw_off%js = js
                 domain%list(m)%recv_nw_off%je = je
             else
                 domain%list(m)%recv_nw%overlap = .TRUE.
                 domain%list(m)%recv_nw%is = is
                 domain%list(m)%recv_nw%ie = ie
                 domain%list(m)%recv_nw%js = js
                 domain%list(m)%recv_nw%je = je
             endif
         else
             domain%list(m)%recv_nw%overlap = .FALSE.
             domain%list(m)%recv_nw_off%overlap = .FALSE.
         end if
!recv_n
         isd = domain%x%compute%begin; ied = domain%x%compute%end
         jsd = domain%y%compute%end+1; jed = domain%y%data%end
         is=isc; ie=iec; js=jsc; je=jec
         domain%list(m)%recv_n%folded = .FALSE.
         if( jsd.GT.jeg )then
             if( domain%y%cyclic )then !try cyclic offset
                 js = js+joff; je = je+joff
             else if( BTEST(domain%fold,NORTH) )then
                 i=is; is = isg+ieg-ie; ie = isg+ieg-i
                 j=js; js = 2*jeg-je+1; je = 2*jeg-j+1
                 domain%list(m)%recv_n%folded = .TRUE.
                 if( BTEST(grid_offset_type,NORTH) )then
                     js = js - 1; je = je - 1
                 end if
             end if
         end if
         is = max(isd,is); ie = min(ied,ie)
         js = max(jsd,js); je = min(jed,je)
         if( ie.GE.is .AND. je.GE.js )then
             domain%list(m)%overlap = .TRUE.
             if( grid_offset_type.NE.AGRID )then
                 domain%list(m)%recv_n_off%overlap = .TRUE.
                 domain%list(m)%recv_n_off%is = is
                 domain%list(m)%recv_n_off%ie = ie
                 domain%list(m)%recv_n_off%js = js
                 domain%list(m)%recv_n_off%je = je
             else
                 domain%list(m)%recv_n%overlap = .TRUE.
                 domain%list(m)%recv_n%is = is
                 domain%list(m)%recv_n%ie = ie
                 domain%list(m)%recv_n%js = js
                 domain%list(m)%recv_n%je = je
             end if
         else
             domain%list(m)%recv_n%overlap = .FALSE.
             domain%list(m)%recv_n_off%overlap = .FALSE.
         end if
!recv_ne
         isd = domain%x%compute%end+1; ied = domain%x%data%end
         jsd = domain%y%compute%end+1; jed = domain%y%data%end
         is=isc; ie=iec; js=jsc; je=jec
         domain%list(m)%recv_ne%folded = .FALSE.
         if( jsd.GT.jeg )then
             if( domain%y%cyclic )then !try cyclic offset
                 js = js+joff; je = je+joff
             else if( BTEST(domain%fold,NORTH) )then
                 i=is; is = isg+ieg-ie; ie = isg+ieg-i
                 j=js; js = 2*jeg-je+1; je = 2*jeg-j+1
                 domain%list(m)%recv_ne%folded = .TRUE.
                 if( BTEST(grid_offset_type,NORTH) )then
                     js = js - 1; je = je - 1
                 end if
             end if
         end if
         if( isd.GT.ieg )then
             if( domain%x%cyclic )then !try cyclic offset
                 is = is+ioff; ie = ie+ioff
             else if( BTEST(domain%fold,EAST) )then
                 i=is; is = 2*ieg-ie+1; ie = 2*ieg-i+1
                 j=js; js = jsg+jeg-je; je = jsg+jeg-j
                 domain%list(m)%recv_ne%folded = .TRUE.
                 if( BTEST(grid_offset_type,EAST) )then
                     is = is - 1; ie = ie - 1
                 end if
             end if
         end if
         is = max(isd,is); ie = min(ied,ie)
         js = max(jsd,js); je = min(jed,je)
         if( ie.GE.is .AND. je.GE.js )then
             domain%list(m)%overlap = .TRUE.
             if( grid_offset_type.NE.AGRID )then
                 domain%list(m)%recv_ne_off%overlap = .TRUE.
                 domain%list(m)%recv_ne_off%is = is
                 domain%list(m)%recv_ne_off%ie = ie
                 domain%list(m)%recv_ne_off%js = js
                 domain%list(m)%recv_ne_off%je = je
             else
                 domain%list(m)%recv_ne%overlap = .TRUE.
                 domain%list(m)%recv_ne%is = is
                 domain%list(m)%recv_ne%ie = ie
                 domain%list(m)%recv_ne%js = js
                 domain%list(m)%recv_ne%je = je
             end if
         else
             domain%list(m)%recv_ne%overlap = .FALSE.
             domain%list(m)%recv_ne_off%overlap = .FALSE.
         end if
      end do
      if( grid_offset_type.EQ.AGRID )domain%remote_domains_initialized = .TRUE.
      if( grid_offset_type.NE.AGRID )domain%remote_off_domains_initialized = .TRUE.
      return
    end subroutine compute_overlaps
! <SUBROUTINE NAME="mpp_define_layout2D" INTERFACE="mpp_define_layout">
!  <IN NAME="global_indices" TYPE="integer" DIM="(4)"></IN>
!  <IN NAME="ndivs" TYPE="integer"></IN>
!  <OUT NAME="layout" TYPE="integer" DIM="(2)"></OUT>
! </SUBROUTINE>
    subroutine mpp_define_layout2D( global_indices, ndivs, layout )
      integer, intent(in) :: global_indices(4) !(/ isg, ieg, jsg, jeg /)
      integer, intent(in) :: ndivs !number of divisions to divide global domain
      integer, intent(out) :: layout(2)

      integer :: isg, ieg, jsg, jeg, isz, jsz, idiv, jdiv

      isg = global_indices(1)
      ieg = global_indices(2)
      jsg = global_indices(3)
      jeg = global_indices(4)

      isz = ieg - isg + 1
      jsz = jeg - jsg + 1
!first try to divide ndivs in the domain aspect ratio: if imperfect aspect, reduce idiv till it divides ndivs
      idiv = nint( sqrt(float(ndivs*isz)/jsz) )
      idiv = max(idiv,1) !for isz=1 line above can give 0
      do while( mod(ndivs,idiv).NE.0 )
         idiv = idiv - 1
      end do                 !will terminate at idiv=1 if not before
      jdiv = ndivs/idiv

      layout = (/ idiv, jdiv /)
      return
    end subroutine mpp_define_layout2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!     MPP_GET and SET routiness: retrieve various components of domains       !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine mpp_get_compute_domain1D( domain, begin, end, size, max_size, is_global )
      type(domain1D), intent(in) :: domain
      integer, intent(out), optional :: begin, end, size, max_size
      logical, intent(out), optional :: is_global

      if( PRESENT(begin)     )begin     = domain%compute%begin
      if( PRESENT(end)       )end       = domain%compute%end
      if( PRESENT(size)      )size      = domain%compute%size
      if( PRESENT(max_size)  )max_size  = domain%compute%max_size
      if( PRESENT(is_global) )is_global = domain%compute%is_global
      return
    end subroutine mpp_get_compute_domain1D

    subroutine mpp_get_data_domain1D( domain, begin, end, size, max_size, is_global )
      type(domain1D), intent(in) :: domain
      integer, intent(out), optional :: begin, end, size, max_size
      logical, intent(out), optional :: is_global

      if( PRESENT(begin)     )begin     = domain%data%begin
      if( PRESENT(end)       )end       = domain%data%end
      if( PRESENT(size)      )size      = domain%data%size
      if( PRESENT(max_size)  )max_size  = domain%data%max_size
      if( PRESENT(is_global) )is_global = domain%data%is_global
      return
    end subroutine mpp_get_data_domain1D

    subroutine mpp_get_global_domain1D( domain, begin, end, size, max_size )
      type(domain1D), intent(in) :: domain
      integer, intent(out), optional :: begin, end, size, max_size

      if( PRESENT(begin)    )begin    = domain%global%begin
      if( PRESENT(end)      )end      = domain%global%end
      if( PRESENT(size)     )size     = domain%global%size
      if( PRESENT(max_size) )max_size = domain%global%max_size
      return
    end subroutine mpp_get_global_domain1D

    subroutine mpp_get_compute_domain2D( domain, xbegin, xend, ybegin, yend, xsize, xmax_size, ysize, ymax_size, &
         x_is_global, y_is_global )
      type(domain2D), intent(in) :: domain
      integer, intent(out), optional :: xbegin, xend, ybegin, yend, xsize, xmax_size, ysize, ymax_size
      logical, intent(out), optional :: x_is_global, y_is_global
      call mpp_get_compute_domain( domain%x, xbegin, xend, xsize, xmax_size, x_is_global )
      call mpp_get_compute_domain( domain%y, ybegin, yend, ysize, ymax_size, y_is_global )
      return
    end subroutine mpp_get_compute_domain2D

    subroutine mpp_get_data_domain2D( domain, xbegin, xend, ybegin, yend, xsize, xmax_size, ysize, ymax_size, &
         x_is_global, y_is_global )
      type(domain2D), intent(in) :: domain
      integer, intent(out), optional :: xbegin, xend, ybegin, yend, xsize, xmax_size, ysize, ymax_size
      logical, intent(out), optional :: x_is_global, y_is_global
      call mpp_get_data_domain( domain%x, xbegin, xend, xsize, xmax_size, x_is_global )
      call mpp_get_data_domain( domain%y, ybegin, yend, ysize, ymax_size, y_is_global )
      return
    end subroutine mpp_get_data_domain2D

    subroutine mpp_get_global_domain2D( domain, xbegin, xend, ybegin, yend, xsize, xmax_size, ysize, ymax_size )
      type(domain2D), intent(in) :: domain
      integer, intent(out), optional :: xbegin, xend, ybegin, yend, xsize, xmax_size, ysize, ymax_size
      call mpp_get_global_domain( domain%x, xbegin, xend, xsize, xmax_size )
      call mpp_get_global_domain( domain%y, ybegin, yend, ysize, ymax_size )
      return
    end subroutine mpp_get_global_domain2D

! <SUBROUTINE NAME="mpp_get_domain_components">
!  <OVERVIEW>
!    Retrieve 1D components of 2D decomposition.
!  </OVERVIEW>
!  <DESCRIPTION>
!    It is sometime necessary to have direct recourse to the domain1D types
!    that compose a domain2D object. This call retrieves them.
!  </DESCRIPTION>
!  <TEMPLATE>
!    call mpp_get_domain_components( domain, x, y )
!  </TEMPLATE>
!  <IN NAME="domain" TYPE="type(domain2D)"></IN>
!  <OUT NAME="x,y"  TYPE="type(domain1D)"></OUT>
! </SUBROUTINE>
    subroutine mpp_get_domain_components( domain, x, y )
      type(domain2D), intent(in) :: domain
      type(domain1D), intent(out), optional :: x, y
      if( PRESENT(x) )x = domain%x
      if( PRESENT(y) )y = domain%y
      return
    end subroutine mpp_get_domain_components

    subroutine mpp_get_compute_domains1D( domain, begin, end, size )
      type(domain1D), intent(in) :: domain
      integer, intent(out), optional, dimension(:) :: begin, end, size 

      if( .NOT.module_is_initialized ) &
           call mpp_error( FATAL, 'MPP_GET_COMPUTE_DOMAINS: must first call mpp_domains_init.' )
!we use shape instead of size for error checks because size is used as an argument
      if( PRESENT(begin) )then
          if( any(shape(begin).NE.shape(domain%list)) ) &
               call mpp_error( FATAL, 'MPP_GET_COMPUTE_DOMAINS: begin array size does not match domain.' )
          begin(:) = domain%list(:)%compute%begin
      end if
      if( PRESENT(end) )then
          if( any(shape(end).NE.shape(domain%list)) ) &
               call mpp_error( FATAL, 'MPP_GET_COMPUTE_DOMAINS: end array size does not match domain.' )
          end(:) = domain%list(:)%compute%end
      end if
      if( PRESENT(size) )then
          if( any(shape(size).NE.shape(domain%list)) ) &
               call mpp_error( FATAL, 'MPP_GET_COMPUTE_DOMAINS: size array size does not match domain.' )
          size(:) = domain%list(:)%compute%size
      end if
      return
    end subroutine mpp_get_compute_domains1D

    subroutine mpp_get_compute_domains2D( domain, xbegin, xend, xsize, ybegin, yend, ysize )
      type(domain2D), intent(in) :: domain
      integer, intent(out), optional, dimension(:) :: xbegin, xend, xsize, ybegin, yend, ysize

      if( .NOT.module_is_initialized ) &
           call mpp_error( FATAL, 'MPP_GET_COMPUTE_DOMAINS: must first call mpp_domains_init.' )

      if( PRESENT(xbegin) )then
          if( size(xbegin).NE.size(domain%list) ) &
               call mpp_error( FATAL, 'MPP_GET_COMPUTE_DOMAINS: xbegin array size does not match domain.' )
          xbegin(:) = domain%list(:)%x%compute%begin
      end if
      if( PRESENT(xend) )then
          if( size(xend).NE.size(domain%list) ) &
               call mpp_error( FATAL, 'MPP_GET_COMPUTE_DOMAINS: xend array size does not match domain.' )
          xend(:) = domain%list(:)%x%compute%end
      end if
      if( PRESENT(xsize) )then
          if( size(xsize).NE.size(domain%list) ) &
               call mpp_error( FATAL, 'MPP_GET_COMPUTE_DOMAINS: xsize array size does not match domain.' )
          xsize(:) = domain%list(:)%x%compute%size
      end if
      if( PRESENT(ybegin) )then
          if( size(ybegin).NE.size(domain%list) ) &
               call mpp_error( FATAL, 'MPP_GET_COMPUTE_DOMAINS: ybegin array size does not match domain.' )
          ybegin(:) = domain%list(:)%y%compute%begin
      end if
      if( PRESENT(yend) )then
          if( size(yend).NE.size(domain%list) ) &
               call mpp_error( FATAL, 'MPP_GET_COMPUTE_DOMAINS: yend array size does not match domain.' )
          yend(:) = domain%list(:)%y%compute%end
      end if
      if( PRESENT(ysize) )then
          if( size(ysize).NE.size(domain%list) ) &
               call mpp_error( FATAL, 'MPP_GET_COMPUTE_DOMAINS: ysize array size does not match domain.' )
          ysize(:) = domain%list(:)%y%compute%size
      end if
      return
    end subroutine mpp_get_compute_domains2D

! <SUBROUTINE NAME="mpp_get_pelist1D" INTERFACE="mpp_get_pelist">
!  <IN NAME="domain" TYPE="type(domain1D)"></IN>
!  <OUT NAME="pelist" TYPE="integer" DIM="(:)"></OUT>
!  <OUT NAME="pos" TYPE="integer"></OUT>
! </SUBROUTINE>
    subroutine mpp_get_pelist1D( domain, pelist, pos )
      type(domain1D), intent(in) :: domain
      integer, intent(out) :: pelist(:)
      integer, intent(out), optional :: pos
      integer :: ndivs

      if( .NOT.module_is_initialized ) &
           call mpp_error( FATAL, 'MPP_GET_PELIST: must first call mpp_domains_init.' )
      ndivs = size(domain%list)
      
      if( size(pelist).NE.ndivs ) &
           call mpp_error( FATAL, 'MPP_GET_PELIST: pelist array size does not match domain.' )

      pelist(:) = domain%list(0:ndivs-1)%pe
      if( PRESENT(pos) )pos = domain%pos
      return
    end subroutine mpp_get_pelist1D

! <SUBROUTINE NAME="mpp_get_pelist2D" INTERFACE="mpp_get_pelist">
!  <IN NAME="domain" TYPE="type(domain2D)"></IN>
!  <OUT NAME="pelist" TYPE="integer" DIM="(:)"></OUT>
!  <OUT NAME="pos" TYPE="integer"></OUT>
! </SUBROUTINE>
    subroutine mpp_get_pelist2D( domain, pelist, pos )
      type(domain2D), intent(in) :: domain
      integer, intent(out) :: pelist(:)
      integer, intent(out), optional :: pos

      if( .NOT.module_is_initialized ) &
           call mpp_error( FATAL, 'MPP_GET_PELIST: must first call mpp_domains_init.' )
      if( size(pelist).NE.size(domain%list) ) &
           call mpp_error( FATAL, 'MPP_GET_PELIST: pelist array size does not match domain.' )

      pelist(:) = domain%list(:)%pe
      if( PRESENT(pos) )pos = domain%pos
      return
    end subroutine mpp_get_pelist2D
! <SUBROUTINE NAME="mpp_get_layout1D" INTERFACE="mpp_get_layout">
!  <IN NAME="domain" TYPE="type(domain1D)"></IN>
!  <OUT NAME="layout" TYPE="integer"></OUT>
! </SUBROUTINE>
    subroutine mpp_get_layout1D( domain, layout )
      type(domain1D), intent(in) :: domain
      integer, intent(out) :: layout
      
      if( .NOT.module_is_initialized ) &
           call mpp_error( FATAL, 'MPP_GET_LAYOUT: must first call mpp_domains_init.' )

      layout = size(domain%list)
      return
    end subroutine mpp_get_layout1D

! <SUBROUTINE NAME="mpp_get_layout2D" INTERFACE="mpp_get_layout">
!  <IN NAME="domain" TYPE="type(domain2D)"></IN>
!  <OUT NAME="layout" TYPE="integer" DIM="(2)"></OUT>
! </SUBROUTINE>
    subroutine mpp_get_layout2D( domain, layout )
      type(domain2D), intent(in) :: domain
      integer, intent(out) :: layout(2)
      
      if( .NOT.module_is_initialized ) &
           call mpp_error( FATAL, 'MPP_GET_LAYOUT: must first call mpp_domains_init.' )

      layout(1) = size(domain%x%list)
      layout(2) = size(domain%y%list)
      return
    end subroutine mpp_get_layout2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!              MPP_UPDATE_DOMAINS: fill halos for 2D decomposition            !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define VECTOR_FIELD_
#define MPP_TYPE_ real(DOUBLE_KIND)
#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain2D_r8_2D
#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain2D_r8_3D
#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain2D_r8_4D
#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain2D_r8_5D
#ifdef  VECTOR_FIELD_
#define MPP_UPDATE_DOMAINS_2D_V_ mpp_update_domain2D_r8_2Dv
#define MPP_UPDATE_DOMAINS_3D_V_ mpp_update_domain2D_r8_3Dv
#define MPP_UPDATE_DOMAINS_4D_V_ mpp_update_domain2D_r8_4Dv
#define MPP_UPDATE_DOMAINS_5D_V_ mpp_update_domain2D_r8_5Dv
#endif
#define MPP_REDISTRIBUTE_2D_ mpp_redistribute_r8_2D
#define MPP_REDISTRIBUTE_3D_ mpp_redistribute_r8_3D
#define MPP_REDISTRIBUTE_4D_ mpp_redistribute_r8_4D
#define MPP_REDISTRIBUTE_5D_ mpp_redistribute_r8_5D
#include <mpp_update_domains2D.h>
#undef  VECTOR_FIELD_

#define MPP_TYPE_ complex(DOUBLE_KIND)
#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain2D_c8_2D
#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain2D_c8_3D
#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain2D_c8_4D
#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain2D_c8_5D
#define MPP_REDISTRIBUTE_2D_ mpp_redistribute_c8_2D
#define MPP_REDISTRIBUTE_3D_ mpp_redistribute_c8_3D
#define MPP_REDISTRIBUTE_4D_ mpp_redistribute_c8_4D
#define MPP_REDISTRIBUTE_5D_ mpp_redistribute_c8_5D
#include <mpp_update_domains2D.h>

#ifndef no_8byte_integers
#define MPP_TYPE_ integer(LONG_KIND)
#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain2D_i8_2D
#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain2D_i8_3D
#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain2D_i8_4D
#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain2D_i8_5D
#define MPP_REDISTRIBUTE_2D_ mpp_redistribute_i8_2D
#define MPP_REDISTRIBUTE_3D_ mpp_redistribute_i8_3D
#define MPP_REDISTRIBUTE_4D_ mpp_redistribute_i8_4D
#define MPP_REDISTRIBUTE_5D_ mpp_redistribute_i8_5D
#include <mpp_update_domains2D.h>

#define MPP_TYPE_ logical(LONG_KIND)
#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain2D_l8_2D
#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain2D_l8_3D
#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain2D_l8_4D
#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain2D_l8_5D
#define MPP_REDISTRIBUTE_2D_ mpp_redistribute_l8_2D
#define MPP_REDISTRIBUTE_3D_ mpp_redistribute_l8_3D
#define MPP_REDISTRIBUTE_4D_ mpp_redistribute_l8_4D
#define MPP_REDISTRIBUTE_5D_ mpp_redistribute_l8_5D
#include <mpp_update_domains2D.h>
#endif

#ifndef no_4byte_reals
#define VECTOR_FIELD_
#define MPP_TYPE_ real(FLOAT_KIND)
#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain2D_r4_2D
#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain2D_r4_3D
#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain2D_r4_4D
#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain2D_r4_5D
#ifdef  VECTOR_FIELD_
#define MPP_UPDATE_DOMAINS_2D_V_ mpp_update_domain2D_r4_2Dv
#define MPP_UPDATE_DOMAINS_3D_V_ mpp_update_domain2D_r4_3Dv
#define MPP_UPDATE_DOMAINS_4D_V_ mpp_update_domain2D_r4_4Dv
#define MPP_UPDATE_DOMAINS_5D_V_ mpp_update_domain2D_r4_5Dv
#endif
#define MPP_REDISTRIBUTE_2D_ mpp_redistribute_r4_2D
#define MPP_REDISTRIBUTE_3D_ mpp_redistribute_r4_3D
#define MPP_REDISTRIBUTE_4D_ mpp_redistribute_r4_4D
#define MPP_REDISTRIBUTE_5D_ mpp_redistribute_r4_5D
#include <mpp_update_domains2D.h>
#undef  VECTOR_FIELD_

#define MPP_TYPE_ complex(FLOAT_KIND)
#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain2D_c4_2D
#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain2D_c4_3D
#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain2D_c4_4D
#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain2D_c4_5D
#define MPP_REDISTRIBUTE_2D_ mpp_redistribute_c4_2D
#define MPP_REDISTRIBUTE_3D_ mpp_redistribute_c4_3D
#define MPP_REDISTRIBUTE_4D_ mpp_redistribute_c4_4D
#define MPP_REDISTRIBUTE_5D_ mpp_redistribute_c4_5D
#include <mpp_update_domains2D.h>
#endif

#define MPP_TYPE_ integer(INT_KIND)
#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain2D_i4_2D
#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain2D_i4_3D
#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain2D_i4_4D
#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain2D_i4_5D
#define MPP_REDISTRIBUTE_2D_ mpp_redistribute_i4_2D
#define MPP_REDISTRIBUTE_3D_ mpp_redistribute_i4_3D
#define MPP_REDISTRIBUTE_4D_ mpp_redistribute_i4_4D
#define MPP_REDISTRIBUTE_5D_ mpp_redistribute_i4_5D
#include <mpp_update_domains2D.h>

#define MPP_TYPE_ logical(INT_KIND)
#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain2D_l4_2D
#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain2D_l4_3D
#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain2D_l4_4D
#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain2D_l4_5D
#define MPP_REDISTRIBUTE_2D_ mpp_redistribute_l4_2D
#define MPP_REDISTRIBUTE_3D_ mpp_redistribute_l4_3D
#define MPP_REDISTRIBUTE_4D_ mpp_redistribute_l4_4D
#define MPP_REDISTRIBUTE_5D_ mpp_redistribute_l4_5D
#include <mpp_update_domains2D.h>

!#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain1D_r8_2D
!#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain1D_r8_3D
!#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain1D_r8_4D
!#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain1D_r8_5D
!#define MPP_TYPE_ real(DOUBLE_KIND)
!#include <mpp_update_domains1D.h>
!
!#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain1D_c8_2D
!#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain1D_c8_3D
!#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain1D_c8_4D
!#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain1D_c8_5D
!#define MPP_TYPE_ complex(DOUBLE_KIND)
!#include <mpp_update_domains1D.h>
!
!#ifndef no_8byte_integers
!#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain1D_i8_2D
!#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain1D_i8_3D
!#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain1D_i8_4D
!#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain1D_i8_5D
!#define MPP_TYPE_ integer(LONG_KIND)
!#include <mpp_update_domains1D.h>
!
!#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain1D_l8_2D
!#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain1D_l8_3D
!#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain1D_l8_4D
!#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain1D_l8_5D
!#define MPP_TYPE_ logical(LONG_KIND)
!#include <mpp_update_domains1D.h>
!#endif
!
!#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain1D_r4_2D
!#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain1D_r4_3D
!#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain1D_r4_4D
!#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain1D_r4_5D
!#define MPP_TYPE_ real(FLOAT_KIND)
!#include <mpp_update_domains1D.h>
!
!#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain1D_c4_2D
!#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain1D_c4_3D
!#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain1D_c4_4D
!#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain1D_c4_5D
!#define MPP_TYPE_ complex(FLOAT_KIND)
!#include <mpp_update_domains1D.h>
!
!#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain1D_i4_2D
!#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain1D_i4_3D
!#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain1D_i4_4D
!#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain1D_i4_5D
!#define MPP_TYPE_ integer(INT_KIND)
!#include <mpp_update_domains1D.h>
!
!#define MPP_UPDATE_DOMAINS_2D_ mpp_update_domain1D_l4_2D
!#define MPP_UPDATE_DOMAINS_3D_ mpp_update_domain1D_l4_3D
!#define MPP_UPDATE_DOMAINS_4D_ mpp_update_domain1D_l4_4D
!#define MPP_UPDATE_DOMAINS_5D_ mpp_update_domain1D_l4_5D
!#define MPP_TYPE_ logical(INT_KIND)
!#include <mpp_update_domains1D.h>

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!              MPP_GLOBAL_REDUCE: get global max/min of field                 !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define MPP_GLOBAL_REDUCE_2D_ mpp_global_max_r8_2d
#define MPP_GLOBAL_REDUCE_3D_ mpp_global_max_r8_3d
#define MPP_GLOBAL_REDUCE_4D_ mpp_global_max_r8_4d
#define MPP_GLOBAL_REDUCE_5D_ mpp_global_max_r8_5d
#define MPP_TYPE_ real(DOUBLE_KIND)
#define REDUCE_VAL_ maxval
#define REDUCE_LOC_ maxloc
#define MPP_REDUCE_ mpp_max
#include <mpp_global_reduce.h>

#define MPP_GLOBAL_REDUCE_2D_ mpp_global_min_r8_2d
#define MPP_GLOBAL_REDUCE_3D_ mpp_global_min_r8_3d
#define MPP_GLOBAL_REDUCE_4D_ mpp_global_min_r8_4d
#define MPP_GLOBAL_REDUCE_5D_ mpp_global_min_r8_5d
#define MPP_TYPE_ real(DOUBLE_KIND)
#define REDUCE_VAL_ minval
#define REDUCE_LOC_ minloc
#define MPP_REDUCE_ mpp_min
#include <mpp_global_reduce.h>

#ifndef no_4byte_reals
#define MPP_GLOBAL_REDUCE_2D_ mpp_global_max_r4_2d
#define MPP_GLOBAL_REDUCE_3D_ mpp_global_max_r4_3d
#define MPP_GLOBAL_REDUCE_4D_ mpp_global_max_r4_4d
#define MPP_GLOBAL_REDUCE_5D_ mpp_global_max_r4_5d
#define MPP_TYPE_ real(FLOAT_KIND)
#define REDUCE_VAL_ maxval
#define REDUCE_LOC_ maxloc
#define MPP_REDUCE_ mpp_max
#include <mpp_global_reduce.h>

#define MPP_GLOBAL_REDUCE_2D_ mpp_global_min_r4_2d
#define MPP_GLOBAL_REDUCE_3D_ mpp_global_min_r4_3d
#define MPP_GLOBAL_REDUCE_4D_ mpp_global_min_r4_4d
#define MPP_GLOBAL_REDUCE_5D_ mpp_global_min_r4_5d
#define MPP_TYPE_ real(FLOAT_KIND)
#define REDUCE_VAL_ minval
#define REDUCE_LOC_ minloc
#define MPP_REDUCE_ mpp_min
#include <mpp_global_reduce.h>
#endif

#ifndef no_8byte_integers
#define MPP_GLOBAL_REDUCE_2D_ mpp_global_max_i8_2d
#define MPP_GLOBAL_REDUCE_3D_ mpp_global_max_i8_3d
#define MPP_GLOBAL_REDUCE_4D_ mpp_global_max_i8_4d
#define MPP_GLOBAL_REDUCE_5D_ mpp_global_max_i8_5d
#define MPP_TYPE_ integer(LONG_KIND)
#define REDUCE_VAL_ maxval
#define REDUCE_LOC_ maxloc
#define MPP_REDUCE_ mpp_max
#include <mpp_global_reduce.h>

#define MPP_GLOBAL_REDUCE_2D_ mpp_global_min_i8_2d
#define MPP_GLOBAL_REDUCE_3D_ mpp_global_min_i8_3d
#define MPP_GLOBAL_REDUCE_4D_ mpp_global_min_i8_4d
#define MPP_GLOBAL_REDUCE_5D_ mpp_global_min_i8_5d
#define MPP_TYPE_ integer(LONG_KIND)
#define REDUCE_VAL_ minval
#define REDUCE_LOC_ minloc
#define MPP_REDUCE_ mpp_min
#include <mpp_global_reduce.h>
#endif

#define MPP_GLOBAL_REDUCE_2D_ mpp_global_max_i4_2d
#define MPP_GLOBAL_REDUCE_3D_ mpp_global_max_i4_3d
#define MPP_GLOBAL_REDUCE_4D_ mpp_global_max_i4_4d
#define MPP_GLOBAL_REDUCE_5D_ mpp_global_max_i4_5d
#define MPP_TYPE_ integer(INT_KIND)
#define REDUCE_VAL_ maxval
#define REDUCE_LOC_ maxloc
#define MPP_REDUCE_ mpp_max
#include <mpp_global_reduce.h>

#define MPP_GLOBAL_REDUCE_2D_ mpp_global_min_i4_2d
#define MPP_GLOBAL_REDUCE_3D_ mpp_global_min_i4_3d
#define MPP_GLOBAL_REDUCE_4D_ mpp_global_min_i4_4d
#define MPP_GLOBAL_REDUCE_5D_ mpp_global_min_i4_5d
#define MPP_TYPE_ integer(INT_KIND)
#define REDUCE_VAL_ minval
#define REDUCE_LOC_ minloc
#define MPP_REDUCE_ mpp_min
#include <mpp_global_reduce.h>

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                   MPP_GLOBAL_SUM: global sum of field                       !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define MPP_GLOBAL_SUM_ mpp_global_sum_r8_2d
#define MPP_EXTRA_INDICES_
#define MPP_TYPE_ real(DOUBLE_KIND)
#include <mpp_global_sum.h>

#define MPP_GLOBAL_SUM_ mpp_global_sum_r8_3d
#define MPP_EXTRA_INDICES_ ,:
#define MPP_TYPE_ real(DOUBLE_KIND)
#include <mpp_global_sum.h>

#define MPP_GLOBAL_SUM_ mpp_global_sum_r8_4d
#define MPP_EXTRA_INDICES_ ,:,:
#define MPP_TYPE_ real(DOUBLE_KIND)
#include <mpp_global_sum.h>

#define MPP_GLOBAL_SUM_ mpp_global_sum_r8_5d
#define MPP_EXTRA_INDICES_ ,:,:,:
#define MPP_TYPE_ real(DOUBLE_KIND)
#include <mpp_global_sum.h>

#define MPP_GLOBAL_SUM_ mpp_global_sum_c8_2d
#define MPP_EXTRA_INDICES_
#define MPP_TYPE_ complex(DOUBLE_KIND)
#include <mpp_global_sum.h>

#define MPP_GLOBAL_SUM_ mpp_global_sum_c8_3d
#define MPP_EXTRA_INDICES_ ,:
#define MPP_TYPE_ complex(DOUBLE_KIND)
#include <mpp_global_sum.h>

#define MPP_GLOBAL_SUM_ mpp_global_sum_c8_4d
#define MPP_EXTRA_INDICES_ ,:,:
#define MPP_TYPE_ complex(DOUBLE_KIND)
#include <mpp_global_sum.h>

#define MPP_GLOBAL_SUM_ mpp_global_sum_c8_5d
#define MPP_EXTRA_INDICES_ ,:,:,:
#define MPP_TYPE_ complex(DOUBLE_KIND)
#include <mpp_global_sum.h>

#ifndef no_4byte_reals
#define MPP_GLOBAL_SUM_ mpp_global_sum_r4_2d
#define MPP_EXTRA_INDICES_
#define MPP_TYPE_ real(FLOAT_KIND)
#include <mpp_global_sum.h>

#define MPP_GLOBAL_SUM_ mpp_global_sum_r4_3d
#define MPP_EXTRA_INDICES_ ,:
#define MPP_TYPE_ real(FLOAT_KIND)
#include <mpp_global_sum.h>

#define MPP_GLOBAL_SUM_ mpp_global_sum_r4_4d
#define MPP_EXTRA_INDICES_ ,:,:
#define MPP_TYPE_ real(FLOAT_KIND)
#include <mpp_global_sum.h>

#define MPP_GLOBAL_SUM_ mpp_global_sum_r4_5d
#define MPP_EXTRA_INDICES_ ,:,:,:
#define MPP_TYPE_ real(FLOAT_KIND)
#include <mpp_global_sum.h>

#define MPP_GLOBAL_SUM_ mpp_global_sum_c4_2d
#define MPP_EXTRA_INDICES_
#define MPP_TYPE_ complex(FLOAT_KIND)
#include <mpp_global_sum.h>

#define MPP_GLOBAL_SUM_ mpp_global_sum_c4_3d
#define MPP_EXTRA_INDICES_ ,:
#define MPP_TYPE_ complex(FLOAT_KIND)
#include <mpp_global_sum.h>

#define MPP_GLOBAL_SUM_ mpp_global_sum_c4_4d
#define MPP_EXTRA_INDICES_ ,:,:
#define MPP_TYPE_ complex(FLOAT_KIND)
#include <mpp_global_sum.h>

#define MPP_GLOBAL_SUM_ mpp_global_sum_c4_5d
#define MPP_EXTRA_INDICES_ ,:,:,:
#define MPP_TYPE_ complex(FLOAT_KIND)
#include <mpp_global_sum.h>
#endif

#ifndef no_8byte_integers
#define MPP_GLOBAL_SUM_ mpp_global_sum_i8_2d
#define MPP_EXTRA_INDICES_
#define MPP_TYPE_ integer(LONG_KIND)
#include <mpp_global_sum.h>

#define MPP_GLOBAL_SUM_ mpp_global_sum_i8_3d
#define MPP_EXTRA_INDICES_ ,:
#define MPP_TYPE_ integer(LONG_KIND)
#include <mpp_global_sum.h>

#define MPP_GLOBAL_SUM_ mpp_global_sum_i8_4d
#define MPP_EXTRA_INDICES_ ,:,:
#define MPP_TYPE_ integer(LONG_KIND)
#include <mpp_global_sum.h>

#define MPP_GLOBAL_SUM_ mpp_global_sum_i8_5d
#define MPP_EXTRA_INDICES_ ,:,:,:
#define MPP_TYPE_ integer(LONG_KIND)
#include <mpp_global_sum.h>
#endif

#define MPP_GLOBAL_SUM_ mpp_global_sum_i4_2d
#define MPP_EXTRA_INDICES_
#define MPP_TYPE_ integer(INT_KIND)
#include <mpp_global_sum.h>

#define MPP_GLOBAL_SUM_ mpp_global_sum_i4_3d
#define MPP_EXTRA_INDICES_ ,:
#define MPP_TYPE_ integer(INT_KIND)
#include <mpp_global_sum.h>

#define MPP_GLOBAL_SUM_ mpp_global_sum_i4_4d
#define MPP_EXTRA_INDICES_ ,:,:
#define MPP_TYPE_ integer(INT_KIND)
#include <mpp_global_sum.h>

#define MPP_GLOBAL_SUM_ mpp_global_sum_i4_5d
#define MPP_EXTRA_INDICES_ ,:,:,:
#define MPP_TYPE_ integer(INT_KIND)
#include <mpp_global_sum.h>
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!              MPP_GLOBAL_FIELD: get global field from domain field           !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define MPP_GLOBAL_FIELD_2D_ mpp_global_field2D_r8_2d
#define MPP_GLOBAL_FIELD_3D_ mpp_global_field2D_r8_3d
#define MPP_GLOBAL_FIELD_4D_ mpp_global_field2D_r8_4d
#define MPP_GLOBAL_FIELD_5D_ mpp_global_field2D_r8_5d
#define MPP_GLOBAL1D_FIELD_2D_ mpp_global_field1D_r8_2d
#define MPP_TYPE_ real(DOUBLE_KIND)
#include <mpp_global_field.h>

#define MPP_GLOBAL_FIELD_2D_ mpp_global_field2D_c8_2d
#define MPP_GLOBAL_FIELD_3D_ mpp_global_field2D_c8_3d
#define MPP_GLOBAL_FIELD_4D_ mpp_global_field2D_c8_4d
#define MPP_GLOBAL_FIELD_5D_ mpp_global_field2D_c8_5d
#define MPP_GLOBAL1D_FIELD_2D_ mpp_global_field1D_c8_2d
#define MPP_TYPE_ complex(DOUBLE_KIND)
#include <mpp_global_field.h>

#ifndef no_8byte_integers
#define MPP_GLOBAL_FIELD_2D_ mpp_global_field2D_i8_2d
#define MPP_GLOBAL_FIELD_3D_ mpp_global_field2D_i8_3d
#define MPP_GLOBAL_FIELD_4D_ mpp_global_field2D_i8_4d
#define MPP_GLOBAL_FIELD_5D_ mpp_global_field2D_i8_5d
#define MPP_GLOBAL1D_FIELD_2D_ mpp_global_field1D_i8_2d
#define MPP_TYPE_ integer(LONG_KIND)
#include <mpp_global_field.h>

#define MPP_GLOBAL_FIELD_2D_ mpp_global_field2D_l8_2d
#define MPP_GLOBAL_FIELD_3D_ mpp_global_field2D_l8_3d
#define MPP_GLOBAL_FIELD_4D_ mpp_global_field2D_l8_4d
#define MPP_GLOBAL_FIELD_5D_ mpp_global_field2D_l8_5d
#define MPP_GLOBAL1D_FIELD_2D_ mpp_global_field1D_l8_2d
#define MPP_TYPE_ logical(LONG_KIND)
#include <mpp_global_field.h>
#endif

#ifndef no_4byte_reals
#define MPP_GLOBAL_FIELD_2D_ mpp_global_field2D_r4_2d
#define MPP_GLOBAL_FIELD_3D_ mpp_global_field2D_r4_3d
#define MPP_GLOBAL_FIELD_4D_ mpp_global_field2D_r4_4d
#define MPP_GLOBAL_FIELD_5D_ mpp_global_field2D_r4_5d
#define MPP_GLOBAL1D_FIELD_2D_ mpp_global_field1D_r4_2d
#define MPP_TYPE_ real(FLOAT_KIND)
#include <mpp_global_field.h>

#define MPP_GLOBAL_FIELD_2D_ mpp_global_field2D_c4_2d
#define MPP_GLOBAL_FIELD_3D_ mpp_global_field2D_c4_3d
#define MPP_GLOBAL_FIELD_4D_ mpp_global_field2D_c4_4d
#define MPP_GLOBAL_FIELD_5D_ mpp_global_field2D_c4_5d
#define MPP_GLOBAL1D_FIELD_2D_ mpp_global_field1D_c4_2d
#define MPP_TYPE_ complex(FLOAT_KIND)
#include <mpp_global_field.h>
#endif

#define MPP_GLOBAL_FIELD_2D_ mpp_global_field2D_i4_2d
#define MPP_GLOBAL_FIELD_3D_ mpp_global_field2D_i4_3d
#define MPP_GLOBAL_FIELD_4D_ mpp_global_field2D_i4_4d
#define MPP_GLOBAL_FIELD_5D_ mpp_global_field2D_i4_5d
#define MPP_GLOBAL1D_FIELD_2D_ mpp_global_field1D_i4_2d
#define MPP_TYPE_ integer(INT_KIND)
#include <mpp_global_field.h>

#define MPP_GLOBAL_FIELD_2D_ mpp_global_field2D_l4_2d
#define MPP_GLOBAL_FIELD_3D_ mpp_global_field2D_l4_3d
#define MPP_GLOBAL_FIELD_4D_ mpp_global_field2D_l4_4d
#define MPP_GLOBAL_FIELD_5D_ mpp_global_field2D_l4_5d
#define MPP_GLOBAL1D_FIELD_2D_ mpp_global_field1D_l4_2d
#define MPP_TYPE_ logical(INT_KIND)
#include <mpp_global_field.h>

end module mpp_domains_mod

#ifdef test_mpp_domains
program mpp_domains_test
  use mpp_mod
  use mpp_domains_mod
  implicit none
  integer :: pe, npes
  integer :: nx=128, ny=128, nz=40, halo=2, stackmax=32768
  real, dimension(:,:,:), allocatable :: global, gcheck
  integer :: unit=7
  logical :: debug=.FALSE., opened
  namelist / mpp_domains_nml / nx, ny, nz, halo, stackmax, debug
  integer :: i, j, k
  integer :: layout(2)
  integer :: is, ie, js, je, isd, ied, jsd, jed
  integer :: id

  call mpp_init()

  call mpp_set_warn_level(FATAL)
!possibly open a file called mpp_domains.nml
  do
     inquire( unit=unit, opened=opened )
     if( .NOT.opened )exit
     unit = unit + 1
     if( unit.EQ.100 )call mpp_error( FATAL, 'Unable to locate unit number.' )
  end do
  open( unit=unit, status='OLD', file='mpp_domains.nml', err=10 )
  read( unit,mpp_domains_nml )
  close(unit)
10 continue

  pe = mpp_pe()
  npes = mpp_npes()

  if( debug )then
      call mpp_domains_init(MPP_DEBUG)
  else
      call mpp_domains_init(MPP_DOMAIN_TIME)
  end if
  call mpp_domains_set_stack_size(stackmax)

  if( pe.EQ.mpp_root_pe() )print '(a,5i4)', 'npes, nx, ny, nz, halo=', npes, nx, ny, nz, halo

  allocate( global(1-halo:nx+halo,1-halo:ny+halo,nz) )
  allocate( gcheck(nx,ny,nz) )

!fill in global array: with k.iiijjj
  do k = 1,nz
     do j = 1,ny
        do i = 1,nx
           global(i,j,k) = k + i*1e-3 + j*1e-6
        end do
     end do
  end do

  call test_halo_update( 'Simple' ) !includes global field, global sum tests
  call test_halo_update( 'Cyclic' )
  call test_halo_update( 'Folded' ) !includes vector field test

  call test_redistribute( 'Complete pelist' )
  call test_redistribute( 'Overlap  pelist' )
  call test_redistribute( 'Disjoint pelist' )
  
  call mpp_domains_exit()
  call mpp_exit()

contains

  subroutine test_redistribute( type )
!test redistribute between two domains
    character(len=*), intent(in) :: type
    type(domain2D) :: domainx, domainy
    real, allocatable, dimension(:,:,:) :: x, y
    integer, allocatable :: pelist(:)
    integer :: pemax

    pemax = npes/2              !the partial pelist will run from 0...pemax

!select pelists
    select case(type)
    case( 'Complete pelist' )
!both pelists run from 0...npes-1
        allocate( pelist(0:npes-1) )
        pelist = (/ (i,i=0,npes-1) /)
        call mpp_declare_pelist( pelist )
    case( 'Overlap  pelist' )
!one pelist from 0...pemax, other from 0...npes-1
        allocate( pelist(0:pemax) )
        pelist = (/ (i,i=0,pemax) /)
        call mpp_declare_pelist( pelist )
    case( 'Disjoint pelist' )
!one pelist from 0...pemax, other from pemax+1...npes-1
        if( pemax+1.GE.npes )return
        allocate( pelist(0:pemax) )
        pelist = (/ (i,i=0,pemax) /)
        call mpp_declare_pelist( pelist )
        call mpp_declare_pelist( (/ (i,i=pemax+1,npes-1) /))
    case default
        call mpp_error( FATAL, 'TEST_REDISTRIBUTE: no such test: '//type )
    end select

!set up x and y arrays
    select case(type)
    case( 'Complete pelist' )
!set up x array
        call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
        call mpp_define_domains( (/1,nx,1,ny/), layout, domainx, name=type )
        call mpp_get_compute_domain( domainx, is,  ie,  js,  je  )
        call mpp_get_data_domain   ( domainx, isd, ied, jsd, jed )
        allocate( x(isd:ied,jsd:jed,nz) )
        x = 0.
        x(is:ie,js:je,:) = global(is:ie,js:je,:)
!set up y array
        call mpp_define_domains( (/1,nx,1,ny/), (/npes,1/), domainy, name=type )
        call mpp_get_compute_domain( domainy, is,  ie,  js,  je  )
        call mpp_get_data_domain   ( domainy, isd, ied, jsd, jed )
        allocate( y(isd:ied,jsd:jed,nz) )
        y = 0.
    case( 'Overlap  pelist' )
!one pelist from 0...pemax, other from 0...npes-1
!set up x array
        call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
        call mpp_define_domains( (/1,nx,1,ny/), layout, domainx, name=type )
        call mpp_get_compute_domain( domainx, is,  ie,  js,  je  )
        call mpp_get_data_domain   ( domainx, isd, ied, jsd, jed )
        allocate( x(isd:ied,jsd:jed,nz) )
        x = 0.
        x(is:ie,js:je,:) = global(is:ie,js:je,:)
!set up y array
        if( ANY(pelist.EQ.pe) )then
            call mpp_set_current_pelist(pelist)
            call mpp_define_layout( (/1,nx,1,ny/), mpp_npes(), layout )
            call mpp_define_domains( (/1,nx,1,ny/), layout, domainy, name=type )
            call mpp_get_compute_domain( domainy, is,  ie,  js,  je  )
            call mpp_get_data_domain   ( domainy, isd, ied, jsd, jed )
            allocate( y(isd:ied,jsd:jed,nz) )
            y = 0.
        end if
    case( 'Disjoint pelist' )
!one pelist from 0...pemax, other from pemax+1...npes-1

!set up y array
        if( ANY(pelist.EQ.pe) )then
            call mpp_set_current_pelist(pelist)
            call mpp_define_layout( (/1,nx,1,ny/), mpp_npes(), layout )
            call mpp_define_domains( (/1,nx,1,ny/), layout, domainy, name=type )
            call mpp_get_compute_domain( domainy, is,  ie,  js,  je  )
            call mpp_get_data_domain   ( domainy, isd, ied, jsd, jed )
            allocate( y(isd:ied,jsd:jed,nz) )
            y = 0.
        else
!set up x array
            call mpp_set_current_pelist( (/ (i,i=pemax+1,npes-1) /) )
            call mpp_define_layout( (/1,nx,1,ny/), mpp_npes(), layout )
            call mpp_define_domains( (/1,nx,1,ny/), layout, domainx, name=type )
            call mpp_get_compute_domain( domainx, is,  ie,  js,  je  )
            call mpp_get_data_domain   ( domainx, isd, ied, jsd, jed )
            allocate( x(isd:ied,jsd:jed,nz) )
            x = 0.
            x(is:ie,js:je,:) = global(is:ie,js:je,:)
         end if
    end select

!go global and redistribute
    call mpp_set_current_pelist()
    call mpp_broadcast_domain(domainx)
    call mpp_broadcast_domain(domainy)

    id = mpp_clock_id( type, flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    call mpp_redistribute( domainx, x, domainy, y )
    call mpp_clock_end  (id)

!check answers on pelist
    if( ANY(pelist.EQ.pe) )then
        call mpp_set_current_pelist(pelist)
        call mpp_global_field( domainy, y, gcheck )
        call compare_checksums( global(1:nx,1:ny,:), gcheck, type )
    end if

    call mpp_set_current_pelist()
    call mpp_sync()

    return
  end subroutine test_redistribute

  subroutine test_halo_update( type )
    character(len=*), intent(in) :: type
    real, allocatable, dimension(:,:,:) :: x, y
    type(domain2D) :: domain
    real :: lsum, gsum
    
    call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
    select case(type)
    case( 'Simple' )
        call mpp_define_domains( (/1,nx,1,ny/), layout, domain, xhalo=halo, yhalo=halo, name=type )
        global(1-halo:0,    :,:) = 0
        global(nx+1:nx+halo,:,:) = 0
        global(:,    1-halo:0,:) = 0
        global(:,ny+1:ny+halo,:) = 0
    case( 'Cyclic' )
        call mpp_define_domains( (/1,nx,1,ny/), layout, domain, xhalo=halo, yhalo=halo, &
             xflags=CYCLIC_GLOBAL_DOMAIN, yflags=CYCLIC_GLOBAL_DOMAIN, name=type )
        global(1-halo:0,    1:ny,:) = global(nx-halo+1:nx,1:ny,:)
        global(nx+1:nx+halo,1:ny,:) = global(1:halo,      1:ny,:)
        global(1-halo:nx+halo,    1-halo:0,:) = global(1-halo:nx+halo,ny-halo+1:ny,:)
        global(1-halo:nx+halo,ny+1:ny+halo,:) = global(1-halo:nx+halo,1:halo,      :)
    case( 'Folded' )
        call mpp_define_domains( (/1,nx,1,ny/), layout, domain, xhalo=halo, yhalo=halo, &
             xflags=CYCLIC_GLOBAL_DOMAIN, yflags=FOLD_NORTH_EDGE, name=type )
        global(1-halo:0,    1:ny,:) = global(nx-halo+1:nx,1:ny,:)
        global(nx+1:nx+halo,1:ny,:) = global(1:halo,      1:ny,:)
        global(1-halo:nx+halo,ny+1:ny+halo,:) = global(nx+halo:1-halo:-1,ny:ny-halo+1:-1,:)
        global(1-halo:nx+halo,1-halo:0,:) = 0
    case default
        call mpp_error( FATAL, 'TEST_MPP_DOMAINS: no such test: '//type )
    end select

!set up x array
    call mpp_get_compute_domain( domain, is,  ie,  js,  je  )
    call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )
    allocate( x(isd:ied,jsd:jed,nz) )
    x = 0.
    x(is:ie,js:je,:) = global(is:ie,js:je,:)

!partial update
    id = mpp_clock_id( type//' partial', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    call mpp_update_domains( x, domain, NUPDATE+EUPDATE )
    call mpp_clock_end  (id)
    call compare_checksums( x(is:ied,js:jed,:), global(is:ied,js:jed,:), type//' partial' )

!full update
    id = mpp_clock_id( type, flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    call mpp_update_domains( x, domain )
    call mpp_clock_end  (id)
    call compare_checksums( x, global(isd:ied,jsd:jed,:), type )

    select case(type)           !extra tests
    case( 'Simple' )

!test mpp_global_field
        id = mpp_clock_id( type//' global field', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
        call mpp_clock_begin(id)
        call mpp_global_field( domain, x, gcheck )
        call mpp_clock_end  (id)
!compare checksums between global and x arrays
        call compare_checksums( global(1:nx,1:ny,:), gcheck, 'mpp_global_field' )

!test mpp_global_sum
        gsum = sum( global(1:nx,1:ny,:) )
        id = mpp_clock_id( type//' sum', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
        call mpp_clock_begin(id)
        lsum = mpp_global_sum( domain, x )
        call mpp_clock_end  (id)
        if( pe.EQ.mpp_root_pe() )print '(a,2es15.8,a,es12.4)', 'Fast sum=', lsum, gsum
!test exact mpp_global_sum
        id = mpp_clock_id( type//' exact sum', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
        call mpp_clock_begin(id)
        lsum = mpp_global_sum( domain, x, BITWISE_EXACT_SUM )
        call mpp_clock_end  (id)
        if( pe.EQ.mpp_root_pe() )print '(a,2es15.8,a,es12.4)', 'Bitwise-exact sum=', lsum, gsum
    case( 'Folded' )!test vector update: folded, with sign flip at fold
!fill in folded north edge, cyclic east and west edge
        global(1-halo:0,    1:ny,:) = global(nx-halo+1:nx,1:ny,:)
        global(nx+1:nx+halo,1:ny,:) = global(1:halo,      1:ny,:)
        global(1-halo:nx+halo-1,ny+1:ny+halo,:) = -global(nx+halo-1:1-halo:-1,ny-1:ny-halo:-1,:)
        global(nx+halo,ny+1:ny+halo,:) = -global(nx-halo,ny-1:ny-halo:-1,:)
        global(1-halo:nx+halo,1-halo:0,:) = 0

        x = 0.
        x(is:ie,js:je,:) = global(is:ie,js:je,:)
!set up y array
        allocate( y(isd:ied,jsd:jed,nz) )
        y = x
        id = mpp_clock_id( type//' vector', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
        call mpp_clock_begin(id)
        call mpp_update_domains( x, y, domain, gridtype=BGRID_NE )
        call mpp_clock_end  (id)
!redundant points must be equal and opposite
        global(nx/2,ny,:) = 0.  !pole points must have 0 velocity
        global(nx  ,ny,:) = 0.  !pole points must have 0 velocity
        global(nx/2+1:nx-1, ny,:) = -global(nx/2-1:1:-1, ny,:)
        global(1-halo:0,    ny,:) = -global(nx-halo+1:nx,ny,:)
        global(nx+1:nx+halo,ny,:) = -global(1:halo,      ny,:)
        call compare_checksums( x, global(isd:ied,jsd:jed,:), type//' X' )
        call compare_checksums( y, global(isd:ied,jsd:jed,:), type//' Y' )
    end select
    return
  end subroutine test_halo_update

  subroutine compare_checksums( a, b, string )
    real, intent(in), dimension(:,:,:) :: a, b
    character(len=*), intent(in) :: string
    integer :: i, j

    call mpp_sync()
    i = mpp_chksum( a, (/pe/) )
    j = mpp_chksum( b, (/pe/) )
    if( i.EQ.j )then
        if( pe.EQ.mpp_root_pe() )call mpp_error( NOTE, trim(string)//': OK.' )
    else
        call mpp_error( FATAL, trim(string)//': chksums are not OK.' )
    end if
  end subroutine compare_checksums
end program mpp_domains_test
#endif

! <INFO>

!   <COMPILER NAME="">     
!     Any module or program unit using <TT>mpp_domains_mod</TT>
!     must contain the line

!     <PRE>
!     use mpp_domains_mod
!     </PRE>

!     <TT>mpp_domains_mod</TT> <TT>use</TT>s <LINK
!     SRC="mpp.html"><TT>mpp_mod</TT></LINK>, and therefore is subject to the <LINK
!     SRC="mpp.html#COMPILING AND LINKING SOURCE">compiling and linking requirements of that module.</LINK>
!   </COMPILER>
!   <PRECOMP FLAG="">      
!     <TT>mpp_domains_mod</TT> uses standard f90, and has no special
!     requirements. There are some OS-dependent
!     pre-processor directives that you might need to modify on
!     non-SGI/Cray systems and compilers. The <LINK
!     SRC="mpp.html#PORTABILITY">portability of <TT>mpp_mod</TT></LINK>
!     obviously is a constraint, since this module is built on top of
!     it. Contact me, Balaji, SGI/GFDL, with questions.
!   </PRECOMP> 
!   <LOADER FLAG="">       
!     The <TT>mpp_domains</TT> source consists of the main source file
!     <TT>mpp_domains.F90</TT> and also requires the following include files:
!    <PRE>
!     <TT>os.h</TT>
!     <TT>mpp_update_domains2D.h</TT>
!     <TT>mpp_global_reduce.h</TT>
!     <TT>mpp_global_sum.h</TT>
!     <TT>mpp_global_field.h</TT>
!    </PRE>
!    GFDL users can check it out of the main CVS repository as part of
!    the <TT>mpp</TT> CVS module. The current public tag is <TT>galway</TT>.
!    External users can download the latest <TT>mpp</TT> package <LINK SRC=
!    "ftp://ftp.gfdl.gov/pub/vb/mpp/mpp.tar.Z">here</LINK>. Public access
!    to the GFDL CVS repository will soon be made available.

!   </LOADER>

! </INFO>
