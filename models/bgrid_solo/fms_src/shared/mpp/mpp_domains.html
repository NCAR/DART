<!--

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

-->
<html>
<head>
<META http-equiv="Content-Type" content="text/html; charset=EUC-JP">
<title>Module mpp_domains_mod</title>
<link type="text/css" href="http://www.gfdl.noaa.gov/~fms/style/doc.css" rel="stylesheet">
<STYLE TYPE="text/css">
          .fixed {
            font-size:medium;
            font-family:monospace;
            border-style:none;
            border-width:0.1em;
            padding:0.1em;
            color:#663366;
          }
        </STYLE>
</head>
<body>
<a name="TOP"></a><font class="header" size="1"><a href="#PUBLIC INTERFACE">PUBLIC INTERFACE </a>~
          <a href="#PUBLIC DATA">PUBLIC DATA </a>~
          <a href="#PUBLIC ROUTINES">PUBLIC ROUTINES </a>~
          <a href="#NAMELIST">NAMELIST </a>~
          <a href="#DIAGNOSTIC FIELDS">DIAGNOSTIC FIELDS </a>~
          <a href="#ERROR MESSAGES">ERROR MESSAGES </a>~
          <a href="#REFERENCES">REFERENCES </a>~ 
          <a href="#NOTES">NOTES</a></font>
<hr>
<h2>module mpp_domains_mod</h2>
<a name="HEADER"></a>
<!-- BEGIN HEADER -->
<div>
<b>Contact:&nbsp;</b><a href="mailto:vb@gfdl.noaa.gov">   V. Balaji </a>
<br>
<b>Reviewers:&nbsp;</b>
<br>
<b>Change History:&nbsp;</b><a href="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/">WebCVS Log</a>
<br>
<b>RCS Log:&nbsp;</b><a href="http://www.gfdl.noaa.gov/~vb/changes_models/bgrid_solo/fms_src/shared/mpp/mpp_domains.html">RCS Log</a>
<br>
<b>Last Modified:&nbsp;</b>2002/07/19 18:10:06<br>
<br>
</div>
<!-- END HEADER -->
<a name="OVERVIEW"></a>
<hr>
<h4>OVERVIEW</h4>
<!-- BEGIN OVERVIEW -->
<p class="text"> 
<tt>mpp_domains_mod</tt>   is a set of simple calls for domain
   decomposition and domain updates on rectilinear grids. It requires the
   module <a href="models/bgrid_solo/fms_src/shared/mpp/mpp.html"></a>, upon which it is built. </p>
<!-- END OVERVIEW -->
<a name="DESCRIPTION"></a>
<!-- BEGIN DESCRIPTION -->
<div>   Scalable implementations of finite-difference codes are generally
   based on decomposing the model domain into subdomains that are
   distributed among processors. These domains will then be obliged to
   exchange data at their boundaries if data dependencies are merely
   nearest-neighbour, or may need to acquire information from the global
   domain if there are extended data dependencies, as in the spectral
   transform. The domain decomposition is a key operation in the
   development of parallel codes.
   <br>
<br> 
<tt>mpp_domains_mod</tt>   provides a domain decomposition and domain
   update API for <i>rectilinear</i>   grids, built on top of the <a href="models/bgrid_solo/fms_src/shared/mpp/mpp.html"></a>   API for message passing. Features
   of <tt>mpp_domains_mod</tt>   include:
   <br>
<br>
   Simple, minimal API, with free access to underlying API for more complicated stuff.
   <br>
<br>
   Design toward typical use in climate/weather CFD codes.
   <br>
<br> 
<h4>Domains</h4>   I have assumed that domain decomposition will mainly be in 2
   horizontal dimensions, which will in general be the two
   fastest-varying indices. There is a separate implementation of 1D
   decomposition on the fastest-varying index, and 1D decomposition on
   the second index, treated as a special case of 2D decomposition, is
   also possible. We define <i>domain</i>   as the grid associated with a <i>task</i>.
   We define the <i>compute domain</i>   as the set of gridpoints that are
   computed by a task, and the <i>data domain</i>   as the set of points
   that are required by the task for the calculation. There can in
   general be more than 1 task per PE, though often
   the number of domains is the same as the processor count. We define
   the <i>global domain</i>   as the global computational domain of the
   entire model (i.e, the same as the computational domain if run on a
   single processor). 2D domains are defined using a derived type <tt>domain2D</tt>,
   constructed as follows (see comments in code for more details):
   <br>
<br> 
<pre>     type, public :: domain_axis_spec
        private
        integer :: begin, end, size, max_size
        logical :: is_global
     end type domain_axis_spec
     type, public :: domain1D
        private
        type(domain_axis_spec) :: compute, data, global, active
        logical :: mustputb, mustgetb, mustputf, mustgetf, folded
        type(domain1D), pointer, dimension(:) :: list
        integer :: pe              !PE to which this domain is assigned
        integer :: pos
     end type domain1D
domaintypes of higher rank can be constructed from type domain1D
typically we only need 1 and 2D, but could need higher (e.g 3D LES)
some elements are repeated below if they are needed once per domain
     type, public :: domain2D
        private
        type(domain1D) :: x
        type(domain1D) :: y
        type(domain2D), pointer, dimension(:) :: list
        integer :: pe              !PE to which this domain is assigned
        integer :: pos
     end type domain2D
     type(domain1D), public :: NULL_DOMAIN1D
     type(domain2D), public :: NULL_DOMAIN2D</pre>   The <tt>domain2D</tt>   type contains all the necessary information to
   define the global, compute and data domains of each task, as well as the PE
   associated with the task. The PEs from which remote data may be
   acquired to update the data domain are also contained in a linked list
   of neighbours. </div>
<br>
<!-- END DESCRIPTION -->
<a name="OTHER MODULES USED"></a>
<hr>
<h4>OTHER MODULES USED</h4>
<!-- BEGIN OTHER MODULES USED -->
<div>
<pre>mpp_mod</pre>
</div>
<!-- END OTHER MODULES USED -->
<!-- BEGIN PUBLIC INTERFACE -->
<a name="PUBLIC INTERFACE"></a>
<hr>
<h4>PUBLIC INTERFACE</h4>
<div>
<pre>use mpp_domains_mod [, only:  mpp_define_domains,<br>                              mpp_update_domains,<br>                              mpp_redistribute,<br>                              mpp_global_field,<br>                              mpp_global_max,<br>                              mpp_global_sum,<br>                              operator,<br>                              mpp_get_compute_domain,<br>                              mpp_get_compute_domains,<br>                              mpp_get_data_domain,<br>                              mpp_get_global_domain,<br>                              mpp_define_layout,<br>                              mpp_get_pelist,<br>                              mpp_get_layout,<br>                              mpp_domains_init,<br>                              mpp_domains_set_stack_size,<br>                              mpp_domains_exit,<br>                              mpp_get_domain_components ]</pre>
<dl>
<dt>
<a href="#mpp_define_domains">mpp_define_domains</a>:</dt>
<dd>   Set up a domain decomposition. </dd>
<dt>
<a href="#mpp_update_domains">mpp_update_domains</a>:</dt>
<dd>   Halo updates. </dd>
<dt>
<a href="#mpp_redistribute">mpp_redistribute</a>:</dt>
<dd>   Reorganization of distributed global arrays. </dd>
<dt>
<a href="#mpp_global_field">mpp_global_field</a>:</dt>
<dd>   Fill in a global array from domain-decomposed arrays. </dd>
<dt>
<a href="#mpp_global_max">mpp_global_max</a>:</dt>
<dd>   Global max/min of domain-decomposed arrays. </dd>
<dt>
<a href="#mpp_global_sum">mpp_global_sum</a>:</dt>
<dd>   Global sum of domain-decomposed arrays. </dd>
<dt>
<a href="#operator">operator</a>:</dt>
<dd>   Equality/inequality operators for domaintypes. </dd>
<dt>
<a href="#mpp_get_compute_domain">mpp_get_compute_domain</a>:</dt>
<dd>   These routines retrieve the axis specifications associated with the compute domains. </dd>
<dt>
<a href="#mpp_get_compute_domains">mpp_get_compute_domains</a>:</dt>
<dd>   Retrieve the entire array of compute domain extents associated with a decomposition. </dd>
<dt>
<a href="#mpp_get_data_domain">mpp_get_data_domain</a>:</dt>
<dd>   These routines retrieve the axis specifications associated with the data domains. </dd>
<dt>
<a href="#mpp_get_global_domain">mpp_get_global_domain</a>:</dt>
<dd>   These routines retrieve the axis specifications associated with the global domains. </dd>
<dt>
<a href="#mpp_define_layout">mpp_define_layout</a>:</dt>
<dd>   Retrieve layout associated with a domain decomposition. </dd>
<dt>
<a href="#mpp_get_pelist">mpp_get_pelist</a>:</dt>
<dd>   Retrieve list of PEs associated with a domain decomposition. </dd>
<dt>
<a href="#mpp_get_layout">mpp_get_layout</a>:</dt>
<dd>   Retrieve layout associated with a domain decomposition. </dd>
<dt>
<a href="#mpp_domains_init">mpp_domains_init</a>:</dt>
<dd>   Initialize domain decomp package. </dd>
<dt>
<a href="#mpp_domains_set_stack_size">mpp_domains_set_stack_size</a>:</dt>
<dd>   Set user stack size. </dd>
<dt>
<a href="#mpp_domains_exit">mpp_domains_exit</a>:</dt>
<dd>   Exit <tt>mpp_domains_mod</tt>. </dd>
<dt>
<a href="#mpp_get_domain_components">mpp_get_domain_components</a>:</dt>
<dd>   Retrieve 1D components of 2D decomposition. </dd>
</dl>
</div>
<br>
<!-- END PUBLIC INTERFACE -->
<a name="PUBLIC DATA"></a>
<hr>
<h4>PUBLIC DATA</h4>
<!-- BEGIN PUBLIC DATA -->
<div>None.<br>
<br>
</div>
<!-- END PUBLIC DATA -->
<a name="PUBLIC ROUTINES"></a>
<hr>
<h4>PUBLIC ROUTINES</h4>
<!-- BEGIN PUBLIC ROUTINES -->
<ol type="a">
<li>
<a name="mpp_define_domains"></a>
<h4>mpp_define_domains</h4>
<pre>
<b>call mpp_define_domains </b>( global_indices, ndivs, domain, &amp; pelist, flags, halo, extent, maskmap )</pre>
<pre>
<b>call mpp_define_domains </b>( global_indices, layout, domain, pelist, &amp; xflags, yflags, xhalo, yhalo, &amp; xextent, yextent, maskmap, name )</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   There are two forms for the <tt>mpp_define_domains</tt>   call. The 2D
   version is generally to be used but is built by repeated calls to the
   1D version, also provided. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>global_indices&nbsp;&nbsp;&nbsp;</tt></td><td>   Defines the global domain. <br>&nbsp;&nbsp;&nbsp;<span class="type">[integer, dimension(2)]</span>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer, dimension(4)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>ndivs&nbsp;&nbsp;&nbsp;</tt></td><td>   Is the number of domain divisions required. <br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>pelist&nbsp;&nbsp;&nbsp;</tt></td><td>   List of PEs to which the domains are to be assigned. <br>&nbsp;&nbsp;&nbsp;<span class="type">[integer, dimension(0:)]</span>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer, dimension(0:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>flags&nbsp;&nbsp;&nbsp;</tt></td><td>   An optional flag to pass additional information
   about the desired domain topology. Useful flags in a 1D decomposition
   include <tt>GLOBAL_DATA_DOMAIN</tt>   and <tt>CYCLIC_GLOBAL_DOMAIN</tt>. Flags are integers: multiple flags may
   be added together. The flag values are public parameters available by
   use association. <br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>halo&nbsp;&nbsp;&nbsp;</tt></td><td>   Width of the halo. <br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>extent&nbsp;&nbsp;&nbsp;</tt></td><td>   Normally <tt>mpp_define_domains</tt>   attempts
   an even division of the global domain across <tt>ndivs</tt>   domains. The <tt>extent</tt>   array can be used by the user to pass a
   custom domain division. The <tt>extent</tt>   array has <tt>ndivs</tt>   elements and holds the compute domain widths, which should add up to
   cover the global domain exactly. <br>&nbsp;&nbsp;&nbsp;<span class="type">[integer, dimension(0:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>maskmap&nbsp;&nbsp;&nbsp;</tt></td><td>   Some divisions may be masked
   (<tt>maskmap=.FALSE.</tt>) to exclude them from the computation (e.g
   for ocean model domains that are all land). The <tt>maskmap</tt>   array
   is dimensioned <tt>ndivs</tt>   and contains <tt>.TRUE.</tt>   values for
   any domain that must be <i>included</i>   in the computation (default
   all). The <tt>pelist</tt>   array length should match the number of
   domains included in the computation. <br>&nbsp;&nbsp;&nbsp;<span class="type">[logical, dimension(0:)]</span>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[logical, dimension(:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>layout&nbsp;&nbsp;&nbsp;</tt></td><td>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer, dimension(2)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>xflags, yflags&nbsp;&nbsp;&nbsp;</tt></td><td>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>xhalo, yhalo&nbsp;&nbsp;&nbsp;</tt></td><td>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>xextent, yextent&nbsp;&nbsp;&nbsp;</tt></td><td>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer, dimension(0:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>name&nbsp;&nbsp;&nbsp;</tt></td><td>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[character(len=*)]</span></td>
</tr>
</table>
</dd>
<br>
<dt>
<b>INPUT/OUTPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>domain&nbsp;&nbsp;&nbsp;</tt></td><td>   Holds the resulting domain decomposition. <br>&nbsp;&nbsp;&nbsp;<span class="type">[type(domain1D)]</span>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[type(domain2D)]</span></td>
</tr>
</table>
</dd>
<br>
<dt>
<b>NOTE</b>
</dt>
<dd>   For example:
   <br>
<br> 
<pre>    call mpp_define_domains( (/1,100/), 10, domain, &amp;
         flags=GLOBAL_DATA_DOMAIN+CYCLIC_GLOBAL_DOMAIN, halo=2 )</pre>   defines 10 compute domains spanning the range [1,100] of the global
   domain. The compute domains are non-overlapping blocks of 10. All the data
   domains are global, and with a halo of 2 span the range [-1:102]. And
   since the global domain has been declared to be cyclic, <tt>domain(9)%next =&gt; domain(0)</tt>   and <tt>domain(0)%prev =&gt;
   domain(9)</tt>. A field is allocated on the data domain, and computations proceed on
   the compute domain. A call to <a href="#mpp_update_domains"></a>   would fill in
   the values in the halo region: <pre>    call mpp_get_data_domain( domain, isd, ied ) !returns -1 and 102
    call mpp_get_compute_domain( domain, is, ie ) !returns (1,10) on PE 0 ...
    allocate( a(isd:ied) )
    do i = is,ie
       a(i) = &lt;perform computations&gt;
    end do
    call mpp_update_domains( a, domain )</pre>   The call to <tt>mpp_update_domains</tt>   fills in the regions outside
   the compute domain. Since the global domain is cyclic, the values at <tt>i=(-1,0)</tt>   are the same as at <tt>i=(99,100)</tt>; and <tt>i=(101,102)</tt>   are the same as <tt>i=(1,2)</tt>.
   <br>
<br>
   The 2D version is just an extension of this syntax to two
   dimensions.
   <br>
<br>
   The 2D version of the above should generally be used in
   codes, including 1D-decomposed ones, if there is a possibility of
   future evolution toward 2D decomposition. The arguments are similar to
   the 1D case, except that now we have optional arguments <tt>flags</tt>, <tt>halo</tt>, <tt>extent</tt>   and <tt>maskmap</tt>   along two axes.
   <br>
<br> 
<tt>flags</tt>   can now take an additional possible value to fold
   one or more edges. This is done by using flags <tt>FOLD_WEST_EDGE</tt>, <tt>FOLD_EAST_EDGE</tt>, <tt>FOLD_SOUTH_EDGE</tt>   or <tt>FOLD_NORTH_EDGE</tt>. When a fold
   exists (e.g cylindrical domain), vector fields reverse sign upon
   crossing the fold. This parity reversal is performed only in the
   vector version of <a href="#mpp_update_domains"></a>. In
   addition, shift operations may need to be applied to vector fields on
   staggered grids, also described in the vector interface to <tt>mpp_update_domains</tt>.
   <br>
<br> 
<tt>name</tt>   is the name associated with the decomposition,
   e.g <tt>'Ocean model'</tt>. If this argument is present, <tt>mpp_define_domains</tt>   will print the domain decomposition
   generated to <tt>stdlog</tt>.
   <br>
<br>
   Examples:
   <br>
<br> 
<pre>    call mpp_define_domains( (/1,100,1,100/), (/2,2/), domain, xhalo=1 )</pre>   will create the following domain layout: <pre>                   |---------|-----------|-----------|-------------|
                   |domain(1)|domain(2)  |domain(3)  |domain(4)    |
    |--------------|---------|-----------|-----------|-------------|
    |Compute domain|1,50,1,50|51,100,1,50|1,50,51,100|51,100,51,100|
    |--------------|---------|-----------|-----------|-------------|
    |Data domain   |0,51,1,50|50,101,1,50|0,51,51,100|50,101,51,100|
    |--------------|---------|-----------|-----------|-------------|</pre>   Again, we allocate arrays on the data domain, perform computations
   on the compute domain, and call <tt>mpp_update_domains</tt>   to update
   the halo region.
   <br>
<br>
   If we wished to perfom a 1D decomposition along <tt>Y</tt>   on the same global domain, we could use: <pre>    call mpp_define_domains( (/1,100,1,100/), layout=(/4,1/), domain, xhalo=1 )</pre>   This will create the following domain layout: <pre>                   |----------|-----------|-----------|------------|
                   |domain(1) |domain(2)  |domain(3)  |domain(4)   |
    |--------------|----------|-----------|-----------|------------|
    |Compute domain|1,100,1,25|1,100,26,50|1,100,51,75|1,100,76,100|
    |--------------|----------|-----------|-----------|------------|
    |Data domain   |0,101,1,25|0,101,26,50|0,101,51,75|1,101,76,100|
    |--------------|----------|-----------|-----------|------------|</pre> 
</dd>
<br>
<br>
</dl>
</li>
<li>
<a name="mpp_update_domains"></a>
<h4>mpp_update_domains</h4>
<pre>
<b>call mpp_update_domains </b>( field, domain, flags )</pre>
<pre>
<b>call mpp_update_domains </b>( fieldx, fieldy, domain, flags, gridtype )</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd> 
<tt>mpp_update_domains</tt>   is used to perform a halo update of a
   domain-decomposed array on each PE. <tt>MPP_TYPE_</tt>   can be of type <tt>complex</tt>, <tt>integer</tt>, <tt>logical</tt>   or <tt>real</tt>;
   of 4-byte or 8-byte kind; of rank up to 5. The vector version (with
   two input data fields) is only present for <tt>real</tt>   types.
   <br>
<br>
   For 2D domain updates, if there are halos present along both <tt>x</tt>   and <tt>y</tt>, we can choose to update one only, by
   specifying <tt>flags=XUPDATE</tt>   or <tt>flags=YUPDATE</tt>. In
   addition, one-sided updates can be performed by setting <tt>flags</tt>   to any combination of <tt>WUPDATE</tt>, <tt>EUPDATE</tt>, <tt>SUPDATE</tt>   and <tt>NUPDATE</tt>, to update the west, east, north
   and south halos respectively. Any combination of halos may be used by
   adding the requisite flags, e.g: <tt>flags=XUPDATE+SUPDATE</tt>   or <tt>flags=EUPDATE+WUPDATE+SUPDATE</tt>   will update the east, west and
   south halos.
   <br>
<br>
   If a call to <tt>mpp_update_domains</tt>   involves at least one E-W
   halo and one N-S halo, the corners involved will also be updated, i.e,
   in the example above, the SE and SW corners will be updated.
   <br>
<br>
   If <tt>flags</tt>   is not supplied, that is
   equivalent to <tt>flags=XUPDATE+YUPDATE</tt>.
   <br>
<br>
   The vector version is passed the <tt>x</tt>   and <tt>y</tt>   components of a vector field in tandem, and both are updated upon
   return. They are passed together to treat parity issues on various
   grids. For example, on a cubic sphere projection, the <tt>x</tt>   and <tt>y</tt>   components may be interchanged when passing from an
   equatorial cube face to a polar face. For grids with folds, vector
   components change sign on crossing the fold.
   <br>
<br>
   Special treatment at boundaries such as folds is also required for
   staggered grids. The following types of staggered grids are
   recognized:
   <br>
<br>
   1) <tt>AGRID</tt>: values are at grid centers.<br>   2) <tt>BGRID_NE</tt>: vector fields are at the NE vertex of a grid
   cell, i.e: the array elements <tt>u(i,j)</tt>   and <tt>v(i,j)</tt>   are
   actually at (i+&frac12;,j+&frac12;) with respect to the grid centers.<br>   3) <tt>BGRID_SW</tt>: vector fields are at the SW vertex of a grid
   cell, i.e: the array elements <tt>u(i,j)</tt>   and <tt>v(i,j)</tt>   are
   actually at (i-&frac12;,j-&frac12;) with respect to the grid centers.<br>   4) <tt>CGRID_NE</tt>: vector fields are at the N and E faces of a
   grid cell, i.e: the array elements <tt>u(i,j)</tt>   and <tt>v(i,j)</tt>   are actually at (i+&frac12;,j) and (i,j+&frac12;) with respect to the
   grid centers.<br>   5) <tt>CGRID_SW</tt>: vector fields are at the S and W faces of a
   grid cell, i.e: the array elements <tt>u(i,j)</tt>   and <tt>v(i,j)</tt>   are actually at (i-&frac12;,j) and (i,j-&frac12;) with respect to the
   grid centers.
   <br>
<br>
   The gridtypes listed above are all available by use association as
   integer parameters. The scalar version of <tt>mpp_update_domains</tt>   assumes that the values of a scalar field are always at <tt>AGRID</tt>   locations, and no special boundary treatment is required. If vector
   fields are at staggered locations, the optional argument <tt>gridtype</tt>   must be appropriately set for correct treatment at
   boundaries.
   <br>
<br>
   It is safe to apply vector field updates to the appropriate arrays
   irrespective of the domain topology: if the topology requires no
   special treatment of vector fields, specifying <tt>gridtype</tt>   will
   do no harm.
   <br>
<br> 
<tt>mpp_update_domains</tt>   internally buffers the date being sent
   and received into single messages for efficiency. A turnable internal
   buffer area in memory is provided for this purpose by <tt>mpp_domains_mod</tt>. The size of this buffer area can be set by
   the user by calling <a href="models/bgrid_solo/fms_src/shared/mpp/mpp_domains.html#mpp_domains_set_stack_size"> </a>. </dd>
<br>
<br>
</dl>
</li>
<li>
<a name="mpp_redistribute"></a>
<h4>mpp_redistribute</h4>
<pre>
<b>call mpp_redistribute </b>( domain_in, field_in, domain_out, field_out )</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd> 
<tt>mpp_redistribute</tt>   is used to reorganize a distributed
   array. <tt>MPP_TYPE_</tt>   can be of type <tt>integer</tt>, <tt>complex</tt>, or <tt>real</tt>; of 4-byte or 8-byte kind; of rank
   up to 5. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>field_in&nbsp;&nbsp;&nbsp;</tt></td><td> <tt>field_in</tt>   is dimensioned on the data domain of <tt>domain_in</tt>. </td>
</tr>
</table>
</dd>
<br>
<dt>
<b>OUTPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>field_out&nbsp;&nbsp;&nbsp;</tt></td><td> <tt>field_out</tt>   on the data domain of <tt>domain_out</tt>. </td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="mpp_global_field"></a>
<h4>mpp_global_field</h4>
<pre>
<b>call mpp_global_field </b>( domain, local, global, flags )</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd> 
<tt>mpp_global_field</tt>   is used to get an entire
   domain-decomposed array on each PE. <tt>MPP_TYPE_</tt>   can be of type <tt>complex</tt>, <tt>integer</tt>, <tt>logical</tt>   or <tt>real</tt>;
   of 4-byte or 8-byte kind; of rank up to 5.
   <br>
<br>
   All PEs in a domain decomposition must call <tt>mpp_global_field</tt>, and each will have a complete global field
   at the end. Please note that a global array of rank 3 or higher could
   occupy a lot of memory. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>domain&nbsp;&nbsp;&nbsp;</tt></td><td></td>
</tr>
<tr>
<td valign="top" align="left"><tt>local&nbsp;&nbsp;&nbsp;</tt></td><td> <tt>local</tt>   is dimensioned on either the compute domain or the
   data domain of <tt>domain</tt>. </td>
</tr>
<tr>
<td valign="top" align="left"><tt>flags&nbsp;&nbsp;&nbsp;</tt></td><td> <tt>flags</tt>   can be given the value <tt>XONLY</tt>   or <tt>YONLY</tt>, to specify a globalization on one axis only. </td>
</tr>
</table>
</dd>
<br>
<dt>
<b>OUTPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>global&nbsp;&nbsp;&nbsp;</tt></td><td> <tt>global</tt>   is dimensioned on the corresponding global domain. </td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="mpp_global_max"></a>
<h4>mpp_global_max</h4>
<pre> 
<b>mpp_global_max</b> ( domain, field, locus )</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd> 
<tt>mpp_global_max</tt>   is used to get the maximum value of a
   domain-decomposed array on each PE. <tt>MPP_TYPE_</tt>   can be of type <tt>integer</tt>   or <tt>real</tt>; of 4-byte or 8-byte kind; of rank
   up to 5. The dimension of <tt>locus</tt>   must equal the rank of <tt>field</tt>.
   <br>
<br>
   All PEs in a domain decomposition must call <tt>mpp_global_max</tt>, and each will have the result upon exit.
   <br>
<br>
   The function <tt>mpp_global_min</tt>, with an identical syntax. is
   also available. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>domain&nbsp;&nbsp;&nbsp;</tt></td><td></td>
</tr>
<tr>
<td valign="top" align="left"><tt>field&nbsp;&nbsp;&nbsp;</tt></td><td> <tt>field</tt>   is dimensioned on either the compute domain or the
   data domain of <tt>domain</tt>. </td>
</tr>
</table>
</dd>
<br>
<dt>
<b>OUTPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>locus&nbsp;&nbsp;&nbsp;</tt></td><td> <tt>locus</tt>, if present, can be used to retrieve the location of
   the maximum (as in the <tt>MAXLOC</tt>   intrinsic of f90). </td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="mpp_global_sum"></a>
<h4>mpp_global_sum</h4>
<pre>
<b>call mpp_global_sum </b>( domain, field, flags )</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd> 
<tt>mpp_global_sum</tt>   is used to get the sum of a
   domain-decomposed array on each PE. <tt>MPP_TYPE_</tt>   can be of type <tt>integer</tt>, <tt>complex</tt>, or <tt>real</tt>; of 4-byte or
   8-byte kind; of rank up to 5. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>domain&nbsp;&nbsp;&nbsp;</tt></td><td></td>
</tr>
<tr>
<td valign="top" align="left"><tt>field&nbsp;&nbsp;&nbsp;</tt></td><td> <tt>field</tt>   is dimensioned on either the compute domain or the
   data domain of <tt>domain</tt>. </td>
</tr>
<tr>
<td valign="top" align="left"><tt>flags&nbsp;&nbsp;&nbsp;</tt></td><td> <tt>flags</tt>, if present, must have the value <tt>BITWISE_EXACT_SUM</tt>. This produces a sum that is guaranteed to
   produce the identical result irrespective of how the domain is
   decomposed. This method does the sum first along the ranks beyond 2,
   and then calls <a href="#mpp_global_field"></a>   to produce a
   global 2D array which is then summed. The default method, which is
   considerably faster, does a local sum followed by <a href="models/bgrid_solo/fms_src/shared/mpp/mpp.html#mpp_sum"></a>   across the domain
   decomposition. </td>
</tr>
</table>
</dd>
<br>
<dt>
<b>NOTE</b>
</dt>
<dd>   All PEs in a domain decomposition must call <tt>mpp_global_sum</tt>, and each will have the result upon exit. </dd>
<br>
<br>
</dl>
</li>
<li>
<a name="operator"></a>
<h4>operator</h4>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   The module provides public operators to check for
   equality/inequality of domaintypes, e.g:
   <br>
<br> 
<pre>    type(domain1D) :: a, b
    type(domain2D) :: c, d
    ...
    if( a.NE.b )then
        ...
    end if
    if( c==d )then
        ...
    end if</pre>   Domains are considered equal if and only if the start and end
   indices of each of their component global, data and compute domains
   are equal. </dd>
<br>
<br>
</dl>
</li>
<li>
<a name="mpp_get_compute_domain"></a>
<h4>mpp_get_compute_domain</h4>
<pre>
<b>call mpp_get_compute_domain </b>
</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   The domain is a derived type with private elements. These routines 
   retrieve the axis specifications associated with the compute domains
   The 2D version of these is a simple extension of 1D. </dd>
<br>
<br>
</dl>
</li>
<li>
<a name="mpp_get_compute_domains"></a>
<h4>mpp_get_compute_domains</h4>
<pre>
<b>call mpp_get_compute_domains </b>( domain, xbegin, xend, xsize, &amp; ybegin, yend, ysize )</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   Retrieve the entire array of compute domain extents associated with a decomposition. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>domain&nbsp;&nbsp;&nbsp;</tt></td><td></td>
</tr>
</table>
</dd>
<br>
<dt>
<b>OUTPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>xbegin,ybegin&nbsp;&nbsp;&nbsp;</tt></td><td></td>
</tr>
<tr>
<td valign="top" align="left"><tt>xend,yend&nbsp;&nbsp;&nbsp;</tt></td><td></td>
</tr>
<tr>
<td valign="top" align="left"><tt>xsize,ysize&nbsp;&nbsp;&nbsp;</tt></td><td></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="mpp_get_data_domain"></a>
<h4>mpp_get_data_domain</h4>
<pre>
<b>call mpp_get_data_domain </b>
</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   The domain is a derived type with private elements. These routines 
   retrieve the axis specifications associated with the data domains.
   The 2D version of these is a simple extension of 1D. </dd>
<br>
<br>
</dl>
</li>
<li>
<a name="mpp_get_global_domain"></a>
<h4>mpp_get_global_domain</h4>
<pre>
<b>call mpp_get_global_domain </b>
</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   The domain is a derived type with private elements. These routines 
   retrieve the axis specifications associated with the global domains.
   The 2D version of these is a simple extension of 1D. </dd>
<br>
<br>
</dl>
</li>
<li>
<a name="mpp_define_layout"></a>
<h4>mpp_define_layout</h4>
<pre>
<b>call mpp_define_layout </b>( global_indices, ndivs, layout )</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   Given a global 2D domain and the number of divisions in the
   decomposition (<tt>ndivs</tt>: usually the PE count unless some
   domains are masked) this calls returns a 2D domain layout.
   <br>
<br>
   By default, <tt>mpp_define_layout</tt>   will attempt to divide the
   2D index space into domains that maintain the aspect ratio of the
   global domain. If this cannot be done, the algorithm favours domains
   that are longer in <tt>x</tt>   than <tt>y</tt>, a preference that could
   improve vector performance. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>global_indices&nbsp;&nbsp;&nbsp;</tt></td><td>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer, dimension(4)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>ndivs&nbsp;&nbsp;&nbsp;</tt></td><td>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
</table>
</dd>
<br>
<dt>
<b>OUTPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>layout&nbsp;&nbsp;&nbsp;</tt></td><td>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer, dimension(2)]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="mpp_get_pelist"></a>
<h4>mpp_get_pelist</h4>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   The 1D version of this call returns an array of the PEs assigned to this 1D domain
   decomposition. In addition the optional argument <tt>pos</tt>   may be
   used to retrieve the 0-based position of the domain local to the
   calling PE, i.e <tt>domain%list(pos)%pe</tt>   is the local PE,
   as returned by <a href="models/bgrid_solo/fms_src/shared/mpp/mpp.html#mpp_pe"></a>.
   The 2D version of this call is identical to 1D version. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>domain&nbsp;&nbsp;&nbsp;</tt></td><td>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[type(domain1D)]</span>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[type(domain2D)]</span></td>
</tr>
</table>
</dd>
<br>
<dt>
<b>OUTPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>pelist&nbsp;&nbsp;&nbsp;</tt></td><td>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer, dimension(:)]</span>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer, dimension(:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>pos&nbsp;&nbsp;&nbsp;</tt></td><td>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="mpp_get_layout"></a>
<h4>mpp_get_layout</h4>
<pre>
<b>call mpp_get_layout </b>( domain, layout )</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   The 1D version of this call returns the number of divisions that was assigned to this
   decomposition axis. The 2D version of this call returns an array of
   dimension 2 holding the results on two axes. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>domain&nbsp;&nbsp;&nbsp;</tt></td><td>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[type(domain1D)]</span>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[type(domain2D)]</span></td>
</tr>
</table>
</dd>
<br>
<dt>
<b>OUTPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>layout&nbsp;&nbsp;&nbsp;</tt></td><td>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer, dimension(2)]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="mpp_domains_init"></a>
<h4>mpp_domains_init</h4>
<pre>
<b>call mpp_domains_init </b>(flags)</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   Called to initialize the <tt>mpp_domains_mod</tt>   package.
   <br>
<br> 
<tt>flags</tt>   can be set to <tt>MPP_VERBOSE</tt>   to have <tt>mpp_domains_mod</tt>   keep you informed of what it's up
   to. <tt>MPP_DEBUG</tt>   returns even more information for debugging.
   <br>
<br> 
<tt>mpp_domains_init</tt>   will call <tt>mpp_init</tt>, to make sure <a href="models/bgrid_solo/fms_src/shared/mpp/mpp.html"></a>   is initialized. (Repeated
   calls to <tt>mpp_init</tt>   do no harm, so don't worry if you already
   called it). </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>flags&nbsp;&nbsp;&nbsp;</tt></td><td>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="mpp_domains_set_stack_size"></a>
<h4>mpp_domains_set_stack_size</h4>
<pre>
<b>call mpp_domains_set_stack_size </b>(n)</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   This sets the size of an array that is used for internal storage by <tt>mpp_domains</tt>. This array is used, for instance, to buffer the
   data sent and received in halo updates.
   <br>
<br>
   This call has implied global synchronization. It should be
   placed somewhere where all PEs can call it. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>n&nbsp;&nbsp;&nbsp;</tt></td><td>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="mpp_domains_exit"></a>
<h4>mpp_domains_exit</h4>
<pre>
<b>call mpp_domains_exit </b>()</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   Serves no particular purpose, but is provided should you require to
   re-initialize <tt>mpp_domains_mod</tt>, for some odd reason. </dd>
<br>
<br>
</dl>
</li>
<li>
<a name="mpp_get_domain_components"></a>
<h4>mpp_get_domain_components</h4>
<pre>
<b>call mpp_get_domain_components </b>( domain, x, y )</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   It is sometime necessary to have direct recourse to the domain1D types
   that compose a domain2D object. This call retrieves them. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>domain&nbsp;&nbsp;&nbsp;</tt></td><td>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[type(domain2D)]</span></td>
</tr>
</table>
</dd>
<br>
<dt>
<b>OUTPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>x,y&nbsp;&nbsp;&nbsp;</tt></td><td>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[type(domain1D)]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
</ol>
<!-- END PUBLIC ROUTINES -->
<a name="NAMELIST"></a>
<!-- BEGIN NAMELIST -->
<!-- END NAMELIST --><a name="DIAGNOSTIC FIELDS"></a>
<!-- BEGIN DIAGNOSTIC FIELDS -->
<!-- END DIAGNOSTIC FIELDS --><a name="DATA SETS"></a>
<!-- BEGIN DATA SETS -->
<hr>
<h4>DATA SETS</h4>
<div>None.<br>
<br>
</div>
<!-- END DATA SETS -->
<a name="ERROR MESSAGES"></a>
<!-- BEGIN ERROR MESSAGES -->
<hr>
<h4>ERROR MESSAGES</h4>
<div>None.<br>
<br>
</div>
<!-- END ERROR MESSAGES -->
<a name="REFERENCES"></a>
<hr>
<h4>REFERENCES</h4>
<!-- BEGIN REFERENCES -->
<div>
        None.
      </div>
<br>
<!-- END REFERENCES -->
<a name="COMPILER SPECIFICS"></a>
<hr>
<h4>COMPILER SPECIFICS</h4>
<!-- BEGIN COMPILER SPECIFICS -->
<div>
<dl>
<dt></dt>
<dd>   Any module or program unit using <tt>mpp_domains_mod</tt>   must contain the line <pre>     use mpp_domains_mod</pre> 
<tt>mpp_domains_mod</tt> <tt>use</tt>s <a href="models/bgrid_solo/fms_src/shared/mpp/mpp.html"></a>, and therefore is subject to the <a href="models/bgrid_solo/fms_src/shared/mpp/mpp.html#COMPILING AND LINKING SOURCE">compiling and linking requirements of that module.</a> 
</dd>
</dl>
</div>
<br>
<!-- END COMPILER SPECIFICS -->
<a name="PRECOMPILER OPTIONS"></a>
<hr>
<h4>PRECOMPILER OPTIONS</h4>
<!-- BEGIN PRECOMPILER OPTIONS -->
<div>
<dl>
<dt></dt>
<dd> 
<tt>mpp_domains_mod</tt>   uses standard f90, and has no special
   requirements. There are some OS-dependent
   pre-processor directives that you might need to modify on
   non-SGI/Cray systems and compilers. The <a href="models/bgrid_solo/fms_src/shared/mpp/mpp.html#PORTABILITY">portability of </a>   obviously is a constraint, since this module is built on top of
   it. Contact me, Balaji, SGI/GFDL, with questions. </dd>
</dl>
</div>
<br>
<!-- END PRECOMPILER OPTIONS -->
<a name="LOADER OPTIONS"></a>
<hr>
<h4>LOADER OPTIONS</h4>
<!-- BEGIN LOADER -->
<div>
<p>   The    source consists of the main source file    and also requires the following include files:    GFDL users can check it out of the main CVS repository as part of
   the    CVS module. The current public tag is .
   External users can download the latest    package . Public access
   to the GFDL CVS repository will soon be made available. </p>
<pre>        
</pre>
</div>
<!-- END LOADER OPTIONS -->
<a name="TEST PROGRAM"></a>
<hr>
<h4>TEST PROGRAM</h4>
<!-- BEGIN TEST PROGRAM -->
<div>None.<br>
</div>
<br>
<!-- END TEST PROGRAM -->
<a name="KNOWN BUGS"></a>
<hr>
<h4>KNOWN BUGS</h4>
<!-- BEGIN KNOWN BUGS -->
<div>
        None.
      </div>
<br>
<!-- END KNOWN BUGS -->
<a name="NOTES"></a>
<hr>
<h4>NOTES</h4>
<!-- BEGIN NOTES -->
<div>None.<br>
</div>
<br>
<!-- END NOTES -->
<a name="FUTURE PLANS"></a>
<hr>
<h4>FUTURE PLANS</h4>
<!-- BEGIN FUTURE PLANS -->
<div>
        None.
      </div>
<br>
<!-- END FUTURE PLANS -->
<hr>
<div align="right">
<font size="-1"><a href="#TOP">top</a></font>
</div>
</body>
</html>
