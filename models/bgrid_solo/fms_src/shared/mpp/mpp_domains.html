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
<!DOCTYPE html>
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<meta http-equiv="Content-Type" content="text/html; charset=utf-8">
<title>Module mpp_domains_mod</title>
<link type="text/css" href=
"http://www.gfdl.noaa.gov/~fms/style/doc.css" rel="stylesheet">
<style type="text/css">
          .fixed {
            font-size:medium;
            font-family:monospace;
            border-style:none;
            border-width:0.1em;
            padding:0.1em;
            color:#663366;
          }
</style>
</head>
<body>
<a name="TOP" id="TOP"></a><font class="header" size="1"><a href=
"#PUBLIC%20INTERFACE">PUBLIC INTERFACE</a> ~ <a href=
"#PUBLIC%20DATA">PUBLIC DATA</a> ~ <a href=
"#PUBLIC%20ROUTINES">PUBLIC ROUTINES</a> ~ <a href=
"#NAMELIST">NAMELIST</a> ~ <a href=
"#DIAGNOSTIC%20FIELDS">DIAGNOSTIC FIELDS</a> ~ <a href=
"#ERROR%20MESSAGES">ERROR MESSAGES</a> ~ <a href=
"#REFERENCES">REFERENCES</a> ~ <a href="#NOTES">NOTES</a></font>
<hr>
<h1>module mpp_domains_mod</h1>
<a name="HEADER" id="HEADER"></a> <a name="OVERVIEW" id=
"OVERVIEW"></a>
<hr>
<h2>OVERVIEW</h2>
<!-- BEGIN OVERVIEW -->
<p class="text"><code>mpp_domains_mod</code> is a set of simple
calls for domain decomposition and domain updates on rectilinear
grids. It requires the module <a href=
"./mpp.html"></a>, upon which it
is built.</p>
<!-- END OVERVIEW -->
<a name="DESCRIPTION" id="DESCRIPTION"></a> 
<!-- BEGIN DESCRIPTION -->
<p>Scalable implementations of finite-difference codes are
generally based on decomposing the model domain into subdomains
that are distributed among processors. These domains will then be
obliged to exchange data at their boundaries if data dependencies
are merely nearest-neighbour, or may need to acquire information
from the global domain if there are extended data dependencies, as
in the spectral transform. The domain decomposition is a key
operation in the development of parallel codes.</p>

<p><code>mpp_domains_mod</code> provides a domain decomposition and
domain update API for <i>rectilinear</i> grids, built on top of the
<a href="./mpp.html"></a> API
for message passing. Features of <code>mpp_domains_mod</code>
include:</p>

<p>Simple, minimal API, with free access to underlying API for more
complicated stuff.</p>

<p>Design toward typical use in climate/weather CFD codes.</p>

<h2>Domains</h2>

<p>I have assumed that domain decomposition will mainly be in 2
horizontal dimensions, which will in general be the two
fastest-varying indices. There is a separate implementation of 1D
decomposition on the fastest-varying index, and 1D decomposition on
the second index, treated as a special case of 2D decomposition, is
also possible. We define <i>domain</i> as the grid associated with
a <i>task</i>. We define the <i>compute domain</i> as the set of
gridpoints that are computed by a task, and the <i>data domain</i>
as the set of points that are required by the task for the
calculation. There can in general be more than 1 task per PE,
though often the number of domains is the same as the processor
count. We define the <i>global domain</i> as the global
computational domain of the entire model (i.e, the same as the
computational domain if run on a single processor). 2D domains are
defined using a derived type <code>domain2D</code>, constructed as
follows (see comments in code for more details):</p>

<pre>
<code>
type, public :: domain_axis_spec
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
    type(domain2D), public :: NULL_DOMAIN2D
</code>
</pre>

<p>The <code>domain2D</code> type contains all the necessary
information to define the global, compute and data domains of each
task, as well as the PE associated with the task. The PEs from
which remote data may be acquired to update the data domain are
also contained in a linked list of neighbours.</p>

<!-- END DESCRIPTION -->
<a name="OTHER MODULES USED"></a>
<hr>
<h2>OTHER MODULES USED</h2>
<!-- BEGIN OTHER MODULES USED -->
<div>
<pre>mpp_mod</pre></div>
<!-- END OTHER MODULES USED -->
<!-- BEGIN PUBLIC INTERFACE -->
<a name="PUBLIC INTERFACE"></a>
<hr>
<h2>PUBLIC INTERFACE</h2>
<div>
<pre>
use mpp_domains_mod [, only:  mpp_define_domains,<br>                              mpp_update_domains,<br>                              mpp_redistribute,<br>                              mpp_global_field,<br>                              mpp_global_max,<br>                              mpp_global_sum,<br>                              operator,<br>                              mpp_get_compute_domain,<br>                              mpp_get_compute_domains,<br>                              mpp_get_data_domain,<br>                              mpp_get_global_domain,<br>                              mpp_define_layout,<br>                              mpp_get_pelist,<br>                              mpp_get_layout,<br>                              mpp_domains_init,<br>                              mpp_domains_set_stack_size,<br>                              mpp_domains_exit,<br>                              mpp_get_domain_components ]</pre>
<dl>
<dt><a href="#mpp_define_domains">mpp_define_domains</a>:</dt>
<dd>Set up a domain decomposition.</dd>
<dt><a href="#mpp_update_domains">mpp_update_domains</a>:</dt>
<dd>Halo updates.</dd>
<dt><a href="#mpp_redistribute">mpp_redistribute</a>:</dt>
<dd>Reorganization of distributed global arrays.</dd>
<dt><a href="#mpp_global_field">mpp_global_field</a>:</dt>
<dd>Fill in a global array from domain-decomposed arrays.</dd>
<dt><a href="#mpp_global_max">mpp_global_max</a>:</dt>
<dd>Global max/min of domain-decomposed arrays.</dd>
<dt><a href="#mpp_global_sum">mpp_global_sum</a>:</dt>
<dd>Global sum of domain-decomposed arrays.</dd>
<dt><a href="#operator">operator</a>:</dt>
<dd>Equality/inequality operators for domaintypes.</dd>
<dt><a href=
"#mpp_get_compute_domain">mpp_get_compute_domain</a>:</dt>
<dd>These routines retrieve the axis specifications associated with
the compute domains.</dd>
<dt><a href=
"#mpp_get_compute_domains">mpp_get_compute_domains</a>:</dt>
<dd>Retrieve the entire array of compute domain extents associated
with a decomposition.</dd>
<dt><a href="#mpp_get_data_domain">mpp_get_data_domain</a>:</dt>
<dd>These routines retrieve the axis specifications associated with
the data domains.</dd>
<dt><a href=
"#mpp_get_global_domain">mpp_get_global_domain</a>:</dt>
<dd>These routines retrieve the axis specifications associated with
the global domains.</dd>
<dt><a href="#mpp_define_layout">mpp_define_layout</a>:</dt>
<dd>Retrieve layout associated with a domain decomposition.</dd>
<dt><a href="#mpp_get_pelist">mpp_get_pelist</a>:</dt>
<dd>Retrieve list of PEs associated with a domain
decomposition.</dd>
<dt><a href="#mpp_get_layout">mpp_get_layout</a>:</dt>
<dd>Retrieve layout associated with a domain decomposition.</dd>
<dt><a href="#mpp_domains_init">mpp_domains_init</a>:</dt>
<dd>Initialize domain decomp package.</dd>
<dt><a href=
"#mpp_domains_set_stack_size">mpp_domains_set_stack_size</a>:</dt>
<dd>Set user stack size.</dd>
<dt><a href="#mpp_domains_exit">mpp_domains_exit</a>:</dt>
<dd>Exit <code>mpp_domains_mod</code>.</dd>
<dt><a href=
"#mpp_get_domain_components">mpp_get_domain_components</a>:</dt>
<dd>Retrieve 1D components of 2D decomposition.</dd>
</dl>
</div>
<br>
<!-- END PUBLIC INTERFACE -->
<a name="PUBLIC DATA"></a>
<hr>
<h2>PUBLIC DATA</h2>
<!-- BEGIN PUBLIC DATA -->
<div>None.<br>
<br></div>
<!-- END PUBLIC DATA -->
<a name="PUBLIC ROUTINES"></a>
<hr>
<h2>PUBLIC ROUTINES</h2>
<!-- BEGIN PUBLIC ROUTINES -->
<ol type="a">
<li><a name="mpp_define_domains" id="mpp_define_domains"></a>
<h2>mpp_define_domains</h2>
<pre>
<b>call mpp_define_domains </b>( global_indices, ndivs, domain, & pelist, flags, halo, extent, maskmap )</pre>
<pre>
<b>call mpp_define_domains </b>( global_indices, layout, domain, pelist, & xflags, yflags, xhalo, yhalo, & xextent, yextent, maskmap, name )</pre>
<dl>
<dt><b>DESCRIPTION</b></dt>
<dd>There are two forms for the <code>mpp_define_domains</code>
call. The 2D version is generally to be used but is built by
repeated calls to the 1D version, also provided.</dd>
<dd><br>
<br></dd>
<dt><b>INPUT</b></dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left">
<code>global_indices   </code></td>
<td>Defines the global domain.<br>
   <span class="type">[integer,
dimension(2)]</span><br>
   <span class="type">[integer,
dimension(4)]</span></td>
</tr>
<tr>
<td valign="top" align="left">
<code>ndivs   </code></td>
<td>Is the number of domain divisions required.<br>
   <span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left">
<code>pelist   </code></td>
<td>List of PEs to which the domains are to be assigned.<br>
   <span class="type">[integer,
dimension(0:)]</span><br>
   <span class="type">[integer,
dimension(0:)]</span></td>
</tr>
<tr>
<td valign="top" align="left">
<code>flags   </code></td>
<td>An optional flag to pass additional information about the
desired domain topology. Useful flags in a 1D decomposition include
<code>GLOBAL_DATA_DOMAIN</code> and
<code>CYCLIC_GLOBAL_DOMAIN</code>. Flags are integers: multiple
flags may be added together. The flag values are public parameters
available by use association.<br>
   <span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left">
<code>halo   </code></td>
<td>Width of the halo.<br>
   <span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left">
<code>extent   </code></td>
<td>Normally <code>mpp_define_domains</code> attempts an even
division of the global domain across <code>ndivs</code> domains.
The <code>extent</code> array can be used by the user to pass a
custom domain division. The <code>extent</code> array has
<code>ndivs</code> elements and holds the compute domain widths,
which should add up to cover the global domain exactly.<br>
   <span class="type">[integer,
dimension(0:)]</span></td>
</tr>
<tr>
<td valign="top" align="left">
<code>maskmap   </code></td>
<td>Some divisions may be masked (<code>maskmap=.FALSE.</code>) to
exclude them from the computation (e.g for ocean model domains that
are all land). The <code>maskmap</code> array is dimensioned
<code>ndivs</code> and contains <code>.TRUE.</code> values for any
domain that must be <i>included</i> in the computation (default
all). The <code>pelist</code> array length should match the number
of domains included in the computation.<br>
   <span class="type">[logical,
dimension(0:)]</span><br>
   <span class="type">[logical,
dimension(:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left">
<code>layout   </code></td>
<td><br>
   <span class="type">[integer,
dimension(2)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><code>xflags,
yflags   </code></td>
<td><br>
   <span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><code>xhalo,
yhalo   </code></td>
<td><br>
   <span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><code>xextent,
yextent   </code></td>
<td><br>
   <span class="type">[integer,
dimension(0:)]</span></td>
</tr>
<tr>
<td valign="top" align="left">
<code>name   </code></td>
<td><br>
   <span class="type">[character(len=*)]</span></td>
</tr>
</table>
</dd>
<dd><br></dd>
<dt><b>INPUT/OUTPUT</b></dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left">
<code>domain   </code></td>
<td>Holds the resulting domain decomposition.<br>
   <span class="type">[type(domain1D)]</span><br>
   <span class="type">[type(domain2D)]</span></td>
</tr>
</table>
</dd>
<dd><br></dd>
<dt><b>NOTE</b></dt>
<dd>For example:<br>
<br>
<pre>    call mpp_define_domains( (/1,100/), 10, domain, &
         flags=GLOBAL_DATA_DOMAIN+CYCLIC_GLOBAL_DOMAIN, halo=2 )</pre>
defines 10 compute domains spanning the range [1,100] of the global
domain. The compute domains are non-overlapping blocks of 10. All
the data domains are global, and with a halo of 2 span the range
[-1:102]. And since the global domain has been declared to be
cyclic, <code>domain(9)%next =&gt; domain(0)</code> and
<code>domain(0)%prev =&gt; domain(9)</code>. A field is allocated
on the data domain, and computations proceed on the compute domain.
A call to <a href="#mpp_update_domains"></a> would fill in the
values in the halo region:
<pre>
call mpp_get_data_domain( domain, isd, ied ) !returns -1 and 102
    call mpp_get_compute_domain( domain, is, ie ) !returns (1,10) on PE 0 ...
    allocate( a(isd:ied) )
    do i = is,ie
       a(i) = &lt;perform computations&gt;
    end do
    call mpp_update_domains( a, domain )</pre>
The call to <code>mpp_update_domains</code> fills in the regions
outside the compute domain. Since the global domain is cyclic, the
values at <code>i=(-1,0)</code> are the same as at
<code>i=(99,100)</code>; and <code>i=(101,102)</code> are the same
as <code>i=(1,2)</code>.<br>
<br>
The 2D version is just an extension of this syntax to two
dimensions.<br>
<br>
The 2D version of the above should generally be used in codes,
including 1D-decomposed ones, if there is a possibility of future
evolution toward 2D decomposition. The arguments are similar to the
1D case, except that now we have optional arguments
<code>flags</code>, <code>halo</code>, <code>extent</code> and
<code>maskmap</code> along two axes.<br>
<br>
<code>flags</code> can now take an additional possible value to
fold one or more edges. This is done by using flags
<code>FOLD_WEST_EDGE</code>, <code>FOLD_EAST_EDGE</code>,
<code>FOLD_SOUTH_EDGE</code> or <code>FOLD_NORTH_EDGE</code>. When
a fold exists (e.g cylindrical domain), vector fields reverse sign
upon crossing the fold. This parity reversal is performed only in
the vector version of <a href="#mpp_update_domains"></a>. In
addition, shift operations may need to be applied to vector fields
on staggered grids, also described in the vector interface to
<code>mpp_update_domains</code>.<br>
<br>
<code>name</code> is the name associated with the decomposition,
e.g <code>'Ocean model'</code>. If this argument is present,
<code>mpp_define_domains</code> will print the domain decomposition
generated to <code>stdlog</code>.<br>
<br>
Examples:<br>
<br>
<pre>
<code>
call mpp_define_domains( (/1,100,1,100/), (/2,2/), domain, xhalo=1 )
</code>
</pre>
will create the following domain layout:
<pre>
<code>  
               |---------|-----------|-----------|-------------|
               |domain(1)|domain(2)  |domain(3)  |domain(4)    |
|--------------|---------|-----------|-----------|-------------|
|Compute domain|1,50,1,50|51,100,1,50|1,50,51,100|51,100,51,100|
|--------------|---------|-----------|-----------|-------------|
|Data domain   |0,51,1,50|50,101,1,50|0,51,51,100|50,101,51,100|
|--------------|---------|-----------|-----------|-------------|
</code>
</pre>
<p>Again, we allocate arrays on the data domain, perform
computations on the compute domain, and call
<code>mpp_update_domains</code> to update the halo region.</p>
<p>If we wished to perfom a 1D decomposition along <code>Y</code>
on the same global domain, we could use:</p>
<pre>
call mpp_define_domains( (/1,100,1,100/), layout=(/4,1/), domain, xhalo=1 )</pre>
<p>This will create the following domain layout:</p>
<pre>
<code>
               |----------|-----------|-----------|------------|
               |domain(1) |domain(2)  |domain(3)  |domain(4)   |
|--------------|----------|-----------|-----------|------------|
|Compute domain|1,100,1,25|1,100,26,50|1,100,51,75|1,100,76,100|
|--------------|----------|-----------|-----------|------------|
|Data domain   |0,101,1,25|0,101,26,50|0,101,51,75|1,101,76,100|
|--------------|----------|-----------|-----------|------------|
</code>
</pre></dd>
<dd><br>
<br></dd>
</dl>
</li>
<li><a name="mpp_update_domains" id="mpp_update_domains"></a>
<h2>mpp_update_domains</h2>
<pre>
<b>call mpp_update_domains </b>( field, domain, flags )</pre>
<pre>
<b>call mpp_update_domains </b>( fieldx, fieldy, domain, flags, gridtype )</pre>
<dl>
<dt><b>DESCRIPTION</b></dt>
<dd><code>mpp_update_domains</code> is used to perform a halo
update of a domain-decomposed array on each PE.
<code>MPP_TYPE_</code> can be of type <code>complex</code>,
<code>integer</code>, <code>logical</code> or <code>real</code>; of
4-byte or 8-byte kind; of rank up to 5. The vector version (with
two input data fields) is only present for <code>real</code>
types.<br>
<br>
For 2D domain updates, if there are halos present along both
<code>x</code> and <code>y</code>, we can choose to update one
only, by specifying <code>flags=XUPDATE</code> or
<code>flags=YUPDATE</code>. In addition, one-sided updates can be
performed by setting <code>flags</code> to any combination of
<code>WUPDATE</code>, <code>EUPDATE</code>, <code>SUPDATE</code>
and <code>NUPDATE</code>, to update the west, east, north and south
halos respectively. Any combination of halos may be used by adding
the requisite flags, e.g: <code>flags=XUPDATE+SUPDATE</code> or
<code>flags=EUPDATE+WUPDATE+SUPDATE</code> will update the east,
west and south halos.<br>
<br>
If a call to <code>mpp_update_domains</code> involves at least one
E-W halo and one N-S halo, the corners involved will also be
updated, i.e, in the example above, the SE and SW corners will be
updated.<br>
<br>
If <code>flags</code> is not supplied, that is equivalent to
<code>flags=XUPDATE+YUPDATE</code>.<br>
<br>
The vector version is passed the <code>x</code> and <code>y</code>
components of a vector field in tandem, and both are updated upon
return. They are passed together to treat parity issues on various
grids. For example, on a cubic sphere projection, the
<code>x</code> and <code>y</code> components may be interchanged
when passing from an equatorial cube face to a polar face. For
grids with folds, vector components change sign on crossing the
fold.<br>
<br>
Special treatment at boundaries such as folds is also required for
staggered grids. The following types of staggered grids are
recognized:<br>
<br>
1) <code>AGRID</code>: values are at grid centers.<br>
2) <code>BGRID_NE</code>: vector fields are at the NE vertex of a
grid cell, i.e: the array elements <code>u(i,j)</code> and
<code>v(i,j)</code> are actually at (i+½,j+½) with respect to the
grid centers.<br>
3) <code>BGRID_SW</code>: vector fields are at the SW vertex of a
grid cell, i.e: the array elements <code>u(i,j)</code> and
<code>v(i,j)</code> are actually at (i-½,j-½) with respect to the
grid centers.<br>
4) <code>CGRID_NE</code>: vector fields are at the N and E faces of
a grid cell, i.e: the array elements <code>u(i,j)</code> and
<code>v(i,j)</code> are actually at (i+½,j) and (i,j+½) with
respect to the grid centers.<br>
5) <code>CGRID_SW</code>: vector fields are at the S and W faces of
a grid cell, i.e: the array elements <code>u(i,j)</code> and
<code>v(i,j)</code> are actually at (i-½,j) and (i,j-½) with
respect to the grid centers.<br>
<br>
The gridtypes listed above are all available by use association as
integer parameters. The scalar version of
<code>mpp_update_domains</code> assumes that the values of a scalar
field are always at <code>AGRID</code> locations, and no special
boundary treatment is required. If vector fields are at staggered
locations, the optional argument <code>gridtype</code> must be
appropriately set for correct treatment at boundaries.<br>
<br>
It is safe to apply vector field updates to the appropriate arrays
irrespective of the domain topology: if the topology requires no
special treatment of vector fields, specifying
<code>gridtype</code> will do no harm.<br>
<br>
<code>mpp_update_domains</code> internally buffers the date being
sent and received into single messages for efficiency. A turnable
internal buffer area in memory is provided for this purpose by
<code>mpp_domains_mod</code>. The size of this buffer area can be
set by the user by calling mpp_domains_set_stack_size in <a href=
"./mpp_domains.html">mpp_domains</a>.</dd>
<dd><br>
<br></dd>
</dl>
</li>
<li><a name="mpp_redistribute" id="mpp_redistribute"></a>
<h2>mpp_redistribute</h2>
<pre>
<b>call mpp_redistribute </b>( domain_in, field_in, domain_out, field_out )</pre>
<dl>
<dt><b>DESCRIPTION</b></dt>
<dd><code>mpp_redistribute</code> is used to reorganize a
distributed array. <code>MPP_TYPE_</code> can be of type
<code>integer</code>, <code>complex</code>, or <code>real</code>;
of 4-byte or 8-byte kind; of rank up to 5.</dd>
<dd><br>
<br></dd>
<dt><b>INPUT</b></dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left">
<code>field_in   </code></td>
<td><code>field_in</code> is dimensioned on the data domain of
<code>domain_in</code>.</td>
</tr>
</table>
</dd>
<dd><br></dd>
<dt><b>OUTPUT</b></dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left">
<code>field_out   </code></td>
<td><code>field_out</code> on the data domain of
<code>domain_out</code>.</td>
</tr>
</table>
</dd>
<dd><br></dd>
</dl>
</li>
<li><a name="mpp_global_field" id="mpp_global_field"></a>
<h2>mpp_global_field</h2>
<pre>
<b>call mpp_global_field </b>( domain, local, global, flags )</pre>
<dl>
<dt><b>DESCRIPTION</b></dt>
<dd><code>mpp_global_field</code> is used to get an entire
domain-decomposed array on each PE. <code>MPP_TYPE_</code> can be
of type <code>complex</code>, <code>integer</code>,
<code>logical</code> or <code>real</code>; of 4-byte or 8-byte
kind; of rank up to 5.<br>
<br>
All PEs in a domain decomposition must call
<code>mpp_global_field</code>, and each will have a complete global
field at the end. Please note that a global array of rank 3 or
higher could occupy a lot of memory.</dd>
<dd><br>
<br></dd>
<dt><b>INPUT</b></dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left">
<code>domain   </code></td>
<td></td>
</tr>
<tr>
<td valign="top" align="left">
<code>local   </code></td>
<td><code>local</code> is dimensioned on either the compute domain
or the data domain of <code>domain</code>.</td>
</tr>
<tr>
<td valign="top" align="left">
<code>flags   </code></td>
<td><code>flags</code> can be given the value <code>XONLY</code> or
<code>YONLY</code>, to specify a globalization on one axis
only.</td>
</tr>
</table>
</dd>
<dd><br></dd>
<dt><b>OUTPUT</b></dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left">
<code>global   </code></td>
<td><code>global</code> is dimensioned on the corresponding global
domain.</td>
</tr>
</table>
</dd>
<dd><br></dd>
</dl>
</li>
<li><a name="mpp_global_max" id="mpp_global_max"></a>
<h2>mpp_global_max</h2>
<pre> 
<b>mpp_global_max</b> ( domain, field, locus )</pre>
<dl>
<dt><b>DESCRIPTION</b></dt>
<dd><code>mpp_global_max</code> is used to get the maximum value of
a domain-decomposed array on each PE. <code>MPP_TYPE_</code> can be
of type <code>integer</code> or <code>real</code>; of 4-byte or
8-byte kind; of rank up to 5. The dimension of <code>locus</code>
must equal the rank of <code>field</code>.<br>
<br>
All PEs in a domain decomposition must call
<code>mpp_global_max</code>, and each will have the result upon
exit.<br>
<br>
The function <code>mpp_global_min</code>, with an identical syntax.
is also available.</dd>
<dd><br>
<br></dd>
<dt><b>INPUT</b></dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left">
<code>domain   </code></td>
<td></td>
</tr>
<tr>
<td valign="top" align="left">
<code>field   </code></td>
<td><code>field</code> is dimensioned on either the compute domain
or the data domain of <code>domain</code>.</td>
</tr>
</table>
</dd>
<dd><br></dd>
<dt><b>OUTPUT</b></dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left">
<code>locus   </code></td>
<td><code>locus</code>, if present, can be used to retrieve the
location of the maximum (as in the <code>MAXLOC</code> intrinsic of
f90).</td>
</tr>
</table>
</dd>
<dd><br></dd>
</dl>
</li>
<li><a name="mpp_global_sum" id="mpp_global_sum"></a>
<h2>mpp_global_sum</h2>
<pre>
<b>call mpp_global_sum </b>( domain, field, flags )</pre>
<dl>
<dt><b>DESCRIPTION</b></dt>
<dd><code>mpp_global_sum</code> is used to get the sum of a
domain-decomposed array on each PE. <code>MPP_TYPE_</code> can be
of type <code>integer</code>, <code>complex</code>, or
<code>real</code>; of 4-byte or 8-byte kind; of rank up to 5.</dd>
<dd><br>
<br></dd>
<dt><b>INPUT</b></dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left">
<code>domain   </code></td>
<td></td>
</tr>
<tr>
<td valign="top" align="left">
<code>field   </code></td>
<td><code>field</code> is dimensioned on either the compute domain
or the data domain of <code>domain</code>.</td>
</tr>
<tr>
<td valign="top" align="left">
<code>flags   </code></td>
<td><code>flags</code>, if present, must have the value
<code>BITWISE_EXACT_SUM</code>. This produces a sum that is
guaranteed to produce the identical result irrespective of how the
domain is decomposed. This method does the sum first along the
ranks beyond 2, and then calls <a href="#mpp_global_field"></a> to
produce a global 2D array which is then summed. The default method,
which is considerably faster, does a local sum followed by mpp_sum in <a href=
"./mpp.html">mpp</a> across
the domain decomposition.</td>
</tr>
</table>
</dd>
<dd><br></dd>
<dt><b>NOTE</b></dt>
<dd>All PEs in a domain decomposition must call
<code>mpp_global_sum</code>, and each will have the result upon
exit.</dd>
<dd><br>
<br></dd>
</dl>
</li>
<li><a name="operator" id="operator"></a>
<h2>operator</h2>
<dl>
<dt><b>DESCRIPTION</b></dt>
<dd>The module provides public operators to check for
equality/inequality of domaintypes, e.g:<br>
<br>
<pre>    type(domain1D) :: a, b
    type(domain2D) :: c, d
    ...
    if( a.NE.b )then
        ...
    end if
    if( c==d )then
        ...
    end if</pre>
Domains are considered equal if and only if the start and end
indices of each of their component global, data and compute domains
are equal.</dd>
<dd><br>
<br></dd>
</dl>
</li>
<li><a name="mpp_get_compute_domain" id=
"mpp_get_compute_domain"></a>
<h2>mpp_get_compute_domain</h2>
<pre>
<b>call mpp_get_compute_domain </b>
</pre>
<dl>
<dt><b>DESCRIPTION</b></dt>
<dd>The domain is a derived type with private elements. These
routines retrieve the axis specifications associated with the
compute domains The 2D version of these is a simple extension of
1D.</dd>
<dd><br>
<br></dd>
</dl>
</li>
<li><a name="mpp_get_compute_domains" id=
"mpp_get_compute_domains"></a>
<h2>mpp_get_compute_domains</h2>
<pre>
<b>call mpp_get_compute_domains </b>( domain, xbegin, xend, xsize, & ybegin, yend, ysize )</pre>
<dl>
<dt><b>DESCRIPTION</b></dt>
<dd>Retrieve the entire array of compute domain extents associated
with a decomposition.</dd>
<dd><br>
<br></dd>
<dt><b>INPUT</b></dt>
<dd>
<code>domain   </code>
</dd>
<dd><br></dd>
<dt><b>OUTPUT</b></dt>
<dd>
<code>xbegin,ybegin   </code>, 
<code>xend,yend   </code>, 
<code>xsize,ysize   </code>
</dd>
<dd><br></dd>
</dl>
</li>
<li><a name="mpp_get_data_domain" id="mpp_get_data_domain"></a>
<h2>mpp_get_data_domain</h2>
<pre>
<b>call mpp_get_data_domain </b>
</pre>
<dl>
<dt><b>DESCRIPTION</b></dt>
<dd>The domain is a derived type with private elements. These
routines retrieve the axis specifications associated with the data
domains. The 2D version of these is a simple extension of 1D.</dd>
<dd><br>
<br></dd>
</dl>
</li>
<li><a name="mpp_get_global_domain" id="mpp_get_global_domain"></a>
<h2>mpp_get_global_domain</h2>
<pre>
<b>call mpp_get_global_domain </b>
</pre>
<dl>
<dt><b>DESCRIPTION</b></dt>
<dd>The domain is a derived type with private elements. These
routines retrieve the axis specifications associated with the
global domains. The 2D version of these is a simple extension of
1D.</dd>
<dd><br>
<br></dd>
</dl>
</li>
<li><a name="mpp_define_layout" id="mpp_define_layout"></a>
<h2>mpp_define_layout</h2>
<pre>
<b>call mpp_define_layout </b>( global_indices, ndivs, layout )</pre>
<dl>
<dt><b>DESCRIPTION</b></dt>
<dd>Given a global 2D domain and the number of divisions in the
decomposition (<code>ndivs</code>: usually the PE count unless some
domains are masked) this calls returns a 2D domain layout.<br>
<br>
By default, <code>mpp_define_layout</code> will attempt to divide
the 2D index space into domains that maintain the aspect ratio of
the global domain. If this cannot be done, the algorithm favours
domains that are longer in <code>x</code> than <code>y</code>, a
preference that could improve vector performance.</dd>
<dd><br>
<br></dd>
<dt><b>INPUT</b></dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left">
<code>global_indices   </code></td>
<td><br>
   <span class="type">[integer,
dimension(4)]</span></td>
</tr>
<tr>
<td valign="top" align="left">
<code>ndivs   </code></td>
<td><br>
   <span class="type">[integer]</span></td>
</tr>
</table>
</dd>
<dd><br></dd>
<dt><b>OUTPUT</b></dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left">
<code>layout   </code></td>
<td><br>
   <span class="type">[integer,
dimension(2)]</span></td>
</tr>
</table>
</dd>
<dd><br></dd>
</dl>
</li>
<li><a name="mpp_get_pelist" id="mpp_get_pelist"></a>
<h2>mpp_get_pelist</h2>
<dl>
<dt><b>DESCRIPTION</b></dt>
<dd>The 1D version of this call returns an array of the PEs
assigned to this 1D domain decomposition. In addition the optional
argument <code>pos</code> may be used to retrieve the 0-based
position of the domain local to the calling PE, i.e
<code>domain%list(pos)%pe</code> is the local PE, as returned by mpp_pe in
<a href=
"./mpp.html">mpp</a>. The 2D
version of this call is identical to 1D version.</dd>
<dd><br>
<br></dd>
<dt><b>INPUT</b></dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left">
<code>domain   </code></td>
<td><br>
   <span class="type">[type(domain1D)]</span><br>
   <span class="type">[type(domain2D)]</span></td>
</tr>
</table>
</dd>
<dd><br></dd>
<dt><b>OUTPUT</b></dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left">
<code>pelist   </code></td>
<td><br>
   <span class="type">[integer,
dimension(:)]</span><br>
   <span class="type">[integer,
dimension(:)]</span></td>
</tr>
<tr>
<td valign="top" align="left">
<code>pos   </code></td>
<td><br>
   <span class="type">[integer]</span><br>
   <span class="type">[integer]</span></td>
</tr>
</table>
</dd>
<dd><br></dd>
</dl>
</li>
<li><a name="mpp_get_layout" id="mpp_get_layout"></a>
<h2>mpp_get_layout</h2>
<pre>
<b>call mpp_get_layout </b>( domain, layout )</pre>
<dl>
<dt><b>DESCRIPTION</b></dt>
<dd>The 1D version of this call returns the number of divisions
that was assigned to this decomposition axis. The 2D version of
this call returns an array of dimension 2 holding the results on
two axes.</dd>
<dd><br>
<br></dd>
<dt><b>INPUT</b></dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left">
<code>domain   </code></td>
<td><br>
   <span class="type">[type(domain1D)]</span><br>
   <span class="type">[type(domain2D)]</span></td>
</tr>
</table>
</dd>
<dd><br></dd>
<dt><b>OUTPUT</b></dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left">
<code>layout   </code></td>
<td><br>
   <span class="type">[integer]</span><br>
   <span class="type">[integer,
dimension(2)]</span></td>
</tr>
</table>
</dd>
<dd><br></dd>
</dl>
</li>
<li><a name="mpp_domains_init" id="mpp_domains_init"></a>
<h2>mpp_domains_init</h2>
<pre>
<b>call mpp_domains_init </b>(flags)</pre>
<dl>
<dt><b>DESCRIPTION</b></dt>
<dd>Called to initialize the <code>mpp_domains_mod</code>
package.<br>
<br>
<code>flags</code> can be set to <code>MPP_VERBOSE</code> to have
<code>mpp_domains_mod</code> keep you informed of what it's up to.
<code>MPP_DEBUG</code> returns even more information for
debugging.<br>
<br>
<code>mpp_domains_init</code> will call <code>mpp_init</code>, to
make sure <a href=
"./mpp.html"></a> is
initialized. (Repeated calls to <code>mpp_init</code> do no harm,
so don't worry if you already called it).</dd>
<dd><br>
<br></dd>
<dt><b>INPUT</b></dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left">
<code>flags   </code></td>
<td><br>
   <span class="type">[integer]</span></td>
</tr>
</table>
</dd>
<dd><br></dd>
</dl>
</li>
<li><a name="mpp_domains_set_stack_size" id=
"mpp_domains_set_stack_size"></a>
<h2>mpp_domains_set_stack_size</h2>
<pre>
<b>call mpp_domains_set_stack_size </b>(n)</pre>
<dl>
<dt><b>DESCRIPTION</b></dt>
<dd>This sets the size of an array that is used for internal
storage by <code>mpp_domains</code>. This array is used, for
instance, to buffer the data sent and received in halo updates.<br>
<br>
This call has implied global synchronization. It should be placed
somewhere where all PEs can call it.</dd>
<dd><br>
<br></dd>
<dt><b>INPUT</b></dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><code>n   </code></td>
<td><br>
   <span class="type">[integer]</span></td>
</tr>
</table>
</dd>
<dd><br></dd>
</dl>
</li>
<li><a name="mpp_domains_exit" id="mpp_domains_exit"></a>
<h2>mpp_domains_exit</h2>
<pre>
<b>call mpp_domains_exit </b>()</pre>
<dl>
<dt><b>DESCRIPTION</b></dt>
<dd>Serves no particular purpose, but is provided should you
require to re-initialize <code>mpp_domains_mod</code>, for some odd
reason.</dd>
<dd><br>
<br></dd>
</dl>
</li>
<li><a name="mpp_get_domain_components" id=
"mpp_get_domain_components"></a>
<h2>mpp_get_domain_components</h2>
<pre>
<b>call mpp_get_domain_components </b>( domain, x, y )</pre>
<dl>
<dt><b>DESCRIPTION</b></dt>
<dd>It is sometime necessary to have direct recourse to the
domain1D types that compose a domain2D object. This call retrieves
them.</dd>
<dd><br>
<br></dd>
<dt><b>INPUT</b></dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left">
<code>domain   </code></td>
<td><br>
   <span class="type">[type(domain2D)]</span></td>
</tr>
</table>
</dd>
<dd><br></dd>
<dt><b>OUTPUT</b></dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left">
<code>x,y   </code></td>
<td><br>
   <span class="type">[type(domain1D)]</span></td>
</tr>
</table>
</dd>
<dd><br></dd>
</dl>
</li>
</ol>
<!-- END PUBLIC ROUTINES -->
<a name="NAMELIST" id="NAMELIST"></a> <!-- BEGIN NAMELIST -->
<!-- END NAMELIST --><a name="DIAGNOSTIC FIELDS"></a> 
<!-- BEGIN DIAGNOSTIC FIELDS -->
<!-- END DIAGNOSTIC FIELDS --><a name="DATA SETS"></a> 
<!-- BEGIN DATA SETS -->
<hr>
<h2>DATA SETS</h2>
<div>None.<br>
<br></div>
<!-- END DATA SETS -->
<a name="ERROR MESSAGES"></a> <!-- BEGIN ERROR MESSAGES -->
<hr>
<h2>ERROR MESSAGES</h2>
<div>None.<br>
<br></div>
<!-- END ERROR MESSAGES -->
<a name="REFERENCES" id="REFERENCES"></a>
<hr>
<h2>REFERENCES</h2>
<!-- BEGIN REFERENCES -->
<div>None.</div>
<br>
<!-- END REFERENCES -->
<a name="COMPILER SPECIFICS"></a>
<hr>
<h2>COMPILER SPECIFICS</h2>
<!-- BEGIN COMPILER SPECIFICS -->
<div>
<dl>
<dd>Any module or program unit using <code>mpp_domains_mod</code>
must contain the line
<pre>     use mpp_domains_mod</pre>
<code>mpp_domains_mod</code> <code>use</code>s <a href=
"./mpp.html"></a>, and therefore
is subject to the compiling and linking requirements of <a href=
"./mpp.html">mpp</a>.</dd>
</dl>
</div>
<br>
<!-- END COMPILER SPECIFICS -->
<a name="PRECOMPILER OPTIONS"></a>
<hr>
<h2>PRECOMPILER OPTIONS</h2>
<!-- BEGIN PRECOMPILER OPTIONS -->
<div>
<dl>
<dd><code>mpp_domains_mod</code> uses standard f90, and has no
special requirements. There are some OS-dependent pre-processor
directives that you might need to modify on non-SGI/Cray systems
and compilers. The portability, as described in <a href="./mpp.html">mpp</a>
obviously is a constraint, since this module is built on top of it. Contact me,
Balaji, SGI/GFDL, with questions.</dd>
</dl>
</div>
<br>
<!-- END PRECOMPILER OPTIONS -->
<a name="LOADER OPTIONS"></a>
<hr>
<h2>LOADER OPTIONS</h2>
<!-- BEGIN LOADER -->
<p>The source consists of the main source file and also requires
the following include files: GFDL users can check it out of the
main CVS repository as part of the CVS module. The current public
tag is . External users can download the latest package . Public
access to the GFDL CVS repository will soon be made available.</p>
<!-- END LOADER OPTIONS -->
<a name="TEST PROGRAM"></a>
<hr>
<h2>TEST PROGRAM</h2>
<!-- BEGIN TEST PROGRAM -->
<div>None.<br></div>
<br>
<!-- END TEST PROGRAM -->
<a name="KNOWN BUGS"></a>
<hr>
<h2>KNOWN BUGS</h2>
<!-- BEGIN KNOWN BUGS -->
<div>None.</div>
<br>
<!-- END KNOWN BUGS -->
<a name="NOTES" id="NOTES"></a>
<hr>
<h2>NOTES</h2>
<!-- BEGIN NOTES -->
<div>None.<br></div>
<br>
<!-- END NOTES -->
<a name="FUTURE PLANS"></a>
<hr>
<h2>FUTURE PLANS</h2>
<!-- BEGIN FUTURE PLANS -->
<div>None.</div>
<br>
<!-- END FUTURE PLANS -->
<hr>
<div align="right"><a href="#TOP"><font size=
"-1">top</font></a></div>
</body>
</html>
