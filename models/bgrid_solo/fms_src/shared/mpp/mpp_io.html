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
<title>Module mpp_io_mod</title>
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
<h2>module mpp_io_mod</h2>
<a name="HEADER"></a>
<!-- BEGIN HEADER -->
<div>
<b>Contact:&nbsp;</b><a href="mailto:vb@gfdl.noaa.gov">   V. Balaji </a>
<br>
<b>Reviewers:&nbsp;</b>
<br>
<b>Change History:&nbsp;</b><a href="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/">WebCVS Log</a>
<br>
<b>RCS Log:&nbsp;</b><a href="http://www.gfdl.noaa.gov/~vb/changes_models/bgrid_solo/fms_src/shared/mpp/mpp_io.html">RCS Log</a>
<br>
<b>Last Modified:&nbsp;</b>2002/07/19 18:10:07<br>
<br>
</div>
<!-- END HEADER -->
<a name="OVERVIEW"></a>
<hr>
<h4>OVERVIEW</h4>
<!-- BEGIN OVERVIEW -->
<p class="text"> 
<tt>mpp_io_mod</tt>, is a set of simple calls for parallel I/O on
   distributed systems. It is geared toward the writing of data in netCDF
   format. It requires the modules <a href="models/bgrid_solo/fms_src/shared/mpp/mpp_domains.html"></a>   and <a href="models/bgrid_solo/fms_src/shared/mpp/mpp.html"></a>, upon which it is built. </p>
<!-- END OVERVIEW -->
<a name="DESCRIPTION"></a>
<!-- BEGIN DESCRIPTION -->
<div>   In massively parallel environments, an often difficult problem is
   the reading and writing of data to files on disk. MPI-IO and MPI-2 IO
   are moving toward providing this capability, but are currently not
   widely implemented. Further, it is a rather abstruse
   API. <tt>mpp_io_mod</tt>   is an attempt at a simple API encompassing a
   certain variety of the I/O tasks that will be required. It does not
   attempt to be an all-encompassing standard such as MPI, however, it
   can be implemented in MPI if so desired. It is equally simple to add
   parallel I/O capability to <tt>mpp_io_mod</tt>   based on vendor-specific
   APIs while providing a layer of insulation for user codes.
   <br>
<br>
   The <tt>mpp_io_mod</tt>   parallel I/O API built on top of the <a href="models/bgrid_solo/fms_src/shared/mpp/mpp_domains.html"></a>   and <a href="models/bgrid_solo/fms_src/shared/mpp/mpp.html"></a>   API for domain decomposition and
   message passing. Features of <tt>mpp_io_mod</tt>   include:
   <br>
<br>
   1) Simple, minimal API, with free access to underlying API for more
   complicated stuff.<br>   2) Self-describing files: comprehensive header information
   (metadata) in the file itself.<br>   3) Strong focus on performance of parallel write: the climate models
   for which it is designed typically read a minimal amount of data
   (typically at the beginning of the run); but on the other hand, tend
   to write copious amounts of data during the run. An interface for
   reading is also supplied, but its performance has not yet been optimized.<br>   4) Integrated netCDF capability: <a href="http://www.unidata.ucar.edu/packages/netcdf/">netCDF</a>   is a
   data format widely used in the climate/weather modeling
   community. netCDF is considered the principal medium of data storage
   for <tt>mpp_io_mod</tt>. But I provide a raw unformatted
   fortran I/O capability in case netCDF is not an option, either due to
   unavailability, inappropriateness, or poor performance.<br>   5) May require off-line post-processing: a tool for this purpose, <tt>mppnccombine</tt>, is available. GFDL users may use <tt>~hnv/pub/mppnccombine</tt>. Outside users may obtain the
   source <a href="ftp://ftp.gfdl.gov/perm/hnv/mpp/mppnccombine.c">here</a>.  It
   can be compiled on any C compiler and linked with the netCDF
   library. The program is free and is covered by the <a href="ftp://ftp.gfdl.gov/perm/hnv/mpp/LICENSE">GPL license</a>.
   <br>
<br>
   The internal representation of the data being written out is
   assumed be the default real type, which can be 4 or 8-byte. Time data
   is always written as 8-bytes to avoid overflow on climatic time scales
   in units of seconds.
   <br>
<br> 
<a href="modes"></a>
<h4>I/O modes in <tt>mpp_io_mod</tt>
</h4>   The I/O activity critical to performance in the models for which <tt>mpp_io_mod</tt>   is designed is typically the writing of large
   datasets on a model grid volume produced at intervals during
   a run. Consider a 3D grid volume, where model arrays are stored as <tt>(i,j,k)</tt>. The domain decomposition is typically along <tt>i</tt>   or <tt>j</tt>: thus to store data to disk as a global
   volume, the distributed chunks of data have to be seen as
   non-contiguous. If we attempt to have all PEs write this data into a
   single file, performance can be seriously compromised because of the
   data reordering that will be required. Possible options are to have
   one PE acquire all the data and write it out, or to have all the PEs
   write independent files, which are recombined offline. These three
   modes of operation are described in the <tt>mpp_io_mod</tt>   terminology
   in terms of two parameters, <i>threading</i>   and <i>fileset</i>,
   as follows:
   <br>
<br> 
<i>Single-threaded I/O:</i>   a single PE acquires all the data
   and writes it out.<br> 
<i>Multi-threaded, single-fileset I/O:</i>   many PEs write to a
   single file.<br> 
<i>Multi-threaded, multi-fileset I/O:</i>   many PEs write to
   independent files. This is also called <i>distributed I/O</i>.
   <br>
<br>
   The middle option is the most difficult to achieve performance. The
   choice of one of these modes is made when a file is opened for I/O, in <a href="#mpp_open"></a>.
   <br>
<br> 
<a href=""></a>
<h4>Metadata in <tt>mpp_io_mod</tt>
</h4>   A requirement of the design of <tt>mpp_io_mod</tt>   is that the file must
   be entirely self-describing: comprehensive header information
   describing its contents is present in the header of every file. The
   header information follows the model of netCDF. Variables in the file
   are divided into <i>axes</i>   and <i>fields</i>. An axis describes a
   co-ordinate variable, e.g <tt>x,y,z,t</tt>. A field consists of data in
   the space described by the axes. An axis is described in <tt>mpp_io_mod</tt>   using the defined type <tt>axistype</tt>:
   <br>
<br> 
<pre>   type, public :: axistype
      sequence
      character(len=128) :: name
      character(len=128) :: units
      character(len=256) :: longname
      character(len=8) :: cartesian
      integer :: len
      integer :: sense           !+/-1, depth or height?
      type(domain1D), pointer :: domain
      real, dimension(:), pointer :: data
      integer :: id, did
      integer :: type  ! external NetCDF type format for axis data
      integer :: natt
      type(atttype), pointer :: Att(:) ! axis attributes
   end type axistype</pre>   A field is described using the type <tt>fieldtype</tt>:
   <br>
<br> 
<pre>   type, public :: fieldtype
      sequence
      character(len=128) :: name
      character(len=128) :: units
      character(len=256) :: longname
      real :: min, max, missing, fill, scale, add
      integer :: pack
      type(axistype), dimension(:), pointer :: axes
      integer, dimension(:), pointer :: size
      integer :: time_axis_index
      integer :: id
      integer :: type ! external NetCDF format for field data
      integer :: natt, ndim
      type(atttype), pointer :: Att(:) ! field metadata
   end type fieldtype</pre>   An attribute (global, field or axis) is described using the <tt>atttype</tt>:
   <br>
<br> 
<pre>   type, public :: atttype
      sequence
      integer :: type, len
      character(len=128) :: name
      character(len=256)  :: catt
      real(FLOAT_KIND), pointer :: fatt(:)
   end type atttype</pre> 
<a href=""></a>This default set of field attributes corresponds
   closely to various conventions established for netCDF files. The <tt>pack</tt>   attribute of a field defines whether or not a
   field is to be packed on output. Allowed values of <tt>pack</tt>   are 1,2,4 and 8. The value of <tt>pack</tt>   is the number of variables written into 8
   bytes. In typical use, we write 4-byte reals to netCDF output; thus
   the default value of <tt>pack</tt>   is 2. For <tt>pack</tt>   = 4 or 8, packing uses a simple-minded linear
   scaling scheme using the <tt>scale</tt>   and <tt>add</tt>   attributes. There is thus likely to be a significant loss of dynamic
   range with packing. When a field is declared to be packed, the <tt>missing</tt>   and <tt>fill</tt>   attributes, if
   supplied, are packed also.
   <br>
<br>
   Please note that the pack values are the same even if the default
   real is 4 bytes, i.e <tt>PACK=1</tt>   still follows the definition
   above and writes out 8 bytes.
   <br>
<br>
   A set of <i>attributes</i>   for each variable is also available. The
   variable definitions and attribute information is written/read by calling <a href="#mpp_write_meta"></a>   or <a href="#mpp_read_meta"></a>. A typical calling
   sequence for writing data might be:
   <br>
<br> 
<pre>   ...
     type(domain2D), dimension(:), allocatable, target :: domain
     type(fieldtype) :: field
     type(axistype) :: x, y, z, t
   ...
     call mpp_define_domains( (/1,nx,1,ny/), domain )
     allocate( a(domain(pe)%x%data%start_index:domain(pe)%x%data%end_index, &amp;
                 domain(pe)%y%data%start_index:domain(pe)%y%data%end_index,nz) )
   ...
     call mpp_write_meta( unit, x, 'X', 'km', 'X distance', &amp;
          domain=domain(pe)%x, data=(/(float(i),i=1,nx)/) )
     call mpp_write_meta( unit, y, 'Y', 'km', 'Y distance', &amp;
          domain=domain(pe)%y, data=(/(float(i),i=1,ny)/) )
     call mpp_write_meta( unit, z, 'Z', 'km', 'Z distance', &amp;
          data=(/(float(i),i=1,nz)/) )
     call mpp_write_meta( unit, t, 'Time', 'second', 'Time' )
   
     call mpp_write_meta( unit, field, (/x,y,z,t/), 'a', '(m/s)', AAA', &amp;
          missing=-1e36 )
   ...
     call mpp_write( unit, x )
     call mpp_write( unit, y )
     call mpp_write( unit, z )
   ...</pre>   In this example, <tt>x</tt>   and <tt>y</tt>   have been
   declared as distributed axes, since a domain decomposition has been
   associated. <tt>z</tt>   and <tt>t</tt>   are undistributed
   axes. <tt>t</tt>   is known to be a <i>record</i>   axis (netCDF
   terminology) since we do not allocate the <tt>data</tt>   element
   of the <tt>axistype</tt>. <i>Only one record axis may be
   associated with a file.</i>   The call to <a href="#mpp_write_meta"></a>   initializes
   the axes, and associates a unique variable ID with each axis. The call
   to <tt>mpp_write_meta</tt>   with argument <tt>field</tt>   declared <tt>field</tt>   to be a 4D variable that is a function
   of <tt>(x,y,z,t)</tt>, and a unique variable ID is associated
   with it. A 3D field will be written at each call to <tt>mpp_write(field)</tt>.
   <br>
<br>
   The data to any variable, including axes, is written by <tt>mpp_write</tt>.
   <br>
<br>
   Any additional attributes of variables can be added through
   subsequent <tt>mpp_write_meta</tt>   calls, using the variable ID as a
   handle. <i>Global</i>   attributes, associated with the dataset as a
   whole, can also be written thus. See the <a href="#mpp_write_meta"></a>   call syntax below
   for further details.
   <br>
<br>
   You cannot interleave calls to <tt>mpp_write</tt>   and <tt>mpp_write_meta</tt>: the first call to <tt>mpp_write</tt>   implies that metadata specification is
   complete.
   <br>
<br>
   A typical calling sequence for reading data might be:
   <br>
<br> 
<pre>   ...
     integer :: unit, natt, nvar, ntime
     type(domain2D), dimension(:), allocatable, target :: domain
     type(fieldtype), allocatable, dimension(:) :: fields
     type(atttype), allocatable, dimension(:) :: global_atts
     real, allocatable, dimension(:) :: times
   ...
     call mpp_define_domains( (/1,nx,1,ny/), domain )
   
     call mpp_read_meta(unit)
     call mpp_get_info(unit,natt,nvar,ntime)
     allocate(global_atts(natt))
     call mpp_get_atts(unit,global_atts)
     allocate(fields(nvar))
     call mpp_get_vars(unit, fields)
     allocate(times(ntime))
     call mpp_get_times(unit, times)
   
     allocate( a(domain(pe)%x%data%start_index:domain(pe)%x%data%end_index, &amp;
                 domain(pe)%y%data%start_index:domain(pe)%y%data%end_index,nz) )
   ...
     do i=1, nvar
       if (fields(i)%name == 'a')  call mpp_read(unit,fields(i),domain(pe), a,
                                                 tindex)
     enddo
   ...</pre>   In this example, the data are distributed as in the previous
   example. The call to <a href="#mpp_read_meta"></a>   initializes
   all of the metadata associated with the file, including global
   attributes, variable attributes and non-record dimension data. The
   call to <tt>mpp_get_info</tt>   returns the number of global
   attributes (<tt>natt</tt>), variables (<tt>nvar</tt>) and
   time levels (<tt>ntime</tt>) associated with the file
   identified by a unique ID (<tt>unit</tt>). <tt>mpp_get_atts</tt>   returns all global attributes for
   the file in the derived type <tt>atttype(natt)</tt>. <tt>mpp_get_vars</tt>   returns variable types
   (<tt>fieldtype(nvar)</tt>).  Since the record dimension data are not allocated for calls to <a href="#mpp_write"></a>, a separate call to <tt>mpp_get_times</tt>   is required to access record dimension data.  Subsequent calls to <tt>mpp_read</tt>   return the field data arrays corresponding to
   the fieldtype.  The <tt>domain</tt>   type is an optional
   argument.  If <tt>domain</tt>   is omitted, the incoming field
   array should be dimensioned for the global domain, otherwise, the
   field data is assigned to the computational domain of a local array.
   <br>
<br> 
<i>Multi-fileset</i>   reads are not supported with <tt>mpp_read</tt>. </div>
<br>
<!-- END DESCRIPTION -->
<a name="OTHER MODULES USED"></a>
<hr>
<h4>OTHER MODULES USED</h4>
<!-- BEGIN OTHER MODULES USED -->
<div>
<pre>        mpp_mod<br>mpp_domains_mod</pre>
</div>
<!-- END OTHER MODULES USED -->
<!-- BEGIN PUBLIC INTERFACE -->
<a name="PUBLIC INTERFACE"></a>
<hr>
<h4>PUBLIC INTERFACE</h4>
<div>
<pre>use mpp_io_mod [, only:  mpp_write_meta,<br>                         mpp_write,<br>                         mpp_read,<br>                         mpp_get_atts,<br>                         mpp_io_init,<br>                         mpp_io_exit,<br>                         mpp_open,<br>                         mpp_close,<br>                         mpp_read_meta,<br>                         mpp_get_info,<br>                         mpp_get_times,<br>                         mpp_flush,<br>                         mpp_get_ncid ]</pre>
<dl>
<dt>
<a href="#mpp_write_meta">mpp_write_meta</a>:</dt>
<dd>   Write metadata. </dd>
<dt>
<a href="#mpp_write">mpp_write</a>:</dt>
<dd>   Write to an open file. </dd>
<dt>
<a href="#mpp_read">mpp_read</a>:</dt>
<dd>   Read from an open file. </dd>
<dt>
<a href="#mpp_get_atts">mpp_get_atts</a>:</dt>
<dd>   Get file global metdata. </dd>
<dt>
<a href="#mpp_io_init">mpp_io_init</a>:</dt>
<dd>   Initialize <tt>mpp_io_mod</tt>. </dd>
<dt>
<a href="#mpp_io_exit">mpp_io_exit</a>:</dt>
<dd>   Exit <tt>mpp_io_mod</tt>. </dd>
<dt>
<a href="#mpp_open">mpp_open</a>:</dt>
<dd>   Open a file for parallel I/O. </dd>
<dt>
<a href="#mpp_close">mpp_close</a>:</dt>
<dd>   Close an open file. </dd>
<dt>
<a href="#mpp_read_meta">mpp_read_meta</a>:</dt>
<dd>   Read metadata. </dd>
<dt>
<a href="#mpp_get_info">mpp_get_info</a>:</dt>
<dd>   Get some general information about a file. </dd>
<dt>
<a href="#mpp_get_times">mpp_get_times</a>:</dt>
<dd>   Get file time data. </dd>
<dt>
<a href="#mpp_flush">mpp_flush</a>:</dt>
<dd>   Flush I/O buffers to disk. </dd>
<dt>
<a href="#mpp_get_ncid">mpp_get_ncid</a>:</dt>
<dd>   Get netCDF ID of an open file. </dd>
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
<a name="mpp_write_meta"></a>
<h4>mpp_write_meta</h4>
<pre>
<b>call mpp_write_meta </b>( unit, axis, name, units, longname, cartesian, sense, domain, data )</pre>
<pre>
<b>call mpp_write_meta </b>( unit, field, axes, name, units, longname, min, max, missing, fill, scale, add, pack )</pre>
<pre>
<b>call mpp_write_meta </b>( unit, id, name, rval=rval, pack=pack )</pre>
<pre>
<b>call mpp_write_meta </b>( unit, id, name, ival=ival )</pre>
<pre>
<b>call mpp_write_meta </b>( unit, id, name, cval=cval )</pre>
<pre>
<b>call mpp_write_meta </b>( unit, name, rval=rval, pack=pack )</pre>
<pre>
<b>call mpp_write_meta </b>( unit, name, ival=ival )</pre>
<pre>
<b>call mpp_write_meta </b>( unit, name, cval=cval )</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   This routine is used to write the <a href="#metadata">metadata</a>   describing the contents of a file being written. Each file can contain
   any number of fields, which are functions of 0-3 space axes and 0-1
   time axes. (Only one time axis can be defined per file). The basic
   metadata defined <a href="#metadata">above</a>   for <tt>axistype</tt>   and <tt>fieldtype</tt>   are written in the first two forms of the call
   shown below. These calls will associate a unique variable ID with each
   variable (axis or field). These can be used to attach any other real,
   integer or character attribute to a variable. The last form is used to
   define a <i>global</i>   real, integer or character attribute that
   applies to the dataset as a whole. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>unit&nbsp;&nbsp;&nbsp;</tt></td><td></td>
</tr>
<tr>
<td valign="top" align="left"><tt>name&nbsp;&nbsp;&nbsp;</tt></td><td></td>
</tr>
<tr>
<td valign="top" align="left"><tt>units&nbsp;&nbsp;&nbsp;</tt></td><td></td>
</tr>
<tr>
<td valign="top" align="left"><tt>longname&nbsp;&nbsp;&nbsp;</tt></td><td></td>
</tr>
<tr>
<td valign="top" align="left"><tt>cartesian&nbsp;&nbsp;&nbsp;</tt></td><td></td>
</tr>
<tr>
<td valign="top" align="left"><tt>sense&nbsp;&nbsp;&nbsp;</tt></td><td></td>
</tr>
<tr>
<td valign="top" align="left"><tt>domain&nbsp;&nbsp;&nbsp;</tt></td><td></td>
</tr>
<tr>
<td valign="top" align="left"><tt>data&nbsp;&nbsp;&nbsp;</tt></td><td></td>
</tr>
<tr>
<td valign="top" align="left"><tt>min, max&nbsp;&nbsp;&nbsp;</tt></td><td></td>
</tr>
<tr>
<td valign="top" align="left"><tt>missing&nbsp;&nbsp;&nbsp;</tt></td><td></td>
</tr>
<tr>
<td valign="top" align="left"><tt>fill&nbsp;&nbsp;&nbsp;</tt></td><td></td>
</tr>
<tr>
<td valign="top" align="left"><tt>scale&nbsp;&nbsp;&nbsp;</tt></td><td></td>
</tr>
<tr>
<td valign="top" align="left"><tt>add&nbsp;&nbsp;&nbsp;</tt></td><td></td>
</tr>
<tr>
<td valign="top" align="left"><tt>pack&nbsp;&nbsp;&nbsp;</tt></td><td></td>
</tr>
<tr>
<td valign="top" align="left"><tt>id&nbsp;&nbsp;&nbsp;</tt></td><td></td>
</tr>
<tr>
<td valign="top" align="left"><tt>cval&nbsp;&nbsp;&nbsp;</tt></td><td></td>
</tr>
<tr>
<td valign="top" align="left"><tt>ival&nbsp;&nbsp;&nbsp;</tt></td><td></td>
</tr>
<tr>
<td valign="top" align="left"><tt>rval&nbsp;&nbsp;&nbsp;</tt></td><td></td>
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
<td valign="top" align="left"><tt>axis&nbsp;&nbsp;&nbsp;</tt></td><td></td>
</tr>
<tr>
<td valign="top" align="left"><tt>field&nbsp;&nbsp;&nbsp;</tt></td><td></td>
</tr>
</table>
</dd>
<br>
<dt>
<b>NOTE</b>
</dt>
<dd>   The first form defines a time or space axis. Metadata corresponding to the type
   above are written to the file on &lt;unit&gt;. A unique ID for subsequen
   references to this axis is returned in axis%id. If the &lt;domain&gt;
   element is present, this is recognized as a distributed data axis
   and domain decomposition information is also written if required (the
   domain decomposition info is required for multi-fileset multi-threaded
   I/O). If the &lt;data&gt; element is allocated, it is considered to be a
   space axis, otherwise it is a time axis with an unlimited dimension. Only
   one time axis is allowed per file. <br>
<br>   The second form defines a field. Metadata corresponding to the type
   above are written to the file on &lt;unit&gt;. A unique ID for subsequen
   references to this field is returned in field%id. At least one axis
   must be associated, 0D variables are not considered. mpp_write_meta
   must previously have been called on all axes associated with this
   field. <br>
<br>   The third form (3 - 5) defines metadata associated with a previously defined
   axis or field, identified to mpp_write_meta by its unique ID &lt;id&gt;.
   The attribute is named &lt;name&gt; and can take on a real, integer
   or character value. &lt;rval&gt; and &lt;ival&gt; can be scalar or 1D arrays.
   This need not be called for attributes already contained in
   the type. <br>
<br>   The last form (6 - 8) defines global metadata associated with the file as a
   whole. The attribute is named &lt;name&gt; and can take on a real, integer
   or character value. &lt;rval&gt; and &lt;ival&gt; can be scalar or 1D arrays. <br>
<br>   Note that <tt>mpp_write_meta</tt>   is expecting axis data on the <i>global</i>   domain even if it is a domain-decomposed axis.
   <br>
<br>
   You cannot interleave calls to <tt>mpp_write</tt>   and <tt>mpp_write_meta</tt>: the first call to <tt>mpp_write</tt>   implies that metadata specification is complete. </dd>
<br>
<br>
</dl>
</li>
<li>
<a name="mpp_write"></a>
<h4>mpp_write</h4>
<pre> 
<b>mpp_write</b> ( unit, axis )</pre>
<pre> 
<b>mpp_write</b> ( unit, field, data, tstamp )</pre>
<pre> 
<b>mpp_write</b> ( unit, field, domain, data, tstamp )</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd> 
<tt>mpp_write</tt>   is used to write data to the file on an I/O unit
   using the file parameters supplied by <a href="#mpp_open"></a>. Axis and field definitions must
   have previously been written to the file using <a href="#mpp_write_meta"></a>.  There are three
   forms of <tt>mpp_write</tt>, one to write axis data, one to write
   distributed field data, and one to write non-distributed field
   data. <i>Distributed</i>   data refer to arrays whose two
   fastest-varying indices are domain-decomposed. Distributed data must
   be 2D or 3D (in space). Non-distributed data can be 0-3D.
   <br>
<br>
   The <tt>data</tt>   argument for distributed data is expected by <tt>mpp_write</tt>   to contain data specified on the <i>data</i>   domain,
   and will write the data belonging to the <i>compute</i>   domain,
   fetching or sending data as required by the parallel I/O <a href="#modes">mode</a>   specified in the <tt>mpp_open</tt>   call. This
   is consistent with our definition of <a href="http:models/bgrid_solo/fms_src/shared/mpp/mpp_domains.html#domains">domains</a>, where all arrays are
   expected to be dimensioned on the data domain, and all operations
   performed on the compute domain.
   <br>
<br>
   The type of the <tt>data</tt>   argument must be a <i>default
   real</i>, which can be 4 or 8 byte. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>tstamp&nbsp;&nbsp;&nbsp;</tt></td><td> <tt>tstamp</tt>   is an optional argument. It is to
   be omitted if the field was defined not to be a function of time.
   Results are unpredictable if the argument is supplied for a time-
   independent field, or omitted for a time-dependent field. Repeated
   writes of a time-independent field are also not recommended. One
   time level of one field is written per call. tstamp must be an 8-byte
   real, even if the default real type is 4-byte. </td>
</tr>
</table>
</dd>
<br>
<dt>
<b>NOTE</b>
</dt>
<dd>   The type of write performed by <tt>mpp_write</tt>   depends on the file
   characteristics on the I/O unit specified at the <a href="#mpp_open"></a>   call. Specifically, the format of
   the output data (e.g netCDF or IEEE), the <tt>threading</tt>   and <tt>fileset</tt>   flags, etc., can be changed there, and require no
   changes to the <tt>mpp_write</tt>   calls.
   <br>
<br>
   Packing is currently not implemented for non-netCDF files, and the <tt>pack</tt>   attribute is ignored. On netCDF files, <tt>NF_DOUBLE</tt>s (8-byte IEEE floating point numbers) are
   written for <tt>pack</tt>=1 and <tt>NF_FLOAT</tt>s for <tt>pack</tt>=2. (<tt>pack</tt>=2 gives the customary
   and default behaviour). We write <tt>NF_SHORT</tt>s (2-byte
   integers) for <tt>pack=4</tt>, or <tt>NF_BYTE</tt>s
   (1-byte integers) for <tt>pack=8</tt>. Integer scaling is done
   using the <tt>scale</tt>   and <tt>add</tt>   attributes at <tt>pack</tt>=4 or 8, satisfying the relation
   <br>
<br> 
<pre>    data = packed_data*scale + add</pre> 
<tt>NOTE: mpp_write</tt>   does not check to see if the scaled
   data in fact fits into the dynamic range implied by the specified
   packing. It is incumbent on the user to supply correct scaling
   attributes.
   <br>
<br>
   You cannot interleave calls to <tt>mpp_write</tt>   and <tt>mpp_write_meta</tt>: the first call to <tt>mpp_write</tt>   implies that metadata specification is
   complete. </dd>
<br>
<br>
</dl>
</li>
<li>
<a name="mpp_read"></a>
<h4>mpp_read</h4>
<pre>
<b>call mpp_read </b>( unit, field, data, time_index )</pre>
<pre>
<b>call mpp_read </b>( unit, field, domain, data, time_index )</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd> 
<tt>mpp_read</tt>   is used to read data to the file on an I/O unit
   using the file parameters supplied by <a href="#mpp_open"></a>. There are two
   forms of <tt>mpp_read</tt>, one to read
   distributed field data, and one to read non-distributed field
   data. <i>Distributed</i>   data refer to arrays whose two
   fastest-varying indices are domain-decomposed. Distributed data must
   be 2D or 3D (in space). Non-distributed data can be 0-3D.
   <br>
<br>
   The <tt>data</tt>   argument for distributed data is expected by <tt>mpp_read</tt>   to contain data specified on the <i>data</i>   domain,
   and will read the data belonging to the <i>compute</i>   domain,
   fetching data as required by the parallel I/O <a href="#modes">mode</a>   specified in the <tt>mpp_open</tt>   call. This
   is consistent with our definition of <a href="http:models/bgrid_solo/fms_src/shared/mpp/mpp_domains.html#domains">domains</a>, where all arrays are
   expected to be dimensioned on the data domain, and all operations
   performed on the compute domain. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>unit&nbsp;&nbsp;&nbsp;</tt></td><td>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>field&nbsp;&nbsp;&nbsp;</tt></td><td>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[type(fieldtype)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>domain&nbsp;&nbsp;&nbsp;</tt></td><td></td>
</tr>
<tr>
<td valign="top" align="left"><tt>time_index&nbsp;&nbsp;&nbsp;</tt></td><td>   time_index is an optional argument. It is to be omitted if the
   field was defined not to be a function of time. Results are
   unpredictable if the argument is supplied for a time- independent
   field, or omitted for a time-dependent field. </td>
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
<td valign="top" align="left"><tt>data&nbsp;&nbsp;&nbsp;</tt></td><td>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:,:,:)]</span></td>
</tr>
</table>
</dd>
<br>
<dt>
<b>NOTE</b>
</dt>
<dd>   The type of read performed by <tt>mpp_read</tt>   depends on
   the file characteristics on the I/O unit specified at the <a href="#mpp_open"></a>   call. Specifically, the
   format of the input data (e.g netCDF or IEEE) and the <tt>threading</tt>   flags, etc., can be changed there, and
   require no changes to the <tt>mpp_read</tt>   calls. (<tt>fileset</tt>   = MPP_MULTI is not supported by <tt>mpp_read</tt>; IEEE is currently not supported).
   <br>
<br>
   Packed variables are unpacked using the <tt>scale</tt>   and <tt>add</tt>   attributes.
   <br>
<br> 
<tt>mpp_read_meta</tt>   must be called prior to calling <tt>mpp_read.</tt> 
</dd>
<br>
<br>
</dl>
</li>
<li>
<a name="mpp_get_atts"></a>
<h4>mpp_get_atts</h4>
<pre>
<b>call mpp_get_atts </b>( unit, global_atts)</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   Get file global metdata. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>unit&nbsp;&nbsp;&nbsp;</tt></td><td>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>global_atts&nbsp;&nbsp;&nbsp;</tt></td><td>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[atttype, dimension(:)]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="mpp_io_init"></a>
<h4>mpp_io_init</h4>
<pre>
<b>call mpp_io_init </b>( flags, maxunit )</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   Called to initialize the <tt>mpp_io_mod</tt>   package. Sets the range
   of valid fortran units and initializes the <tt>mpp_file</tt>   array of <tt>type(filetype)</tt>. <tt>mpp_io_init</tt>   will call <tt>mpp_init</tt>   and <tt>mpp_domains_init</tt>, to make sure its parent modules have been
   initialized. (Repeated calls to the <tt>init</tt>   routines do no harm,
   so don't worry if you already called it). </dd>
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
<tr>
<td valign="top" align="left"><tt>maxunit&nbsp;&nbsp;&nbsp;</tt></td><td>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="mpp_io_exit"></a>
<h4>mpp_io_exit</h4>
<pre>
<b>call mpp_io_exit </b>()</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   It is recommended, though not at present required, that you call this
   near the end of a run. This will close all open files that were opened
   with <a href="#mpp_open"></a>. Files opened otherwise
   are not affected. </dd>
<br>
<br>
</dl>
</li>
<li>
<a name="mpp_open"></a>
<h4>mpp_open</h4>
<pre>
<b>call mpp_open </b>( unit, file, action, form, access, threading, fileset, iospec, nohdrs, recl, pelist )</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   Open a file for parallel I/O. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>file&nbsp;&nbsp;&nbsp;</tt></td><td>   file is the filename: REQUIRED
   we append .nc to filename if it is a netCDF file
   we append .&lt;pppp&gt; to filename if fileset is private (pppp is PE number) <br>&nbsp;&nbsp;&nbsp;<span class="type">[character(len=*)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>action&nbsp;&nbsp;&nbsp;</tt></td><td>   action is one of MPP_RDONLY, MPP_APPEND, MPP_WRONLY or MPP_OVERWR. <br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>form&nbsp;&nbsp;&nbsp;</tt></td><td>   form is one of MPP_ASCII:  formatted read/write
   MPP_NATIVE: unformatted read/write with no conversion
   MPP_IEEE32: unformatted read/write with conversion to IEEE32
   MPP_NETCDF: unformatted read/write with conversion to netCDF <br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>access&nbsp;&nbsp;&nbsp;</tt></td><td>   access is one of MPP_SEQUENTIAL or MPP_DIRECT (ignored for netCDF).
   RECL argument is REQUIRED for direct access IO. <br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>threading&nbsp;&nbsp;&nbsp;</tt></td><td>   threading is one of MPP_SINGLE or MPP_MULTI
   single-threaded IO in a multi-PE run is done by PE0. <br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>fileset&nbsp;&nbsp;&nbsp;</tt></td><td>   fileset is one of MPP_MULTI and MPP_SINGLE
   fileset is only used for multi-threaded I/O
   if all I/O PEs in &lt;pelist&gt; use a single fileset, they write to the same file
   if all I/O PEs in &lt;pelist&gt; use a multi  fileset, they each write an independent file <br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>pelist&nbsp;&nbsp;&nbsp;</tt></td><td>   pelist is the list of I/O PEs (currently ALL). <br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>recl&nbsp;&nbsp;&nbsp;</tt></td><td>   recl is the record length in bytes. <br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>iospec&nbsp;&nbsp;&nbsp;</tt></td><td>   iospec is a system hint for I/O organization, e.g assign(1) on SGI/Cray systems. <br>&nbsp;&nbsp;&nbsp;<span class="type">[character(len=*)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>nohdrs&nbsp;&nbsp;&nbsp;</tt></td><td>   nohdrs has no effect when action=MPP_RDONLY|MPP_APPEND or when form=MPP_NETCDF <br>&nbsp;&nbsp;&nbsp;<span class="type">[logical]</span></td>
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
<td valign="top" align="left"><tt>unit&nbsp;&nbsp;&nbsp;</tt></td><td>   unit is intent(OUT): always _returned_by_ mpp_open(). <br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
</table>
</dd>
<br>
<dt>
<b>NOTE</b>
</dt>
<dd>   The integer parameters to be passed as flags (<tt>MPP_RDONLY</tt>,
   etc) are all made available by use association. The <tt>unit</tt>   returned by <tt>mpp_open</tt>   is guaranteed unique. For non-netCDF I/O
   it is a valid fortran unit number and fortran I/O can be directly called
   on the file.
   <br>
<br> 
<tt>MPP_WRONLY</tt>   will guarantee that existing files named <tt>file</tt>   will not be clobbered. <tt>MPP_OVERWR</tt>   allows overwriting of files.
   <br>
<br>
   Files opened read-only by many processors will give each processor
   an independent pointer into the file, i.e:
   <br>
<br> 
<pre>      namelist / nml / ...
   ...
      call mpp_open( unit, 'input.nml', action=MPP_RDONLY )
      read(unit,nml)</pre>   will result in each PE independently reading the same namelist.
   <br>
<br>
   Metadata identifying the file and the version of <tt>mpp_io_mod</tt>   are written to a file that is opened <tt>MPP_WRONLY</tt>   or <tt>MPP_OVERWR</tt>. If this is a
   multi-file set, and an additional global attribute <tt>NumFilesInSet</tt>   is written to be used by post-processing
   software.
   <br>
<br>
   If <tt>nohdrs=.TRUE.</tt>   all calls to write attributes will
   return successfully <i>without</i>   performing any writes to the
   file. The default is <tt>.FALSE.</tt>.
   <br>
<br>
   For netCDF files, headers are always written even if <tt>nohdrs=.TRUE.</tt>   The string <tt>iospec</tt>   is passed to the OS to
   characterize the I/O to be performed on the file opened on <tt>unit</tt>. This is typically used for I/O optimization. For
   example, the FFIO layer on SGI/Cray systems can be used for
   controlling synchronicity of reads and writes, buffering of data
   between user space and disk for I/O optimization, striping across
   multiple disk partitions, automatic data conversion and the like
   (<tt>man intro_ffio</tt>). All these actions are controlled through
   the <tt>assign</tt>   command. For example, to specify asynchronous
   caching of data going to a file open on <tt>unit</tt>, one would do:
   <br>
<br> 
<pre>   call mpp_open( unit, ... iospec='-F cachea' )</pre>   on an SGI/Cray system, which would pass the supplied <tt>iospec</tt>   to the <tt>assign(3F)</tt>   system call.
   <br>
<br>
   Currently <tt>iospec </tt>performs no action on non-SGI/Cray
   systems. The interface is still provided, however: users are cordially
   invited to add the requisite system calls for other systems. </dd>
<br>
<br>
</dl>
</li>
<li>
<a name="mpp_close"></a>
<h4>mpp_close</h4>
<pre>
<b>call mpp_close </b>( unit, action )</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   Closes the open file on <tt>unit</tt>. Clears the <tt>type(filetype)</tt>   object <tt>mpp_file(unit)</tt>   making it
   available for reuse. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>unit&nbsp;&nbsp;&nbsp;</tt></td><td> 
<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>action&nbsp;&nbsp;&nbsp;</tt></td><td> 
<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="mpp_read_meta"></a>
<h4>mpp_read_meta</h4>
<pre>
<b>call mpp_read_meta </b>(unit)</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   This routine is used to read the <a href="#metadata">metadata</a>   describing the contents of a file. Each file can contain any number of
   fields, which are functions of 0-3 space axes and 0-1 time axes. (Only
   one time axis can be defined per file). The basic metadata defined <a href="#metadata">above</a>   for <tt>axistype</tt>   and <tt>fieldtype</tt>   are stored in <tt>mpp_io_mod</tt>   and
   can be accessed outside of <tt>mpp_io_mod</tt>   using calls to <tt>mpp_get_info</tt>, <tt>mpp_get_atts</tt>, <tt>mpp_get_vars</tt>   and <tt>mpp_get_times</tt>. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>unit&nbsp;&nbsp;&nbsp;</tt></td><td> 
<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
</table>
</dd>
<br>
<dt>
<b>NOTE</b>
</dt>
<dd> 
<tt>mpp_read_meta</tt>   must be called prior to <tt>mpp_read</tt>. </dd>
<br>
<br>
</dl>
</li>
<li>
<a name="mpp_get_info"></a>
<h4>mpp_get_info</h4>
<pre>
<b>call mpp_get_info </b>( unit, ndim, nvar, natt, ntime )</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   Get some general information about a file. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>unit&nbsp;&nbsp;&nbsp;</tt></td><td> 
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
<td valign="top" align="left"><tt>ndim&nbsp;&nbsp;&nbsp;</tt></td><td> 
<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>nvar&nbsp;&nbsp;&nbsp;</tt></td><td> 
<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>natt&nbsp;&nbsp;&nbsp;</tt></td><td> 
<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>ntime&nbsp;&nbsp;&nbsp;</tt></td><td> 
<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="mpp_get_times"></a>
<h4>mpp_get_times</h4>
<pre>
<b>call mpp_get_times </b>( unit, time_values )</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   Get file time data. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>unit&nbsp;&nbsp;&nbsp;</tt></td><td> 
<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
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
<td valign="top" align="left"><tt>time_values&nbsp;&nbsp;&nbsp;</tt></td><td> 
<br>&nbsp;&nbsp;&nbsp;<span class="type">[real(DOUBLE_KIND), dimension(:)]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="mpp_flush"></a>
<h4>mpp_flush</h4>
<pre>
<b>call mpp_flush </b>(unit)</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   Flushes the open file on <tt>unit</tt>   to disk. Any outstanding
   asynchronous writes will be completed. Any buffer layers between the
   user and the disk (e.g the FFIO layer on SGI/Cray systems) will be
   flushed. Calling <tt>mpp_flush</tt>   on a unit opened with the <tt>MPP_RDONLY</tt>   attribute is likely to lead to erroneous behaviour. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>unit&nbsp;&nbsp;&nbsp;</tt></td><td> 
<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="mpp_get_ncid"></a>
<h4>mpp_get_ncid</h4>
<pre> 
<b>mpp_get_ncid</b> (unit)</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   This returns the <tt>ncid</tt>   associated with the open file on <tt>unit</tt>. It is used in the instance that the user desires to
   perform netCDF calls upon the file that are not provided by the <tt>mpp_io_mod</tt>   API itself. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>unit&nbsp;&nbsp;&nbsp;</tt></td><td> 
<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
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
<dd>   Any module or program unit using <tt>mpp_io_mod</tt>   must contain the line <pre>     use mpp_io_mod</pre>   If netCDF output is desired, the cpp flag <tt>-Duse_netCDF</tt>   must be turned on. The loader step requires an explicit link to the
   netCDF library (typically something like <tt>-L/usr/local/lib
   -lnetcdf</tt>, depending on the path to the netCDF library). <a href="http://www.unidata.ucar.edu/packages/netcdf/guidef">netCDF
   release 3 for fortran</a>   is required.
   <br>
<br>
   Please also consider the compiling and linking requirements of <a href="models/bgrid_solo/fms_src/shared/mpp/mpp_domains.html#linking"></a>   and <a href="models/bgrid_solo/fms_src/shared/mpp/mpp.html#linking"></a>, which are <tt>use</tt>d by this module. </dd>
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
<tt>mpp_io_mod</tt>   uses standard f90. On SGI/Cray systems, certain I/O
   characteristics are specified using <tt>assign(3F)</tt>. On other
   systems, the user may have to provide similar capability if required.
   <br>
<br>
   There are some OS-dependent
   pre-processor directives that you might need to modify on
   non-SGI/Cray systems and compilers. </dd>
</dl>
</div>
<br>
<!-- END PRECOMPILER OPTIONS -->
<a name="LOADER OPTIONS"></a>
<hr>
<h4>LOADER OPTIONS</h4>
<!-- BEGIN LOADER -->
<div>
<p>    The    source consists of the main source file    and also requires the following include files:    (when compiled with )       GFDL users can check it out of the main CVS repository as part of
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
