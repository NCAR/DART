<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>module location_mod (threed_cartesian)</title>
<link rel="stylesheet" type="text/css" href=
"../../../docs/html/doc.css">
<link href="../../../docs/images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>MODULE location_mod (threed_cartesian)</h1>
<table border="0" summary="" cellpadding="5">
<tr>
<td valign="middle"><img src="../../../docs/images/Dartboard7.png"
alt="DART project logo" height="70"></td>
<td>Jump to <a href="../../../docs/index.html">DART Documentation
Main Index</a></td>
</tr>
</table>
<a href="#Namelist">NAMELIST</a> / <a href=
"#Interface">INTERFACES</a> / <a href="#FilesUsed">FILES</a> /
<a href="#References">REFERENCES</a> / <a href="#Errors">ERRORS</a>
/ <a href="#FuturePlans">PLANS</a> / <a href=
"#PrivateComponents">PRIVATE COMPONENTS</a> / <a href=
"#Legalese">TERMS OF USE</a>
<h2>Overview</h2>
<p>The DART framework needs to be able to compute distances between
locations, to pass location information to and from the model
interface code (model_mod.f90), and to be able to read and write
location information to files. DART isolates all this location
information into separate modules so that the main algorithms can
operate with the same code independent of whether the model uses
latitude/longitude/height, 1D unit cartesian coordinates,
cylindrical coordinates, etc. DART provides about half a dozen
possible coordinate systems, and others can be added. The most
common one for geophysical models is the <a href=
"../threed_sphere/location_mod.html">threed_sphere</a> version.
This document describes an alternative 3D cartesian coordinate
system.</p>
<p><strong>Note that only one location module can be compiled into
any single DART executable, and most earth observational data is
generated in [latitude, longitude, vertical pressure or height]
coordinates - the threed_sphere option. The cartesian and 3D sphere
locations cannot be mixed or used together.</strong></p>
<p>This location module provides a representation of a physical
location in an [X, Y, Z] cartesian coordinate space. A type that
abstracts the location is provided along with operators to set,
get, read, write, and compute distances between locations. This is
a member of a class of similar location modules that provide the
same abstraction for different represenations of physical
space.</p>
<h4>Location-independent code</h4>
<p>All types of location modules define the same module name
<em class="code">location_mod</em>. Therefore, the DART framework
and any user code should include a Fortran 90 <em class=
"code">use</em> statement of <em class="code">location_mod</em>.
The selection of which location module will be compiled into the
program is controlled by which source file name is specified in the
<em class="file">path_names_xxx</em> file, which is used by the
<em class="file">mkmf_xxx</em> scripts.</p>
<p>All types of location modules define the same Fortran 90 derived
type <em class="code">location_type</em>. Programs that need to
pass location information to subroutines but do not need to
interpret the contents can declare, receive, and pass this derived
type around in their code independent of which location module is
specified at compile time. Model and location-independent utilities
should be written in this way. However, as soon as the contents of
the location type needs to be accessed by user code then it becomes
dependent on the exact type of location module that it is compiled
with.</p>
<h4>Usage of distance routines</h4>
<a name="Distance" id="Distance"></a>
<p>Regardless of the fact that the distance subroutine names
include the string 'obs', there is nothing specific to observations
in these routines. They work to compute distances between any set
of locations. The most frequent use of these routines in the filter
code is to compute the distance between a single observation and
items in the state vector, and also between a single observation
and other nearby observations. However, any source for locations is
supported.</p>
<p>In simpler location modules (like the <em class="file">oned</em>
version) there is no need for anything other than a brute force
search between the base location and all available state vector
locations. However in the case of large geophysical models which
typically use the <em class="file">threed_cartesian</em> locations
code, the brute-force search time is prohibitive. The location code
pre-processes all locations into a set of <em>bins</em> and then
only needs to search the lists of locations in nearby bins when
looking for locations that are within a specified distance.</p>
<p>The expected calling sequence of the <em class=
"code">get_close</em> routines is as follows:</p>
<pre>
<code>
call get_close_maxdist_init()  ! is called before get_close_obs_init()
call get_close_obs_init()

call get_close_obs()           ! called many, many times

call get_close_obs_destroy()
</code>
</pre>
<p>In the <em class="file">threed_cartesian</em> implementation the
first routine initializes some data structures, the second one bins
up the list of locations, and then the third one is called multiple
times to find all locations within a given radius of some reference
location, and to optionally compute the exact separation distance
from the reference location. The last routine deallocates the
space. See the documentation below for the specific details for
each routine.</p>
<p>All 4 of these routines must be present in every location module
but in most other versions all but <em class=
"code">get_close_obs()</em> are stubs. In this <em class=
"file">threed_cartesian</em> version of the locations module all
are fully implemented.</p>
<h4>Interaction with model_mod.f90 code</h4>
<a name="ModelMod" id="ModelMod"></a>
<p>The filter and other DART programs could call the <em class=
"code">get_close</em> routines directly, but typically do not. They
declare them (in a <em class="code">use</em> statement) to be in
the <em class="code">model_mod</em> module, and all model interface
modules are required to supply them. However in many cases the
model_mod only needs to contain another <em class="code">use</em>
statement declaring them to come from the <em class=
"code">location_mod</em> module. Thus they 'pass through' the
model_mod but the user does not need to provide a subroutine or any
code for them.</p>
<p>However, if the model interface code wants to intercept and
alter the default behavior of the get_close routines, it is able
to. Typically the model_mod still calls the location_mod routines
and then adjusts the results before passing them back to the
calling code. To do that, the model_mod must be able to call the
routines in the location_mod which have the same names as the
subroutines it is providing. To allow the compiler to distinguish
which routine is to be called where, we use the Fortran 90 feature
which allows a module routine to be renamed in the use statement.
For example, a common case is for the model_mod to want to supply
additions to the get_close_obs() routine only. At the top of the
model_mod code it would declare:</p>
<pre>
<em class="code">
use location_mod, only :: location_get_close_obs =&gt; get_close_obs,    &amp;
                          get_close_maxdist_init, get_close_obs_init, &amp;
                          get_close_obs_destroy
</em>
</pre>
<p>That makes calls to the maxdist_init, init, and destroy routines
simply pass through to the code in the location_mod, but the
model_mod must supply a get_close_obs() subroutine. When it wants
to call the code in the location_mod it calls <em class=
"code">location_get_close_obs()</em>.</p>
<p>One use pattern is for the model_mod to call the location
get_close_obs() routine without the <em class="code">dist</em>
argument. This returns a list of any potentially close locations
without computing the exact distance from the base location. At
this point the list of locations is a copy and the model_mod
routine is free to alter the list in any way it chooses: for
example, it can change the locations to make certain types of
locations appear closer or further away from the base location.
Then typically the model_mod code loops over the list calling the
<em class="code">get_dist()</em> routine to get the actual
distances to be returned to the calling code.</p>
<h4>Horizontal Distance Only</h4>
<p>This option is not supported for the threed_cartesian
option.</p>
<h4>Precomputation for Run-time Search Efficiency</h4>
<p>For search efficiency all locations are pre-binned. For the
non-octree option, the total list of locations is divided up into
<em>nx</em> by <em>ny</em> by <em>nz</em> boxes and the index
numbers of all items (both state vector entries and observations)
are stored in the appropriate box. To locate all points close to a
given location, only the locations listed in the boxes within the
search radius must be checked. This speeds up the computations, for
example, when localization controls which state vector items are
impacted by any given observation. The search radius is the
localization distance and only those state vector items in boxes
closer than the radius to the observation location are
processed.</p>
<p>The default values have given good performance on many of our
existing model runs, but for tuning purposes the box counts have
been added to the namelist to allow adjustment. By default the code
prints some summary information about how full the average box is,
how many are empty, and how many items were in the box with the
largest count. The namelist value <em>output_box_info</em> can be
set to .true. to get even more information about the box
statistics. The best performance will be obtained somewhere between
two extremes; the worst extreme is all the points are located in
just a few boxes. This degenerates into a (slow) linear search
through the index list. The other extreme is a large number of
empty or sparsely filled boxes. The overhead of creating, managing,
and searching a long list of boxes will impact performance. The
best performance lies somewhere in the middle, where each box
contains a reasonable number of values, more or less evenly
distributed across boxes. The absolute numbers for best performance
will certainly vary from case to case.</p>
<!--=====================================================================-->
<!--===================== DESCRIPTION OF A NAMELIST =====================-->
<!--=====================================================================-->
<a name="Namelist" id="Namelist"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>NAMELIST</h2>
<p>This namelist is read from the file <em class=
"file">input.nml</em>. Namelists start with an ampersand '&amp;'
and terminate with a slash '/'. Character strings that contain a
'/' must be enclosed in quotes to prevent them from prematurely
terminating the namelist.</p>
<div class="namelist">
<pre>
&amp;location_nml
   nx                  = 10
   ny                  = 10
   nz                  = 10
   x_is_periodic       = .false.
   min_x_for_periodic  = -888888.0
   max_x_for_periodic  = -888888.0
   y_is_periodic       = .false.
   min_y_for_periodic  = -888888.0
   max_y_for_periodic  = -888888.0
   z_is_periodic       = .false.
   min_z_for_periodic  = -888888.0
   max_z_for_periodic  = -888888.0
   compare_to_correct  = .false.
   output_box_info     = .false.
   print_box_level     = 0
   debug               = 0
  /
</pre></div>
<br>
<br>
<p>Items in this namelist either control the way in which distances
are computed and/or influence the code performance.</p>
<div>
<table border="0" cellpadding="10" width="100%" summary=
'namelist description'>
<thead align="left">
<tr>
<th>Item</th>
<th>Type</th>
<th>Description</th>
</tr>
</thead>
<tbody>
<tr>
<td>nx, ny, nz</td>
<td>integer</td>
<td>The number of boxes in each dimension to use to speed the
searches. This is <strong>not</strong> the number of
gridcells.</td>
</tr>
<tr>
<td>x_is_periodic, y_is_periodic, z_is_periodic</td>
<td>logical</td>
<td>If .true., the domain wraps in the coordinate.</td>
</tr>
<tr>
<td>min_x_for_periodic, max_x_for_periodic</td>
<td>real(r8)</td>
<td>The minimum and maximum values that are considered to be
identical locations if <em class=
"code">x_is_periodic = .true.</em></td>
</tr>
<tr>
<td>min_y_for_periodic, max_y_for_periodic</td>
<td>real(r8)</td>
<td>The minimum and maximum values that are considered to be
identical locations if <em class=
"code">y_is_periodic = .true.</em></td>
</tr>
<tr>
<td>min_z_for_periodic, max_z_for_periodic</td>
<td>real(r8)</td>
<td>The minimum and maximum values that are considered to be
identical locations if <em class=
"code">z_is_periodic = .true.</em></td>
</tr>
<tr>
<td>compare_to_correct</td>
<td>logical</td>
<td>If true, do an exhaustive search for the closest point. Only
useful for debugging because the performance cost is
prohibitive.</td>
</tr>
<tr>
<td>output_box_info</td>
<td>logical</td>
<td>Print out debugging info.</td>
</tr>
<tr>
<td>print_box_level</td>
<td>logical</td>
<td>If output_box_info is true, how detailed should the output
be.</td>
</tr>
<tr>
<td>debug</td>
<td>integer</td>
<td>The higher the number, the more verbose the run-time output. 0
(zero) is the minimum run-time output.</td>
</tr>
</tbody>
</table>
</div>
<br>
<br>
<!--==================================================================-->
 <a name="Interface" id="Interface"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>OTHER MODULES USED</h2>
<pre>
types_mod
utilities_mod
random_seq_mod
</pre>
<!--==================================================================-->
<!-- Declare all public entities ...                                  -->
<!-- duplicate public routines template as many times as necessary    -->
<!-- make sure you replace all yyyroutine?? strings                   -->
<!--==================================================================-->
<!--Note to authors. The first row of the table is different.         -->
<!--==================================================================-->
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>PUBLIC INTERFACES</h2>
<table>
<tr>
<td><em class="code">use location_mod, only :</em></td>
<td><a href="#location_type">location_type</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_close_type">get_close_type</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_location">get_location</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#set_location">set_location</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#write_location">write_location</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#read_location">read_location</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#interactive_location">interactive_location</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#set_location_missing">set_location_missing</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#query_location">query_location</a></td>
</tr>
<tr>
<td> </td>
<td><a href=
"#get_close_maxdist_init">get_close_maxdist_init</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_close_obs_init">get_close_obs_init</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_close_obs">get_close_obs</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_close_obs_destroy">get_close_obs_destroy</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_dist">get_dist</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#LocationDims">LocationDims</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#LocationName">LocationName</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#LocationLName">LocationLName</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#horiz_dist_only">horiz_dist_only</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#vert_is_undef">vert_is_undef</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#vert_is_surface">vert_is_surface</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#vert_is_pressure">vert_is_pressure</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#vert_is_scale_height">vert_is_scale_height</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#vert_is_level">vert_is_level</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#vert_is_height">vert_is_height</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#equal">operator(==)</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#not_equal">operator(/=)</a></td>
</tr>
</table>
<p>Namelist interface <a href="#Namelist"><em class=
"code">&amp;location_nml</em></a> must be read from file <em class=
"file">input.nml</em>.</p>
<p>A note about documentation style. Optional arguments are
enclosed in brackets <em class="optionalcode">[like this]</em>.</p>
<!--===================== DESCRIPTION OF A LOCAL TYPE =====================-->
<a name="location_type" id="location_type"></a><br>
<div class="type"><em class="call">type location_type</em>
<pre>
   private
   real(r8) :: x, y, z
end type location_type
</pre></div>
<div class="indent1"><!-- Description -->
<p>Provides an abstract representation of physical location in a 3D
cartesian space.</p>
<table border="0" cellpadding="3" width="100%">
<tr>
<th align="left">Component</th>
<th align="left">Description</th>
</tr>
<tr>
<td valign="top">x, y, z</td>
<td>location in each dimension</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A LOCAL TYPE =====================-->
 <a name="get_close_type" id="get_close_type"></a><br>
<div class="type"><em class="call">type get_close_type</em>
<pre>
   private
   integer, pointer  :: loc_box(:)           ! (nloc); List of loc indices in boxes
   integer, pointer  :: count(:, :, :)       ! (nx, ny, nz); # of locs in each box
   integer, pointer  :: start(:, :, :)       ! (nx, ny, nz); Start of list of locs in this box
   real(r8)          :: bot_x, top_x         ! extents in x, y, z
   real(r8)          :: bot_y, top_y
   real(r8)          :: bot_z, top_z
   real(r8)          :: x_width, y_width, z_width    ! widths of boxes in x,y,z
   real(r8)          :: nboxes_x, nboxes_y, nboxes_z ! based on maxdist how far to search
end type get_close_type
</pre></div>
<div class="indent1"><!-- Description -->
<p>Provides a structure for doing efficient computation of close
locations.</p>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_location" id="get_location"></a><br>
<div class="routine"><em class="call">var = get_location(loc)</em>
<pre>
real(r8), dimension(3)          :: <em class=
"code">get_location</em>
type(location_type), intent(in) :: <em class="code">loc</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Extracts the x, y, z locations from a location type and returns
in a 3 element real array.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">get_location</em></td>
<td>The x,y,z values</td>
</tr>
<tr>
<td valign="top"><em class="code">loc</em></td>
<td>A location type</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="set_location" id="set_location"></a><br>
<div class="routine"><em class="call">var = set_location(x, y,
z)</em><br>
<em class="call">var = set_location(lon, lat, height, radius)</em>
<pre>
type(location_type) :: <em class="code">set_location</em>
real(r8), intent(in)    :: <em class="code">x</em>
real(r8), intent(in)    :: <em class="code">y</em>
real(r8), intent(in)    :: <em class="code">z</em>
</pre>
<p>or</p>
<pre>
type(location_type) :: <em class="code">set_location</em>
real(r8), intent(in)    :: <em class="code">lon</em>
real(r8), intent(in)    :: <em class="code">lat</em>
real(r8), intent(in)    :: <em class="code">height</em>
real(r8), intent(in)    :: <em class="code">radius</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns a location type with the input [x,y,z] or allows the
input to be specified as locations on the surface of a sphere with
a specified radius and height above the surface.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">set_location</em></td>
<td>A location type</td>
</tr>
<tr>
<td valign="top"><em class="code">x, y, z</em></td>
<td>Coordinates along each axis</td>
</tr>
<tr>
<td valign="top"><em class="code">lon, lat</em></td>
<td>Longitude, Latitude in degrees</td>
</tr>
<tr>
<td valign="top"><em class="code">height</em></td>
<td>Vertical location in same units as radius (e.g. meters)</td>
</tr>
<tr>
<td valign="top"><em class="code">radius</em></td>
<td>The radius of the sphere in same units as height (e.g.
meters)</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="write_location" id="write_location"></a><br>
<div class="routine"><em class="call">call write_location(locfile,
loc <em class="optionalcode">[, fform, charstring]</em>)</em>
<pre>
integer,               intent(in)       :: <em class=
"code"> locfile </em>
type(location_type),   intent(in)       :: <em class=
"code"> loc </em>
character(len=*), optional, intent(in)  :: <em class=
"optionalcode"> fform </em>
character(len=*), optional, intent(out) :: <em class=
"optionalcode"> charstring </em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Given an integer IO channel of an open file and a location,
writes the location to this file. The <em class=
"optionalcode">fform</em> argument controls whether write is
"FORMATTED" or "UNFORMATTED" with default being formatted. If the
final <em class="optionalcode">charstring</em> argument is
specified, the formatted location information is written to the
character string only, and the <em class="code">locfile</em>
argument is ignored.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">locfile</em></td>
<td>the unit number of an open file.</td>
</tr>
<tr>
<td valign="top"><em class="code">loc</em></td>
<td>location type to be written.</td>
</tr>
<tr>
<td valign="top"><em class="optionalcode">fform</em></td>
<td>Format specifier ("FORMATTED" or "UNFORMATTED"). Default is
"FORMATTED" if not specified.</td>
</tr>
<tr>
<td valign="top"><em class="optionalcode">charstring</em></td>
<td>Character buffer where formatted location string is written if
present, and no output is written to the file unit.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="read_location" id="read_location"></a><br>
<div class="routine"><em class="call">var = read_location(locfile
<em class="optionalcode">[, fform]</em>)</em>
<pre>
type(location_type)                    :: <em class=
"code">read_location</em>
integer, intent(in)                    :: <em class=
"code">locfile</em>
character(len=*), optional, intent(in) :: <em class=
"optionalcode">fform</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Reads a location_type from a file open on channel locfile using
format <em class="optionalcode">fform</em> (default is
formatted).</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">read_location</em></td>
<td>Returned location type read from file</td>
</tr>
<tr>
<td valign="top"><em class="code">locfile</em></td>
<td>Integer channel opened to a file to be read</td>
</tr>
<tr>
<td valign="top"><em class="optionalcode">fform</em></td>
<td>Optional format specifier ("FORMATTED" or "UNFORMATTED").
Default "FORMATTED".</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="interactive_location" id="interactive_location"></a><br>
<div class="routine"><em class="call">call
interactive_location(location <em class="optionalcode">[,
set_to_default]</em>)</em>
<pre>
type(location_type), intent(out) :: <em class="code">location</em>
logical, optional, intent(in)    :: <em class=
"optionalcode">set_to_default</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Use standard input to define a location type. With
set_to_default true get one with all elements set to 0.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">location</em></td>
<td>Location created from standard input</td>
</tr>
<tr>
<td valign="top"><em class="optionalcode">set_to_default</em></td>
<td>If true, sets all elements of location type to 0</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="query_location" id="query_location"></a><br>
<div class="routine"><em class="call">var = query_location(loc
<em class="optionalcode">[, attr]</em>)</em>
<pre>
real(r8)                               :: <em class=
"code">query_location</em>
type(location_type), intent(in)        :: <em class="code">loc</em>
character(len=*), optional, intent(in) :: <em class=
"optionalcode">attr</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns the value of x, y, z depending on the attribute
specification. If attr is not present, returns 'x'.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">query_location</em></td>
<td>Returns x, y, or z.</td>
</tr>
<tr>
<td valign="top"><em class="code">loc</em></td>
<td>A location type</td>
</tr>
<tr>
<td valign="top"><em class="optionalcode">attr</em></td>
<td>Selects 'X', 'Y', 'Z'. If not specified, 'X' is returned.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="set_location_missing" id="set_location_missing"></a><br>
<div class="routine"><em class="call">var =
set_location_missing()</em>
<pre>
type(location_type) :: <em class="code">set_location_missing</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns a location with all elements set to missing values
defined in types module.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">set_location_missing</em></td>
<td>A location with all elements set to missing values</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_close_maxdist_init" id=
"get_close_maxdist_init"></a><br>
<div class="routine"><em class="call">call
get_close_maxdist_init(gc,maxdist, <em class=
"optionalcode">[maxdist_list]</em>)</em>
<pre>
type(get_close_type), intent(inout) :: <em class="code">gc</em>
real(r8), intent(in)                :: <em class=
"code">maxdist</em>
real(r8), intent(in), optional      :: <em class=
"optionalcode">maxdist_list(:)</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Sets the threshhold distance. <em class="code">maxdist</em> is
in units of radians. Anything closer than this is deemed to be
close. This routine must be called first, before the other
<em class="code">get_close</em> routines. It allocates space so it
is necessary to call <em class="code">get_close_obs_destroy</em>
when completely done with getting distances between locations.</p>
<p>If the last optional argument is not specified, maxdist applies
to all locations. If the last argument is specified, it must be a
list of exactly the length of the number of specific types in the
obs_kind_mod.f90 file. This length can be queried with the <a href=
"../../modules/observations/obs_kind_mod.html#get_num_types_of_obs">
get_num_types_of_obs()</a> function to get count of obs types. It
allows a different maximum distance to be set per base type when
get_close() is called.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">gc</em></td>
<td>Data for efficiently finding close locations.</td>
</tr>
<tr>
<td valign="top"><em class="code">maxdist</em></td>
<td>Anything closer than this number of radians is a close
location.</td>
</tr>
<tr>
<td valign="top"><em class="optionalcode">maxdist</em></td>
<td>If specified, must be a list of real values. The length of the
list must be exactly the same length as the number of observation
types defined in the obs_def_kind.f90 file. (See <a href=
"../../modules/observations/obs_kind_mod.html#get_num_types_of_obs">
get_num_types_of_obs()</a> to get count of obs types.) The values
in this list are used for the obs types as the close distance
instead of the maxdist argument.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_close_obs_init" id="get_close_obs_init"></a><br>
<div class="routine"><em class="call">call get_close_obs_init(gc,
num, obs)</em>
<pre>
type(get_close_type),             intent(inout) :: <em class=
"code">gc</em>
integer,                          intent(in)    :: <em class=
"code">num</em>
type(location_type), dimension(:) intent(in)    :: <em class=
"code">obs</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Initialize storage for efficient identification of locations
close to a given location. Allocates storage for keeping track of
which 'box' each location in the list is in. Must be called after
<em class="code">get_close_maxdist_init</em>, and the list of
locations here must be the same as the list of locations passed
into <em class="code">get_close_obs()</em>. If the list changes,
<em class="code">get_close_obs_destroy()</em> must be called, and
both the initialization routines must be called again. It allocates
space so it is necessary to call <em class=
"code">get_close_obs_destroy</em> when completely done with getting
distances between locations.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">gc</em></td>
<td>Structure that contains data to efficiently find locations
close to a given location.</td>
</tr>
<tr>
<td valign="top"><em class="code">num</em></td>
<td>The number of locations in the list.</td>
</tr>
<tr>
<td valign="top"><em class="code">obs</em></td>
<td>The locations of each element in the list, not used in 1D
implementation.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_close_obs" id="get_close_obs"></a><br>
<div class="routine"><em class="call">call get_close_obs(gc,
base_obs_loc, base_obs_type, obs, obs_kind, num_close, close_ind,
dist)</em>
<pre>
type(get_close_type),              intent(in)  :: <em class=
"code">gc</em>
type(location_type),               intent(in)  :: <em class=
"code">base_obs_loc</em>
integer,                           intent(in)  :: <em class=
"code">base_obs_type</em>
type(location_type), dimension(:), intent(in)  :: <em class=
"code">obs</em>
integer,             dimension(:), intent(in)  :: <em class=
"code">obs_kind</em>
integer,                           intent(out) :: <em class=
"code">num_close</em>
integer,             dimension(:), intent(out) :: <em class=
"code">close_ind</em>
real(r8), optional,  dimension(:), intent(out) :: <em class=
"optionalcode">dist</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Given a single location and a list of other locations, returns
the indices of all the locations close to the single one along with
the number of these and the distances for the close ones. The list
of locations passed in via the <em class="code">obs</em> argument
must be identical to the list of <em class="code">obs</em> passed
into the most recent call to <em class=
"code">get_close_obs_init()</em>. If the list of locations of
interest changes <em class="code">get_close_obs_destroy()</em> must
be called and then the two initialization routines must be called
before using <em class="code">get_close_obs()</em> again.</p>
<p>Note that the base location is passed with the specific type
associated with that location. The list of potential close
locations is matched with a list of generic kinds. This is because
in the current usage in the DART system the base location is always
associated with an actual observation, which has both a specific
type and generic kind. The list of potentially close locations is
used both for other observation locations but also for state
variable locations which only have a generic kind.</p>
<p>If called without the optional <em class=
"optionalcode">dist</em> argument, all locations that are
potentially close are returned, which is likely a superset of the
locations that are within the threshold distance specified in the
<em class="code">get_close_maxdist_init()</em> call.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">gc</em></td>
<td>Structure to allow efficient identification of locations close
to a given location.</td>
</tr>
<tr>
<td valign="top"><em class="code">base_obs_loc</em></td>
<td>Single given location.</td>
</tr>
<tr>
<td valign="top"><em class="code">base_obs_type</em></td>
<td>Specific type of the single location.</td>
</tr>
<tr>
<td valign="top"><em class="code">obs</em></td>
<td>List of locations from which close ones are to be found.</td>
</tr>
<tr>
<td valign="top"><em class="code">obs_kind</em></td>
<td>Generic kind associated with locations in obs list.</td>
</tr>
<tr>
<td valign="top"><em class="code">num_close</em></td>
<td>Number of locations close to the given location.</td>
</tr>
<tr>
<td valign="top"><em class="code">close_ind</em></td>
<td>Indices of those locations that are close.</td>
</tr>
<tr>
<td valign="top"><em class="optionalcode">dist</em></td>
<td>Distance between given location and the close ones identified
in close_ind.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_close_obs_destroy" id=
"get_close_obs_destroy"></a><br>
<div class="routine"><em class="call">call
get_close_obs_destroy(gc)</em>
<pre>
type(get_close_type), intent(inout) :: <em class="code">gc</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Releases memory associated with the <em class="code">gc</em>
derived type. Must be called whenever the list of locations
changes, and then <em class="code">get_close_maxdist_init</em> and
<em class="code">get_close_obs_init</em> must be called again with
the new locations list.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">gc</em></td>
<td>Data for efficiently finding close locations.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_dist" id="get_dist"></a><br>
<div class="routine"><em class="call">var = get_dist(loc1, loc2,
<em class="optionalcode">[, type1, kind2, no_vert]</em>)</em>
<pre>
real(r8)                        :: get_dist
type(location_type), intent(in) :: <em class="code">loc1</em>
type(location_type), intent(in) :: <em class="code">loc2</em>
integer, optional,   intent(in) :: <em class=
"optionalcode">type1</em>
integer, optional,   intent(in) :: <em class=
"optionalcode">kind2</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns the distance between two locations.</p>
<p>The type and kind arguments are not used by the default location
code, but are available to any user-supplied distance routines
which want to do specialized calculations based on the types/kinds
associated with each of the two locations.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">loc1</em></td>
<td>First of two locations to compute distance between.</td>
</tr>
<tr>
<td valign="top"><em class="code">loc2</em></td>
<td>Second of two locations to compute distance between.</td>
</tr>
<tr>
<td valign="top"><em class="optionalcode">type1</em></td>
<td>DART specific type associated with location 1.</td>
</tr>
<tr>
<td valign="top"><em class="optionalcode">kind2</em></td>
<td>DART generic kind associated with location 2.</td>
</tr>
<tr>
<td valign="top"><em class="code">var</em></td>
<td>distance between loc1 and loc2.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="vert_is_undef" id="vert_is_undef"></a><br>
<div class="routine"><em class="call">var = vert_is_undef(loc)</em>
<pre>
logical                         :: <em class=
"code">vert_is_undef</em>
type(location_type), intent(in) :: <em class="code">loc</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Always returns .false.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">vert_is_undef</em></td>
<td>Always returns .false.</td>
<td></td>
</tr>
<tr>
<td valign="top"><em class="code">loc</em></td>
<td>A location type</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="vert_is_surface" id="vert_is_surface"></a><br>
<div class="routine"><em class="call">var =
vert_is_surface(loc)</em>
<pre>
logical                         :: <em class=
"code">vert_is_surface</em>
type(location_type), intent(in) :: <em class="code">loc</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Always returns .false.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">vert_is_surface</em></td>
<td>Always returns .false.</td>
<td></td>
</tr>
<tr>
<td valign="top"><em class="code">loc</em></td>
<td>A location type</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="vert_is_pressure" id="vert_is_pressure"></a><br>
<div class="routine"><em class="call">var =
vert_is_pressure(loc)</em>
<pre>
logical                         :: <em class=
"code">vert_is_pressure</em>
type(location_type), intent(in) :: <em class="code">loc</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Always returns .false.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">vert_is_pressure</em></td>
<td>Always returns .false.</td>
<td></td>
</tr>
<tr>
<td valign="top"><em class="code">loc</em></td>
<td>A location type</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="vert_is_scale_height" id="vert_is_scale_height"></a><br>
<div class="routine"><em class="call">var =
vert_is_scale_height(loc)</em>
<pre>
logical                         :: <em class=
"code">vert_is_scale_height</em>
type(location_type), intent(in) :: <em class="code">loc</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Always returns .false.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">vert_is_scale_height</em></td>
<td>Always returns .false.</td>
<td></td>
</tr>
<tr>
<td valign="top"><em class="code">loc</em></td>
<td>A location type</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="vert_is_level" id="vert_is_level"></a><br>
<div class="routine"><em class="call">var = vert_is_level(loc)</em>
<pre>
logical                         :: <em class=
"code">vert_is_level</em>
type(location_type), intent(in) :: <em class="code">loc</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Always returns .false.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">vert_is_level</em></td>
<td>Always returns .false.</td>
<td></td>
</tr>
<tr>
<td valign="top"><em class="code">loc</em></td>
<td>A location type</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="vert_is_height" id="vert_is_height"></a><br>
<div class="routine"><em class="call">var =
vert_is_height(loc)</em>
<pre>
logical                         :: <em class=
"code">vert_is_height</em>
type(location_type), intent(in) :: <em class="code">loc</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Always returns .false.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">vert_is_height</em></td>
<td>Always returns .false.</td>
<td></td>
</tr>
<tr>
<td valign="top"><em class="code">loc</em></td>
<td>A location type</td>
</tr>
</table>
</div>
<br>
<!--=============== DESCRIPTION OF A ROUTINE =================-->
 <a name="has_vertical_localization" id=
"has_vertical_localization"></a><br>
<div class="routine"><em class="call">var =
has_vertical_localization()</em>
<pre>
logical :: has_vertical_localization
</pre></div>
<div class="indent1"><!-- Description -->
<p>Always returns .false.</p>
<p>This routine should perhaps be renamed to something like
'using_vertical_for_distance' or something similar. The current use
for it is in the localization code inside filter, but that doesn't
make this a representative function name. And at least in current
usage, returning the opposite setting of the namelist item makes
the code read more direct (fewer double negatives).</p>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="equal" id="equal"></a><br>
<div class="routine"><em class="call">loc1 == loc2</em>
<pre>
type(location_type), intent(in) :: <em class="code">loc1, loc2</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns true if the two location types have identical values,
else false.</p>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="not_equal" id="not_equal"></a><br>
<div class="routine"><em class="call">loc1 /= loc2</em>
<pre>
type(location_type), intent(in) :: <em class="code">loc1, loc2</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns true if the two location types do NOT have identical
values, else false.</p>
</div>
<br>
<!--=============== DESCRIPTION OF A PUBLIC CONSTANT =================-->
 <a name="LocationDims" id="LocationDims"></a><br>
<div class="routine">
<pre>
<em class="call">integer, parameter :: LocationDims = 3</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>This is a <b>constant</b>. Contains the number of real values in
a location type. Useful for output routines that must deal
transparently with many different location modules.</p>
</div>
<br>
<!--=============== DESCRIPTION OF A PUBLIC CONSTANT =================-->
 <a name="LocationName" id="LocationName"></a><br>
<div class="routine">
<pre>
<em class=
"call">character(len=129), parameter :: LocationName = "loc3Dcartesian"</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>This is a <b>constant</b>. A parameter to identify this location
module in output metadata.</p>
</div>
<br>
<!--============= DESCRIPTION OF A PUBLIC CONSTANT =================-->
 <a name="LocationLName" id="LocationLName"></a><br>
<div class="routine">
<pre>
<em class=
"call">character(len=129), parameter :: LocationLName = <br>
       "threed cartesian locations: x, y, z"</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>This is a <b>constant</b>. A parameter set to "threed cartesian
locations: x, y, z" used to identify this location module in output
long name metadata.</p>
</div>
<br>
<!--==================================================================-->
<!-- Describe the Files Used by this module.                          -->
<!--==================================================================-->
 <a name="FilesUsed" id="FilesUsed"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>FILES</h2>
<table border="0">
<tr>
<th>filename</th>
<th>purpose</th>
</tr>
<tr>
<td>input.nml</td>
<td>to read the location_mod namelist</td>
</tr>
</table>
<!--==================================================================-->
<!-- Cite references, if need be.                                     -->
<!--==================================================================-->
<a name="References" id="References"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>REFERENCES</h2>
<ol>
<li>none</li>
</ol>
<!--==================================================================-->
<!-- Describe all the error conditions and codes.                     -->
<!--==================================================================-->
<a name="Errors" id="Errors"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>ERROR CODES and CONDITIONS</h2>
<div class="errors">
<table border="1" cellspacing="1" cellpadding="10" width="100%">
<tr>
<th>Routine</th>
<th>Message</th>
<th>Comment</th>
</tr>
<tr><!-- routine -->
<td valign="top">nc_write_location</td>
<!-- message -->
<td valign="top">Various NetCDF-f90 interface error messages</td>
<!-- comment -->
<td valign="top">From one of the NetCDF calls in
nc_write_location</td>
</tr>
</table>
</div>
<h2>KNOWN BUGS</h2>
<p>The octree code works fine to store values, but the search for
all points within a given radius of a base point is not supported.
So for this module the 3D box option (use_octree = .false.) should
be used.</p>
<!--==================================================================-->
<!-- Describe Future Plans.                                           -->
<!--==================================================================-->
<a name="FuturePlans" id="FuturePlans"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>FUTURE PLANS</h2>
<p>Want to fix octree code, and make it easier to detect when bad
combinations of tuning parameters are being used.</p>
<p>See the note in the 'has_vertical_localization()' about a better
name for this routine.</p>
<p>The functions of 'get_close_maxdist_init()' and
'get_close_obs_init()' appear to be able to be combined into a
single init routine. This impacts all model_mods, however, since
they can intercept these routines. Doing this will be a
non-backwards compatible change.</p>
<p>The use of 'obs' in all these routine names should probably be
changed to 'loc' since there is no particular dependence that they
be observations. They may need to have an associated DART kind, but
these routines are used for DART state vector entries so it's
misleading to always call them 'obs'.</p>
<!--==================================================================-->
<!-- PrivateComponents                                                -->
<!--==================================================================-->
<a name="PrivateComponents" id="PrivateComponents"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>PRIVATE COMPONENTS</h2>
<p>N/A</p>
<!--==================================================================-->
<!-- Legalese & Metadata                                              -->
<!--==================================================================-->
<a name="Legalese" id="Legalese"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>Terms of Use</h2>
<p>DART software - Copyright UCAR. This open source software is
provided by UCAR, "as is", without charge, subject to all terms of
use at <a href=
"http://www.image.ucar.edu/DAReS/DART/DART_download">http://www.image.ucar.edu/DAReS/DART/DART_download</a></p>
<!--==================================================================-->
</body>
</html>
