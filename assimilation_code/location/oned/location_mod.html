<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>module location_mod (1D)</title>
<link rel="stylesheet" type="text/css" href=
"../../../docs/html/doc.css">
<link href="../../../docs/images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>MODULE (1D) location_mod</h1>
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
interface code (in model_mod.f90), and to be able to read and write
location information to files. DART isolates all this location
information into separate modules so that the main algorithms can
operate with the same code independent of whether the model uses
latitude/longitude/height, one-d unit sphere coordinates,
cylindrical coordinates, etc. DART provides about half a dozen
possible coordinate systems, and others can be added.</p>
<p>This locations module provides a representation of a physical
location on a periodic 1D domain with location values between 0 and
1. A type that abstracts the location is provided along with
operators to set, get, read, write, and compute distances between
locations. This is a member of a class of similar location modules
that provide the same abstraction for different represenations of
physical space.</p>
<p>All possible location modules define the same module name
<em class="code">location_mod</em>. Therefore, the DART framework
and any user code should include a Fortran 90 'use' statement of
'location_mod'. The selection of exactly which location module is
compiled is specified by the source file name in the <em class=
"file">path_names_xxx</em> file, which is read by the <em class=
"file">mkmf_xxx</em> scripts.</p>
<p>The model-specific <em class="file">model_mod.f90</em> files
need to define four <em class="code">get_close</em> routines, but
in most cases they can simply put a <em class="code">use</em>
statement at the top which uses the routines in the locations
module, and they do not have to provide any additional code.</p>
<p>However, if the model interface code wants to intercept and
alter the default behavior of the get_close routines, they are able
to. The correct usage of the <em class="code">get_close</em>
routines is as follows:</p>
<pre>
<em class="code">
call get_close_maxdist_init()  ! must be called before get_close_obs_init()
call get_close_obs_init()
...
call get_close_obs()           ! many, many times
...
call get_close_obs_destroy()
</em>
</pre>
<p>Regardless of the fact that the names include the string 'obs',
they are intended for use with any group of locations in the
system, frequently state vector items or observations, but any
location is acceptable.</p>
<!--==================================================================-->
<!--================== DESCRIPTION OF A NAMELIST =====================-->
<!--==================================================================-->
<a name="Namelist" id="Namelist"></a>
<hr>
<h2>NAMELIST</h2>
<p>This version of the locations module does not have any namelist
input.</p>
<!--==================================================================-->
<a name="Interface" id="Interface"></a>
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
<td><a href="#vert_is_level">vert_is_level</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#vert_is_height">vert_is_height</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#vert_constants">VERTISUNDEF</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#vert_constants">VERTISSURFACE</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#vert_constants">VERTISLEVEL</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#vert_constants">VERTISPRESSURE</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#vert_constants">VERTISHEIGHT</a></td>
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
<p>There is currently no namelist interface for the 1D location
module.</p>
<p>A note about documentation style. Optional arguments are
enclosed in brackets <em class="optionalcode">[like this]</em>.</p>
<!--===================== DESCRIPTION OF A LOCAL TYPE =====================-->
<a name="location_type" id="location_type"></a><br>
<div class="type"><em class="call">type location_type</em>
<pre>
   private
   real(r8) :: x
end type location_type
</pre></div>
<div class="indent1"><!-- Description -->
<p>Provides an abstract representation of physical location on a
one-dimensional periodic domain.</p>
<table border="0" cellpadding="3" width="100%">
<tr>
<th align="left">Component</th>
<th align="left">Description</th>
</tr>
<tr>
<td valign="top">x</td>
<td>Location has range 0 to 1</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A LOCAL TYPE =====================-->
 <a name="get_close_type" id="get_close_type"></a><br>
<div class="type"><em class="call">type get_close_type</em>
<pre>
   private
   integer  :: num
   real(r8) :: maxdist
end type get_close_type
</pre></div>
<div class="indent1"><!-- Description -->
<p>Provides a structure for doing efficient computation of close
locations. Doesn't do anything in the 1D implementation except
provide appropriate stubs.</p>
<table border="0" cellpadding="3" width="100%">
<tr>
<th align="left">Component</th>
<th align="left">Description</th>
</tr>
<tr>
<td valign="top">num</td>
<td>Number of locations in list</td>
</tr>
<tr>
<td valign="top">maxdist</td>
<td>Threshhold distance. Anything closer is close.</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_location" id="get_location"></a><br>
<div class="routine"><em class="call">var = get_location(loc)</em>
<pre>
real(r8)                        :: <em class=
"code">get_location</em>
type(location_type), intent(in) :: <em class="code">loc</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Extracts the real location value, range 0 to 1, from a location
type.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">get_location</em></td>
<td>The real value for a location</td>
</tr>
<tr>
<td valign="top"><em class="code">loc</em></td>
<td>A location derived type</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="set_location" id="set_location"></a><br>
<div class="routine"><em class="call">var = set_location(x)</em>
<pre>
type(location_type)   :: <em class="code">set_location</em>
real(r8), intent(in)  :: <em class="code">x</em>

</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns a location type with the location x.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">set_location</em></td>
<td>A location derived type</td>
</tr>
<tr>
<td valign="top"><em class="code">x</em></td>
<td>Location value in the range 0. to 1.</td>
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
<p>Returns the location value if attr = 'X' or if attr is not
passed.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">query_location</em></td>
<td>Returns value of x.</td>
</tr>
<tr>
<td valign="top"><em class="code">loc</em></td>
<td>A location type</td>
</tr>
<tr>
<td valign="top"><em class="optionalcode">attr</em></td>
<td>Selects 'X'</td>
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
<p>Returns a location with location set to missing value from
types_mod.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">set_location_missing</em></td>
<td>A location set to missing value</td>
</tr>
</table>
</div>
<br>
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_close_maxdist_init" id=
"get_close_maxdist_init"></a><br>
<div class="routine"><em class="call">call
get_close_maxdist_init(gc,maxdist <em class="optionalcode">,
[maxdist_list]</em>)</em>
<pre>
type(get_close_type), intent(inout) :: <em class="code">gc</em>
real(r8), intent(in)                :: <em class=
"code">maxdist</em>
real(r8), intent(in), optional      :: <em class=
"optionalcode">maxdist_list(:)</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Sets the threshhold distance. Anything closer than this is
deemed to be close. This routine must be called first, before the
other <em class="code">get_close</em> routines. It allocates space
so it is necessary to call <em class=
"code">get_close_obs_destroy</em> when completely done with getting
distances between locations.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">gc</em></td>
<td>Data for efficiently finding close locations.</td>
</tr>
<tr>
<td valign="top"><em class="code">maxdist</em></td>
<td>Anything closer than this distance is a close location.</td>
</tr>
<tr>
<td valign="top"><em class="optionalcode">maxdist_list</em></td>
<td>Ignored for this location type.</td>
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
close to a given location. The oned implementation is minimal and
just records the number of locations here. Must be called after
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
base_obs_loc, base_obs_kind, obs, obs_kind, num_close, close_ind,
dist)</em>
<pre>
type(get_close_type),              intent(in)  :: <em class=
"code">gc</em>
type(location_type),               intent(in)  :: <em class=
"code">base_obs_loc</em>
integer,                           intent(in)  :: <em class=
"code">base_obs_kind</em>
type(location_type), dimension(:), intent(in)  :: <em class=
"code">obs</em>
integer, dimension(:),             intent(in)  :: <em class=
"code">obs_kind</em>
integer,                           intent(out) :: <em class=
"code">num_close</em>
integer, dimension(:),             intent(out) :: <em class=
"code">close_ind</em>
real(r8), dimension(:),            intent(out) :: <em class=
"code">dist</em>
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
<td valign="top"><em class="code">base_obs_kind</em></td>
<td>Kind of the single location.</td>
</tr>
<tr>
<td valign="top"><em class="code">obs</em></td>
<td>List of locations from which close ones are to be found.</td>
</tr>
<tr>
<td valign="top"><em class="code">obs_kind</em></td>
<td>Kind associated with locations in obs list.</td>
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
<td valign="top"><em class="code">dist</em></td>
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
<em class="optionalcode">[, kind1, kind2]</em>)</em>
<pre>
real(r8)                        :: get_dist
type(location_type), intent(in) :: <em class="code">loc1</em>
type(location_type), intent(in) :: <em class="code">loc2</em>
integer, optional,   intent(in) :: <em class=
"optionalcode">kind1</em>
integer, optional,   intent(in) :: <em class=
"optionalcode">kind2</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Return the distance between 2 locations. Since this is a
periodic domain, the shortest distance may wrap around.</p>
<p>The kind arguments are not used by the default location code,
but are available to any user-supplied distance routines which want
to do specialized calculations based on the kinds associated with
each of the two locations.</p>
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
<td valign="top"><em class="optionalcode">kind1</em></td>
<td>DART kind associated with location 1.</td>
</tr>
<tr>
<td valign="top"><em class="optionalcode">kind2</em></td>
<td>DART kind associated with location 2.</td>
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
<p>Always returns false; this locations module has no vertical
coordinates. Provided only for compile-time compatibility with
other location modules.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">vert_is_undef</em></td>
<td>Always returns .FALSE.</td>
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
<p>Always returns false; this locations module has no vertical
coordinates. Provided only for compile-time compatibility with
other location modules.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">vert_is_surface</em></td>
<td>Always returns .FALSE.</td>
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
<p>Always returns false; this locations module has no vertical
coordinates. Provided only for compile-time compatibility with
other location modules.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">vert_is_pressure</em></td>
<td>Always returns .FALSE.</td>
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
<p>Always returns false; this locations module has no vertical
coordinates. Provided only for compile-time compatibility with
other location modules.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">vert_is_level</em></td>
<td>Always returns .FALSE.</td>
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
<p>Always returns false; this locations module has no vertical
coordinates. Provided only for compile-time compatibility with
other location modules.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">vert_is_height</em></td>
<td>Always returns .FALSE.</td>
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
<p>Always returns false; this locations module has no vertical
coordinates. Provided only for compile-time compatibility with
other location modules.</p>
<p>See note in threed_sphere locations module about the function
name.</p>
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
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="vert_constants" id="vert_constants"></a><br>
<div class="routine">
<pre>
<em class="call">integer, parameter :: VERTISUNDEF    = -2</em>
<em class="call">integer, parameter :: VERTISSURFACE  = -1</em>
<em class="call">integer, parameter :: VERTISLEVEL    =  1</em>
<em class="call">integer, parameter :: VERTISPRESSURE =  2</em>
<em class="call">integer, parameter :: VERTISHEIGHT   =  3</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>This locations module has no vertical coordinate, but for
compatibility with other location modules, these are defined.</p>
</div>
<br>
<!--=============== DESCRIPTION OF A PUBLIC CONSTANT =================-->
 <a name="LocationDims" id="LocationDims"></a><br>
<div class="routine">
<pre>
<em class="call">integer, parameter :: LocationDims = 1</em>
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
"call">character(len=129), parameter :: LocationName = "loc1d"</em>
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
"call">character(len=129), parameter :: LocationLName = "location on unit circle"</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>This is a <b>constant</b>. A parameter to identify this location
module in output long name metadata.</p>
</div>
<br>
<!--==================================================================-->
<!-- Describe the Files Used by this module.                          -->
<!--==================================================================-->
 <a name="FilesUsed" id="FilesUsed"></a>
<hr>
<h2>FILES</h2>
<p>None.</p>
<!--==================================================================-->
<!-- Cite references, if need be.                                     -->
<!--==================================================================-->
<a name="References" id="References"></a>
<hr>
<h2>REFERENCES</h2>
<ol>
<li>none</li>
</ol>
<!--==================================================================-->
<!-- Describe all the error conditions and codes.                     -->
<!--==================================================================-->
<a name="Errors" id="Errors"></a>
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
<td valign="top">set_location</td>
<!-- message -->
<td valign="top">Value of x is out of 0-&gt;1 range</td>
<!-- comment -->
<td valign="top">x cannot be less than 0 or greater than 1</td>
</tr>
<tr><!-- routine -->
<td valign="top">query_location</td>
<!-- message -->
<td valign="top">Only x is legal attribute to request from
location</td>
<!-- comment -->
<td valign="top">attr must be 'x' or 'X'</td>
</tr>
<tr><!-- routine -->
<td valign="top">read_location</td>
<!-- message -->
<td valign="top">Expected location header "loc1d" in input
file</td>
<!-- comment -->
<td valign="top">Can only read one-dimensional location files</td>
</tr>
</table>
</div>
<h2>KNOWN BUGS</h2>
<p>none at this time</p>
<!--==================================================================-->
<!-- Describe Future Plans.                                           -->
<!--==================================================================-->
<a name="FuturePlans" id="FuturePlans"></a>
<hr>
<h2>FUTURE PLANS</h2>
<p>None.</p>
<!--==================================================================-->
<!-- PrivateComponents                                                -->
<!--==================================================================-->
<a name="PrivateComponents" id="PrivateComponents"></a>
<hr>
<h2>PRIVATE COMPONENTS</h2>
<p>N/A</p>
<!--==================================================================-->
<!-- Legalese & Metadata                                              -->
<!--==================================================================-->
<a name="Legalese" id="Legalese"></a>
<hr>
<h2>Terms of Use</h2>
<p>DART software - Copyright UCAR. This open source software is
provided by UCAR, "as is", without charge, subject to all terms of
use at <a href=
"http://www.image.ucar.edu/DAReS/DART/DART_download">http://www.image.ucar.edu/DAReS/DART/DART_download</a></p>
<!--==================================================================-->
</body>
</html>
