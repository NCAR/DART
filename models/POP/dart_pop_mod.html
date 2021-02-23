<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
    "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0" />
<title>module dart_pop_mod</title>
<link rel="stylesheet" type="text/css" href=
"../../docs/html/doc.css" />
<link href="../../docs/images/dart.ico" rel="shortcut icon" />
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>MODULE dart_pop_mod (POP)</h1>
<table border="0" summary="" cellpadding="5">
<tr>
<td valign="middle"><img src="../../docs/images/Dartboard7.png"
alt="DART project logo" height="70" /></td>
<td>Jump to <a href="../../docs/index.html">DART Documentation Main
Index</a></td>
</tr>
</table>
<a href="#Namelist">NAMELIST</a> / <a href=
"#Interface">INTERFACES</a> / <a href="#FilesUsed">FILES</a> /
<a href="#References">REFERENCES</a> / <a href="#Errors">ERRORS</a>
/ <a href="#FuturePlans">PLANS</a> / <a href=
"#PrivateComponents">PRIVATE COMPONENTS</a> / <a href=
"#Legalese">TERMS OF USE</a>
<h2>Overview</h2>
<p><em class="program">dart_pop_mod</em> provides a consistent
collection of routines that are useful for multiple programs e.g.
<em class="program">dart_to_pop</em>, <em class=
"program">pop_to_dart</em>, etc.</p>
<!--==================================================================-->
<!--=================== DESCRIPTION OF A NAMELIST  ===================-->
<!--==================================================================-->
<a name="Namelist" id="Namelist"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr />
<h2>NAMELIST</h2>
<p>There are no namelists unique to this module. It is necessary
for this module to read some of the POP namelists, and so they are
declared in this module. In one instance, DART will read the
<em class="code">time_manager_nml</em> namelist and <b>write</b> an
updated version to control the length of the integration of POP.
All other information is simply read from the namelists and is used
in the same context as POP itself. The POP documentation should be
consulted. <strong>Only the variables of interest to DART are
described in this document.</strong><br />
<br />
All namelists are read from a file named <em class=
"file">pop_in</em>.</p>
<!--==================================================================-->
<a name="time_manager_nml" id="time_manager_nml"></a>
<div class="namelist">
<pre>
<em class=
"call">namelist /time_manager_nml/ </em> allow_leapyear, stop_count, stop_option
</pre></div>
<div class="indent1"><!-- Description -->
<p><em class="program">dart_to_pop</em> controls the model advance
of LANL/POP by creating a <em class=
"code">&amp;time_manager_nml</em> in <em class=
"file">pop_in.DART</em> <strong>IFF</strong> the DART state being
converted has the 'advance_to_time' record. The <em class=
"file">pop_in.DART</em> must be concatenated with the other
namelists needed by POP into a file called <em class=
"file">pop_in</em> . We have chosen to store the other
namelists (which contain static information) in a file called
<em class="file">pop_in.part2</em>. Initially, the <em class=
"code">time_manager_nml</em> is stored in a companion file called
<em class="file">pop_in.part1</em> and the two files are
concatenated into the expected <em class="file">pop_in</em> - then,
during the course of an assimilation experiment, DART keeps writing
out a new <em class="code">time_manager_nml</em> with new
integration information - which gets appended with the static
information in <em class="file">pop_in.part2</em> </p>
<table border="0" cellpadding="3" width="100%">
<tr>
<th align="left">Contents</th>
<th align="left">Type</th>
<th align="left">Description</th>
</tr>
<tr><!--contents-->
<td valign="top">allow_leapyear   </td>
<!--  type  -->
<td valign="top">logical</td>
<!--descript-->
<td valign="top">DART ignores the setting of this parameter. All
observations must use a Gregorian calendar. There are pathological
cases, but if you are doing data assimilation, just use the
Gregorian calendar. <em class=
"units">[default: .true.]</em></td>
</tr>
<tr><!--contents-->
<td valign="top">stop_count</td>
<!--  type  -->
<td valign="top">integer</td>
<!--descript-->
<td valign="top">the number of model advance steps to take.
<em class="units">[default: 1]</em></td>
</tr>
<tr><!--contents-->
<td valign="top">stop_option</td>
<!--  type  -->
<td valign="top">character(len=64)   </td>
<!--descript-->
<td valign="top">The units for the number of model advance steps
(<em class="code">stop_count</em>) to take. <em class=
"units">[default: 'ndays']</em></td>
</tr>
</table>
</div>
<br />
<!--==================================================================-->
 <a name="io_nml" id="io_nml"></a>
<div class="namelist">
<pre>
<em class=
"call">namelist /io_nml/ </em> luse_pointer_files, pointer_filename
</pre></div>
<div class="indent1"><!-- Description -->
<table border="0" cellpadding="3" width="100%">
<tr>
<th align="left">Contents</th>
<th align="left">Type</th>
<th align="left">Description</th>
</tr>
<tr><!--contents-->
<td valign="top">luse_pointer_files</td>
<!--  type  -->
<td valign="top">logical</td>
<!--descript-->
<td valign="top">switch to indicate the use of pointer files or
not. If <em class="code">.true.</em>, a pointer file is used to
contain the name of the restart file to be used. DART requires this
to be <em class="code">.true</em>. <em class=
"units">[default: .true.]</em></td>
</tr>
<tr><!--contents-->
<td valign="top">pointer_filename</td>
<!--  type  -->
<td valign="top">character(len=100)</td>
<!--descript-->
<td valign="top">The name of the pointer file. All of the DART
scripts presume and require the use of the default. Each ensmeble
member gets its own pointer file. <em class=
"units">[default: rpointer.ocn.[1-N].restart]</em></td>
</tr>
</table>
</div>
<br />
<!--==================================================================-->
 <a name="restart_nml" id="restart_nml"></a><br />
<div class="namelist">
<pre>
<em class=
"call">namelist /restart_nml/ </em> restart_freq_opt, restart_freq
</pre></div>
<div class="indent1"><!-- Description -->
<table border="0" cellpadding="3" width="100%">
<tr>
<th align="left">Contents</th>
<th align="left">Type</th>
<th align="left">Description</th>
</tr>
<tr><!--contents-->
<td valign="top">luse_pointer_files   </td>
<!--  type  -->
<td valign="top">logical</td>
<!--descript-->
<td valign="top">switch to indicate the use of pointer files or
not. If <em class="code">.true.</em>, a pointer file is used to
contain the name of the restart file to be used. DART requires this
to be <em class="code">.true</em>. <em class=
"units">[default: .true.]</em></td>
</tr>
<tr><!--contents-->
<td valign="top">pointer_filename</td>
<!--  type  -->
<td valign="top">character(len=100)   </td>
<!--descript-->
<td valign="top">The name of the pointer file. All of the DART
scripts presume and require the use of the default. Each ensmeble
member gets its own pointer file. <em class=
"units">[default: rpointer.ocn.[1-N].restart]</em></td>
</tr>
</table>
</div>
<br />
<!--==================================================================-->
 <a name="init_ts_nml" id="init_ts_nml"></a><br />
<div class="namelist">
<pre>
<em class=
"call">namelist /init_ts_nml/ </em> init_ts_option, init_ts_file, init_ts_file_fmt
</pre></div>
<div class="indent1"><!-- Description -->
<p>The <em class="program">dart_pop_mod:initialize_module()</em>
routine reads <em class="file">pop_in</em> . There are several
code stubs for future use that may allow for a more fully-supported
POP namelist implementation. This namelist is one of them. Until
further notice, the <em class="code">init_ts_nml</em> is completely
ignored by DART.</p>
<table border="0" cellpadding="3" width="100%">
<tr>
<th align="left">Contents</th>
<th align="left">Type</th>
<th align="left">Description</th>
</tr>
<tr><!--contents-->
<td valign="top">init_ts_option</td>
<!--  type  -->
<td valign="top">character(len=64)</td>
<!--descript-->
<td valign="top">NOT USED by DART. All T,S information comes from a
netCDF restart file named <em class="file">pop.r.nc</em> <em class=
"units">[default: 'restart']</em></td>
</tr>
<tr><!--contents-->
<td valign="top">init_ts_file</td>
<!--  type  -->
<td valign="top">character(len=100)   </td>
<!--descript-->
<td valign="top">NOT USED by DART. All T,S information comes from
<em class="file">pop.r.nc</em> <em class=
"units">[default: 'pop.r']</em></td>
</tr>
<tr><!--contents-->
<td valign="top">init_ts_file_fmt   </td>
<!--  type  -->
<td valign="top">character(len=64)</td>
<!--descript-->
<td valign="top">NOT USED by DART. The file format is <em class=
"code">'nc'</em> <em class="units">[default: 'nc']</em></td>
</tr>
</table>
</div>
<br />
<!--==================================================================-->
 <a name="domain_nml" id="domain_nml"></a><br />
<div class="namelist">
<pre>
<em class="call">namelist /domain_nml/ </em> ew_boundary_type
</pre></div>
<div class="indent1"><!-- Description -->
<p>DART needs to know if the East-West domain is cyclic for spatial
interpolations. Presently, DART has only been tested for the dipole
grid, which is cyclic E-W and closed N-S.</p>
<table border="0" cellpadding="3" width="100%">
<tr>
<th align="left">Contents</th>
<th align="left">Type</th>
<th align="left">Description</th>
</tr>
<tr><!--contents-->
<td valign="top">ew_boundary_type   </td>
<!--  type  -->
<td valign="top">character(len=64)   </td>
<!--descript-->
<td valign="top">switch to indicate whether the East-West domain is
cyclic or not. DART/POP has not been tested in a regional
configuration, so DART requires this to be <em class=
"code">'cyclic'</em>. <em class=
"units">[default: 'cyclic']</em></td>
</tr>
</table>
</div>
<br />
<!--==================================================================-->
 <a name="grid_nml" id="grid_nml"></a><br />
<div class="namelist">
<pre>
<em class=
"call">namelist /grid_nml/ </em> horiz_grid_opt,  vert_grid_opt,  topography_opt, &amp;
                               horiz_grid_file, vert_grid_file, topography_file
</pre></div>
<div class="indent1"><!-- Description -->
<p>The POP grid information comes in several files: horizontal grid
lat/lons in one, the vertical grid spacing in another, and the
topography (lowest valid vertical level) in a third.<br />
<br />
Here is what we can get from the (binary) horizontal grid file:</p>
<pre>
real(r8), dimension(:,:) :: ULAT,  &amp;! latitude  (radians) of U points
real(r8), dimension(:,:) :: ULON,  &amp;! longitude (radians) of U points
real(r8), dimension(:,:) :: HTN ,  &amp;! length (cm) of north edge of T box
real(r8), dimension(:,:) :: HTE ,  &amp;! length (cm) of east  edge of T box
real(r8), dimension(:,:) :: HUS ,  &amp;! length (cm) of south edge of U box
real(r8), dimension(:,:) :: HUW ,  &amp;! length (cm) of west  edge of U box
real(r8), dimension(:,:) :: ANGLE  &amp;! angle
</pre>
<p>The vertical grid file is ascii, with 3 columns/line:</p>
<pre>
cell thickness(in cm)   cell center(in m)   cell bottom(in m)</pre>
<p>Here is what we can get from the topography file:</p>
<pre>
integer, dimension(:,:), :: KMT    &amp;! k index of deepest grid cell on T grid
</pre>
<p>These must be derived or come from someplace else ...</p>
<pre>
KMU               k index of deepest grid cell on U grid
HT                real(r8) value of deepest valid T depth (in cm)
HU                real(r8) value of deepest valid U depth (in cm)
</pre>
<table border="0" cellpadding="3" width="100%">
<tr>
<th align="left">Contents</th>
<th align="left">Type</th>
<th align="left">Description</th>
</tr>
<tr><!--contents-->
<td valign="top">horiz_grid_opt, vert_grid_opt, topography_opt</td>
<!--  type  -->
<td valign="top">character(len=64)</td>
<!--descript-->
<td valign="top">switch to indicate whether or not the grids will
come from an external file or not. DART requires ALL of these to be
<em class="code">'file'</em>. <em class=
"units">[default: 'file']</em></td>
</tr>
<tr><!--contents-->
<td valign="top">horiz_grid_file</td>
<!--  type  -->
<td valign="top">character(len=100)</td>
<!--descript-->
<td valign="top">The name of the binary file containing the values
for the horizontal grid. The <strong>dimensions</strong> of the
grid are read from <em class="file">pop.r.nc</em>. It would have
been nice to include the actual grid information in the netCDF
files. <em class=
"units">[default: 'horiz_grid.gx3v5.r8ieee.le']</em></td>
</tr>
<tr><!--contents-->
<td valign="top">vert_grid_file</td>
<!--  type  -->
<td valign="top">character(len=100)</td>
<!--descript-->
<td valign="top">The name of the ASCII file containing the values
for the vertical grid. The file must contain three columns of data
pertaining to the cell thickness (in cm), the cell center (in
meters), and the cell bottom (in meters). Again, it would have been
nice to include the vertical grid information in the netCDF files.
<em class="units">[default: 'vert_grid.gx3v5']</em></td>
</tr>
<tr><!--contents-->
<td valign="top">topography_grid_file</td>
<!--  type  -->
<td valign="top">character(len=100)</td>
<!--descript-->
<td valign="top">The name of the binary file containing the values
for the topography information. The <strong>dimensions</strong> of
the grid are read from <em class="file">pop.r.nc</em>. <em class=
"units">[default: 'topography.gx3v5.r8ieee.le']</em></td>
</tr>
</table>
</div>
<br />
<!--==================================================================-->
 <a name="OtherModulesUsed" id="OtherModulesUsed"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr />
<h2>OTHER MODULES USED</h2>
<pre>
types_mod
time_manager_mod
utilities_mod
typesizes
netcdf
</pre>
<!--==================================================================-->
<!-- Note to authors. The first row of the table is different.        -->
<!--==================================================================-->
<!-- Declare all public entities ...                                  -->
<!-- duplicate public routines template as many times as necessary    -->
<!-- make sure you replace all yyyroutine?? strings                   -->
<!--==================================================================-->
<a name="Interface" id="Interface"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr />
<h2>PUBLIC INTERFACES</h2>
<p>Only a select number of interfaces used are discussed here. Each
module has its own discussion of their routines.</p>
<h3 class="indent1">Interface Routines</h3>
<table width="80%">
<tr>
<td><em class="call">use dart_pop_mod, only :</em></td>
<td><a href="#get_pop_calendar">get_pop_calendar</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#set_model_time_step">set_model_time_step</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_horiz_grid_dims">get_horiz_grid_dims</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_vert_grid_dim">get_vert_grid_dim</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#read_horiz_grid">read_horiz_grid</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#read_topography">read_topography</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#read_vert_grid">read_vert_grid</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#write_pop_namelist">write_pop_namelist</a></td>
</tr>
<tr>
<td> </td>
<td><a href=
"#get_pop_restart_filename">get_pop_restart_filename</a></td>
</tr>
</table>
<!--==================================================================-->
<h3 class="indent1">Required Interface Routines</h3>
<!--==================================================================-->
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
<a name="get_pop_calendar" id="get_pop_calendar"></a><br />
<div class="routine"><em class="call">call
get_pop_calendar(calstring)</em>
<pre>
character(len=*), intent(out) :: <em class="code">calstring</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns a string containing the type of calendar in use.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">calstring</em></td>
<td>DART/POP uses a 'gregorian' calendar.</td>
</tr>
</table>
</div>
<br />
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="set_model_time_step" id="set_model_time_step"></a><br />
<div class="routine"><em class="call">poptimestep =
set_model_time_step()</em>
<pre>
type(time_type), intent(out) :: <em class="code">poptimestep</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p><em class="code">set_model_time_step</em> returns the model time
step that was set in the <a href=
"#restart_nml">restart_nml</a><em class="code">restart_freq</em>.
This is the minimum amount of time DART thinks the POP model can
advance. Indirectly, this specifies the minimum assimilation
interval.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">poptimestep</em></td>
<td>the minimum assimilation interval</td>
</tr>
</table>
</div>
<br />
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_horiz_grid_dims" id="get_horiz_grid_dims"></a><br />
<div class="routine"><em class="call">call get_horiz_grid_dims(Nx,
Ny)</em>
<pre>
integer, intent(out) :: <em class="code">Nx, Ny</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p><em class="code">get_horiz_grid_dims</em> reads <em class=
"file">pop.r.nc</em> to determine the number of longitudes and
latitudes.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">Nx   </em></td>
<td>the length of the 'i' dimension in the POP restart file. The
number of longitudes in use.</td>
</tr>
<tr>
<td valign="top"><em class="code">Ny   </em></td>
<td>the length of the 'j' dimension in the POP restart file. The
number of latitudes in use.</td>
</tr>
</table>
</div>
<br />
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_vert_grid_dim" id="get_vert_grid_dim"></a><br />
<div class="routine"><em class="call">call get_vert_grid_dim( Nz
)</em>
<pre>
integer, intent(out) :: <em class="code">Nz</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p><em class="code">get_vert_grid_dim</em> reads <em class=
"file">pop.r.nc</em> to determine the number of vertical levels in
use.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">Nz   </em></td>
<td>the length of the 'k' dimension in the POP restart file. The
number of vertical levels in use.</td>
</tr>
</table>
</div>
<br />
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="read_horiz_grid" id="read_horiz_grid"></a><br />
<div class="routine"><em class="call">call read_horiz_grid(nx, ny,
ULAT, ULON, TLAT, TLON)</em>
<pre>
integer,                    intent(in)  :: <em class=
"code">nx, ny</em>
real(r8), dimension(nx,ny), intent(out) :: <em class=
"code">ULAT, ULON, TLAT, TLON</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p><em class="code">read_horiz_grid</em> reads the direct access
binary files containing the POP grid information. <strong>The first
record is REQUIRED to be 'ULAT', the second record is REQUIRED to
be 'ULON'.</strong></p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">nx   </em></td>
<td>The number of longitudes in the grid.</td>
</tr>
<tr>
<td valign="top"><em class="code">ny   </em></td>
<td>The number of latitudes in the grid.</td>
</tr>
<tr>
<td valign="top"><em class="code">ULAT   </em></td>
<td>The matrix of latitudes for the UVEL and VVEL variables. Units
are degrees [-90,90].</td>
</tr>
<tr>
<td valign="top"><em class="code">ULON   </em></td>
<td>The matrix of longitudes for the UVEL and VVEL variables. Units
are degrees. [0,360]</td>
</tr>
<tr>
<td valign="top"><em class="code">TLAT   </em></td>
<td>The matrix of latitudes for the SALT and TEMP variables. Units
are degrees [-90,90].</td>
</tr>
<tr>
<td valign="top"><em class="code">TLON   </em></td>
<td>The matrix of longitudes for the SALT and TEMP variables. Units
are degrees. [0,360]</td>
</tr>
</table>
</div>
<br />
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="read_topography" id="read_topography"></a><br />
<div class="routine"><em class="call">call read_topography(nx, ny,
KMT, KMU)</em>
<pre>
integer,                   intent(in)  :: <em class=
"code">nx, ny</em>
integer, dimension(nx,ny), intent(out) :: <em class=
"code">KMT, KMU</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p><em class="code">read_topography</em> reads the direct access
binary files containing the POP topography information. <strong>The
first record is REQUIRED to be 'KMT'.</strong> 'KMU' is calculated
from 'KMT'.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">nx   </em></td>
<td>The number of longitudes in the grid.</td>
</tr>
<tr>
<td valign="top"><em class="code">ny   </em></td>
<td>The number of latitudes in the grid.</td>
</tr>
<tr>
<td valign="top"><em class="code">KMT   </em></td>
<td>The matrix containing the lowest valid depth index at grid
centroids.</td>
</tr>
<tr>
<td valign="top"><em class="code">KMU   </em></td>
<td>The matrix containing the lowest valid depth index at grid
corners.</td>
</tr>
</table>
</div>
<br />
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="read_vert_grid" id="read_vert_grid"></a><br />
<div class="routine"><em class="call">call read_vert_grid(nz, ZC,
ZG)</em>
<pre>
integer,                 intent(in)  :: <em class="code">nz</em>
real(r8), dimension(nz), intent(out) :: <em class=
"code">ZC, ZG</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p><em class="code">read_vert_grid</em> reads the ASCII file
containing the information about the vertical levels. The file must
contain three columns of data pertaining to; 1) the cell thickness
(in cm),<br />
2) the cell center (in meters),<br />
and 3) the cell bottom (in meters).</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">nz   </em></td>
<td>The number of vertical levels.</td>
</tr>
<tr>
<td valign="top"><em class="code">ZC   </em></td>
<td>The depth (in meters) at the grid centers.</td>
</tr>
<tr>
<td valign="top"><em class="code">ZG   </em></td>
<td>The depth (in meters) at the grid edges.</td>
</tr>
</table>
</div>
<br />
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="write_pop_namelist" id="write_pop_namelist"></a><br />
<div class="routine"><em class="call">call
write_pop_namelist(model_time, adv_to_time)</em>
<pre>
type(time_type), intent(in)  :: <em class="code">model_time</em>
type(time_type), intent(in)  :: <em class="code">adv_to_time</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p><em class="code">write_pop_namelist</em> writes the POP namelist
<em class="code">time_manager_nml</em> with the information
necessary to advance POP to the next assimilation time. The
namelist is written to a file called <em class=
"code">pop_in.DART</em>. Presently, DART is configured to minimally
advance POP for 86400 seconds - i.e. 1 day. The forecast length
(the difference between 'model_time' and 'adv_to_time') must be an
integer number of days with the current setup. An error will result
if it is not.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class=
"code">model_time   </em></td>
<td>The 'valid' time of the current model state.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">adv_to_time   </em></td>
<td>The time of the next assimilation.</td>
</tr>
</table>
</div>
<br />
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_pop_restart_filename" id=
"get_pop_restart_filename"></a><br />
<div class="routine"><em class="call">call
get_pop_restart_filename( filename )</em>
<pre>
character(len=*), intent(out) :: <em class="code">filename</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p><em class="code">get_pop_restart_filename</em> returns the
filename containing the POP restart information. At this point the
filename is <strong>hardwired</strong> to <em class=
"file">pop.r.nc</em>, but may become more flexible in future
versions. The filename may be derived from the <em class=
"code">restart_nml</em> but is currently ignored.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class=
"code">filename   </em></td>
<td>The name of the POP restart file.</td>
</tr>
</table>
</div>
<br />
<!--==================================================================-->
<!-- Describe the Files Used by this module.                          -->
<!--==================================================================-->
 <a name="FilesUsed" id="FilesUsed"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr />
<h2>FILES</h2>
<table border="0" width="100%">
<tr>
<th align="left">filename</th>
<th align="left">purpose</th>
</tr>
<tr>
<td>pop_in</td>
<td>to read the POP namelists</td>
</tr>
<tr>
<td>pop.r.nc</td>
<td>provides grid dimensions and 'valid_time' of the model
state</td>
</tr>
<tr>
<td><em class="code">&amp;grid_nml</em> "horiz_grid_file"</td>
<td>contains the values of the horizontal grid</td>
</tr>
<tr>
<td><em class="code">&amp;grid_nml</em> "vert_grid_file"</td>
<td>contains the number and values of the vertical levels</td>
</tr>
<tr>
<td><em class="code">&amp;grid_nml</em> "topography_grid_file"</td>
<td>contains the indices of the wet/dry cells</td>
</tr>
<tr>
<td>pop_in.DART</td>
<td>to control the integration of the POP model advance</td>
</tr>
</table>
<br />
<!--==================================================================-->
<!-- Cite references, if need be.                                     -->
<!--==================================================================-->
 <a name="References" id="References"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr />
<h2>REFERENCES</h2>
<ul>
<li>none</li>
</ul>
<!--==================================================================-->
<!-- Describe all the error conditions and codes.                     -->
<!--==================================================================-->
<a name="Errors" id="Errors"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr />
<h2>ERROR CODES and CONDITIONS</h2>
<div class="errors">
<table border="1" cellspacing="1" cellpadding="10" width="100%">
<tr>
<th>Routine</th>
<th width="50%">Message</th>
<th>Comment</th>
</tr>
<tr><!-- routine -->
<td valign="top">initialize_module</td>
<!-- message -->
<td valign="top">pop_in:init_ts_file pop.r.nc not found'</td>
<!-- comment -->
<td valign="top">The POP restart file MUST be called 'pop.r.nc'.
Make a soft link if necessary.</td>
</tr>
<tr><!-- routine -->
<td valign="top">get_horiz_grid_dims</td>
<!-- message -->
<td valign="top">unable to find either 'i' or 'nlon' in file
pop.r.nc</td>
<!-- comment -->
<td valign="top">The POP restart file must contain dimensions named
either 'i' or 'nlon'.</td>
</tr>
<tr><!-- routine -->
<td valign="top">get_horiz_grid_dims</td>
<!-- message -->
<td valign="top">unable to find either 'j' or 'nlat' in file
pop.r.nc</td>
<!-- comment -->
<td valign="top">The POP restart file must contain dimensions named
either 'j' or 'nlat'.</td>
</tr>
<tr><!-- routine -->
<td valign="top">set_model_time_step</td>
<!-- message -->
<td valign="top">restart_freq_opt must be nday</td>
<!-- comment -->
<td valign="top">Pretty self-explanatory. The POP namelist must
specify the forecast length as a multiple of 'days'.</td>
</tr>
<tr><!-- routine -->
<td valign="top">write_pop_namelist</td>
<!-- message -->
<td valign="top">adv_to_time has seconds == xxx must be zero'</td>
<!-- comment -->
<td valign="top">DART is asking POP to advance to a time that is a
fraction of a day away. This should not be possible. Contact the
DART developers.</td>
</tr>
<tr><!-- routine -->
<td valign="top">write_pop_namelist</td>
<!-- message -->
<td valign="top">stop_option must be "nday"</td>
<!-- comment -->
<td valign="top">the POP <em class=
"code">time_manager_nml:stop_option</em> is not set to 'nday'. This
is required by DART.</td>
</tr>
<tr><!-- routine -->
<td valign="top">read_horiz_grid</td>
<!-- message -->
<td valign="top">pop_in:horiz_grid_file 'XYZ' not found</td>
<!-- comment -->
<td valign="top">The horizontal grid filename specified in
<em class="file">pop_in</em><em class="code">grid_nml</em> cannot
be found.</td>
</tr>
<tr><!-- routine -->
<td valign="top">calc_tpoints</td>
<!-- message -->
<td valign="top">pop_in&amp;domain_nml:ew_boundary_type 'X'
unknown</td>
<!-- comment -->
<td valign="top">The <em class="code">ew_boundary_type</em> must be
'cyclic' - until DART/POP gets tested with non-cyclic domains.</td>
</tr>
<tr><!-- routine -->
<td valign="top">read_topography</td>
<!-- message -->
<td valign="top">pop_in:topography_file 'XYZ' not found</td>
<!-- comment -->
<td valign="top">The topography file specified in <em class=
"file">pop_in</em><em class="code">grid_nml</em> cannot be
found.</td>
</tr>
<tr><!-- routine -->
<td valign="top">read_vert_grid</td>
<!-- message -->
<td valign="top">pop_in:vert_grid_file 'XYZ' not found</td>
<!-- comment -->
<td valign="top">The vertical grid file specified in <em class=
"file">pop_in</em><em class="code">grid_nml</em> cannot be
found.</td>
</tr>
<tr><!-- routine -->
<td valign="top">read_vert_grid</td>
<!-- message -->
<td valign="top">error reading depths, line 'X'</td>
<!-- comment -->
<td valign="top">The vertical grid file is corrupt or does not have
the expected three pieces of information per line.</td>
</tr>
</table>
</div>
<h2>KNOWN BUGS</h2>
<p>There are no known bugs, but there sure is a lot of dependence
on assimilating on daily boundaries - and the pop.r.nc file.</p>
<!--==================================================================-->
<!-- Describe Future Plans.                                           -->
<!--==================================================================-->
<a name="FuturePlans" id="FuturePlans"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr />
<h2>FUTURE PLANS</h2>
<p>none at this time</p>
<!--==================================================================-->
<!-- PrivateComponents                                                -->
<!--==================================================================-->
<a name="PrivateComponents" id="PrivateComponents"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr />
<h2>PRIVATE COMPONENTS</h2>
<p>N/A</p>
<!--==================================================================-->
<!-- Legalese & Metadata                                              -->
<!--==================================================================-->
<a name="Legalese" id="Legalese"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr />
<h2>Terms of Use</h2>
<p>DART software - Copyright UCAR. This open source software is
provided by UCAR, "as is", without charge, subject to all terms of
use at <a href=
"http://www.image.ucar.edu/DAReS/DART/DART_download">http://www.image.ucar.edu/DAReS/DART/DART_download</a></p>
<!--==================================================================-->
</body>
</html>
