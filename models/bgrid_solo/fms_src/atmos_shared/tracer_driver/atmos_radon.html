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
<title>Module atmos_radon_mod</title>
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
<h1>module atmos_radon_mod</h1>
<a name="HEADER" id="HEADER"></a> <a name="OVERVIEW" id=
"OVERVIEW"></a>
<hr>
<h2>OVERVIEW</h2>
<!-- BEGIN OVERVIEW -->
<p class="text">This code allows the implementation of an extremely
simplified radon tracer in the FMS framework.<br>
<br>
It should be taken as the implementation of a very simple tracer
which bears some characteristics of radon.</p>
<!-- END OVERVIEW -->
<a name="DESCRIPTION" id="DESCRIPTION"></a> 
<!-- BEGIN DESCRIPTION -->
<div>This module presents an implementation of a tracer. It should
be taken as representing radon only in a rudimentary manner.</div>
<br>
<!-- END DESCRIPTION -->
<a name="OTHER MODULES USED"></a>
<hr>
<h2>OTHER MODULES USED</h2>
<!-- BEGIN OTHER MODULES USED -->
<div>
<pre>
fms_mod<br>    time_manager_mod<br>    diag_manager_mod<br>  tracer_manager_mod<br>   field_manager_mod<br>tracer_utilities_mod</pre></div>
<!-- END OTHER MODULES USED -->
<!-- BEGIN PUBLIC INTERFACE -->
<a name="PUBLIC INTERFACE"></a>
<hr>
<h2>PUBLIC INTERFACE</h2>
<div>
<pre>
use atmos_radon_mod [, only:  atmos_radon_sourcesink,<br>                              atmos_radon_init,<br>                              atmos_radon_end ]</pre>
<dl>
<dt><a href=
"#atmos_radon_sourcesink">atmos_radon_sourcesink</a>:</dt>
<dd>The routine that calculate the sources and sinks of radon.</dd>
<dt><a href="#atmos_radon_init">atmos_radon_init</a>:</dt>
<dd>The constructor routine for the radon module.</dd>
<dt><a href="#atmos_radon_end">atmos_radon_end</a>:</dt>
<dd>The destructor routine for the radon module.</dd>
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
<li><a name="atmos_radon_sourcesink" id=
"atmos_radon_sourcesink"></a>
<h2>atmos_radon_sourcesink</h2>
<pre>
<b>call atmos_radon_sourcesink </b>(lon, lat, land, pwt, radon, radon_dt, Time, kbot)</pre>
<dl>
<dt><b>DESCRIPTION</b></dt>
<dd>This is a very rudimentary implementation of radon.<br>
<br>
It is assumed that the Rn222 flux is 3.69e-21 kg/m*m/sec over land
for latitudes &lt; 60N<br>
<br>
Between 60N and 70N the source = source * .5<br>
<br>
Rn222 has a half-life time of 3.83 days, which corresponds to an
e-folding time of 5.52 days.<br>
<br></dd>
<dd><br>
<br></dd>
<dt><b>INPUT</b></dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>lon   </tt></td>
<td>Longitude of the centre of the model gridcells<br>
   <span class="type">[real,
dimension(:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>lat   </tt></td>
<td>Latitude of the centre of the model gridcells<br>
   <span class="type">[real,
dimension(:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>land   </tt></td>
<td>Land/sea mask.<br>
   <span class="type">[real,
dimension(:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>pwt   </tt></td>
<td>The pressure weighting array. = dP/grav<br>
   <span class="type">[real,
dimension(:,:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>radon   </tt></td>
<td>The array of the radon mixing ratio.<br>
   <span class="type">[real,
dimension(:,:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>Time   </tt></td>
<td>Model time.<br>
   <span class="type">[type(time_type)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>kbot   </tt></td>
<td>Integer array describing which model layer intercepts the
surface.<br>
   <span class="type">[integer, optional,
dimension(:,:)]</span></td>
</tr>
</table>
</dd>
<dd><br></dd>
<dt><b>OUTPUT</b></dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left">
<tt>radon_dt   </tt></td>
<td>The array of the tendency of the radon mixing ratio.<br>
   <span class="type">[real,
dimension(:,:,:)]</span></td>
</tr>
</table>
</dd>
<dd><br></dd>
</dl>
</li>
<li><a name="atmos_radon_init" id="atmos_radon_init"></a>
<h2>atmos_radon_init</h2>
<pre>
<b>call atmos_radon_init </b>
</pre>
<dl>
<dt><b>DESCRIPTION</b></dt>
<dd>A routine to initialize the radon module.</dd>
<dd><br>
<br></dd>
<dt><b>INPUT</b></dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>mask   </tt></td>
<td>optional mask (0. or 1.) that designates which grid points are
above (=1.) or below (=0.) the ground dimensioned as
(nlon,nlat,nlev).<br>
   <span class="type">[real, optional,
dimension(:,:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>Time   </tt></td>
<td>Model time.<br>
   <span class="type">[type(time_type)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>axes   </tt></td>
<td>The axes relating to the tracer array dimensioned as (nlon,
nlat, nlev, ntime)<br>
   <span class="type">[integer,
dimension(4)]</span></td>
</tr>
</table>
</dd>
<dd><br></dd>
<dt><b>INPUT/OUTPUT</b></dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>r   </tt></td>
<td>Tracer fields dimensioned as (nlon,nlat,nlev,ntrace).<br>
   <span class="type">[real,
dimension(:,:,:,:)]</span></td>
</tr>
</table>
</dd>
<dd><br></dd>
</dl>
</li>
<li><a name="atmos_radon_end" id="atmos_radon_end"></a>
<h2>atmos_radon_end</h2>
<pre>
<b>call atmos_radon_end </b>
</pre>
<dl>
<dt><b>DESCRIPTION</b></dt>
<dd>This subroutine writes the version name to logfile and
exits.</dd>
<dd><br>
<br></dd>
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
<hr>
<div align="right"><a href="#TOP"><font size=
"-1">top</font></a></div>
</body>
</html>
