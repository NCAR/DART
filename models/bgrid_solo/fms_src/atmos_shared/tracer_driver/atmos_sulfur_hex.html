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
<title>Module atmos_sulfur_hex_mod</title>
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
<h1>module atmos_sulfur_hex_mod</h1>
<a name="HEADER" id="HEADER"></a> <a name="OVERVIEW" id=
"OVERVIEW"></a>
<hr>
<h2>OVERVIEW</h2>
<!-- BEGIN OVERVIEW -->
<p class="text">This code allows the implementation of sulfur
hexafluoride tracer in the FMS framework.</p>
<!-- END OVERVIEW -->
<a name="DESCRIPTION" id="DESCRIPTION"></a> 
<!-- BEGIN DESCRIPTION -->
<br>
<!-- END DESCRIPTION -->
<a name="OTHER MODULES USED"></a>
<hr>
<h2>OTHER MODULES USED</h2>
<!-- BEGIN OTHER MODULES USED -->
<div>
<pre>
fms_mod<br>    time_manager_mod<br>    diag_manager_mod<br>  tracer_manager_mod<br>   field_manager_mod<br>tracer_utilities_mod<br>       constants_mod</pre></div>
<!-- END OTHER MODULES USED -->
<!-- BEGIN PUBLIC INTERFACE -->
<a name="PUBLIC INTERFACE"></a>
<hr>
<h2>PUBLIC INTERFACE</h2>
<div>
<pre>
use atmos_sulfur_hex_mod [, only:  atmos_sf6_sourcesink,<br>                                   atmos_sulfur_hex_init,<br>                                   sulfur_hex_end ]</pre>
<dl>
<dt><a href="#atmos_sf6_sourcesink">atmos_sf6_sourcesink</a>:</dt>
<dd>A routine to calculate the sources and sinks of sulfur
hexafluoride.</dd>
<dt><a href=
"#atmos_sulfur_hex_init">atmos_sulfur_hex_init</a>:</dt>
<dd>The constructor routine for the sulfur hexafluoride
module.</dd>
<dt><a href="#sulfur_hex_end">sulfur_hex_end</a>:</dt>
<dd>The destructor routine for the sulfur hexafluoride module.</dd>
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
<li><a name="atmos_sf6_sourcesink" id="atmos_sf6_sourcesink"></a>
<h2>atmos_sf6_sourcesink</h2>
<pre>
<b>call atmos_sf6_sourcesink </b>(lon, lat, land, pwt, sf6, sf6_dt, Time, is, ie, js, je, kbot)</pre>
<dl>
<dt><b>DESCRIPTION</b></dt>
<dd>A routine to calculate the sources and sinks of sulfur
hexafluoride.</dd>
<dd><br>
<br></dd>
<dt><b>INPUT</b></dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>lon      </tt></td>
<td>Longitude of the centre of the model gridcells.<br>
      <span class="type">[real,
dimension(:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>lat      </tt></td>
<td>Latitude of the centre of the model gridcells.<br>
      <span class="type">[real,
dimension(:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>land      </tt></td>
<td>Land/sea mask.<br>
      <span class="type">[real,
dimension(:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>pwt      </tt></td>
<td>The pressure weighting array. = dP/grav<br>
      <span class="type">[real,
dimension(:,:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>sf6      </tt></td>
<td>The array of the sulfur hexafluoride mixing ratio.<br>
      <span class="type">[real,
dimension(:,:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>Time      </tt></td>
<td>Model time.<br>
      <span class="type">[type(time_type)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>is, ie, js,
je      </tt></td>
<td>Local domain boundaries.<br>
      <span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>kbot      </tt></td>
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
<tt>sf6_dt      </tt></td>
<td>The array of the tendency of the sulfur hexafluoride mixing
ratio.<br>
      <span class="type">[real,
dimension(:,:,:)]</span></td>
</tr>
</table>
</dd>
<dd><br></dd>
</dl>
</li>
<li><a name="atmos_sulfur_hex_init" id="atmos_sulfur_hex_init"></a>
<h2>atmos_sulfur_hex_init</h2>
<pre>
<b>call atmos_sulfur_hex_init </b>(lonb, latb, r, axes, Time, mask)</pre>
<dl>
<dt><b>DESCRIPTION</b></dt>
<dd>A routine to initialize the sulfur hexafluoride module.</dd>
<dd><br>
<br></dd>
<dt><b>INPUT</b></dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>lonb      </tt></td>
<td>The longitudes for the local domain.<br>
      <span class="type">[real,
dimension(:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>latb      </tt></td>
<td>The latitudes for the local domain.<br>
      <span class="type">[real,
dimension(:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>mask      </tt></td>
<td>optional mask (0. or 1.) that designates which grid points are
above (=1.) or below (=0.) the ground dimensioned as
(nlon,nlat,nlev).<br>
      <span class="type">[real, optional,
dimension(:,:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>Time      </tt></td>
<td>Model time.<br>
      <span class="type">[type(time_type)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>axes      </tt></td>
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
<td valign="top" align="left"><tt>r      </tt></td>
<td>Tracer fields dimensioned as (nlon,nlat,nlev,ntrace).<br>
      <span class="type">[real,
dimension(:,:,:,:)]</span></td>
</tr>
</table>
</dd>
<dd><br></dd>
</dl>
</li>
<li><a name="sulfur_hex_end" id="sulfur_hex_end"></a>
<h2>sulfur_hex_end</h2>
<pre>
<b>call sulfur_hex_end </b>
</pre>
<dl>
<dt><b>DESCRIPTION</b></dt>
<dd>This subroutine is the exit routine for the sulfur hexafluoride
module.</dd>
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
<div>
<dl>
<dt>Sulfur hexaflouride emissions</dt>
<dd>Monthly.emissions contains the estimated global emission rate
of SF6 in Gg/yr for 62 months between December 1988 and January
1994, inclusive. These are based on the annual estimates of Levin
and Hesshaimer (submitted), and have been linearly interpolated to
monthly values. The last half of 1993 has been extrapolated using
the trend for the previous 12 months.<br>
<br>
The dataset can be obtained from the contact person above.</dd>
</dl>
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
<div>
<ol>
<li>Levin, I. and V. Hessahimer: Refining of atmospheric transport
model entries by the globally observed passive tracer distributions
of 85Krypton and Sulfur Hexafluoride (SF6). Submitted to the
Journal of Geophysical Research.</li>
</ol>
</div>
<br>
<!-- END REFERENCES -->
<a name="COMPILER SPECIFICS"></a>
<hr>
<h2>COMPILER SPECIFICS</h2>
<!-- BEGIN COMPILER SPECIFICS -->
<div>None.</div>
<br>
<!-- END COMPILER SPECIFICS -->
<a name="PRECOMPILER OPTIONS"></a>
<hr>
<h2>PRECOMPILER OPTIONS</h2>
<!-- BEGIN PRECOMPILER OPTIONS -->
<div>None.</div>
<br>
<!-- END PRECOMPILER OPTIONS -->
<a name="LOADER OPTIONS"></a>
<hr>
<h2>LOADER OPTIONS</h2>
<!-- BEGIN LOADER -->
<div>None.<br>
<br></div>
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
