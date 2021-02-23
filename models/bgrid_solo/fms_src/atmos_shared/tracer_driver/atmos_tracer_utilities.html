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
<title>Module atmos_tracer_utilities_mod</title>
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
<h1>module atmos_tracer_utilities_mod</h1>
<a name="HEADER" id="HEADER"></a> <a name="OVERVIEW" id=
"OVERVIEW"></a>
<hr>
<h2>OVERVIEW</h2>
<!-- BEGIN OVERVIEW -->
<p class="text">This code provides some utility routines for
atmospheric tracers in the FMS framework.</p>
<!-- END OVERVIEW -->
<a name="DESCRIPTION" id="DESCRIPTION"></a> 
<!-- BEGIN DESCRIPTION -->
<div>This module gives utility routines which can be used to
provide consistent removal mechanisms for atmospheric tracers.<br>
<br>
In particular it provides schemes for wet and dry deposiiton that
can be easily utilized.<br>
<br></div>
<br>
<!-- END DESCRIPTION -->
<a name="OTHER MODULES USED"></a>
<hr>
<h2>OTHER MODULES USED</h2>
<!-- BEGIN OTHER MODULES USED -->
<div>
<pre>
fms_mod<br>  time_manager_mod<br>  diag_manager_mod<br>tracer_manager_mod<br> field_manager_mod<br>     constants_mod<br>  horiz_interp_mod</pre></div>
<!-- END OTHER MODULES USED -->
<!-- BEGIN PUBLIC INTERFACE -->
<a name="PUBLIC INTERFACE"></a>
<hr>
<h2>PUBLIC INTERFACE</h2>
<div>
<pre>
use atmos_tracer_utilities_mod [, only:  atmos_tracer_utilities_init,<br>                                         dry_deposition,<br>                                         wet_deposition,<br>                                         interp_emiss,<br>                                         tracer_utilities_end ]</pre>
<dl>
<dt><a href=
"#atmos_tracer_utilities_init">atmos_tracer_utilities_init</a>:</dt>
<dd>This is a routine to create and register the dry and wet
deposition fields of the tracers.</dd>
<dt><a href="#dry_deposition">dry_deposition</a>:</dt>
<dd>Routine to calculate the fraction of tracer to be removed by
dry deposition.</dd>
<dt><a href="#wet_deposition">wet_deposition</a>:</dt>
<dd>Routine to calculate the fraction of tracer removed by wet
deposition</dd>
<dt><a href="#interp_emiss">interp_emiss</a>:</dt>
<dd>A routine to interpolate emission fields of arbitrary
resolution onto the resolution of the model.</dd>
<dt><a href="#tracer_utilities_end">tracer_utilities_end</a>:</dt>
<dd>The destructor routine for the tracer utilities module.</dd>
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
<li><a name="atmos_tracer_utilities_init" id=
"atmos_tracer_utilities_init"></a>
<h2>atmos_tracer_utilities_init</h2>
<pre>
<b>call atmos_tracer_utilities_init </b>(lonb,latb, mass_axes, Time)</pre>
<dl>
<dt><b>DESCRIPTION</b></dt>
<dd>This routine creates diagnostic names for dry and wet
deposition fields of the tracers. It takes the tracer name and
appends "ddep" for the dry deposition field and "wdep" for the wet
deposition field. This names can then be entered in the diag_table
for diagnostic output of the tracer dry and wet deposition. The
module name associated with these fields in "tracers". The units of
the deposition fields are assumed to be kg/m2/s.</dd>
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
<td valign="top" align="left">
<tt>mass_axes      </tt></td>
<td>The axes relating to the tracer array.<br>
      <span class="type">[integer,
dimension(3)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>Time      </tt></td>
<td>Model time.<br>
      <span class="type">[type(time_type)]</span></td>
</tr>
</table>
</dd>
<dd><br></dd>
</dl>
</li>
<li><a name="dry_deposition" id="dry_deposition"></a>
<h2>dry_deposition</h2>
<pre>
<b>call dry_deposition </b>( n, is, js, u, v, T, pwt, pfull, u_star, landmask, dsinku, tracer, Time)</pre>
<dl>
<dt><b>DESCRIPTION</b></dt>
<dd>There are two types of dry deposition coded.<br>
<br>
1) Wind driven derived dry deposition velocity.<br>
<br>
2) Fixed dry deposition velocity.<br>
<br>
The theory behind the wind driven dry deposition velocity
calculation assumes that the deposition can be modeled as a
parallel resistance type problem.<br>
<br>
Total resistance to HNO3-type dry deposition,
<pre>       R = Ra + Rb
  resisa = aerodynamic resistance
  resisb = surface resistance (laminar layer + uptake)
         = 5/u*  [s/cm]        for neutral stability
      Vd = 1/R</pre>
For the fixed dry deposition velocity, there is no change in the
deposition velocity but the variation of the depth of the surface
layer implies that there is variation in the amount deposited.<br>
<br>
To utilize this section of code add one of the following lines as a
method for the tracer of interest in the field table.
<pre> "dry_deposition","wind_driven","surfr=XXX"
     where XXX is the total resistance defined above.

 "dry_deposition","fixed","land=XXX, sea=YYY"
     where XXX is the dry deposition velocity (m/s) over land
       and YYY is the dry deposition velocity (m/s) over sea.</pre></dd>
<dd><br>
<br></dd>
<dt><b>INPUT</b></dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>n      </tt></td>
<td>The tracer number.<br>
      <span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>is,
js      </tt></td>
<td>Start indices for array (computational indices).<br>
      <span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>u      </tt></td>
<td>U wind field.<br>
      <span class="type">[real,
dimension(:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>v      </tt></td>
<td>V wind field.<br>
      <span class="type">[real,
dimension(:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>T      </tt></td>
<td>Temperature.<br>
      <span class="type">[real,
dimension(:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>pwt      </tt></td>
<td>Pressure differential of half levels.<br>
      <span class="type">[real,
dimension(:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>pfull      </tt></td>
<td>Full pressure levels.<br>
      <span class="type">[real,
dimension(:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left">
<tt>u_star      </tt></td>
<td>Friction velocity.<br>
      <span class="type">[real,
dimension(:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left">
<tt>landmask      </tt></td>
<td>Land - sea mask.<br>
      <span class="type">[logical]</span></td>
</tr>
</table>
</dd>
<dd><br></dd>
<dt><b>OUTPUT</b></dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left">
<tt>dsinku      </tt></td>
<td>The amount of tracer in the surface layer which is dry
deposited per second.<br>
      <span class="type">[real,
dimension(:,:)]</span></td>
</tr>
</table>
</dd>
<dd><br></dd>
</dl>
</li>
<li><a name="wet_deposition" id="wet_deposition"></a>
<h2>wet_deposition</h2>
<pre>
<b>call wet_deposition </b>(n, T, pfull, phalf, rain, snow, qdt, tracer, tracer_dt, Time, cloud_param, is, js)</pre>
<dl>
<dt><b>DESCRIPTION</b></dt>
<dd>Schemes allowed here are<br>
<br>
1) Deposition removed in the same fractional amount as the modeled
precipitation rate is to a standardized precipitation rate.
Basically this scheme assumes that a fractional area of the gridbox
is affected by precipitation and that this precipitation rate is
due to a cloud of standardized cloud liquid water content. Removal
is constant throughout the column where precipitation is
occuring.<br>
<br>
2) Removal according to Henry's Law. This law states that the ratio
of the concentation in cloud water and the partial pressure in the
interstitial air is a constant. In this instance, the units for
Henry's constant are kg/L/Pa (normally it is M/L/Pa) Parameters for
a large number of species can be found at
http://www.mpch-mainz.mpg.de/~sander/res/henry.html To utilize this
section of code add one of the following lines as a method for the
tracer of interest in the field table.
<pre> "wet_deposition","henry","henry=XXX, dependence=YYY"
     where XXX is the Henry's constant for the tracer in question
       and YYY is the temperature dependence of the Henry's Law constant.

 "wet_deposition","fraction","lslwc=XXX, convlwc=YYY"
     where XXX is the liquid water content of a standard large scale cloud
       and YYY is the liquid water content of a standard convective cloud.</pre></dd>
<dd><br>
<br></dd>
<dt><b>INPUT</b></dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>n      </tt></td>
<td>Tracer number<br>
      <span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>is,
js      </tt></td>
<td>start indices for array (computational indices)<br>
      <span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>T      </tt></td>
<td>Temperature<br>
      <span class="type">[real,
dimension(:,:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>pfull      </tt></td>
<td>Full level pressure field<br>
      <span class="type">[real,
dimension(:,:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>phalf      </tt></td>
<td>Half level pressure field<br>
      <span class="type">[real,
dimension(:,:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>rain      </tt></td>
<td>Precipitation in the form of rain<br>
      <span class="type">[real,
dimension(:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>snow      </tt></td>
<td>Precipitation in the form of snow<br>
      <span class="type">[real,
dimension(:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>qdt      </tt></td>
<td>The tendency of the specific humidity due to the cloud
parametrization<br>
      <span class="type">[real,
dimension(:,:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left">
<tt>tracer      </tt></td>
<td>The tracer field<br>
      <span class="type">[real,
dimension(:,:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>Time      </tt></td>
<td>The time structure for submitting wet deposition as a
diagnostic<br>
      <span class="type">[type(time_type)]</span></td>
</tr>
<tr>
<td valign="top" align="left">
<tt>cloud_param      </tt></td>
<td>Is this a convective (convect) or large scale (lscale) cloud
parametrization?<br>
      <span class="type">[character]</span></td>
</tr>
</table>
</dd>
<dd><br></dd>
<dt><b>OUTPUT</b></dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left">
<tt>tracer_dt      </tt></td>
<td>The tendency of the tracer field due to wet deposition.<br>
      <span class="type">[real,
dimension(:,:,:)]</span></td>
</tr>
</table>
</dd>
<dd><br></dd>
</dl>
</li>
<li><a name="interp_emiss" id="interp_emiss"></a>
<h2>interp_emiss</h2>
<pre>
<b>call interp_emiss </b>(global_source, start_lon, start_lat, & lon_resol, lat_resol, data_out)</pre>
<dl>
<dt><b>DESCRIPTION</b></dt>
<dd>Routine to interpolate emission fields (or any 2D field) to the
model resolution. The local section of the global field is returned
to the local processor.</dd>
<dd><br>
<br></dd>
<dt><b>INPUT</b></dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left">
<tt>global_source      </tt></td>
<td>Global emission field.<br>
      <span class="type">[real,
dimension(:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left">
<tt>start_lon      </tt></td>
<td>Longitude of starting point of emission field (in radians).
This is the westernmost boundary of the global field.<br>
      <span class="type">[real]</span></td>
</tr>
<tr>
<td valign="top" align="left">
<tt>start_lat      </tt></td>
<td>Latitude of starting point of emission field (in radians). This
is the southern boundary of the global field.<br>
      <span class="type">[real]</span></td>
</tr>
<tr>
<td valign="top" align="left">
<tt>lon_resol      </tt></td>
<td>Longitudinal resolution of the emission data (in radians).<br>
      <span class="type">[real]</span></td>
</tr>
<tr>
<td valign="top" align="left">
<tt>lat_resol      </tt></td>
<td>Latitudinal resolution of the emission data (in radians).<br>
      <span class="type">[real]</span></td>
</tr>
</table>
</dd>
<dd><br></dd>
<dt><b>OUTPUT</b></dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left">
<tt>data_out      </tt></td>
<td>Interpolated emission field on the local PE.<br>
      <span class="type">[real,
dimension(:,:)]</span></td>
</tr>
</table>
</dd>
<dd><br></dd>
</dl>
</li>
<li><a name="tracer_utilities_end" id="tracer_utilities_end"></a>
<h2>tracer_utilities_end</h2>
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
