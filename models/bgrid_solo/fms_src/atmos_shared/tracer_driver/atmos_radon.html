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
<title>Module atmos_radon_mod</title>
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
<h2>module atmos_radon_mod</h2>
<a name="HEADER"></a>
<!-- BEGIN HEADER -->
<div>
<b>Contact:</b>&nbsp;<a href="mailto:wfc@gfdl.noaa.gov">
   William Cooke</a>
<br>
<b>Reviewers:</b>&nbsp;<a href="mailto:lwh@gfdl.noaa.gov">
   Larry Horowitz</a>
<br>
<b>Change History: </b>&nbsp;<a href="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/">WebCVS Log</a>
<br>
<b>Last Modified:</b>&nbsp;2002/06/07 20:01:39<br>
<br>
</div>
<!-- END HEADER -->
<a name="OVERVIEW"></a>
<hr>
<h4>OVERVIEW</h4>
<!-- BEGIN OVERVIEW -->
<p class="text">
   This code allows the implementation of an extremely simplified 
   radon tracer in the FMS framework.
   <br>
<br>
   It should be taken as the implementation of a very simple tracer 
   which bears some characteristics of radon.</p>
<!-- END OVERVIEW -->
<a name="DESCRIPTION"></a>
<!-- BEGIN DESCRIPTION -->
<div>
   This module presents an implementation of a tracer.
   It should be taken as representing radon only in a rudimentary manner.</div>
<br>
<!-- END DESCRIPTION -->
<a name="OTHER MODULES USED"></a>
<hr>
<h4>OTHER MODULES USED</h4>
<!-- BEGIN OTHER MODULES USED -->
<div>
<pre>             fms_mod<br>    time_manager_mod<br>    diag_manager_mod<br>  tracer_manager_mod<br>   field_manager_mod<br>tracer_utilities_mod</pre>
</div>
<!-- END OTHER MODULES USED -->
<!-- BEGIN PUBLIC INTERFACE -->
<a name="PUBLIC INTERFACE"></a>
<hr>
<h4>PUBLIC INTERFACE</h4>
<div>
<pre>use atmos_radon_mod [, only:  atmos_radon_sourcesink,<br>                              atmos_radon_init,<br>                              atmos_radon_end ]</pre>
<dl>
<dt>
<a href="#atmos_radon_sourcesink">atmos_radon_sourcesink</a>:</dt>
<dd>
   The routine that calculate the sources and sinks of radon.</dd>
<dt>
<a href="#atmos_radon_init">atmos_radon_init</a>:</dt>
<dd>
   The constructor routine for the radon module.</dd>
<dt>
<a href="#atmos_radon_end">atmos_radon_end</a>:</dt>
<dd>
   The destructor routine for the radon module.</dd>
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
<a name="atmos_radon_sourcesink"></a>
<h4>atmos_radon_sourcesink</h4>
<pre>
<b>call atmos_radon_sourcesink </b>(lon, lat, land, pwt, radon, radon_dt, Time, kbot)</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>
   This is a very rudimentary implementation of radon.
   <br>
<br>
   It is assumed that the Rn222 flux is 3.69e-21 kg/m*m/sec over land 
   for latitudes &lt; 60N
   <br>
<br>
   Between 60N and 70N the source  = source * .5
   <br>
<br>
   Rn222 has a half-life time of 3.83 days, which corresponds to an 
   e-folding time of 5.52 days.
   <br>
<br>
</dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>lon&nbsp;&nbsp;&nbsp;</tt></td><td>Longitude of the centre of the model gridcells<br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>lat&nbsp;&nbsp;&nbsp;</tt></td><td>Latitude of the centre of the model gridcells<br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>land&nbsp;&nbsp;&nbsp;</tt></td><td>Land/sea mask.<br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>pwt&nbsp;&nbsp;&nbsp;</tt></td><td>The pressure weighting array. = dP/grav<br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:,:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>radon&nbsp;&nbsp;&nbsp;</tt></td><td>The array of the radon mixing ratio.<br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:,:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>Time&nbsp;&nbsp;&nbsp;</tt></td><td>Model time.<br>&nbsp;&nbsp;&nbsp;<span class="type">[type(time_type)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>kbot&nbsp;&nbsp;&nbsp;</tt></td><td>Integer array describing which model layer intercepts the surface.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer, optional, dimension(:,:)]</span></td>
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
<td valign="top" align="left"><tt>radon_dt&nbsp;&nbsp;&nbsp;</tt></td><td>The array of the tendency of the radon mixing ratio.<br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:,:,:)]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="atmos_radon_init"></a>
<h4>atmos_radon_init</h4>
<pre>
<b>call atmos_radon_init </b>
</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>
   A routine to initialize the radon module.</dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>mask&nbsp;&nbsp;&nbsp;</tt></td><td>optional mask (0. or 1.) that designates which grid points
          are above (=1.) or below (=0.) the ground dimensioned as
          (nlon,nlat,nlev).<br>&nbsp;&nbsp;&nbsp;<span class="type">[real, optional, dimension(:,:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>Time&nbsp;&nbsp;&nbsp;</tt></td><td>Model time.<br>&nbsp;&nbsp;&nbsp;<span class="type">[type(time_type)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>axes&nbsp;&nbsp;&nbsp;</tt></td><td>The axes relating to the tracer array dimensioned as
          (nlon, nlat, nlev, ntime)<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer, dimension(4)]</span></td>
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
<td valign="top" align="left"><tt>r&nbsp;&nbsp;&nbsp;</tt></td><td>Tracer fields dimensioned as (nlon,nlat,nlev,ntrace).<br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:,:,:,:)]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="atmos_radon_end"></a>
<h4>atmos_radon_end</h4>
<pre>
<b>call atmos_radon_end </b>
</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>
   This subroutine writes the version name to logfile and exits.</dd>
<br>
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
<hr>
<div align="right">
<font size="-1"><a href="#TOP">top</a></font>
</div>
</body>
</html>
