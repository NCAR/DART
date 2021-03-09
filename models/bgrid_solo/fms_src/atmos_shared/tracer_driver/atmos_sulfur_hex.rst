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
<title>Module atmos_sulfur_hex_mod</title>
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
<h2>module atmos_sulfur_hex_mod</h2>
<a name="HEADER"></a>
<!-- BEGIN HEADER -->
<div>
<b>Contact:</b>&nbsp;<a href="mailto:jbg@gfdl.noaa.gov">
   Jeff Greenblatt</a>
<br>
<b>Reviewers:</b>&nbsp;<a href="mailto:wfc@gfdl.noaa.gov">
   William Cooke</a>
<br>
<b>Change History: </b>&nbsp;<a href="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/">WebCVS Log</a>
<br>
<b>Last Modified:</b>&nbsp;2002/06/14 16:58:21<br>
<br>
</div>
<!-- END HEADER -->
<a name="OVERVIEW"></a>
<hr>
<h4>OVERVIEW</h4>
<!-- BEGIN OVERVIEW -->
<p class="text">
   This code allows the implementation of sulfur hexafluoride
   tracer in the FMS framework.</p>
<!-- END OVERVIEW -->
<a name="DESCRIPTION"></a>
<!-- BEGIN DESCRIPTION -->
<div>

</div>
<br>
<!-- END DESCRIPTION -->
<a name="OTHER MODULES USED"></a>
<hr>
<h4>OTHER MODULES USED</h4>
<!-- BEGIN OTHER MODULES USED -->
<div>
<pre>             fms_mod<br>    time_manager_mod<br>    diag_manager_mod<br>  tracer_manager_mod<br>   field_manager_mod<br>tracer_utilities_mod<br>       constants_mod</pre>
</div>
<!-- END OTHER MODULES USED -->
<!-- BEGIN PUBLIC INTERFACE -->
<a name="PUBLIC INTERFACE"></a>
<hr>
<h4>PUBLIC INTERFACE</h4>
<div>
<pre>use atmos_sulfur_hex_mod [, only:  atmos_sf6_sourcesink,<br>                                   atmos_sulfur_hex_init,<br>                                   sulfur_hex_end ]</pre>
<dl>
<dt>
<a href="#atmos_sf6_sourcesink">atmos_sf6_sourcesink</a>:</dt>
<dd>
   A routine to calculate the sources and sinks of sulfur hexafluoride.</dd>
<dt>
<a href="#atmos_sulfur_hex_init">atmos_sulfur_hex_init</a>:</dt>
<dd>
   The constructor routine for the sulfur hexafluoride module.</dd>
<dt>
<a href="#sulfur_hex_end">sulfur_hex_end</a>:</dt>
<dd>
   The destructor routine for the sulfur hexafluoride module.</dd>
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
<a name="atmos_sf6_sourcesink"></a>
<h4>atmos_sf6_sourcesink</h4>
<pre>
<b>call atmos_sf6_sourcesink </b>(lon, lat, land, pwt, sf6, sf6_dt, Time, is, ie, js, je, kbot)</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>
   A routine to calculate the sources and sinks of sulfur hexafluoride.</dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>lon&nbsp;&nbsp;&nbsp;</tt></td><td>Longitude of the centre of the model gridcells.<br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>lat&nbsp;&nbsp;&nbsp;</tt></td><td>Latitude of the centre of the model gridcells.<br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>land&nbsp;&nbsp;&nbsp;</tt></td><td>Land/sea mask.<br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>pwt&nbsp;&nbsp;&nbsp;</tt></td><td>The pressure weighting array. = dP/grav<br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:,:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>sf6&nbsp;&nbsp;&nbsp;</tt></td><td>The array of the sulfur hexafluoride mixing ratio.<br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:,:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>Time&nbsp;&nbsp;&nbsp;</tt></td><td>Model time.<br>&nbsp;&nbsp;&nbsp;<span class="type">[type(time_type)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>is, ie, js, je&nbsp;&nbsp;&nbsp;</tt></td><td>Local domain boundaries.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
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
<td valign="top" align="left"><tt>sf6_dt&nbsp;&nbsp;&nbsp;</tt></td><td>The array of the tendency of the sulfur hexafluoride mixing ratio.<br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:,:,:)]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="atmos_sulfur_hex_init"></a>
<h4>atmos_sulfur_hex_init</h4>
<pre>
<b>call atmos_sulfur_hex_init </b>(lonb, latb, r, axes, Time, mask)</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>
   A routine to initialize the sulfur hexafluoride module.</dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>lonb&nbsp;&nbsp;&nbsp;</tt></td><td>The longitudes for the local domain.<br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>latb&nbsp;&nbsp;&nbsp;</tt></td><td>The latitudes for the local domain.<br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:)]</span></td>
</tr>
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
<a name="sulfur_hex_end"></a>
<h4>sulfur_hex_end</h4>
<pre>
<b>call sulfur_hex_end </b>
</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>
   This subroutine is the exit routine for the sulfur hexafluoride module.</dd>
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
<div>
<dl>
<dt>Sulfur hexaflouride emissions</dt>
<dd>
   Monthly.emissions contains the estimated global emission rate of SF6 in
   Gg/yr for 62 months between December 1988 and January 1994, inclusive.
   These are based on the annual estimates of Levin and Hesshaimer
   (submitted), and have been linearly interpolated to monthly values. The
   last half of 1993 has been extrapolated using the trend for the previous 12
   months. 
   <br>
<br>
   The dataset can be obtained from the contact person above.</dd>
</dl>
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
<ol>
<li>
   Levin, I. and V. Hessahimer: Refining of atmospheric
   transport model entries by the globally observed passive tracer
   distributions of 85Krypton and Sulfur Hexafluoride (SF6). Submitted to the
   Journal of Geophysical Research.</li>
</ol>
</div>
<br>
<!-- END REFERENCES -->
<a name="COMPILER SPECIFICS"></a>
<hr>
<h4>COMPILER SPECIFICS</h4>
<!-- BEGIN COMPILER SPECIFICS -->
<div>
        None.
      </div>
<br>
<!-- END COMPILER SPECIFICS -->
<a name="PRECOMPILER OPTIONS"></a>
<hr>
<h4>PRECOMPILER OPTIONS</h4>
<!-- BEGIN PRECOMPILER OPTIONS -->
<div>
        None.
      </div>
<br>
<!-- END PRECOMPILER OPTIONS -->
<a name="LOADER OPTIONS"></a>
<hr>
<h4>LOADER OPTIONS</h4>
<!-- BEGIN LOADER -->
<div>None.<br>
<br>
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
