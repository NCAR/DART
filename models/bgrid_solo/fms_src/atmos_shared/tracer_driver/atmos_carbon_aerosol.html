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
<title>Module atmos_carbon_aerosol_mod</title>
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
<h1>module atmos_carbon_aerosol_mod</h1>
<a name="HEADER" id="HEADER"></a> <!-- END HEADER -->
<a name="OVERVIEW" id="OVERVIEW"></a>
<hr>
<h2>OVERVIEW</h2>
<!-- BEGIN OVERVIEW -->
<p class="text">This code allows the implementation of black and
organic carbon tracers in the FMS framework.</p>
<!-- END OVERVIEW -->
<a name="DESCRIPTION" id="DESCRIPTION"></a> 
<!-- BEGIN DESCRIPTION -->
<div>This module presents the method of Cooke et al. (1999, 2002)
In its present implementation the black and organic carbon tracers
are from the combustion of fossil fuel. While the code here should
provide insights into the carbonaceous aerosol cycle it is provided
here more as an example of how to implement a tracer module in the
FMS infrastructure. The parameters of the model should be checked
and set to values corresponding to previous works if a user wishes
to try to reproduce those works.</div>
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
use atmos_carbon_aerosol_mod [, only:  atmos_blackc_sourcesink,<br>                                       atmos_organic_sourcesink,<br>                                       atmos_carbon_aerosol_init,<br>                                       atmos_carbon_aerosol_end ]</pre>
<dl>
<dt><a href=
"#atmos_blackc_sourcesink">atmos_blackc_sourcesink</a>:</dt>
<dd>A subroutine to calculate the source and sinks of black carbon
aerosol.</dd>
<dt><a href=
"#atmos_organic_sourcesink">atmos_organic_sourcesink</a>:</dt>
<dd>A subroutine to calculate the source and sinks of organic
carbon aerosol.</dd>
<dt><a href=
"#atmos_carbon_aerosol_init">atmos_carbon_aerosol_init</a>:</dt>
<dd>Subroutine to initialize the carbon aerosol module.</dd>
<dt><a href=
"#atmos_carbon_aerosol_end">atmos_carbon_aerosol_end</a>:</dt>
<dd>The destructor routine for the carbon aerosol module.</dd>
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
<li><a name="atmos_blackc_sourcesink" id=
"atmos_blackc_sourcesink"></a>
<h2>atmos_blackc_sourcesink</h2>
<pre>
<b>call atmos_blackc_sourcesink </b>(lon, lat, land, pwt, & black_cphob, black_cphob_dt, & black_cphil, black_cphil_dt, & Time, is, ie, js, je, kbot)</pre>
<dl>
<dt><b>DESCRIPTION</b></dt>
<dd>This routine calculates the source and sink terms for black
carbon. Simply put, the hydrophobic aerosol has sources from
emissions and sinks from dry deposition and transformation into
hydrophilic aerosol. The hydrophilic aerosol also has emission
sources and has sinks of wet and dry deposition.<br>
<br>
The following schematic shows how the black carbon scheme is
implemented.
<pre> +------------+  Trans-   +------------+
 |  Hydro-    | formation |  Hydro-    |
 |  phobic    |           |  philic    |
 |  black     |----------&gt;|  black     |
 |  carbon    |           |  carbon    |
 |            |           |            |
 +------------+           +------------+
    ^      |                ^    |   |
    |      |                |    |   |
    |      =                |    =   =
  Source  Dry            Source Dry Wet
          Dep.                  Dep Dep</pre>
The transformation time used here is 1 day, which corresponds to an
e-folding time of 1.44 days. This can be varied as necessary.</dd>
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
<td valign="top" align="left">
<tt>black_cphob   </tt></td>
<td>The array of the hydrophobic black carbon aerosol mixing
ratio<br>
   <span class="type">[real,
dimension(:,:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left">
<tt>black_cphil   </tt></td>
<td>The array of the hydrophilic black carbon aerosol mixing
ratio<br>
   <span class="type">[real,
dimension(:,:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>Time   </tt></td>
<td>Model time.<br>
   <span class="type">[type(time_type)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>is, ie, js,
je   </tt></td>
<td>Local domain boundaries.<br>
   <span class="type">[integer]</span></td>
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
<tt>black_cphob_dt   </tt></td>
<td>The array of the tendency of the hydrophobic black carbon
aerosol mixing ratio.<br>
   <span class="type">[real,
dimension(:,:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left">
<tt>black_cphil_dt   </tt></td>
<td>The array of the tendency of the hydrophilic black carbon
aerosol mixing ratio.<br>
   <span class="type">[real,
dimension(:,:,:)]</span></td>
</tr>
</table>
</dd>
<dd><br></dd>
</dl>
</li>
<li><a name="atmos_organic_sourcesink" id=
"atmos_organic_sourcesink"></a>
<h2>atmos_organic_sourcesink</h2>
<pre>
<b>call atmos_organic_sourcesink </b>(lon, lat, land, pwt, organic_carbon, organic_carbon_dt, & Time, is, ie, js, je, kbot)</pre>
<dl>
<dt><b>DESCRIPTION</b></dt>
<dd>This routine calculates the source and sink terms for organic
carbon. Simply put, the hydrophobic aerosol has sources from
emissions and sinks from dry deposition and transformation into
hydrophilic aerosol. The hydrophilic aerosol also has emission
sources and has sinks of wet and dry deposition.<br>
<br>
The following schematic shows how the organic carbon scheme is
implemented.
<pre> +------------+  Trans-   +------------+
 |  Hydro-    | formation |  Hydro-    |
 |  phobic    |           |  philic    |
 |  organic   |----------&gt;|  organic   |
 |  carbon    |           |  carbon    |
 |            |           |            |
 +------------+           +------------+
    ^      |                ^    |   |
    |      |                |    |   |
    |      =                |    =   =
  Source  Dry            Source Dry Wet
          Dep.                  Dep Dep</pre>
The transformation time used here is 2 days, which corresponds to
an e-folding time of 2.88 days. This can be varied as
necessary.<br>
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
<td valign="top" align="left">
<tt>organic_carbon   </tt></td>
<td>The array of the organic carbon aerosol mixing ratio<br>
   <span class="type">[real,
dimension(:,:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>Time   </tt></td>
<td>Model time.<br>
   <span class="type">[type(time_type)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>is, ie, js,
je   </tt></td>
<td>Local domain boundaries.<br>
   <span class="type">[integer]</span></td>
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
<tt>organic_carbon_dt   </tt></td>
<td>The array of the tendency of the organic carbon aerosol mixing
ratio.<br>
   <span class="type">[real,
dimension(:,:,:)]</span></td>
</tr>
</table>
</dd>
<dd><br></dd>
</dl>
</li>
<li><a name="atmos_carbon_aerosol_init" id=
"atmos_carbon_aerosol_init"></a>
<h2>atmos_carbon_aerosol_init</h2>
<pre>
<b>call atmos_carbon_aerosol_init </b>(lonb, latb, r, axes, Time, mask)</pre>
<dl>
<dt><b>DESCRIPTION</b></dt>
<dd>This subroutine querys the tracer manager to find the indices
for the various carbonaceous aerosol tracers. It also registers the
emission fields for diagnostic purposes.<br>
<br></dd>
<dd><br>
<br></dd>
<dt><b>INPUT</b></dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>lonb   </tt></td>
<td>The longitudes for the local domain.<br>
   <span class="type">[real,
dimension(:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>latb   </tt></td>
<td>The latitudes for the local domain.<br>
   <span class="type">[real,
dimension(:)]</span></td>
</tr>
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
<li><a name="atmos_carbon_aerosol_end" id=
"atmos_carbon_aerosol_end"></a>
<h2>atmos_carbon_aerosol_end</h2>
<pre>
<b>call atmos_carbon_aerosol_end </b>
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
<div>
<dl>
<dt>Black carbon emissions</dt>
<dd>The black carbon emission dataset is that derived in Cooke et
al. (1999) The dataset can be obtained from the contact person
above.</dd>
<dt>Organic carbon emissions</dt>
<dd>The organic carbon emission dataset is that derived in Cooke et
al. (1999) The dataset can be obtained from the contact person
above.</dd>
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
<li>Cooke, W. F. and J. J. N. Wilson, A global black carbon aerosol
model, J. Geophys. Res., 101, 19395-19409, 1996.</li>
<li>Cooke, W. F., C. Liousse, H. Cachier and J. Feichter,
Construction of a 1 x 1 fossil fuel emission dataset for
carbonaceous aerosol and implementation and radiative impact in the
ECHAM-4 model, J. Geophys. Res., 104, 22137-22162, 1999</li>
<li>Cooke, W.F., V. Ramaswamy and P. Kasibathla, A GCM study of the
global carbonaceous aerosol distribution. J. Geophys. Res., 107,
accepted, 2002</li>
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
