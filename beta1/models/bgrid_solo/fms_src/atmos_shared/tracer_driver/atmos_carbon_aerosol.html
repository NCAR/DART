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
<title>Module atmos_carbon_aerosol_mod</title>
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
<h2>module atmos_carbon_aerosol_mod</h2>
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
<b>Last Modified:</b>&nbsp;2002/06/14 16:02:12<br>
<br>
</div>
<!-- END HEADER -->
<a name="OVERVIEW"></a>
<hr>
<h4>OVERVIEW</h4>
<!-- BEGIN OVERVIEW -->
<p class="text">
   This code allows the implementation of black and organic carbon 
   tracers in the FMS framework.</p>
<!-- END OVERVIEW -->
<a name="DESCRIPTION"></a>
<!-- BEGIN DESCRIPTION -->
<div>
   This module presents the method of Cooke et al. (1999, 2002) 
   In its present implementation the black and organic carbon tracers 
   are from the combustion of fossil fuel.
   While the code here should provide insights into the carbonaceous 
   aerosol cycle it is provided here more as an example of how to implement 
   a tracer module in the FMS infrastructure. The parameters of the model 
   should be checked and set to values corresponding to previous works if 
   a user wishes to try to reproduce those works.</div>
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
<pre>use atmos_carbon_aerosol_mod [, only:  atmos_blackc_sourcesink,<br>                                       atmos_organic_sourcesink,<br>                                       atmos_carbon_aerosol_init,<br>                                       atmos_carbon_aerosol_end ]</pre>
<dl>
<dt>
<a href="#atmos_blackc_sourcesink">atmos_blackc_sourcesink</a>:</dt>
<dd>
   A subroutine to calculate the source and sinks of black carbon aerosol.</dd>
<dt>
<a href="#atmos_organic_sourcesink">atmos_organic_sourcesink</a>:</dt>
<dd>
   A subroutine to calculate the source and sinks of organic carbon aerosol.</dd>
<dt>
<a href="#atmos_carbon_aerosol_init">atmos_carbon_aerosol_init</a>:</dt>
<dd>
   Subroutine to initialize the carbon aerosol module.</dd>
<dt>
<a href="#atmos_carbon_aerosol_end">atmos_carbon_aerosol_end</a>:</dt>
<dd>
   The destructor routine for the carbon aerosol module.</dd>
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
<a name="atmos_blackc_sourcesink"></a>
<h4>atmos_blackc_sourcesink</h4>
<pre>
<b>call atmos_blackc_sourcesink </b>(lon, lat, land, pwt, &amp; black_cphob, black_cphob_dt, &amp; black_cphil, black_cphil_dt, &amp; Time, is, ie, js, je, kbot)</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>
   This routine calculates the source and sink terms for black carbon.
   Simply put, the hydrophobic aerosol has sources from emissions and 
   sinks from dry deposition and transformation into hydrophilic aerosol.
   The hydrophilic aerosol also has emission sources and has sinks of wet 
   and dry deposition.
   <br>
<br>
   The following schematic shows how the black carbon scheme 
   is implemented.<pre> +------------+  Trans-   +------------+
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
<td valign="top" align="left"><tt>black_cphob&nbsp;&nbsp;&nbsp;</tt></td><td>The array of the hydrophobic black carbon aerosol mixing ratio<br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:,:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>black_cphil&nbsp;&nbsp;&nbsp;</tt></td><td>The array of the hydrophilic black carbon aerosol mixing ratio<br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:,:,:)]</span></td>
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
<td valign="top" align="left"><tt>black_cphob_dt&nbsp;&nbsp;&nbsp;</tt></td><td>The array of the tendency of the hydrophobic black carbon aerosol mixing ratio.<br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:,:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>black_cphil_dt&nbsp;&nbsp;&nbsp;</tt></td><td>The array of the tendency of the hydrophilic black carbon aerosol mixing ratio.<br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:,:,:)]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="atmos_organic_sourcesink"></a>
<h4>atmos_organic_sourcesink</h4>
<pre>
<b>call atmos_organic_sourcesink </b>(lon, lat, land, pwt, organic_carbon, organic_carbon_dt, &amp; Time, is, ie, js, je, kbot)</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>
   This routine calculates the source and sink terms for organic carbon.
   Simply put, the hydrophobic aerosol has sources from emissions and 
   sinks from dry deposition and transformation into hydrophilic aerosol.
   The hydrophilic aerosol also has emission sources and has sinks of wet 
   and dry deposition.
   <br>
<br>
   The following schematic shows how the organic carbon scheme 
   is implemented.<pre> +------------+  Trans-   +------------+
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
   The transformation time used here is 2 days, which corresponds to an
   e-folding time of 2.88 days. This can be varied as necessary.
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
<td valign="top" align="left"><tt>organic_carbon&nbsp;&nbsp;&nbsp;</tt></td><td>The array of the organic carbon aerosol mixing ratio<br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:,:,:)]</span></td>
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
<td valign="top" align="left"><tt>organic_carbon_dt&nbsp;&nbsp;&nbsp;</tt></td><td>The array of the tendency of the organic carbon aerosol mixing ratio.<br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:,:,:)]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="atmos_carbon_aerosol_init"></a>
<h4>atmos_carbon_aerosol_init</h4>
<pre>
<b>call atmos_carbon_aerosol_init </b>(lonb, latb, r, axes, Time, mask)</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>
   This subroutine querys the tracer manager to find the indices for the 
   various carbonaceous aerosol tracers. It also registers the emission 
   fields for diagnostic purposes.
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
<a name="atmos_carbon_aerosol_end"></a>
<h4>atmos_carbon_aerosol_end</h4>
<pre>
<b>call atmos_carbon_aerosol_end </b>
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
<div>
<dl>
<dt>Black carbon emissions</dt>
<dd>
   The black carbon emission dataset is that derived in Cooke et al. (1999)
   The dataset can be obtained from the contact person above.</dd>
<dt>Organic carbon emissions</dt>
<dd>
   The organic carbon emission dataset is that derived in Cooke et al. (1999)
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
   Cooke, W. F. and J. J. N. Wilson,  A global black carbon aerosol model,
   J. Geophys. Res., 101, 19395-19409, 1996.</li>
<li>
   Cooke, W. F., C. Liousse, H. Cachier and J. Feichter, 
   Construction of a 1 x 1 fossil fuel emission dataset for carbonaceous
   aerosol and implementation and radiative impact in the ECHAM-4 model, 
   J. Geophys. Res., 104, 22137-22162, 1999</li>
<li>
   Cooke, W.F., V. Ramaswamy and P. Kasibathla, 
   A GCM study of the global carbonaceous aerosol distribution.
   J. Geophys. Res., 107, accepted, 2002</li>
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
