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
<title>Module atmos_tracer_driver_mod</title>
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
<h2>module atmos_tracer_driver_mod</h2>
<a name="HEADER"></a>
<!-- BEGIN HEADER -->
<div>
<b>Contact:</b>&nbsp;<a href="mailto:wfc@gfdl.noaa.gov">
   William Cooke</a>
<br>
<b>Reviewers:</b>&nbsp;<a href="mailto:mjh@gfdl.noaa.gov">
   Matt Harrison</a>,&nbsp;
    <a href="mailto:bw@gfdl.noaa.gov">
   Bruce Wyman</a>
<br>
<b>Change History: </b>&nbsp;<a href="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/">WebCVS Log</a>
<br>
<b>Last Modified:</b>&nbsp;2002/06/14 18:32:06<br>
<br>
</div>
<!-- END HEADER -->
<a name="OVERVIEW"></a>
<hr>
<h4>OVERVIEW</h4>
<!-- BEGIN OVERVIEW -->
<p class="text">
   This code allows the user to easily add tracers to the FMS framework.</p>
<!-- END OVERVIEW -->
<a name="DESCRIPTION"></a>
<!-- BEGIN DESCRIPTION -->
<div>
   This code allows a user to easily implement tracer code in the FMS
   framework.  The tracer and tracer tendency arrays are supplied along
   with longtitude,  latitude, wind, temperature, and pressure
   information which allows a  user to implement sources and sinks of the
   tracer which depend on these parameters.
   <br>
<br>
   In the following example, radon being implemented in the atmosphere
   will be used as an example of how to implement a tracer in the FMS
   framework.
   <br>
<br>
   Within the global scope of tracer_driver_mod a use
   statement should be inserted for each tracer to be added.<pre>      use radon_mod, only : radon_sourcesink, radon_init, radon_end</pre>
   An integer parameter, which  will be used as an identifier for the
   tracer, should be assigned.<pre>      integer :: nradon</pre>
   Within tracer_driver_init a call to the tracer manager is needed in
   order  to identify which tracer it has set the tracer as.<pre>      nradon = get_tracer_index(MODEL_ATMOS,'radon')</pre>
   Here MODEL_ATMOS is a parameter defined in field_manager. 
   'radon' is the name of the tracer within the field_table.
   <br>
<br>
   If the tracer exists then the integer returned will be positive and it
   can be used to call the initialization routines for the individual
   tracers.<pre>      if (nradon &gt; 0) then
           call radon_init(Argument list)
      endif</pre>
   Within tracer_driver the user can also use the identifier to surround
   calls to the source-sink routines for the tracer of interest.
   <br>
<br>
<pre>      if (nradon &gt; 0 .and. nradon &lt;= nt) then
          call radon_sourcesink (Argument list)
          rdt(:,:,:,nradon)=rdt(:,:,:,nradon)+rtnd(:,:,:)
      endif</pre>
   It is the users responsibility to add the tendency generated by the
   sourcesink routine.
   <br>
<br>
   Within tracer_driver_end the user can add calls to the
   terminators for the appropriate source sink routines.
   <br>
<br>
<pre>      call radon_end</pre>
   This may simply be a deallocation statement or a routine to send
   output to the logfile stating that the termination routine has been
   called.</div>
<br>
<!-- END DESCRIPTION -->
<a name="OTHER MODULES USED"></a>
<hr>
<h4>OTHER MODULES USED</h4>
<!-- BEGIN OTHER MODULES USED -->
<div>
<pre>                 fms_mod<br>        time_manager_mod<br>      tracer_manager_mod<br>       field_manager_mod<br>    tracer_utilities_mod<br>           constants_mod<br>         atmos_radon_mod<br>atmos_carbon_aerosol_mod<br>    atmos_sulfur_hex_mod</pre>
</div>
<!-- END OTHER MODULES USED -->
<!-- BEGIN PUBLIC INTERFACE -->
<a name="PUBLIC INTERFACE"></a>
<hr>
<h4>PUBLIC INTERFACE</h4>
<div>
<pre>use atmos_tracer_driver_mod [, only:  atmos_tracer_driver,<br>                                      atmos_tracer_driver_init,<br>                                      atmos_tracer_driver_end ]</pre>
<dl>
<dt>
<a href="#atmos_tracer_driver">atmos_tracer_driver</a>:</dt>
<dd>
   A routine which allows tracer code to be called.</dd>
<dt>
<a href="#atmos_tracer_driver_init">atmos_tracer_driver_init</a>:</dt>
<dd>
   Subroutine to initialize the tracer driver module.</dd>
<dt>
<a href="#atmos_tracer_driver_end">atmos_tracer_driver_end</a>:</dt>
<dd>
   Subroutine to terminate the tracer driver module.</dd>
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
<a name="atmos_tracer_driver"></a>
<h4>atmos_tracer_driver</h4>
<pre>
<b>call atmos_tracer_driver </b>(is, ie, js, je, Time, lon, lat, land, phalf, pfull, r, &amp; u, v, t, q, u_star, rdt, rm, rdiag, kbot)</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>
   This subroutine calls the source sink routines for atmospheric
   tracers. This is the interface between the dynamical core of the 
   model and the tracer code. It should supply all the necessary 
   information to a user that they need in order to calculate the 
   tendency of that tracer with respect to emissions or chemical losses.
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
<td valign="top" align="left"><tt>is, ie, js, je&nbsp;&nbsp;&nbsp;</tt></td><td>Local domain boundaries.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>Time&nbsp;&nbsp;&nbsp;</tt></td><td>Model time.<br>&nbsp;&nbsp;&nbsp;<span class="type">[type(time_type)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>lon&nbsp;&nbsp;&nbsp;</tt></td><td>Longitude of the centre of the model gridcells<br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>lat&nbsp;&nbsp;&nbsp;</tt></td><td>Latitude of the centre of the model gridcells<br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>land&nbsp;&nbsp;&nbsp;</tt></td><td>Land/sea mask.<br>&nbsp;&nbsp;&nbsp;<span class="type">[logical, dimension(:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>phalf&nbsp;&nbsp;&nbsp;</tt></td><td>Pressures on the model half levels.<br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:,:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>pfull&nbsp;&nbsp;&nbsp;</tt></td><td>Pressures on the model full levels.<br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:,:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>r&nbsp;&nbsp;&nbsp;</tt></td><td>The tracer array in the component model.<br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:,:,:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>u&nbsp;&nbsp;&nbsp;</tt></td><td>Zonal wind speed.<br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:,:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>v&nbsp;&nbsp;&nbsp;</tt></td><td>Meridonal wind speed.<br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:,:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>t&nbsp;&nbsp;&nbsp;</tt></td><td>Temperature.<br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:,:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>q&nbsp;&nbsp;&nbsp;</tt></td><td>Specific humidity. This may also be accessible as a
       portion of the tracer array.<br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:,:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>u_star&nbsp;&nbsp;&nbsp;</tt></td><td>Friction velocity :: 
            The magnitude of the wind stress is density*(ustar**2)
            The drag coefficient for momentum is u_star**2/(u**2+v**2)<br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>rm&nbsp;&nbsp;&nbsp;</tt></td><td>The tracer array in the component model for the previous timestep.<br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:,:,:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>kbot&nbsp;&nbsp;&nbsp;</tt></td><td>Integer array describing which model layer intercepts the surface.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer, optional, dimension(:,:)]</span></td>
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
<td valign="top" align="left"><tt>rdt&nbsp;&nbsp;&nbsp;</tt></td><td>The tendency of the tracer array in the compenent
         model. The tendency due to sources and sinks computed
         in the individual tracer routines should be added to
         this array before exiting tracer_driver.<br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:,:,:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>rdiag&nbsp;&nbsp;&nbsp;</tt></td><td>The array of diagnostic tracers. As these may be changed within the
           tracer routines for diagnostic purposes, they need to be writable.<br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:,:,:,:)]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="atmos_tracer_driver_init"></a>
<h4>atmos_tracer_driver_init</h4>
<pre>
<b>call atmos_tracer_driver_init </b>(lonb,latb, r, mask, axes, Time)</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>
   The purpose of the arguments here are for passing on to the individual
   tracer code. The user may wish to provide initial values which can be
   implemented in the initialization part of the tracer code. Remember that
   the tracer manager will provide a simple fixed or exponential profile if
   the user provides data for this within the field table. However if a more
   complicated profile is required then it should be set up in the
   initialization section of the user tracer code.</dd>
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
<a name="atmos_tracer_driver_end"></a>
<h4>atmos_tracer_driver_end</h4>
<pre>
<b>call atmos_tracer_driver_end </b>
</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>
   Termination routine for tracer_driver. It should also call
   the destructors for the individual tracer routines.</dd>
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
<div>
<dl>
<dt>
<b>FATAL in atmos_tracer_driver</b>
</dt>
<dd>
<span class="errmsg">tracer_driver_init must be called first.</span>
</dd>
<dd>
   Tracer_driver_init needs to be called before tracer_driver.</dd>
</dl>
<br>
</div>
<!-- END ERROR MESSAGES -->
<hr>
<div align="right">
<font size="-1"><a href="#TOP">top</a></font>
</div>
</body>
</html>
