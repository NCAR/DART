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
<title>Module constants_mod</title>
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
<h2>module constants_mod</h2>
<a name="HEADER"></a>
<!-- BEGIN HEADER -->
<div>
<b>Contact:</b>&nbsp;<a href="mailto:bw@gfdl.noaa.gov">   Bruce Wyman </a>
<br>
<b>Reviewers:</b>&nbsp;<br>
<b>Change History: </b>&nbsp;<a href="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/">WebCVS Log</a>
<br>
<b>Last Modified:</b>&nbsp;2002/03/01 23:43:42<br>
<br>
</div>
<!-- END HEADER -->
<a name="OVERVIEW"></a>
<hr>
<h4>OVERVIEW</h4>
<!-- BEGIN OVERVIEW -->
<p class="text">   Defines useful constants for Earth in mks units. </p>
<!-- END OVERVIEW -->
<a name="DESCRIPTION"></a>
<!-- BEGIN DESCRIPTION -->
<div>   Constants are defined as real parameters.They are accessed through
   the "use" statement. While the local name of constant may be changed,
   their values can not be redefined. </div>
<br>
<!-- END DESCRIPTION -->
<a name="OTHER MODULES USED"></a>
<hr>
<h4>OTHER MODULES USED</h4>
<!-- BEGIN OTHER MODULES USED -->
<div>
<pre>fms_mod</pre>
</div>
<!-- END OTHER MODULES USED -->
<!-- BEGIN PUBLIC INTERFACE -->
<a name="PUBLIC INTERFACE"></a>
<hr>
<h4>PUBLIC INTERFACE</h4>
<div>
<pre>use constants_mod [, only:  constants_init ]</pre>
<dl>
<dt>
<a href="#constants_init">constants_init</a>:</dt>
<dd>   A optional initialization routine. The only purpose of this routine 
   is to write the version and tag name information to the log file. </dd>
</dl>
</div>
<br>
<!-- END PUBLIC INTERFACE -->
<a name="PUBLIC DATA"></a>
<hr>
<h4>PUBLIC DATA</h4>
<!-- BEGIN PUBLIC DATA -->
<div>
<table align="center" cellspacing="2" CELLPADDING="2" BORDER="2">
<tr>
<th> Name  </th><th> Type  </th><th> Value  </th><th> Units  </th><th> Description  </th>
</tr>
<tr>
<td> RADIUS  </td><td> real  </td><td> 6376.e3  </td><td> meters  </td><td>    radius of the earth   </td>
</tr>
<tr>
<td> OMEGA  </td><td> real  </td><td> 7.292e-5  </td><td> 1/sec  </td><td>    rotation rate of the planet (earth)   </td>
</tr>
<tr>
<td> GRAV  </td><td> real  </td><td> 9.80  </td><td> m/s2  </td><td>    acceleration due to gravity   </td>
</tr>
<tr>
<td> RDGAS  </td><td> real  </td><td> 287.04  </td><td> J/kg/deg  </td><td>    gas constant for dry air   </td>
</tr>
<tr>
<td>  KAPPA  </td><td> real  </td><td> 2./7.  </td><td>   </td><td>    RDGAS / CP   </td>
</tr>
<tr>
<td> CP  </td><td> real  </td><td> RDGAS/KAPPA  </td><td> J/kg/deg  </td><td>    specific heat capacity of dry air at constant pressure   </td>
</tr>
<tr>
<td> RVGAS  </td><td> real  </td><td> 461.50  </td><td> J/Kg/deg  </td><td>    gas constant for water vapor   </td>
</tr>
<tr>
<td> DENS_H2O  </td><td> real  </td><td> 1000.  </td><td> Kg/m3  </td><td>    density of liquid water   </td>
</tr>
<tr>
<td> HLV  </td><td> real  </td><td> 2.500e6  </td><td> J/Kg  </td><td>    latent heat of evaporation   </td>
</tr>
<tr>
<td> HLF  </td><td> real  </td><td> 3.34e5  </td><td> J/kg  </td><td>    latent heat of fusion   </td>
</tr>
<tr>
<td> HLS  </td><td> real  </td><td> 2.834e6  </td><td> J/Kg  </td><td>    latent heat of sublimation   </td>
</tr>
<tr>
<td> TFREEZE  </td><td> real  </td><td> 273.16  </td><td> deg K  </td><td>    temp where fresh water freezes   </td>
</tr>
<tr>
<td> STEFAN  </td><td> real  </td><td> 5.6734e-8  </td><td> (W/m2/deg4  </td><td>    Stefan-Boltzmann constant   </td>
</tr>
<tr>
<td> VONKARM  </td><td> real  </td><td> 0.40  </td><td> ---  </td><td>    Von Karman constant   </td>
</tr>
<tr>
<td> PI  </td><td> real  </td><td> 3.14159265358979323846  </td><td> ---  </td><td>    is it enough?   </td>
</tr>
</table>
<br>
</div>
<!-- END PUBLIC DATA -->
<a name="PUBLIC ROUTINES"></a>
<hr>
<h4>PUBLIC ROUTINES</h4>
<!-- BEGIN PUBLIC ROUTINES -->
<ol type="a">
<li>
<a name="constants_init"></a>
<h4>constants_init</h4>
<pre>
<b>call constants_init </b>
</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   The only purpose of this routine is to write the version and
   tag name information to the log file.  This routine does not
   have to be called.  If it is called more than once or called
   from other than the root PE it will return silently.
   There are no arguments. </dd>
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
<a name="REFERENCES"></a>
<hr>
<h4>REFERENCES</h4>
<!-- BEGIN REFERENCES -->
<div>
        None.
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
<div>   &lt;B&gt;NOTES ON USAGE:&lt;/B&gt;
   <br>
<br>
   All constants have been declared as type REAL, PARAMETER.
   <br>
<br>
   The value a constant can not be changed in a users program.
   New constants can be defined in terms of values from the
   constants module using a parameter statement.&lt;br&gt;&lt;br&gt;
   <br>
<br>
   The name given to a particular constant may be changed.&lt;br&gt;&lt;br&gt;
   <br>
<br>
   Constants can be used on the right side on an assignment statement
   (their value can not be reassigned). 
   <br>
<br>
   &lt;B&gt;EXAMPLES:&lt;/B&gt; <pre>     use constants_mod, only:  TFREEZE, grav_new =&gt; GRAV
     real, parameter :: grav_inv = 1.0 / grav_new
     tempc(:,:,:) = tempk(:,:,:) - TFREEZE
     geopotential(:,:) = height(:,:) * grav_new</pre> 
</div>
<br>
<!-- END NOTES -->
<a name="FUTURE PLANS"></a>
<hr>
<h4>FUTURE PLANS</h4>
<!-- BEGIN FUTURE PLANS -->
<div>
<p>   1.  Renaming of constants. </p>
<p>   2.  Additional constants. </p>
</div>
<br>
<!-- END FUTURE PLANS -->
<hr>
<div align="right">
<font size="-1"><a href="#TOP">top</a></font>
</div>
</body>
</html>
