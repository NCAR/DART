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
<title>Module gaussian_topog_mod</title>
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
<h2>module gaussian_topog_mod</h2>
<a name="HEADER"></a>
<!-- BEGIN HEADER -->
<div>
<b>Contact:</b>&nbsp;<a href="mailto:bw@gfdl.noaa.gov">   Bruce Wyman </a>
<br>
<b>Reviewers:</b>&nbsp;<br>
<b>Change History: </b>&nbsp;<a href="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/">WebCVS Log</a>
<br>
<b>Last Modified:</b>&nbsp;2002/03/22 01:42:43<br>
<br>
</div>
<!-- END HEADER -->
<a name="OVERVIEW"></a>
<hr>
<h4>OVERVIEW</h4>
<!-- BEGIN OVERVIEW -->
<p class="text">   Routines for creating Gaussian-shaped land surface topography
   for latitude-longitude grids. </p>
<!-- END OVERVIEW -->
<a name="DESCRIPTION"></a>
<!-- BEGIN DESCRIPTION -->
<div>   Interfaces generate simple Gaussian-shaped mountains from
   parameters specified by either argument list or namelist input.
   The mountain shapes are controlled by the height, half-width,
   and ridge-width parameters. </div>
<br>
<!-- END DESCRIPTION -->
<a name="OTHER MODULES USED"></a>
<hr>
<h4>OTHER MODULES USED</h4>
<!-- BEGIN OTHER MODULES USED -->
<div>
<pre>      fms_mod<br>constants_mod</pre>
</div>
<!-- END OTHER MODULES USED -->
<!-- BEGIN PUBLIC INTERFACE -->
<a name="PUBLIC INTERFACE"></a>
<hr>
<h4>PUBLIC INTERFACE</h4>
<div>
<pre>use gaussian_topog_mod [, only:  gaussian_topog_init,<br>                                 get_gaussian_topog ]</pre>
<dl>
<dt>
<a href="#gaussian_topog_init">gaussian_topog_init</a>:</dt>
<dd>   Returns a surface height field that consists
   of the sum of one or more Gaussian-shaped mountains. </dd>
<dt>
<a href="#get_gaussian_topog">get_gaussian_topog</a>:</dt>
<dd>   Returns a simple surface height field that consists of a single
   Gaussian-shaped mountain. </dd>
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
<a name="gaussian_topog_init"></a>
<h4>gaussian_topog_init</h4>
<pre>&lt;B&gt;call <b>gaussian_topog_init</b> &lt;/B&gt; ( lon, lat, zsurf )</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   Returns a land surface topography that consists of a "set" of
   simple Gaussian-shaped mountains.  The height, position,
   width, and elongation of the mountains can be controlled
   by variables in namelist <a href="#NAMELIST">&amp;gaussian_topog_nml</a>. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>lon&nbsp;&nbsp;&nbsp;</tt></td><td>   The mean grid box longitude in radians. <br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>lat&nbsp;&nbsp;&nbsp;</tt></td><td>   The mean grid box latitude in radians. <br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:)]</span></td>
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
<td valign="top" align="left"><tt>zsurf&nbsp;&nbsp;&nbsp;</tt></td><td>   The surface height (in meters).
   The size of this field must be size(lon) by size(lat). <br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:,:)]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="get_gaussian_topog"></a>
<h4>get_gaussian_topog</h4>
<pre>zsurf = &lt;B&gt; <b>get_gaussian_topog</b> &lt;/B&gt; ( lon, lat, height [, olond, olatd, wlond, wlatd, rlond, rlatd ] )</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   Returns a single Gaussian-shaped mountain.
   The height, position, width, and elongation of the mountain
   is controlled by optional arguments. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>lon&nbsp;&nbsp;&nbsp;</tt></td><td>   The mean grid box longitude in radians. <br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>lat&nbsp;&nbsp;&nbsp;</tt></td><td>   The mean grid box latitude in radians. <br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>height&nbsp;&nbsp;&nbsp;</tt></td><td>   Maximum surface height in meters. <br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(scalar)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>olond, olatd&nbsp;&nbsp;&nbsp;</tt></td><td>   Position/origin of mountain in degrees longitude and latitude.
   This is the location of the maximum height. <br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(scalar)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>wlond, wlatd&nbsp;&nbsp;&nbsp;</tt></td><td>   Gaussian half-width of mountain in degrees longitude and latitude. <br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(scalar)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>rlond, rlatd&nbsp;&nbsp;&nbsp;</tt></td><td>   Ridge half-width of mountain in degrees longitude and latitude.
   This is the elongation of the maximum height. <br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(scalar)]</span></td>
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
<td valign="top" align="left"><tt>zsurf&nbsp;&nbsp;&nbsp;</tt></td><td>   The surface height (in meters).
   The size of the returned field is size(lon) by size(lat). <br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:,:)]</span></td>
</tr>
</table>
</dd>
<br>
<dt>
<b>NOTE</b>
</dt>
<dd>   Mountains do not wrap around the poles. </dd>
<br>
<br>
</dl>
</li>
</ol>
<!-- END PUBLIC ROUTINES -->
<a name="NAMELIST"></a>
<!-- BEGIN NAMELIST -->
<hr>
<h4>NAMELIST</h4>
<div>
<b>&amp;gaussian_topog_nml</b>
<br>
<br>
<dl>
<dt>
<tt>height</tt>
</dt>
<dl>   Height in meters of the Gaussian mountains. <br>
<span class="type">[real, dimension(mxmtns), units: meter, default: 0.]</span>
</dl>
<dt>
<tt>olon, olat</tt>
</dt>
<dl>   The longitude and latitude of mountain origins (in degrees). <br>
<span class="type">[real, dimension(mxmtns), units: degree, default: 0.]</span>
</dl>
<dt>
<tt>wlon, wlat</tt>
</dt>
<dl>   The longitude and latitude half-width of mountain tails (in degrees). <br>
<span class="type">[real, dimension(mxmtns), units: degree, default: 0.]</span>
</dl>
<dt>
<tt>rlon, rlat</tt>
</dt>
<dl>   The longitude and latitude half-width of mountain ridges (in degrees).  For a
   "standard" Gaussian mountain set rlon=rlat=0. <br>
<span class="type">[real, dimension(mxmtns), units: degree, default: 0.]</span>
</dl>
<dt>
<tt>NOTE</tt>
</dt>
<dl>   The variables in this namelist are only used when routine
   &lt;TT&gt;gaussian_topog_init&lt;/TT&gt; is called.  The namelist variables
   are dimensioned (by 10), so that multiple mountains can be generated.
   <br>
<br>
   Internal parameter mxmtns = 10. By default no mountains are generated. <br>
<span class="type">[]</span>
</dl>
</dl>
</div>
<br>
<!-- END NAMELIST -->
<a name="DIAGNOSTIC FIELDS"></a>
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
<b>FATAL in get_gaussian_topog</b>
</dt>
<dd>
<span class="errmsg">shape(zsurf) is not equal to (/size(lon),size(lat)/)</span>
</dd>
<dd>   Check the input grid size and output field size.
   The input grid is defined at the midpoint of grid boxes. </dd>
</dl>
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
<div>   NAMELIST FOR GENERATING GAUSSIAN MOUNTAINS
   <br>
<br>
   * multiple mountains can be generated
   * the final mountains are the sum of all
   <br>
<br>
   height = height in meters
   olon, olat = longitude,latitude origin              (degrees)
   rlon, rlat = longitude,latitude half-width of ridge (degrees)
   wlon, wlat = longitude,latitude half-width of tail  (degrees)
   <br>
<br>
   Note: For the standard gaussian mountain
   set rlon = rlat = 0 .
   <br>
<br> 
<pre>       height --&gt;   ___________________________
                   /                           \
                  /              |              \
    gaussian     /               |               \
      sides --&gt; /                |                \
               /               olon                \
         _____/                olat                 \______

              |    |             |
              |&lt;--&gt;|&lt;-----------&gt;|
              |wlon|    rlon     |
               wlat     rlat</pre>   See the <a href="topography.html#TEST PROGRAM">topography </a>module documentation for a test program. </div>
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
