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
<title>Module topography_mod</title>
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
<h2>module topography_mod</h2>
<a name="HEADER"></a>
<!-- BEGIN HEADER -->
<div>
<b>Contact:</b>&nbsp;<a href="mailto:bw@gfdl.noaa.gov">   Bruce Wyman </a>
<br>
<b>Reviewers:</b>&nbsp;<br>
<b>Change History: </b>&nbsp;<a href="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/">WebCVS Log</a>
<br>
<b>Last Modified:</b>&nbsp;2002/03/22 01:45:05<br>
<br>
</div>
<!-- END HEADER -->
<a name="OVERVIEW"></a>
<hr>
<h4>OVERVIEW</h4>
<!-- BEGIN OVERVIEW -->
<p class="text">   Routines for creating land surface topography fields and land-water masks
   for latitude-longitude grids. </p>
<!-- END OVERVIEW -->
<a name="DESCRIPTION"></a>
<!-- BEGIN DESCRIPTION -->
<div>   This module generates realistic mountains and land-water masks
   on a specified latitude-longitude grid by interpolating from the
   1/6 degree Navy mean topography and percent water data sets.
   The fields that can be generated are mean and standard deviation
   of topography within the specified grid boxes; and land-ocean (or
   water) mask and land-ocean (or water) fractional area.
   <br>
<br>
   The interpolation scheme conserves the area-weighted average
   of the input data by using module horiz_interp.
   <br>
<br>
   The interfaces get_gaussian_topog and gaussian_topog_init are documented in <a href="gaussian_topog.html">gaussian_topog_mod</a>. </div>
<br>
<!-- END DESCRIPTION -->
<a name="OTHER MODULES USED"></a>
<hr>
<h4>OTHER MODULES USED</h4>
<!-- BEGIN OTHER MODULES USED -->
<div>
<pre>gaussian_topog_mod<br>  horiz_interp_mod<br>           fms_mod</pre>
</div>
<!-- END OTHER MODULES USED -->
<!-- BEGIN PUBLIC INTERFACE -->
<a name="PUBLIC INTERFACE"></a>
<hr>
<h4>PUBLIC INTERFACE</h4>
<div>
<pre>use topography_mod [, only:  get_topog_mean,<br>                             get_topog_stdev,<br>                             get_ocean_frac,<br>                             get_ocean_mask,<br>                             get_water_frac,<br>                             get_water_mask ]</pre>
<dl>
<dt>
<a href="#get_topog_mean">get_topog_mean</a>:</dt>
<dd>   Returns a "realistic" mean surface height field. </dd>
<dt>
<a href="#get_topog_stdev">get_topog_stdev</a>:</dt>
<dd>   Returns a standard deviation of higher resolution topography with 
   the given model grid boxes. </dd>
<dt>
<a href="#get_ocean_frac">get_ocean_frac</a>:</dt>
<dd>   Returns fractional area covered by ocean in a grid box. </dd>
<dt>
<a href="#get_ocean_mask">get_ocean_mask</a>:</dt>
<dd>   Returns a land-ocean mask in a grid box. </dd>
<dt>
<a href="#get_water_frac">get_water_frac</a>:</dt>
<dd>   Returns fractional area covered by water. </dd>
<dt>
<a href="#get_water_mask">get_water_mask</a>:</dt>
<dd>   Returns a land-water mask in a grid box. </dd>
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
<a name="get_topog_mean"></a>
<h4>get_topog_mean</h4>
<pre>flag = &lt;B&gt; <b>get_topog_mean</b> &lt;/B&gt; ( blon, blat, zmean )</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   Returns realistic mountains on a latitude-longtude grid.
   The returned field is the mean topography for the given grid boxes.
   Computed using a conserving area-weighted interpolation.
   The current input data set is the 1/6 degree Navy mean topography. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>blon&nbsp;&nbsp;&nbsp;</tt></td><td>   The longitude (in radians) at grid box boundaries. <br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>blat&nbsp;&nbsp;&nbsp;</tt></td><td>   The latitude (in radians) at grid box boundaries. <br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:)]</span></td>
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
<td valign="top" align="left"><tt>zmean&nbsp;&nbsp;&nbsp;</tt></td><td>   The mean surface height (meters).
   The size of this field must be size(blon)-1 by size(blat)-1. <br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>get_topog_mean&nbsp;&nbsp;&nbsp;</tt></td><td>   A logical value of TRUE is returned if the surface height field
   was successfully created. A value of FALSE may be returned if the
   input topography data set was not readable. <br>&nbsp;&nbsp;&nbsp;<span class="type">[logical]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="get_topog_stdev"></a>
<h4>get_topog_stdev</h4>
<pre>flag = &lt;B&gt; <b>get_topog_stdev</b> &lt;/B&gt; ( blon, blat, stdev )</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   Returns the standard deviation of the "finer" input topography data set,
   currently the Navy 1/6 degree mean topography data, within the
   boundaries of the given input grid. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>blon&nbsp;&nbsp;&nbsp;</tt></td><td>   The longitude (in radians) at grid box boundaries. <br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>blat&nbsp;&nbsp;&nbsp;</tt></td><td>   The latitude (in radians) at grid box boundaries. <br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:)]</span></td>
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
<td valign="top" align="left"><tt>stdev&nbsp;&nbsp;&nbsp;</tt></td><td>   The standard deviation of surface height (in meters) within
   given input model grid boxes.
   The size of this field must be size(blon)-1 by size(blat)-1. <br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>get_topog_stdev&nbsp;&nbsp;&nbsp;</tt></td><td>   A logical value of TRUE is returned if the output field was
   successfully created. A value of FALSE may be returned if the
   input topography data set was not readable. <br>&nbsp;&nbsp;&nbsp;<span class="type">[logical]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="get_ocean_frac"></a>
<h4>get_ocean_frac</h4>
<pre>flag = &lt;B&gt; <b>get_ocean_frac</b> &lt;/B&gt; ( blon, blat, ocean_frac )</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   Returns fractional area covered by ocean in the given model grid boxes. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>blon&nbsp;&nbsp;&nbsp;</tt></td><td>   The longitude (in radians) at grid box boundaries. <br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>blat&nbsp;&nbsp;&nbsp;</tt></td><td>   The latitude (in radians) at grid box boundaries. <br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:)]</span></td>
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
<td valign="top" align="left"><tt>ocean_frac&nbsp;&nbsp;&nbsp;</tt></td><td>   The fractional amount (0 to 1) of ocean in a grid box.
   The size of this field must be size(blon)-1 by size(blat)-1. <br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>get_ocean_frac&nbsp;&nbsp;&nbsp;</tt></td><td>   A logical value of TRUE is returned if the output field
   was successfully created. A value of FALSE may be returned
   if the Navy 1/6 degree percent water data set was not readable. <br>&nbsp;&nbsp;&nbsp;<span class="type">[logical]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="get_ocean_mask"></a>
<h4>get_ocean_mask</h4>
<pre>flag = &lt;B&gt; <b>get_ocean_mask</b> &lt;/B&gt; ( blon, blat, ocean_mask )</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   Returns a land-ocean mask in the given model grid boxes. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>blon&nbsp;&nbsp;&nbsp;</tt></td><td>   The longitude (in radians) at grid box boundaries. <br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>blat&nbsp;&nbsp;&nbsp;</tt></td><td>   The latitude (in radians) at grid box boundaries. <br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:)]</span></td>
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
<td valign="top" align="left"><tt>ocean_frac&nbsp;&nbsp;&nbsp;</tt></td><td>   The fractional amount (0 to 1) of ocean in a grid box.
   The size of this field must be size(blon)-1 by size(blat)-1. <br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>get_ocean_mask&nbsp;&nbsp;&nbsp;</tt></td><td>   A logical value of TRUE is returned if the output field
   was successfully created. A value of FALSE may be returned
   if the Navy 1/6 degree percent water data set was not readable. <br>&nbsp;&nbsp;&nbsp;<span class="type">[logical]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="get_water_frac"></a>
<h4>get_water_frac</h4>
<pre>flag = &lt;B&gt; <b>get_water_frac</b> &lt;/B&gt; ( blon, blat, water_frac )</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   Returns the percent of water in a grid box. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>blon&nbsp;&nbsp;&nbsp;</tt></td><td>   The longitude (in radians) at grid box boundaries. <br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>blat&nbsp;&nbsp;&nbsp;</tt></td><td>   The latitude (in radians) at grid box boundaries. <br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:)]</span></td>
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
<td valign="top" align="left"><tt>water_frac&nbsp;&nbsp;&nbsp;</tt></td><td>   The fractional amount (0 to 1) of water in a grid box.
   The size of this field must be size(blon)-1 by size(blat)-1. <br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>get_water_frac&nbsp;&nbsp;&nbsp;</tt></td><td>   A logical value of TRUE is returned if the output field
   was successfully created. A value of FALSE may be returned
   if the Navy 1/6 degree percent water data set was not readable. <br>&nbsp;&nbsp;&nbsp;<span class="type">[logical]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="get_water_mask"></a>
<h4>get_water_mask</h4>
<pre>flag = &lt;B&gt; <b>get_water_mask</b> &lt;/B&gt; ( blon, blat, water_mask )</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   Returns a land-water mask in the given model grid boxes. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>blon&nbsp;&nbsp;&nbsp;</tt></td><td>   The longitude (in radians) at grid box boundaries. <br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>blat&nbsp;&nbsp;&nbsp;</tt></td><td>   The latitude (in radians) at grid box boundaries. <br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:)]</span></td>
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
<td valign="top" align="left"><tt>water_mask&nbsp;&nbsp;&nbsp;</tt></td><td>   A binary mask for water (true) or land (false).
   The size of this field must be size(blon)-1 by size(blat)-1. <br>&nbsp;&nbsp;&nbsp;<span class="type">[real, dimension(:,:)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>get_water_mask&nbsp;&nbsp;&nbsp;</tt></td><td>   A logical value of TRUE is returned if the output field
   was successfully created. A value of FALSE may be returned
   if the Navy 1/6 degree percent water data set was not readable. <br>&nbsp;&nbsp;&nbsp;<span class="type">[logical]</span></td>
</tr>
</table>
</dd>
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
<b>&amp;topography_nml</b>
<br>
<br>
<dl>
<dt>
<tt>topog_file</tt>
</dt>
<dl>   Name of topography file. <br>
<span class="type">[character, default: DATA/navy_topography.data]</span>
</dl>
<dt>
<tt>water_file</tt>
</dt>
<dl>   Name of percent water file. <br>
<span class="type">[character, default: DATA/navy_pctwater.data]</span>
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
<div>
<dl>
<dt></dt>
<dd>   This module uses the 1/6 degree U.S. Navy mean topography
   and percent water data sets.
   <br>
<br>
   These data sets have been re-formatted to separate 32-bit IEEE files.
   The names of these files is specified by the <a href="#NAMELIST">namelist</a>   input.
   <br>
<br>
   The format for both files is as follows: <pre>     record = 1    nlon, nlat
     record = 2    blon, blat
     record = 3    data</pre>   where: <pre>     nlon, nlat = The number of longitude and latitude points
                  in the horizontal grid.  For the 1/6 degree
                  data sets this is 2160 x 1080. [integer]
     blon, blat = The longitude and latitude grid box boundaries in degrees.
                     [real :: blon(nlon+1), blat(nlat+1)]

     data       = The topography or percent water data.
                    [real :: data(nlon,nlat)]</pre> 
</dd>
</dl>
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
<b>FATAL in get_topog_mean</b>
</dt>
<dd>
<span class="errmsg">shape(zmean) is not equal to (/size(blon)-1,size(blat)-1/))</span>
</dd>
<dd>   Check the input grid size and output field size. </dd>
<dt>
<b>FATAL in get_water_frac</b>
</dt>
<dd>
<span class="errmsg">shape(water_frac) is not equal to (/size(blon)-1,size(blat)-1/))</span>
</dd>
<dd>   Check the input grid size and output field size. </dd>
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
<div>
<dl>
<dt></dt>
<dd>   To run this program you will need the topography and percent water
   data sets and use the following namelist (in file input.nml).
   <br>
<br>
   &amp;gaussian_topog_nml
   height = 5000., 3000., 3000., 3000.,
   olon   =   90.,  255.,  285.,    0.,
   olat   =   45.,   45.,  -15.,  -90.,
   wlon   =   15.,   10.,    5.,  180.,
   wlat   =   15.,   25.,   25.,   20., /
   <br>
<br>
   program test
   <br>
<br>
   test program for topography and gaussian_topog modules <pre>  use topography_mod
  implicit none
  
  integer, parameter :: nlon=24, nlat=18
  real :: x(nlon), y(nlat), xb(nlon+1), yb(nlat+1), z(nlon,nlat)
  real :: hpi, rtd
  integer :: i,j
  logical :: a
  
 gaussian mountain parameters
  real, parameter :: ht=4000.
  real, parameter :: x0=90., y0=45. ! origin in degrees
  real, parameter :: xw=15., yw=15. ! half-width in degees
  real, parameter :: xr=30., yr= 0. ! ridge-width in degrees
  
 create lat/lon grid in radians
    hpi = acos(0.0)
    rtd = 90./hpi ! rad to deg
    do i=1,nlon
      xb(i) = 4.*hpi*real(i-1)/real(nlon)
    enddo
      xb(nlon+1) = xb(1)+4.*hpi
      yb(1) = -hpi
    do j=2,nlat
      yb(j) = yb(j-1) + 2.*hpi/real(nlat)
    enddo
      yb(nlat+1) = hpi
 mid-point of grid boxes
    x(1:nlon) = 0.5*(xb(1:nlon)+xb(2:nlon+1))
    y(1:nlat) = 0.5*(yb(1:nlat)+yb(2:nlat+1))
 test topography_mod routines
    a = get_topog_mean(xb,yb,z)
    call printz ('get_topog_mean')
  
    a = get_water_frac(xb,yb,z)
    z = z*100. ! in percent
    call printz ('get_water_frac')
  
    a = get_ocean_frac(xb,yb,z)
    z = z*100. ! in percent
    call printz ('get_ocean_frac')
  
 test gaussian_topog_mod routines
    a = .true.
    z = get_gaussian_topog(x,y,ht,x0,y0,xw,yw,xr,yr)
    call printz ('get_gaussian_topog')
  
    call gaussian_topog_init (x,y,z)
    call printz ('gaussian_topog_init')
  
  contains
  
 simple printout of topog/water array
    subroutine printz (lab)
    character(len=*), intent(in) :: lab
     if (a) then
        print '(/a)', trim(lab)
     else
        print '(/a)', 'no data available: '//trim(lab)
        return
     endif
 print full grid
        print '(3x,25i5)', (nint(x(i)*rtd),i=1,nlon)
      do j=nlat,1,-1
        print '(i3,25i5)', nint(y(j)*rtd), (nint(z(i,j)),i=1,nlon)
      enddo
    end subroutine printz
  
  end program test</pre> 
</dd>
</dl>
</div>
<br>
<!-- END TEST PROGRAM -->
<a name="KNOWN BUGS"></a>
<hr>
<h4>KNOWN BUGS</h4>
<!-- BEGIN KNOWN BUGS -->
<div>
<p>   Water mask produces some possible erroneous water points along
   the coast of Antarctic (at about 90W). </p>
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
<p>Use of netcdf data sets. </p>
<p>Incorporate other topography and ocean data sets. </p>
</div>
<br>
<!-- END FUTURE PLANS -->
<hr>
<div align="right">
<font size="-1"><a href="#TOP">top</a></font>
</div>
</body>
</html>
