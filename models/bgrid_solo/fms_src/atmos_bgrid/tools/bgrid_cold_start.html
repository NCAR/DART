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
<title>module bgrid_cold_start_mod</title>
</head>
<body bgcolor="#AABBCC" text="#332211">
<div align="center"><font size="-2"><a href="#INTERFACE">PUBLIC
INTERFACE</a> / <a href="#ROUTINES">ROUTINES</a> / <a href=
"#NAMELIST">NAMELIST</a> / <a href="#ERRORS">ERRORS</a> / <a href=
"#NOTES">NOTES</a></font><br>
<br></div>
<hr>
<h1>module bgrid_cold_start_mod</h1>
<a name="OVERVIEW" id="OVERVIEW">
<hr>
<h2>OVERVIEW</h2>
<!-- BEGIN OVERVIEW -->
<pre>

     Generates initial conditions for the B-grid dynamics core.

</pre></a><!-- END OVERVIEW -->
<!--------------------------------------------------------------------->
 <a name="DESCRIPTION" id="DESCRIPTION"><!-- BEGIN DESCRIPTION -->
<pre>
     Provides initial data (flow at rest) for all fields that
     are saved on the B-grid core restart file.

     The intended purpose of this module is to create initial data
     when running an idealized version of the model. However,
     there are options to generate (or read) realistic topography,
     but there is no way to generate a land-sea mask.

</pre></a><!-- END DESCRIPTION -->
<!--------------------------------------------------------------------->
 <a name="MODULES_USED" id="MODULES_USED">
<hr>
<h2>OTHER MODULES USED</h2>
<!-- BEGIN MODULES_USED -->
<pre>

    bgrid_horiz_mod
     bgrid_halo_mod
     topography_mod
            fms_mod
      constants_mod

</pre></a><!-- END MODULES_USED -->
<!--------------------------------------------------------------------->
 <a name="INTERFACE" id="INTERFACE">
<hr>
<h2>PUBLIC INTERFACE</h2>
<!-- BEGIN INTERFACE -->
<pre>

   <b>use bgrid_cold_start_mod</b> [ ,only: cold_start_resol, cold_start ]

   <a href="#cold_start_resol">cold_start_resol</a>
        Returns the model resolution (read from namelist bgrid_cold_start_nml)
        to the calling program so that storage for the prognostic variables
        can be allocated.

   <a href="#cold_start">cold_start</a>
        Generates initial data (flow at rest) for all fields saved on
        the B-grid dynamical core restart file.


   <b>NOTES</b>

      * A namelist interface called <a href=
"#NAMELIST">bgrid_cold_start_nml</a> is read
        from file <b>input.nml</b>.

</pre></a><!-- END INTERFACE -->
<!--------------------------------------------------------------------->
 <a name="ROUTINES" id="ROUTINES">
<hr>
<h2>PUBLIC ROUTINES</h2>
<!-- BEGIN ROUTINES -->
<pre>
<a name="cold_start_resol" id="cold_start_resol">

<b>call cold_start_resol</b> ( nx, ny, nz ) OUTPUT nx number of grid boxes along the longitude (x) axis ny number of grid boxes along the latitude (y) axis nz number of vertical model levels NOTES The values returned by this routine are from the namelist &amp;bgrid_cold_start_nml. -------------------------------------------------------------------- </a><a name="cold_start"
id="cold_start">

<b>call cold_start</b> ( Hgrid, eta, peta, fis, res, ps, pssl, u, v, t ) INPUT/OUTPUT Hgrid Horizontal grid constants (see horiz_grid_mod) OUTPUT eta Values of sigma/eta at the interface between model levels. Must be dimensioned by the number of levels + 1. NOTE: This value corresponds to the "B" term when determining pressure (i.e., p = A + B*ps). [real, dimension(:)] peta Values of constant pressure at the interface between model levels. Must be dimensioned by the number of levels + 1. NOTE: This value corresponds to the "A" term when determining pressure (i.e., p = A + B*ps). [real, dimension(:)] fis geopotential height of the surface (m2/s2) [real, dimension(:,:)] res reciprocal of eta at the surface [real, dimension(:,:)] ps surface pressure (pascals) [real, dimension(:,:)] pssl surface pressure at eta=1. (pascals) [real, dimension(:,:)] u zonal wind component (m/s) [real, dimension(:,:,:)] v meridional wind component (m/s) [real, dimension(:,:,:)] t temperature (deg K) [real, dimension(:,:,:)] NOTES The storage for all arguments must be allocated by the calling program. </a></pre></a><!-- END ROUTINES -->
<!--------------------------------------------------------------------->
 <a name="NAMELIST" id="NAMELIST">
<hr>
<h2>NAMELIST</h2>
<!-- BEGIN NAMELIST -->
<pre>

 <b>&amp;bgrid_cold_start_nml</b>

   nlon   = number of grid points along the longitude axis (1st dimension)
              [integer, default: nlon = 0]

   nlat   = number of grid points along the latitude axis (2nd dimension)
              [integer, default: nlat = 0]

   nlev   = number of vertical levels
              [integer, default: nlev = 0]

   pref   = initial surface pressure in pascals
              [real, default: pref = 1000.e2]

   tref   = initial temperature in deg kelvin
              [real, default: tref = 255.]

   equal_vert_spacing  = Should the levels be equally spaced in sigma or
                         spaced using the formula of Smagorinski (1969).
                           [logical, default: equal_vert_spacing = .true.]

 NOTES

   1) If nlon, nlat, or nlev are not specified the program will terminate.
   2) There is no option for the hybrid or eta coordinate. 
      All vertical coordinates use the pure sigma system.


</pre></a><!-- END NAMELIST -->
<!--------------------------------------------------------------------->
 <a name="ERRORS" id="ERRORS">
<hr>
<h2>ERROR MESSAGES</h2>
<!-- BEGIN ERRORS -->
<pre>

<b>FATAL errors in bgrid_cold_start_mod</b>

    <b>resolution not specified</b>
        One of two possible problems may have occurred.
        1) Must specify values for namelist variables: nlon, nlat, nlev.
        2) If you expected to read a restart file, the restart was probably
           not found and the model has inadvertently tried to generate one.

    <b>incorrect resolution in file topography.res</b>
        The resolution read from the topography file does not agree with
        the resolution given by the namelist.

</pre></a><!-- END ERRORS -->
<!--------------------------------------------------------------------->
 <a name="NOTES" id="NOTES">
<hr>
<h2>NOTES</h2>
<!-- BEGIN NOTES -->
<pre>

  Initial data (momentum,temperature,surface pressure) are computed as:

         u  = 0.0
         v  = 0.0
         t  = tref
         ps = pref * exp( -fis / (tref * Rd) )

  where

      fis = geopotential height at surface (m2/s2)
      Rd  = dry gas constant
      tref, pref = namelist parameters

</pre></a><!-- END NOTES -->
<!--------------------------------------------------------------------->
 <a name="PLANS" id="PLANS">
<hr>
<h2>FUTURE PLANS</h2>
<!-- BEGIN PLANS -->
<pre>

   Coordinate development of this module with the "off-line" program
   for generating spin-up initial conditions.

</pre></a><!-- END PLANS -->
<!--------------------------------------------------------------------->
<hr>
</body>
</html>
