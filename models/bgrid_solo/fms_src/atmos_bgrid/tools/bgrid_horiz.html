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
<title>module bgrid_horiz_mod</title>
</head>
<body bgcolor="#AABBCC" text="#000000">
<div align="center"><font size="-2"><a href="#INTERFACE">PUBLIC
INTERFACE</a> / <a href="#DATA_TYPES">DATA</a> / <a href=
"#ROUTINES">ROUTINES</a> / <a href="#ERRORS">ERRORS</a> / <a href=
"#REFERENCES">REFERENCES</a> / <a href=
"#NOTES">NOTES</a></font><br>
<br></div>
<hr>
<h1>module bgrid_horiz_mod</h1>
<a name="OVERVIEW" id="OVERVIEW">
<hr>
<h2>OVERVIEW</h2>
<!-- BEGIN OVERVIEW -->
<pre>

     Initializes horizontal grid constants needed by the B-grid dynamical core
     and determines the domain decomposition for running on distributed
     memory machines. This routine calls the FMS mpp_domains package.
     A derived-type variable (horiz_grid_type) is returned that contains
     the grid constants and domain data types.

</pre></a><!-- END OVERVIEW -->
<!--------------------------------------------------------------------->
 <a name="DESCRIPTION" id="DESCRIPTION"><!-- BEGIN DESCRIPTION -->
<pre>

The B-Grid
----------

The B-grid can be visualized as two overlapping grids, one that
contains momentum (u and v) and the other mass fields (surface
pressure, temperature, and tracers). These separate grids are
rectangular in shape, defined by longitudes along the x-axis and
latitude along the y-axis.  The grids are diagonally shifted from
each other, such that, the center of a momentum grid box is located
at the corners where four mass grid boxes intersect.
Horizontal indexing increases from west to east, and from south to north.
Indexing is set up so that a velocity point with the same i,j is
located to the north and east of the mass (temperature) point.


       T(i-1,j+1)            T(i,j+1)          T(i+1,j+1)         

                  v(i-1,j  )          v(i,j  )

       T(i-1,j  )            T(i,j  )          T(i+1,j  )          

                  v(i-1,j-1)          v(i,j-1)

       T(i-1,j-1)            T(i,j-1)          T(i+1,j-1)     

For computational purposes, extra rows and columns of grid boxes
are carried on the perimeter of the horizontal grid. These extra
points are called halo points. The number of halo points along
the west/east boundaries and south/north boundaries are specified
as separate arguments when initializing the horizontal grid. The
default is for one halo row along all boundaries. Because of the
grid staggering, the global momentum grid does not use it's
northernmost row.

The interior grid (excluding halo points) is called the compute domain.
The global size of the longitude axis (first dimension) can have an
even or odd number of points.  The latitude axis (2nd dimension) should
always be an even number of points.  When there is an even number of
global latitude points, mass points will straddle the equator,
while velocity points lie along the equator.

</pre></a><!-- END DESCRIPTION -->
<!--------------------------------------------------------------------->
 <a name="MODULES_USED" id="MODULES_USED">
<hr>
<h2>OTHER MODULES USED</h2>
<!-- BEGIN MODULES_USED -->
<pre>

     mpp_domains_mod
       constants_mod
             fms_mod

</pre></a><!-- END MODULES_USED -->
<!--------------------------------------------------------------------->
 <a name="INTERFACE" id="INTERFACE">
<hr>
<h2>PUBLIC INTERFACE</h2>
<!-- BEGIN INTERFACE -->
<pre>

     <b>use bgrid_horiz_mod</b> [ ,only: horiz_grid_type,
                                  bgrid_type,
                                  horiz_grid_init,
                                  get_horiz_grid_size,
                                  get_horiz_grid_bound,
                                  update_np,  update_sp,
                                  TGRID, VGRID    ]

     <a href="#DATA_TYPES">bgrid_type</a>,  <a href=
"#DATA_TYPES">horiz_grid_type</a>
          Derived-type variables containing horizontal grid constants needed
          by the B-grid dynamical core (see below).

     <a href="#horiz_grid_init">horiz_grid_init</a>
          Function that initializes/returns a variable of
          type (horiz_grid_type).

     <a href="#get_horiz_grid_size">get_horiz_grid_size</a>
          Subroutine that returns the number of grid boxes along the x and 
          y axes of the temperature grid or velocity grid.  There is an
          option for returning either the global or local compute grid sizes.
          Usually used in conjunction with routine get_horiz_grid_bound.

     <a href="#get_horiz_grid_bound">get_horiz_grid_bound</a>
          Subroutine that returns the longitude and latitude boundaries (edges)
          of temperature or velocity grid boxes. There is an option for
          returning either the global or local compute grid edges.
          Routine get_horiz_grid_size is usually called first to get
          the axis sizes.

     <a href="#update_np">update_np</a>,  <a href=
"#update_sp">update_sp</a>
          Logical functions that determine if the northern (or southern)
          halo rows for the local processor lie within the global halo
          rows (ie.e, need a polar boundary update).

     TGRID, VGRID
          Integer parameters to be used as the "grid" argument with interfaces
          get_horiz_grid_size, get_horiz_grid_bound, update_np, and update_sp.

</pre></a><!-- END INTERFACE -->
<!--------------------------------------------------------------------->
 <a name="DATA_TYPES" id="DATA_TYPES">
<hr>
<h2>PUBLIC DATA</h2>
<!-- BEGIN DATA_TYPES -->
<h2>TYPE bgrid_type</h2>
<pre>
    is, ie        <font size=
"-1">first, last x-axis index in the compute domain <font color=
"#000099">[integer,scalar]</font></font>
    js, je        <font size=
"-1">first, last y-axis index in the compute domain <font color=
"#000099">[integer,scalar]</font></font>
    isg, ieg      <font size=
"-1">first, last x-axis index in the global  domain <font color=
"#000099">[integer,scalar]</font></font>
    jsg, jeg      <font size=
"-1">first, last y-axis index in the global  domain <font color=
"#000099">[integer,scalar]</font></font>
    dx            <font size=
"-1">local grid spacing for x-axis (in meters) <font color=
"#000099">[real,dimension(ilb:iub,jlb:jub)]</font></font>
    rdx           <font size=
"-1">reciprocal of dx (1/m) <font color=
"#000099">[real,dimension(ilb:iub,jlb:jub)]</font></font>
    dy            <font size=
"-1">local grid spacing for y-axis (in meters) <font color=
"#000099">[real,scalar]</font></font>
    rdy           <font size=
"-1">reciprocal of dy (1/m) <font color=
"#000099">[real,scalar]</font></font>
    area          <font size=
"-1">area of a local grid box (in m2) <font color=
"#000099">[real,dimension(ilb:iub,jlb:jub)]</font></font>
    rarea         <font size=
"-1">reciprocal of area (1/m2) <font color=
"#000099">[real,dimension(ilb:iub,jlb:jub)]</font></font>
    tph, tlm      <font size=
"-1">latitude, longitude at the center of a local grid box (in radians)
                        <font color=
"#000099">[real,dimension(ilb:iub,jlb:jub)]</font></font>
    aph, alm      <font size=
"-1">actual latitude, longitude at the center of a local grid box (in radians)
                        <font color=
"#000099">[real,dimension(ilb:iub,jlb:jub)]</font></font>
    blatg         <font size=
"-1">latitude boundaries of grid boxes along the global y-axis (in radians)
                        <font color=
"#000099">[real,dimension(jsg:jeg+1)]</font></font>
    blong         <font size=
"-1">longitude boundaries grid boxes along the global x-axis (in radians)
                        <font color=
"#000099">[real,dimension(isg:ieg+1)]</font></font>
    Domain        <font size=
"-1">domain2D derived-type variable with halosize ihalo,jhalo
                        <font color=
"#000099">[type(domain2d)]</font></font>
    Domain_nohalo <font size=
"-1">domain2D derived-type variable without halos, used for outputting diagnostic fields
                        <font color=
"#000099">[type(domain2d)]</font></font>
</pre>
<h2>TYPE horiz_grid_type</h2>
<pre>
    Tmp          <font size=
"-1">grid constants for the temperature/tracer/mass grid <font color="#000099">[type(bgrid_type)]</font></font>
    Vel          <font size=
"-1">grid constants for the u/v wind component grid <font color=
"#000099">[type(bgrid_type)]</font></font>
    sinphv       <font size="-1">sine of Vel%aph <font size="2"
color="#000099">[real,dimension(ilb:iub,jlb:jub)]</font></font>
    tanphv       <font size="-1">tangent of Vel%aph <font size="2"
color="#000099">[real,dimension(ilb:iub,jlb:jub)]</font></font>
    nlon         <font size=
"-1">number of longitude (x-axis) grid points in the global grid (no halo points)
                      <font size="-1" color=
"#000099">[integer,scalar]</font></font>
    nlat         <font size=
"-1">number of latitude (y-axis) grid points in the global grid (no halo points)
                      <font size="-1" color=
"#000099">[integer,scalar]</font></font>
    isize, jsize <font size=
"-1">size of arrays on the local processor's grid (including halo points)
                     Note: isize=iub-ilb+1, jsize=jub-jlb+1
                      <font size="-1" color=
"#000099">[integer,scalar]</font></font>
    ilb, iub     <font size=
"-1">lower, upper bounds of local x-axis indexing <font size="2"
color="#000099">[integer,scalar]</font></font>
    jlb, jub     <font size=
"-1">lower, upper bounds of local y-axis indexing <font size="2"
color="#000099">[integer,scalar]</font></font>
    ihalo, jhalo <font size=
"-1">number of halo points along the east and west boundaries (ihalo) or
                     south and north boundaries (jhalo) <font color="#000099">[integer,scalar]</font></font>
    dlm          <font size=
"-1">grid spacing for x-axis (in radians) <font color=
"#000099">[real,scalar]</font></font>
    dph          <font size=
"-1">grid spacing for y-axis (in radians) <font color=
"#000099">[real,scalar]</font></font>
    dlmd         <font size=
"-1">grid spacing for x-axis (in degrees of longitude) <font color=
"#000099">[real,scalar]</font></font>
    dphd         <font size=
"-1">grid spacing for y-axis (in degrees of latitude) <font color=
"#000099">[real,scalar]</font></font>
    decompx      <font size=
"-1">indicates if the x-axis has been decomposed across processors <font color="#000099">[logical,scalar]</font></font>
    channel      <font size=
"-1">indicates if the channel model feature has been implemented (NOT RECOMMENDED)
                      <font color=
"#000099">[logical,scalar]</font></font>
</pre>
<pre>
NOTES:

   The local indexing variables: is, ie, js, je, ilb, iub, jlb, jub,
   are consistent with the global index values (isg, ieg, jsg, jeg).

   All local arrays (on either the temperature or velocity grid) have the same
   horizontal size [dimension(ilb:iub,jlb:jub)]. Due to the grid staggering,
   the northernmost velocity domains do not use the last j-row (jub).
</pre></a><!-- END DATA_TYPES -->
<!--------------------------------------------------------------------->
 <a name="ROUTINES" id="ROUTINES">
<hr>
<h2>PUBLIC ROUTINES</h2>
<!-- BEGIN ROUTINES --></a> <a name="horiz_grid_init" id=
"horiz_grid_init">
<pre>
Hgrid = <b>horiz_grid_init</b> ( nlon, nlat 
                          <font color=
"#007700">[, ihalo, jhalo, decomp, channel, tph0d, tlm0d]</font> )

INPUT
   nlon, nlat     The number of global longitude, latitude grid points (respectively)
                  for the mass/temperature grid.
                     <font size="-1" color=
"#000099">[integer,scalar]</font>

OPTIONAL INPUT
   ihalo, jhalo   The number of halo points for the x and y axis
                  (both mass and velocity grids).
                     <font size="-1" color=
"#000099">[integer, default: ihalo=1, jhalo=1]</font>

   decomp         The domain decomposition for the x and y axis, where
                  decomp(1) = x-axis decomposition, decomp(2) = y-axis
                  decomposition. If decomp(1)*decomp(2) does not equal
                  the number of processors the model will fail.
                  If either or both decomp(1), decomp(2) is zero,
                  <b>the rules below apply</b>.
                     <font size="-1" color=
"#000099">[integer,dimension(2), default: decomp=0,0]</font>

   channel        Flag for running in channel model mode.
                  This option has not been recently used and is currently not supported.
                     <font size="-1" color=
"#000099">[logical, default: channel=FALSE]</font>

   tph0d,tlm0d    Latitude/longitude for shifting the position of the poles.
                  Set both to zero for no transformation (the default).
                  This option has not been recently used and is currently
                  not supported.
                     <font size="-1" color=
"#000099">[real, default: tph0d=0., tlm0d=0.]</font>

RETURNED VALUE
                  Derived-type variable that contains all necessary horizontal grid information
                  and domain decomposition variables needed by the model.
                     <font size="-1" color=
"#000099">[type(horiz_grid_type)]</font>

<b>NOTES</b>

   1) If decomp(1)=0 and decomp(2)&gt;0, then decomp(1)=decomp(2)/(# PEs), or vice versa.
      The final decomposition must always satisfy: decomp(1)*decomp(2)=number of processors.
   2) If decomp(1)=decomp(2)=0, then one-dimensional decomposition of
      the y-axis will be used. If there is fewer than two rows per PE, then
      two-dimensional decomposition will be automatically implemented.
</pre>
<!------------------------>
<hr width="50%" align="left"></a> <a name="get_horiz_grid_size" id=
"get_horiz_grid_size">
<pre>
<b>call get_horiz_grid_size</b> ( Hgrid, grid, nlon, nlat <font color="#007700">[, global]</font> )

INPUT
   Hgrid       The derived-type variable returned by a previous call to horiz_grid_init.
                  <font size="-1" color=
"#000099">[type(horiz_grid_type)]</font>

   grid        Specifies which grid (temperature or velocity) the returned grid size
                  will be for. Use the publicly accessible parameters: <b>TGRID</b> or <b>VGRID</b>.
                  <font size="-1" color=
"#000099">[integer,scalar]</font>

OUTPUT
   nlon, nlat  The number of grid points along the x- and y-axis, respectively.
               The returned values will be for either the global grid size or
               the local processor's grid size depending on the value of optional
               argument "global".
                  <font size="-1" color=
"#000099">[integer,scalar]</font>

OPTIONAL INPUT
   global      Flag that determines if the returned values are for the global
               grid (global=.TRUE.) or the local compute grid (global=.FALSE.).
                  <font size="-1" color=
"#000099">[logical,scalar, default: FALSE]</font>
</pre>
<!------------------------>
<hr width="50%" align="left"></a> <a name="get_horiz_grid_bound"
id="get_horiz_grid_bound">
<pre>
<b>call get_horiz_grid_bound</b> ( Hgrid, grid, blon, blat <font color="#007700">[, global]</font> )

INPUT
   Hgrid    The derived-type variable returned by a previous call to horiz_grid_init.
               <font size="-1" color=
"#000099">[type(horiz_grid_type)]</font>

   grid     Specifies which grid (temperature or velocity) the returned grid boundaries
            will be for. Use the publicly accessible parameters: <b>TGRID</b> or <b>VGRID</b>.
               <font size="-1" color=
"#000099">[integer,scalar]</font>

OUTPUT
   blon     Longitude edges of grid boxes for either the global grid or the
            local processor's grid depending on the value of optional argument "global".
               <font size="-1" color=
"#000099">[real,dimension(nlon+1)]</font>

   blat     Latitude edges of grid boxes for either the global grid or the
            local processor's grid depending on the value of optional argument "global".
               <font size="-1" color=
"#000099">[real,dimension(nlat+1)]</font>

OPTIONAL INPUT
   global   Flag that determines if the returned values are for the global
            grid (global=.TRUE.) or the local compute grid (global=.FALSE.).
               <font size="-1" color=
"#000099">[logical,scalar, default: FALSE]</font>
</pre>
<!------------------------>
<hr width="50%" align="left"></a> <a name="update_np" id=
"update_np">
<pre>
answer = <b>update_np</b> ( Hgrid, grid )

INPUT
   Hgrid    The derived-type variable returned by a previous call to horiz_grid_init.
               <font size="-1" color=
"#000099">[type(horiz_grid_type)]</font>

   grid     Specifies which grid (temperature or velocity) the returned value
            will be for. Use the publicly accessible parameters: <b>TGRID</b> or <b>VGRID</b>.
               <font size="-1" color=
"#000099">[integer,scalar]</font>

RETURNED VALUE
   Returns TRUE when the halo rows along the <b>north boundary</b> lie within
   global halo rows, otherwise FALSE is returned. Note that for single
   processor runs the returned value will always be TRUE.
       <font size="-1" color="#000099">[logical]</font>
</pre>
<!------------------------>
<hr width="50%" align="left"></a> <a name="update_sp" id=
"update_sp">
<pre>
answer = <b>update_sp</b> ( Hgrid, grid )

INPUT
   Hgrid    The derived-type variable returned by a previous call to horiz_grid_init.
               <font size="-1" color=
"#000099">[type(horiz_grid_type)]</font>

   grid     Specifies which grid (temperature or velocity) the returned value
               will be for. Use the publicly accessible parameters: <b>TGRID</b> or <b>VGRID</b>.
               <font size="-1" color=
"#000099">[integer,scalar]</font>

RETURNED VALUE
   Returns TRUE when the halo rows along the <b>south boundary</b> lie within
   global halo rows, otherwise FALSE is returned. Note that for single
   processor runs the returned value will always be TRUE.
       <font size="-1" color="#000099">[logical]</font>
</pre></a><!-- END ROUTINES -->
<!--------------------------------------------------------------------->
 <a name="ERRORS" id="ERRORS">
<hr>
<h2>ERROR MESSAGES</h2>
<!-- BEGIN ERRORS -->
<pre>

<b>FATAL errors in horiz_grid_init</b>

    negative halo size
          <font size=
"-1">The halo input arguments must be &gt;= 1.</font>

    number of processor requested not compatible with grid
          <font size=
"-1">The domain decomposition either requested or computed was
            not compatible with the resolution of the grid.
            If you requested a decomposition other than the default,
            check to make sure it is correct.
            Also, the number of processor you are running on may not be
            compatible with the resolution of the model.
            Otherwise there may be a code error.</font>

    j halo size too big for decomposition
          <font size=
"-1">Either decrease the halo size or decrease the decomposition
            of the y-axis. The later can be done by increasing the decomposition
            of the x-axis or decreasing the number of processors.</font>

<b>FATAL errors in get_horiz_grid_bound</b>

    invalid argument dimension for blon
          <font size=
"-1">The size of the blon output argument must equal the
            number of longitude grid boxes plus one.</font>

    invalid argument dimension for blat
          <font size=
"-1">The size of the blat output argument must equal the
            number of latitude grid boxes plus one.</font>

<b>FATAL errors in get_horiz_grid_size, get_horiz_grid_bound,
                 update_sp, update_np</b>

    invalid grid
          <font size=
"-1">The input argument "grid" has an incorrect value.
            Make sure you are using one of public parameters, <b>TGRID</b> or <b>VGRID</b>.</font>

</pre></a><!-- END ERRORS -->
<!--------------------------------------------------------------------->
 <a name="REFERENCES" id="REFERENCES">
<hr>
<h2>REFERENCES</h2>
<!-- BEGIN REFERENCES -->
<pre>

     None.

</pre></a><!-- END REFERENCES -->
<!--------------------------------------------------------------------->
 <a name="BUGS" id="BUGS">
<hr>
<h2>KNOWN BUGS</h2>
<!-- BEGIN BUGS -->
<pre>

   * The grid transformation option (tph0d, tlm0d) has not been
     extensively tested. If want to use this option check with
     the developer.

   * The channel model mode has not been extensively tested.

</pre></a><!-- END BUGS -->
<!--------------------------------------------------------------------->
 <a name="NOTES" id="NOTES">
<hr>
<h2>NOTES</h2>
<!-- BEGIN NOTES -->
<pre>

     None.

</pre></a><!-- END NOTES -->
<!--------------------------------------------------------------------->
 <a name="PLANS" id="PLANS">
<hr>
<h2>FUTURE PLANS</h2>
<!-- BEGIN PLANS -->
<pre>

     None.

</pre></a><!-- END PLANS -->
<!--------------------------------------------------------------------->
<hr>
</body>
</html>
