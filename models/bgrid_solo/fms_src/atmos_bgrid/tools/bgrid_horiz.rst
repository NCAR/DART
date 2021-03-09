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
<HTML>
<TITLE>module bgrid_horiz_mod</TITLE>
<BODY BGCOLOR="#AABBCC" TEXT="#000000" >

<DIV ALIGN="CENTER"> <FONT SIZE=-2>
<A HREF="#INTERFACE">PUBLIC INTERFACE</A> / 
<A HREF="#DATA_TYPES">DATA</A> / 
<A HREF="#ROUTINES">ROUTINES</A> / 
<A HREF="#ERRORS">ERRORS</A> / 
<A HREF="#REFERENCES">REFERENCES</A> / 
<A HREF="#NOTES">NOTES</A> 
</FONT>
<BR><BR></DIV><HR>


<H2>module bgrid_horiz_mod</H2>
<A NAME="HEADER">
<PRE>
     <B>Contact:</B>   B. Wyman
     <B>Reviewers:</B>
     <B>Change history:</B> <A HREF="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS
/atmos/bgrid/tools/bgrid_horiz.f90">WebCVS Log for bgrid_horiz.f90</A>

</PRE>
</A><!-- END HEADER -->
<!--------------------------------------------------------------------->
<A NAME="OVERVIEW">
<HR>
<H4>OVERVIEW</H4>
<!-- BEGIN OVERVIEW -->
<PRE>

     Initializes horizontal grid constants needed by the B-grid dynamical core
     and determines the domain decomposition for running on distributed
     memory machines. This routine calls the FMS mpp_domains package.
     A derived-type variable (horiz_grid_type) is returned that contains
     the grid constants and domain data types.

</PRE>
</A><!-- END OVERVIEW -->
<!--------------------------------------------------------------------->
<A NAME="DESCRIPTION">
<!-- BEGIN DESCRIPTION -->
<PRE>

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

</PRE>
</A><!-- END DESCRIPTION -->
<!--------------------------------------------------------------------->
<A NAME="MODULES_USED">
<HR>
<H4>OTHER MODULES USED</H4>
<!-- BEGIN MODULES_USED -->
<PRE>

     mpp_domains_mod
       constants_mod
             fms_mod

</PRE>
</A><!-- END MODULES_USED -->
<!--------------------------------------------------------------------->
<A NAME="INTERFACE">
<HR>
<H4>PUBLIC INTERFACE</H4>
<!-- BEGIN INTERFACE -->
<PRE>

     <b>use bgrid_horiz_mod</b> [ ,only: horiz_grid_type,
                                  bgrid_type,
                                  horiz_grid_init,
                                  get_horiz_grid_size,
                                  get_horiz_grid_bound,
                                  update_np,  update_sp,
                                  TGRID, VGRID    ]

     <a href="#DATA_TYPES">bgrid_type</a>,  <a href="#DATA_TYPES">horiz_grid_type</a>
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

     <a href="#update_np">update_np</a>,  <a href="#update_sp">update_sp</a>
          Logical functions that determine if the northern (or southern)
          halo rows for the local processor lie within the global halo
          rows (ie.e, need a polar boundary update).

     TGRID, VGRID
          Integer parameters to be used as the "grid" argument with interfaces
          get_horiz_grid_size, get_horiz_grid_bound, update_np, and update_sp.

</PRE>
</A><!-- END INTERFACE -->
<!--------------------------------------------------------------------->
<A NAME="DATA_TYPES">
<HR>
<H4>PUBLIC DATA</H4>
<!-- BEGIN DATA_TYPES -->
<H4>TYPE bgrid_type</H4>
<PRE>
    is, ie        <FONT SIZE=-1>first, last x-axis index in the compute domain <FONT COLOR="#000099">[integer,scalar]</FONT></FONT>
    js, je        <FONT SIZE=-1>first, last y-axis index in the compute domain <FONT COLOR="#000099">[integer,scalar]</FONT></FONT>
    isg, ieg      <FONT SIZE=-1>first, last x-axis index in the global  domain <FONT COLOR="#000099">[integer,scalar]</FONT></FONT>
    jsg, jeg      <FONT SIZE=-1>first, last y-axis index in the global  domain <FONT COLOR="#000099">[integer,scalar]</FONT></FONT>
    dx            <FONT SIZE=-1>local grid spacing for x-axis (in meters) <FONT COLOR="#000099">[real,dimension(ilb:iub,jlb:jub)]</FONT></FONT>
    rdx           <FONT SIZE=-1>reciprocal of dx (1/m) <FONT COLOR="#000099">[real,dimension(ilb:iub,jlb:jub)]</FONT></FONT>
    dy            <FONT SIZE=-1>local grid spacing for y-axis (in meters) <FONT COLOR="#000099">[real,scalar]</FONT></FONT>
    rdy           <FONT SIZE=-1>reciprocal of dy (1/m) <FONT COLOR="#000099">[real,scalar]</FONT></FONT>
    area          <FONT SIZE=-1>area of a local grid box (in m2) <FONT COLOR="#000099">[real,dimension(ilb:iub,jlb:jub)]</FONT></FONT>
    rarea         <FONT SIZE=-1>reciprocal of area (1/m2) <FONT COLOR="#000099">[real,dimension(ilb:iub,jlb:jub)]</FONT></FONT>
    tph, tlm      <FONT SIZE=-1>latitude, longitude at the center of a local grid box (in radians)
                        <FONT COLOR="#000099">[real,dimension(ilb:iub,jlb:jub)]</FONT></FONT>
    aph, alm      <FONT SIZE=-1>actual latitude, longitude at the center of a local grid box (in radians)
                        <FONT COLOR="#000099">[real,dimension(ilb:iub,jlb:jub)]</FONT></FONT>
    blatg         <FONT SIZE=-1>latitude boundaries of grid boxes along the global y-axis (in radians)
                        <FONT COLOR="#000099">[real,dimension(jsg:jeg+1)]</FONT></FONT>
    blong         <FONT SIZE=-1>longitude boundaries grid boxes along the global x-axis (in radians)
                        <FONT COLOR="#000099">[real,dimension(isg:ieg+1)]</FONT></FONT>
    Domain        <FONT SIZE=-1>domain2D derived-type variable with halosize ihalo,jhalo
                        <FONT COLOR="#000099">[type(domain2d)]</FONT></FONT>
    Domain_nohalo <FONT SIZE=-1>domain2D derived-type variable without halos, used for outputting diagnostic fields
                        <FONT COLOR="#000099">[type(domain2d)]</FONT></FONT>
</PRE>

<H4>TYPE horiz_grid_type</H4>
<PRE>
    Tmp          <FONT SIZE=-1>grid constants for the temperature/tracer/mass grid <FONT COLOR="#000099">[type(bgrid_type)]</FONT></FONT>
    Vel          <FONT SIZE=-1>grid constants for the u/v wind component grid <FONT COLOR="#000099">[type(bgrid_type)]</FONT></FONT>
    sinphv       <FONT SIZE=-1>sine of Vel%aph <FONT SIZE=2 COLOR="#000099">[real,dimension(ilb:iub,jlb:jub)]</FONT></FONT>
    tanphv       <FONT SIZE=-1>tangent of Vel%aph <FONT SIZE=2 COLOR="#000099">[real,dimension(ilb:iub,jlb:jub)]</FONT></FONT>
    nlon         <FONT SIZE=-1>number of longitude (x-axis) grid points in the global grid (no halo points)
                      <FONT SIZE=-1 COLOR="#000099">[integer,scalar]</FONT></FONT>
    nlat         <FONT SIZE=-1>number of latitude (y-axis) grid points in the global grid (no halo points)
                      <FONT SIZE=-1 COLOR="#000099">[integer,scalar]</FONT></FONT>
    isize, jsize <FONT SIZE=-1>size of arrays on the local processor's grid (including halo points)
                     Note: isize=iub-ilb+1, jsize=jub-jlb+1
                      <FONT SIZE=-1 COLOR="#000099">[integer,scalar]</FONT></FONT>
    ilb, iub     <FONT SIZE=-1>lower, upper bounds of local x-axis indexing <FONT SIZE=2 COLOR="#000099">[integer,scalar]</FONT></FONT>
    jlb, jub     <FONT SIZE=-1>lower, upper bounds of local y-axis indexing <FONT SIZE=2 COLOR="#000099">[integer,scalar]</FONT></FONT>
    ihalo, jhalo <FONT SIZE=-1>number of halo points along the east and west boundaries (ihalo) or
                     south and north boundaries (jhalo) <FONT COLOR="#000099">[integer,scalar]</FONT></FONT>
    dlm          <FONT SIZE=-1>grid spacing for x-axis (in radians) <FONT COLOR="#000099">[real,scalar]</FONT></FONT>
    dph          <FONT SIZE=-1>grid spacing for y-axis (in radians) <FONT COLOR="#000099">[real,scalar]</FONT></FONT>
    dlmd         <FONT SIZE=-1>grid spacing for x-axis (in degrees of longitude) <FONT COLOR="#000099">[real,scalar]</FONT></FONT>
    dphd         <FONT SIZE=-1>grid spacing for y-axis (in degrees of latitude) <FONT COLOR="#000099">[real,scalar]</FONT></FONT>
    decompx      <FONT SIZE=-1>indicates if the x-axis has been decomposed across processors <FONT COLOR="#000099">[logical,scalar]</FONT></FONT>
    channel      <FONT SIZE=-1>indicates if the channel model feature has been implemented (NOT RECOMMENDED)
                      <FONT COLOR="#000099">[logical,scalar]</FONT></FONT>
</PRE>

<PRE>
NOTES:

   The local indexing variables: is, ie, js, je, ilb, iub, jlb, jub,
   are consistent with the global index values (isg, ieg, jsg, jeg).

   All local arrays (on either the temperature or velocity grid) have the same
   horizontal size [dimension(ilb:iub,jlb:jub)]. Due to the grid staggering,
   the northernmost velocity domains do not use the last j-row (jub).
</PRE>

</A><!-- END DATA_TYPES -->
<!--------------------------------------------------------------------->
<A NAME="ROUTINES">
<HR>
<H4>PUBLIC ROUTINES</H4>
<!-- BEGIN ROUTINES -->
<a name="horiz_grid_init">

<PRE>
Hgrid = <B>horiz_grid_init</B> ( nlon, nlat 
                          <FONT COLOR="#007700">[, ihalo, jhalo, decomp, channel, tph0d, tlm0d]</FONT> )

INPUT
   nlon, nlat     The number of global longitude, latitude grid points (respectively)
                  for the mass/temperature grid.
                     <FONT SIZE=-1 COLOR="#000099">[integer,scalar]</FONT>

OPTIONAL INPUT
   ihalo, jhalo   The number of halo points for the x and y axis
                  (both mass and velocity grids).
                     <FONT SIZE=-1 COLOR="#000099">[integer, default: ihalo=1, jhalo=1]</FONT>

   decomp         The domain decomposition for the x and y axis, where
                  decomp(1) = x-axis decomposition, decomp(2) = y-axis
                  decomposition. If decomp(1)*decomp(2) does not equal
                  the number of processors the model will fail.
                  If either or both decomp(1), decomp(2) is zero,
                  <b>the rules below apply</b>.
                     <FONT SIZE=-1 COLOR="#000099">[integer,dimension(2), default: decomp=0,0]</FONT>

   channel        Flag for running in channel model mode.
                  This option has not been recently used and is currently not supported.
                     <FONT SIZE=-1 COLOR="#000099">[logical, default: channel=FALSE]</FONT>

   tph0d,tlm0d    Latitude/longitude for shifting the position of the poles.
                  Set both to zero for no transformation (the default).
                  This option has not been recently used and is currently
                  not supported.
                     <FONT SIZE=-1 COLOR="#000099">[real, default: tph0d=0., tlm0d=0.]</FONT>

RETURNED VALUE
                  Derived-type variable that contains all necessary horizontal grid information
                  and domain decomposition variables needed by the model.
                     <FONT SIZE=-1 COLOR="#000099">[type(horiz_grid_type)]</FONT>

<B>NOTES</B>

   1) If decomp(1)=0 and decomp(2)&gt;0, then decomp(1)=decomp(2)/(# PEs), or vice versa.
      The final decomposition must always satisfy: decomp(1)*decomp(2)=number of processors.
   2) If decomp(1)=decomp(2)=0, then one-dimensional decomposition of
      the y-axis will be used. If there is fewer than two rows per PE, then
      two-dimensional decomposition will be automatically implemented.
</PRE>

<!------------------------>
<HR WIDTH="50%" ALIGN=LEFT>
<a name="get_horiz_grid_size">

<PRE>
<B>call get_horiz_grid_size</B> ( Hgrid, grid, nlon, nlat <FONT COLOR="#007700">[, global]</FONT> )

INPUT
   Hgrid       The derived-type variable returned by a previous call to horiz_grid_init.
                  <FONT SIZE=-1 COLOR="#000099">[type(horiz_grid_type)]</FONT>

   grid        Specifies which grid (temperature or velocity) the returned grid size
                  will be for. Use the publicly accessible parameters: <B>TGRID</B> or <B>VGRID</B>.
                  <FONT SIZE=-1 COLOR="#000099">[integer,scalar]</FONT>

OUTPUT
   nlon, nlat  The number of grid points along the x- and y-axis, respectively.
               The returned values will be for either the global grid size or
               the local processor's grid size depending on the value of optional
               argument "global".
                  <FONT SIZE=-1 COLOR="#000099">[integer,scalar]</FONT>

OPTIONAL INPUT
   global      Flag that determines if the returned values are for the global
               grid (global=.TRUE.) or the local compute grid (global=.FALSE.).
                  <FONT SIZE=-1 COLOR="#000099">[logical,scalar, default: FALSE]</FONT>
</PRE>

<!------------------------>
<HR WIDTH="50%" ALIGN=LEFT>
<a name="get_horiz_grid_bound">

<PRE>
<B>call get_horiz_grid_bound</B> ( Hgrid, grid, blon, blat <FONT COLOR="#007700">[, global]</FONT> )

INPUT
   Hgrid    The derived-type variable returned by a previous call to horiz_grid_init.
               <FONT SIZE=-1 COLOR="#000099">[type(horiz_grid_type)]</FONT>

   grid     Specifies which grid (temperature or velocity) the returned grid boundaries
            will be for. Use the publicly accessible parameters: <B>TGRID</B> or <B>VGRID</B>.
               <FONT SIZE=-1 COLOR="#000099">[integer,scalar]</FONT>

OUTPUT
   blon     Longitude edges of grid boxes for either the global grid or the
            local processor's grid depending on the value of optional argument "global".
               <FONT SIZE=-1 COLOR="#000099">[real,dimension(nlon+1)]</FONT>

   blat     Latitude edges of grid boxes for either the global grid or the
            local processor's grid depending on the value of optional argument "global".
               <FONT SIZE=-1 COLOR="#000099">[real,dimension(nlat+1)]</FONT>

OPTIONAL INPUT
   global   Flag that determines if the returned values are for the global
            grid (global=.TRUE.) or the local compute grid (global=.FALSE.).
               <FONT SIZE=-1 COLOR="#000099">[logical,scalar, default: FALSE]</FONT>
</PRE>

<!------------------------>
<HR WIDTH="50%" ALIGN=LEFT>
<a name="update_np">

<PRE>
answer = <B>update_np</B> ( Hgrid, grid )

INPUT
   Hgrid    The derived-type variable returned by a previous call to horiz_grid_init.
               <FONT SIZE=-1 COLOR="#000099">[type(horiz_grid_type)]</FONT>

   grid     Specifies which grid (temperature or velocity) the returned value
            will be for. Use the publicly accessible parameters: <B>TGRID</B> or <B>VGRID</B>.
               <FONT SIZE=-1 COLOR="#000099">[integer,scalar]</FONT>

RETURNED VALUE
   Returns TRUE when the halo rows along the <B>north boundary</B> lie within
   global halo rows, otherwise FALSE is returned. Note that for single
   processor runs the returned value will always be TRUE.
       <FONT SIZE=-1 COLOR="#000099">[logical]</FONT>
</PRE>

<!------------------------>
<HR WIDTH="50%" ALIGN=LEFT>
<a name="update_sp">

<PRE>
answer = <B>update_sp</B> ( Hgrid, grid )

INPUT
   Hgrid    The derived-type variable returned by a previous call to horiz_grid_init.
               <FONT SIZE=-1 COLOR="#000099">[type(horiz_grid_type)]</FONT>

   grid     Specifies which grid (temperature or velocity) the returned value
               will be for. Use the publicly accessible parameters: <B>TGRID</B> or <B>VGRID</B>.
               <FONT SIZE=-1 COLOR="#000099">[integer,scalar]</FONT>

RETURNED VALUE
   Returns TRUE when the halo rows along the <B>south boundary</B> lie within
   global halo rows, otherwise FALSE is returned. Note that for single
   processor runs the returned value will always be TRUE.
       <FONT SIZE=-1 COLOR="#000099">[logical]</FONT></FONT>
</PRE>

</A><!-- END ROUTINES -->
<!--------------------------------------------------------------------->
<A NAME="ERRORS">
<HR>
<H4>ERROR MESSAGES</H4>
<!-- BEGIN ERRORS -->
<PRE>

<b>FATAL errors in horiz_grid_init</b>

    negative halo size
          <FONT SIZE=-1>The halo input arguments must be >= 1.</FONT>

    number of processor requested not compatible with grid
          <FONT SIZE=-1>The domain decomposition either requested or computed was
            not compatible with the resolution of the grid.
            If you requested a decomposition other than the default,
            check to make sure it is correct.
            Also, the number of processor you are running on may not be
            compatible with the resolution of the model.
            Otherwise there may be a code error.</FONT>

    j halo size too big for decomposition
          <FONT SIZE=-1>Either decrease the halo size or decrease the decomposition
            of the y-axis. The later can be done by increasing the decomposition
            of the x-axis or decreasing the number of processors.</FONT>

<b>FATAL errors in get_horiz_grid_bound</b>

    invalid argument dimension for blon
          <FONT SIZE=-1>The size of the blon output argument must equal the
            number of longitude grid boxes plus one.</FONT>

    invalid argument dimension for blat
          <FONT SIZE=-1>The size of the blat output argument must equal the
            number of latitude grid boxes plus one.</FONT>

<b>FATAL errors in get_horiz_grid_size, get_horiz_grid_bound,
                 update_sp, update_np</b>

    invalid grid
          <FONT SIZE=-1>The input argument "grid" has an incorrect value.
            Make sure you are using one of public parameters, <b>TGRID</b> or <b>VGRID</b>.</FONT>

</PRE>
</A><!-- END ERRORS -->
<!--------------------------------------------------------------------->
<A NAME="REFERENCES">
<HR>
<H4>REFERENCES</H4>
<!-- BEGIN REFERENCES -->
<PRE>

     None.

</PRE>
</A><!-- END REFERENCES -->
<!--------------------------------------------------------------------->
<A NAME="BUGS">
<HR>
<H4>KNOWN BUGS</H4>
<!-- BEGIN BUGS -->
<PRE>

   * The grid transformation option (tph0d, tlm0d) has not been
     extensively tested. If want to use this option check with
     the developer.

   * The channel model mode has not been extensively tested.

</PRE>
</A><!-- END BUGS -->
<!--------------------------------------------------------------------->
<A NAME="NOTES">
<HR>
<H4>NOTES</H4>
<!-- BEGIN NOTES -->
<PRE>

     None.

</PRE>
</A><!-- END NOTES -->
<!--------------------------------------------------------------------->
<A NAME="PLANS">
<HR>
<H4>FUTURE PLANS</H4>
<!-- BEGIN PLANS -->
<PRE>

     None.

</PRE>
</A><!-- END PLANS -->
<!--------------------------------------------------------------------->

<HR>
</BODY>
</HTML>
