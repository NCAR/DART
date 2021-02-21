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
<title>module bgrid_halo_mod</title>
</head>
<body bgcolor="#AABBCC" text="#332211">
<div align="center"><font size="-2"><a href="#INTERFACE">PUBLIC
INTERFACE</a> / <a href="#ROUTINES">ROUTINES</a> / <a href=
"#ERRORS">ERRORS</a> / <a href="#NOTES">NOTES</a></font><br>
<br></div>
<hr>
<h1>module bgrid_halo_mod</h1>
<a name="OVERVIEW" id="OVERVIEW">
<hr>
<h2>OVERVIEW</h2>
<!-- BEGIN OVERVIEW -->
<pre>

     Provides an interface for updating the B-grid dynamical core
     halo rows and columns, including the polar halo rows.
     See the <a href=
"#NOTES">notes</a> section for a description of the polar 
     boundary condition.

</pre></a><!-- END OVERVIEW -->
<!--------------------------------------------------------------------->
 <a name="DESCRIPTION" id="DESCRIPTION"><!-- BEGIN DESCRIPTION -->
<pre>
     Additional interfaces are provided for updating the halo rows
     of horiz_grid_type derived-type variables and for setting
     fluxes at the sub-polar to zero.

</pre></a><!-- END DESCRIPTION -->
<!--------------------------------------------------------------------->
 <a name="MODULES_USED" id="MODULES_USED">
<hr>
<h2>OTHER MODULES USED</h2>
<!-- BEGIN MODULES_USED -->
<pre>

    horiz_grid_mod
     utilities_mod
   mpp_domains_mod

</pre></a><!-- END MODULES_USED -->
<!--------------------------------------------------------------------->
 <a name="INTERFACE" id="INTERFACE">
<hr>
<h2>PUBLIC INTERFACE</h2>
<!-- BEGIN INTERFACE -->
<pre>

   <b>use bgrid_halo_mod</b> [, only: update_halo,
                               horiz_grid_boundary,
                               vel_flux_boundary,
                               TEMP, UWND, VWND,
                               SOUTH, NORTH, WEST, EAST, NOPOLE, POLEONLY ]

   <a href="#update_halo">update_halo</a>
        Updates all requested rows in the halo region for a requested
        field of type: <b>TEMP</b>, <b>UWND</b>, or <b>VWND</b>.  By default all boundaries
        will be updated.  To update only specific boundaries use the
        public parameters: <b>SOUTH</b>, <b>NORTH</b>, <b>WEST</b>, <b>EAST</b>, <b>NOPOLE</b>, and <b>POLEONLY</b>.

   <a href="#horiz_grid_boundary">horiz_grid_boundary</a>
        Updates the halo region for variables of derived-type horiz_grid_type.
        Should be called after horiz_grid_init.

   <a href="#vel_flux_boundary">vel_flux_boundary</a>
        Zeros the sub-polar row of fields on the velocity grid.
        This routine is needed by several dynamics routines for conservation,
        but should not normally be needed by the user.

   TEMP, UWND, VWND
        Integer parameters to be used as the "field" argument in
        interface update_halo.  The VWND value will result in a sign flip
        beyond the polar rows.  Data at auxilary mass flux points can
        also use these values: U-flux points use TEMP and V-flux points use
        either UWND or VWND depending whether a sign flip is desired.

   SOUTH, NORTH, WEST, EAST, NOPOLE, POLEONLY
        Integer parameters to be used as the optional "flags" argument to
        interface update_halo. NOPOLE and POLEONLY apply only to fields
        UWND and VWND.


</pre></a><!-- END INTERFACE -->
<!--------------------------------------------------------------------->
 <a name="ROUTINES" id="ROUTINES">
<hr>
<h2>PUBLIC ROUTINES</h2>
<!-- BEGIN ROUTINES -->
<pre>
<a name="update_halo" id="update_halo">

<b>call update_halo</b> ( Hgrid, field, data <font color=
"#007700">[,flags]</font> ) INPUT field Specifies which field/grid the data is on. You must use the publicly accessible parameters: <b>TEMP</b>, <b>UWND</b>, or <b>VWND</b>. See the description of these parameters above. <font size="-1"
color=
"#000099">[integer,scalar]</font> INPUT/OUTPUT Hgrid The derived-type variable returned by a previous call to horiz_grid_init. See the module horiz_grid_mod for details. <font size="-1"
color=
"#000099">[type(horiz_grid_type)]</font> data Data array on any valid grid, may have 2, 3, or 4 dimensions. The dimensions correspond to the x, y and z axes, and tracer number. <font size="-1"
color=
"#000099">[real,dimension(:,:) or dimension(:,:,:) or dimension(:,:,:,:)]</font> OPTIONAL INPUT flags Integer flag that describes which halo regions should be updated. By default all halo regions are updated. The value of flags should be some combination of the public parameter values SOUTH, NORTH, WEST, EAST, and NOPOLE. For example, to only update the north and east halo regions set flags=NORTH+EAST. The flag for NOPOLE suppresses the halo update of velocity at the poles. An additional flag POLEONLY may be used independently to update only the north and south polar halo rows without updating the interior halo rows. This kind of update will require no processor to processor communication. <font size="-1"
color=
"#000099">[integer,scalar]</font> ------------------------------------------------------------------------ </a><a name="horiz_grid_boundary"
id="horiz_grid_boundary">

<b>call horiz_grid_boundary</b> ( Hgrid ) INPUT/OUTPUT Hgrid The derived-type variable returned by a previous call to horiz_grid_init. See the module horiz_grid_mod for details. <font size="-1"
color=
"#000099">[type(horiz_grid_type)]</font> ------------------------------------------------------------------------ </a><a name="vel_flux_boundary"
id="vel_flux_boundary">

<b>call vel_flux_boundary</b> ( Hgrid, data ) INPUT Hgrid The derived-type variable returned by a previous call to horiz_grid_init. See the module horiz_grid_mod for details. <font size="-1"
color=
"#000099">[type(horiz_grid_type)]</font> INPUT/OUTPUT data Real data array on velocity grid, may have 2 or 3 dimensions. The dimensions correspond to the x, y and z axes. The data in the sub-polar row will be set to zero. <font size="-1"
color=
"#000099">[real,dimension(:,:) or dimension(:,:,:)]</font> </a></pre></a><!-- END ROUTINES -->
<!--------------------------------------------------------------------->
 <a name="ERRORS" id="ERRORS">
<hr>
<h2>ERROR MESSAGES</h2>
<!-- BEGIN ERRORS -->
<pre>

<b>Fatal errors in update_halo</b>

    i dimension has wrong size
        <font size=
"-1">The 1st (i) dimension of input/output argument data must
          have a size equal to Hgrid % isize (the entire i dimension).</font>

    j dimension has wrong size
        <font size=
"-1">The 2nd (j) dimension of input/output argument data must
          have a size equal to Hgrid % jsize (the entire j dimension).</font>

    invalid value for flags
        <font size=
"-1">The value of optional argument flags was invalid.  This can only 
          happen when the user specifies a value.  The value of flags should
          be set using the public module parameters, and must be in the
          following range: 0 &lt; flags &lt;= SOUTH+NORTH+WEST+EAST+NOPOLE+POLEONLY.</font>

    invalid field
        <font size=
"-1">The input argument "field" has an incorrect value.  Make sure
          you are using one of public parameters: <b>TEMP</b>, <b>UWND</b>, or <b>VWND</b>.</font>

</pre></a><!-- END ERRORS -->
<!--------------------------------------------------------------------->
 <a name="BUGS" id="BUGS">
<hr>
<h2>KNOWN BUGS</h2>
<!-- BEGIN BUGS -->
<pre>

     None.

</pre></a><!-- END BUGS -->
<!--------------------------------------------------------------------->
 <a name="NOTES" id="NOTES">
<hr>
<h2>NOTES</h2>
<!-- BEGIN NOTES -->
<pre>
   At the north-south polar boundaries the following 
   boundary conditions are applied:

      a) Velocities at the poles are equal to zero

             u(i,p) = 0
             v(i,p) = 0

      b) Halo points along the north and south polar boundaries
         are set as follows:

             T(i,p+1/2) =  T(i,p-1/2)
             u(i,p+1)   =  u(i,p-1)
             v(i,p+1)   = -v(i,p-1)
      
         where p + # is a halo row and p - # is a row within 
         the computational domain.

   If there is no decomposition along the x-axis then the east-west
   boundaries are updated for global cyclic continuity.

   All other halo points are updated using MPP_UPDATE_DOMAINS.
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
