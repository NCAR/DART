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
<TITLE>module bgrid_halo_mod</TITLE>
<BODY BGCOLOR="#AABBCC" TEXT="#332211" >

<DIV ALIGN="CENTER"> <FONT SIZE=-2>
<A HREF="#INTERFACE">PUBLIC INTERFACE</A> / 
<A HREF="#ROUTINES">ROUTINES</A> / 
<A HREF="#ERRORS">ERRORS</A> / 
<A HREF="#NOTES">NOTES</A> 
</FONT>
<BR><BR></DIV><HR>


<H2>module bgrid_halo_mod</H2>
<A NAME="HEADER">
<PRE>
     <B>Contact:</B>   B. Wyman
     <B>Reviewers:</B>
     <B>Change history:</B> <A HREF="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS
/atmos/bgrid/tools/bgrid_halo.f90">WebCVS Log for bgrid_halo.f90</A>

</PRE>
</A><!-- END HEADER -->
<!--------------------------------------------------------------------->
<A NAME="OVERVIEW">
<HR>
<H4>OVERVIEW</H4>
<!-- BEGIN OVERVIEW -->
<PRE>

     Provides an interface for updating the B-grid dynamical core
     halo rows and columns, including the polar halo rows.
     See the <A HREF="#NOTES">notes</A> section for a description of the polar 
     boundary condition.

</PRE>
</A><!-- END OVERVIEW -->
<!--------------------------------------------------------------------->
<A NAME="DESCRIPTION">
<!-- BEGIN DESCRIPTION -->
<PRE>
     Additional interfaces are provided for updating the halo rows
     of horiz_grid_type derived-type variables and for setting
     fluxes at the sub-polar to zero.

</PRE>
</A><!-- END DESCRIPTION -->
<!--------------------------------------------------------------------->
<A NAME="MODULES_USED">
<HR>
<H4>OTHER MODULES USED</H4>
<!-- BEGIN MODULES_USED -->
<PRE>

    horiz_grid_mod
     utilities_mod
   mpp_domains_mod

</PRE>
</A><!-- END MODULES_USED -->
<!--------------------------------------------------------------------->
<A NAME="INTERFACE">
<HR>
<H4>PUBLIC INTERFACE</H4>
<!-- BEGIN INTERFACE -->
<PRE>

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


</PRE>
</A><!-- END INTERFACE -->
<!--------------------------------------------------------------------->
<A NAME="ROUTINES">
<HR>
<H4>PUBLIC ROUTINES</H4>
<!-- BEGIN ROUTINES -->
<PRE>
<a name="update_halo">

<b>call update_halo</b> ( Hgrid, field, data <FONT COLOR="#007700">[,flags]</FONT> )

INPUT
   field     Specifies which field/grid the data is on.
             You must use the publicly accessible parameters: <B>TEMP</B>, <B>UWND</B>, or <B>VWND</B>.
             See the description of these parameters above.
                <FONT SIZE=-1 COLOR="#000099">[integer,scalar]</FONT>

INPUT/OUTPUT
   Hgrid     The derived-type variable returned by a previous call to horiz_grid_init.
             See the module horiz_grid_mod for details.
                <FONT SIZE=-1 COLOR="#000099">[type(horiz_grid_type)]</FONT>

   data      Data array on any valid grid, may have 2, 3, or 4 dimensions.
             The dimensions correspond to the x, y and z axes, and tracer number.
                <FONT SIZE=-1 COLOR="#000099">[real,dimension(:,:) or dimension(:,:,:) or dimension(:,:,:,:)]</FONT>

OPTIONAL INPUT
   flags     Integer flag that describes which halo regions should be updated.
             By default all halo regions are updated. The value of flags 
             should be some combination of the public parameter values
             SOUTH, NORTH, WEST, EAST, and NOPOLE. For example, to only
             update the north and east halo regions set flags=NORTH+EAST.
             The flag for NOPOLE suppresses the halo update of velocity at the poles.

             An additional flag POLEONLY may be used independently to 
             update only the north and south polar halo rows without updating
             the interior halo rows. This kind of update will require no
             processor to processor communication.
                <FONT SIZE=-1 COLOR="#000099">[integer,scalar]</FONT>

------------------------------------------------------------------------
<a name="horiz_grid_boundary">

<b>call horiz_grid_boundary</b> ( Hgrid )

INPUT/OUTPUT
   Hgrid     The derived-type variable returned by a previous call to horiz_grid_init.
             See the module horiz_grid_mod for details.
                <FONT SIZE=-1 COLOR="#000099">[type(horiz_grid_type)]</FONT>

------------------------------------------------------------------------
<a name="vel_flux_boundary">

<b>call vel_flux_boundary</b> ( Hgrid, data )

INPUT
   Hgrid     The derived-type variable returned by a previous call to horiz_grid_init.
             See the module horiz_grid_mod for details.
                <FONT SIZE=-1 COLOR="#000099">[type(horiz_grid_type)]</FONT>
INPUT/OUTPUT
   data      Real data array on velocity grid, may have 2 or 3 dimensions.
             The dimensions correspond to the x, y and z axes.
             The data in the sub-polar row will be set to zero.
                <FONT SIZE=-1 COLOR="#000099">[real,dimension(:,:) or dimension(:,:,:)]</FONT>

</PRE>
</A><!-- END ROUTINES -->
<!--------------------------------------------------------------------->
<A NAME="ERRORS">
<HR>
<H4>ERROR MESSAGES</H4>
<!-- BEGIN ERRORS -->
<PRE>

<b>Fatal errors in update_halo</b>

    i dimension has wrong size
        <FONT SIZE=-1>The 1st (i) dimension of input/output argument data must
          have a size equal to Hgrid % isize (the entire i dimension).</FONT>

    j dimension has wrong size
        <FONT SIZE=-1>The 2nd (j) dimension of input/output argument data must
          have a size equal to Hgrid % jsize (the entire j dimension).</FONT>

    invalid value for flags
        <FONT SIZE=-1>The value of optional argument flags was invalid.  This can only 
          happen when the user specifies a value.  The value of flags should
          be set using the public module parameters, and must be in the
          following range: 0 &lt; flags &lt;= SOUTH+NORTH+WEST+EAST+NOPOLE+POLEONLY.</FONT>

    invalid field
        <FONT SIZE=-1>The input argument "field" has an incorrect value.  Make sure
          you are using one of public parameters: <b>TEMP</b>, <b>UWND</b>, or <b>VWND</b>.</FONT>

</PRE>
</A><!-- END ERRORS -->
<!--------------------------------------------------------------------->
<A NAME="BUGS">
<HR>
<H4>KNOWN BUGS</H4>
<!-- BEGIN BUGS -->
<PRE>

     None.

</PRE>
</A><!-- END BUGS -->
<!--------------------------------------------------------------------->
<A NAME="NOTES">
<HR>
<H4>NOTES</H4>
<!-- BEGIN NOTES -->
<PRE>
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
