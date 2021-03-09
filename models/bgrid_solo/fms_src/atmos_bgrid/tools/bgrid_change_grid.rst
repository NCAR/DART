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
<TITLE>module bgrid_change_grid_mod</TITLE>
<BODY BGCOLOR="#AABBCC" TEXT="#332211" >

<DIV ALIGN="CENTER"> <FONT SIZE=-2>
<A HREF="#INTERFACE">PUBLIC INTERFACE</A> / 
<A HREF="#ROUTINES">ROUTINES</A> / 
<A HREF="#ERRORS">ERRORS</A> / 
<A HREF="#NOTES">NOTES</A> 
</FONT>
<BR><BR></DIV><HR>


<H2>module bgrid_change_grid_mod</H2>
<A NAME="HEADER">
<PRE>
     <B>Contact:</B>   Bruce Wyman
     <B>Reviewers:</B>
     <B>Change history:</B> <A HREF="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS
/atmos/bgrid/tools/bgrid_change_grid.f90">WebCVS Log for bgrid_change_grid.f90</A>

</PRE>
</A><!-- END HEADER -->
<!--------------------------------------------------------------------->
<A NAME="OVERVIEW">
<HR>
<H4>OVERVIEW</H4>
<!-- BEGIN OVERVIEW -->
<PRE>

    Provides interfaces to interpolate between the B-grid mass and
    velocity grids.

</PRE>
</A><!-- END OVERVIEW -->
<!--------------------------------------------------------------------->
<A NAME="DESCRIPTION">
<!-- BEGIN DESCRIPTION -->
<PRE>
    The interpolation to a grid box is performed using the (four) 
    grid boxes centered at it's corners.  The interfaces have been
    overloaded for area-weighted interpolation and simple
    equal-weighted (4-point) interpolation.

</PRE>
</A><!-- END DESCRIPTION -->
<!--------------------------------------------------------------------->
<A NAME="MODULES_USED">
<HR>
<H4>OTHER MODULES USED</H4>
<!-- BEGIN MODULES_USED -->
<PRE>

    bgrid_horiz_mod

</PRE>
</A><!-- END MODULES_USED -->
<!--------------------------------------------------------------------->
<A NAME="INTERFACE">
<HR>
<H4>PUBLIC INTERFACE</H4>
<!-- BEGIN INTERFACE -->
<PRE>

   <b>use bgrid_change_grid_mod</b> [,only:  mass_to_vel, vel_to_mass ]

   <a href="#mass_to_vel">mass_to_vel</a>
        Interpolates one field from the mass grid to the velocity grid.

   <a href="#vel_to_mass">vel_to_mass</a>
        Interpolates the velocity components from the velocity grid to
        the mass grid.

</PRE>
</A><!-- END INTERFACE -->
<!--------------------------------------------------------------------->
<A NAME="ROUTINES">
<HR>
<H4>PUBLIC ROUTINES</H4>
<!-- BEGIN ROUTINES -->
<PRE>
<a name="mass_to_vel">

<b>call mass_to_vel</b> (Hgrid, fm, fv)
         OR
<b>call mass_to_vel</b>        (fm, fv)

  INPUT

      Hgrid    horiz_grid_type (see horiz_grid_mod)
               When this variable is present, area weighted averaging
               is performed. When this variable is not present,
               simple four-point averaging is performed.

      fm       2-d or 3-d real array located at mass points

  OUTPUT

      fv       2-d or 3-d real array located at velocity points


  NOTES

      1) No output value is calculated in the east-most and north-most rows.
      2) If the Hgrid interface is used, then the horizontal dimensions of
         the input/output arrays must be consistent with the size of the
         local data domain.
      3) The input and output arrays can be the same since a temporary array
         is used for the result.

--------------------------------------------------------------------
<a name="vel_to_mass">

<b>call vel_to_mass</b> (Hgrid, u, v, um, vm, mask)
        OR
<b>call vel_to_mass</b>        (u, v, um, vm, mask)

  INPUT

      Hgrid    horiz_grid_type (see horiz_grid_mod)
               When this variable is present, area weighted averaging
               is performed. When this variable is not present,
               simple four-point averaging is performed.
   
      u, v     2-d or 3-d real arrays of velocity components array

  OUTPUT

      um, vm   2-d or 3-d real arrays averaged to mass points

  INPUT

      mask     2-d or 3-d topography mask array (real, 0. or 1.) located
               at velocity points for the eta/step-mountain coordinate


  NOTES

      1) No output value is calculated in the west-most and south-most rows.
      2) If the Hgrid interface is used, then the horizontal dimensions of
         the input/output arrays must be consistent with the size of the
         local data domain.
      3) The input and output arrays can be the same since a temporary array
         is used for the result.

</PRE>
</A><!-- END ROUTINES -->
<!--------------------------------------------------------------------->
<A NAME="ERRORS">
<HR>
<H4>ERROR MESSAGES</H4>
<!-- BEGIN ERRORS -->
<PRE>

     None.

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

     The 2-d versions call the 3-d versions.

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
