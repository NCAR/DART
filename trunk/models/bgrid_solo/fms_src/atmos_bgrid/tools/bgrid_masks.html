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
<TITLE>module bgrid_masks_mod</TITLE>
<BODY BGCOLOR="#AABBCC" TEXT="#332211" >

<DIV ALIGN="CENTER"> <FONT SIZE=1>
<A HREF="#INTERFACE">PUBLIC INTERFACE</A> / 
<A HREF="#ROUTINES">ROUTINES</A> / 
<A HREF="#ERRORS">ERRORS</A> / 
<A HREF="#NOTES">NOTES</A> 
</FONT>
<BR><BR></DIV><HR>


<H2>module bgrid_masks_mod</H2>
<A NAME="HEADER">
<PRE>
     <B>Contact:</B>   Bruce Wyman
     <B>Reviewers:</B>
     <B>Change history:</B> <A HREF="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS
/atmos/bgrid/tools/bgrid_masks.f90">WebCVS Log for bgrid_masks.f90</A>

</PRE>
</A><!-- END HEADER -->
<!--------------------------------------------------------------------->
<A NAME="OVERVIEW">
<HR>
<H4>OVERVIEW</H4>
<!-- BEGIN OVERVIEW -->
<PRE>

     Provides a data structure for three-dimensional masks used
     to define the step-mountain/eta coordinate topography.

</PRE>
</A><!-- END OVERVIEW -->
<!--------------------------------------------------------------------->
<A NAME="DESCRIPTION">
<!-- BEGIN DESCRIPTION -->
<PRE>

</PRE>
</A><!-- END DESCRIPTION -->
<!--------------------------------------------------------------------->
<A NAME="MODULES_USED">
<HR>
<H4>OTHER MODULES USED</H4>
<!-- BEGIN MODULES_USED -->
<PRE>

   bgrid_horiz_mod
   bgrid_vert_mod
   bgrid_halo_mod
   fms_mod
   mpp_mod

</PRE>
</A><!-- END MODULES_USED -->
<!--------------------------------------------------------------------->
<A NAME="INTERFACE">
<HR>
<H4>PUBLIC INTERFACE</H4>
<!-- BEGIN INTERFACE -->
<PRE>

   <b>use bgrid_masks_mod</b> [, only: grid_mask_type, mask_type,
                                grid_masks_init ]

   <a href="#DATA TYPES">grid_mask_type, mask_type</a>
        Data structures that contain the 3d step-mountain topography masks
        and 2d indexing for the lowest model level.

   <a href="#grid_masks_init">grid_masks_init</a>
        Initializes data with the grid_mask_type.

</PRE>
</A><!-- END INTERFACE -->
<!--------------------------------------------------------------------->
<A NAME="DATA TYPES">
<HR>
<H4>DATA TYPES</H4>
<!-- BEGIN DATA TYPES -->
<PRE>
   <b>type grid_mask_type</b>
      type(mask_type) :: Tmp, Vel
      logical :: sigma
   end type grid_mask_type

   Tmp = grid masking values for the temperature/mass grid
   Vel = grid masking values for the velocity/momentum grid
   sigma = logical flag that specific whether vertical coordinate is
             the step-mountain (eta) or sigma coordinate

   <b>type mask_type</b>
      real,    pointer, dimension(:,:,:) :: mask 
      integer, pointer, dimension(:,:)   :: kbot 
      integer                            :: kbotmin 
   end type mask_type

   mask  = step-mountain topography mask (0.0 or 1.0) for
             mass (height) grid points
   kbot  = lowest model level above ground

   note:  for the sigma coordinate, mask = 1.0 everywhere, and
          kbot = number of vertical model levels


</PRE>
</A><!-- END DATA TYPES -->
<!--------------------------------------------------------------------->
<A NAME="ROUTINES">
<HR>
<H4>PUBLIC ROUTINES</H4>
<!-- BEGIN ROUTINES -->
<PRE>
<a name="grid_masks_init">

Mask = <b>grid_masks_init</b> ( Hgrid, Vgrid, res )

INPUT
   Vgrid     The derived-type variable returned by a previous call to vert_grid_init.
                <FONT SIZE=-1 COLOR="#000099">[type(vert_grid_type)]</FONT>

   res       Reciporal of eta at the surface (i.e., the model interface that
             coincides with the step-mountain height). Note: for sigma
             coordinate model res=1. everywhere.
                <FONT SIZE=-1 COLOR="#000099">[real, dimension(:,:)]</FONT>

INPUT/OUTPUT
   Hgrid     The derived-type variable returned by a previous call to horiz_grid_init.
                See the module horiz_grid_mod for details.
                <FONT SIZE=-1 COLOR="#000099">[type(horiz_grid_type)]</FONT>

RETURNS
   Mask      The derived-type variable containing grid masking arrays.
                <FONT SIZE=-1 COLOR="#000099">[type(grid_mask_type)]</FONT>

</PRE>
</A><!-- END ROUTINES -->
<!--------------------------------------------------------------------->
<A NAME="ERRORS">
<HR>
<H4>ERROR MESSAGES</H4>
<!-- BEGIN ERRORS -->
<PRE>

 There are no error messages printed by this module.

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

   The interface grid_masks_init prints a message to STDOUT that
   describes the type of vertical coordinate that was initialized.

   If the eta coordinate is detected, a note of caution is printed.
   The eta coordinate is currently not supported.

   These messages are probably more appropiately printed from
   module bgrid_vert_mod (maybe in a future version).


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
