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
<title>module bgrid_masks_mod</title>
</head>
<body bgcolor="#AABBCC" text="#332211">
<div align="center"><font size="1"><a href="#INTERFACE">PUBLIC
INTERFACE</a> / <a href="#ROUTINES">ROUTINES</a> / <a href=
"#ERRORS">ERRORS</a> / <a href="#NOTES">NOTES</a></font><br>
<br></div>
<hr>
<h1>module bgrid_masks_mod</h1>
<!--------------------------------------------------------------------->
<a name="OVERVIEW" id="OVERVIEW">
<hr>
<h2>OVERVIEW</h2>
<!-- BEGIN OVERVIEW -->
<pre>

     Provides a data structure for three-dimensional masks used
     to define the step-mountain/eta coordinate topography.

</pre></a><!-- END OVERVIEW -->
<!--------------------------------------------------------------------->
 <a name="DESCRIPTION" id="DESCRIPTION"><!-- BEGIN DESCRIPTION -->
<pre>

</pre></a><!-- END DESCRIPTION -->
<!--------------------------------------------------------------------->
 <a name="MODULES_USED" id="MODULES_USED">
<hr>
<h2>OTHER MODULES USED</h2>
<!-- BEGIN MODULES_USED -->
<pre>

   bgrid_horiz_mod
   bgrid_vert_mod
   bgrid_halo_mod
   fms_mod
   mpp_mod

</pre></a><!-- END MODULES_USED -->
<!--------------------------------------------------------------------->
 <a name="INTERFACE" id="INTERFACE">
<hr>
<h2>PUBLIC INTERFACE</h2>
<!-- BEGIN INTERFACE -->
<pre>

   <b>use bgrid_masks_mod</b> [, only: grid_mask_type, mask_type,
                                grid_masks_init ]

   <a href="#DATA%20TYPES">grid_mask_type, mask_type</a>
        Data structures that contain the 3d step-mountain topography masks
        and 2d indexing for the lowest model level.

   <a href="#grid_masks_init">grid_masks_init</a>
        Initializes data with the grid_mask_type.

</pre></a><!-- END INTERFACE -->
<!--------------------------------------------------------------------->
 <a name="DATA TYPES">
<hr>
<h2>DATA TYPES</h2>
<!-- BEGIN DATA TYPES -->
<pre>
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


</pre></a><!-- END DATA TYPES -->
<!--------------------------------------------------------------------->
 <a name="ROUTINES" id="ROUTINES">
<hr>
<h2>PUBLIC ROUTINES</h2>
<!-- BEGIN ROUTINES -->
<pre>
<a name="grid_masks_init" id="grid_masks_init">

Mask = <b>grid_masks_init</b> ( Hgrid, Vgrid, res ) INPUT Vgrid The derived-type variable returned by a previous call to vert_grid_init. <font size="-1"
color=
"#000099">[type(vert_grid_type)]</font> res Reciporal of eta at the surface (i.e., the model interface that coincides with the step-mountain height). Note: for sigma coordinate model res=1. everywhere. <font size="-1"
color=
"#000099">[real, dimension(:,:)]</font> INPUT/OUTPUT Hgrid The derived-type variable returned by a previous call to horiz_grid_init. See the module horiz_grid_mod for details. <font size="-1"
color=
"#000099">[type(horiz_grid_type)]</font> RETURNS Mask The derived-type variable containing grid masking arrays. <font size="-1"
color="#000099">[type(grid_mask_type)]</font> </a></pre></a> 
<!-- END ROUTINES -->
<!--------------------------------------------------------------------->
 <a name="ERRORS" id="ERRORS">
<hr>
<h2>ERROR MESSAGES</h2>
<!-- BEGIN ERRORS -->
<pre>

 There are no error messages printed by this module.

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

   The interface grid_masks_init prints a message to STDOUT that
   describes the type of vertical coordinate that was initialized.

   If the eta coordinate is detected, a note of caution is printed.
   The eta coordinate is currently not supported.

   These messages are probably more appropiately printed from
   module bgrid_vert_mod (maybe in a future version).


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
