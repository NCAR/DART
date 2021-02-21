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
<title>module bgrid_vert_mod</title>
</head>
<body bgcolor="#AABBCC" text="#332211">
<div align="center"><font size="-2"><a href="#INTERFACE">PUBLIC
INTERFACE</a> / <a href="#DATA_TYPES">DATA</a> / <a href=
"#ROUTINES">ROUTINES</a> / <a href="#ERRORS">ERRORS</a> / <a href=
"#REFERENCES">REFERENCES</a> / <a href=
"#NOTES">NOTES</a></font><br>
<br></div>
<hr>
<h1>module bgrid_vert_mod</h1>
<a name="OVERVIEW" id="OVERVIEW">
<hr>
<h2>OVERVIEW</h2>
<!-- BEGIN OVERVIEW -->
<pre>

     Initializes the derived-type variable vert_grid_type,
     which contains vertical grid constants needed by the
     B-grid dynamical core and many other B-grid routines.

</pre></a><!-- END OVERVIEW -->
<!--------------------------------------------------------------------->
 <a name="DESCRIPTION" id="DESCRIPTION"><!-- BEGIN DESCRIPTION -->
<pre>
 
     This module initializes the constants needed for a
     hybrid sigma/pressure vertical coordinate that also
     includes a step-mountain (eta) coordinate option.

The B-grid model can be run with either "sigma" (terrain following)
coordinates or "eta" (the step-mountain coordinate).  In the sigma
coordinate system, the coordinate surfaces follow the topography,
there is always the same number of model levels above ground at all
grid boxes.  In the eta coordinate system, the coordinate surfaces
are nearly horizontal, model grid boxes will be beneath the topography.
The step-mountain topography can be thought of as "building blocks".
The model levels below the topography are masked out.

In both coordinate systems, vertical indexing increases from the
top of the atmosphere towards the surface.

In addition, a hybrid coordinate may be specified.  Near the surface
the coordinate must be a sigma/eta coordinate, with a transition
between arbitrary levels to a pressure coordinate in the uppermost
model levels.

     Pressure p is determined by the following formula:

            p(k) = peta(k) + eta(k) * pssl

     where   (k) = the vertical index
             peta = reference pressure
             eta  = sigma/eta values (between 0. and 1.)
             pssl = ps * res
             ps   = surface pressure
             res  = 1./eta_surface (=1 for non-eta case)

             peta, eta      are dimensioned only in the vertical
             pssl, ps, res  are dimensioned only in the horizontal


</pre></a><!-- END DESCRIPTION -->
<!--------------------------------------------------------------------->
 <a name="MODULES_USED" id="MODULES_USED">
<hr>
<h2>OTHER MODULES USED</h2>
<!-- BEGIN MODULES_USED -->
<pre>

     constants_mod
     fms_mod

</pre></a><!-- END MODULES_USED -->
<!--------------------------------------------------------------------->
 <a name="INTERFACE" id="INTERFACE">
<hr>
<h2>PUBLIC INTERFACE</h2>
<!-- BEGIN INTERFACE -->
<pre>

   <b>use bgrid_vert_mod</b> [,only: vert_grid_type      ,
                              vert_grid_init      ,
                              compute_pres_depth  ,
                              compute_pres_full   ,
                              compute_pres_half   ,
                              compute_pres_weights,
                              compute_pressures   ,
                              compute_geop_height ,
                              compute_height      ]

   <a href="#DATA_TYPES">vert_grid_type</a>
        Defined-type variable containing the vertical grid constants
        needed by the B-grid dynamical core (see below).

   <a href="#vert_grid_init">vert_grid_init</a>
        Initializes a vert_grid_type.

   <a href="#compute_pres_depth">compute_pres_depth</a>
        Computes the pressure weight (mass) of a model layer.

   <a href="#compute_pres_full">compute_pres_full</a>
        Computes the pressure at full model levels.

   <a href="#compute_pres_half">compute_pres_half</a>
        Computes the pressure at half model levels.

   <a href="#compute_pres_weights">compute_pres_weights</a>
        Computes weighting terms for averaging data from half levels to
        full model levels.

   <a href="#compute_pressures">compute_pressures</a>
        Computes everything from routines compute_pres_depth,
        compute_pres_full, compute_pres_half, and compute_pres_weights.

   <a href="#compute_geop_height">compute_geop_height</a>
        Computes the geopotential height (in m2/s2) at full model levels
        and (optionally) at half model levels.

   <a href="#compute_height">compute_height</a>
        Computes the height (in meters) at full and half model levels.

   <a href="#compute_height_bottom">compute_height_bottom</a>
        Computes the height (in meters) and pressure at the lowest
        full model level.

</pre></a><!-- END INTERFACE -->
<!--------------------------------------------------------------------->
 <a name="DATA_TYPES" id="DATA_TYPES">
<hr>
<h2>PUBLIC DATA</h2>
<!-- BEGIN DATA_TYPES -->
<pre>

 <b>type (vert_grid_type)</b>

   nlev  = number of vertical levels
             [integer]
   nplev = number of pure pressure levels at the top of the model
             [integer]
   deta  = vertical eta thickness/depth of model layers 
             [real, dimension(nlev)]
   aeta  = eta values at model levels (full levels)
             [real, dimension(nlev)]
   eta   = eta values at model layer interfaces (half levels)
             [real, dimension(nlev+1)]
   dpeta = vertical pressure thickness/depth of model layers
             [real, dimension(nlev)]
   apeta = pressure values at model levels (full levels)
             [real, dimension(nlev)]
   peta  = pressure values at model layer interfaces (half levels)
             [real, dimension(nlev+1)]
   dfl   = reference values of geopotential height at half levels (m2/s2)
             [real, dimension(nlev+1)]
   wta,
     wtb = reference values for the weighting terms at full levels
             [real, dimension(nlev)]

   psmin = minimum allowable surface pressure for the hybrid
           coordinate, below this value model layers will
           have negative pressure weights
             [real]
   hybrid = logical flag that indicates if the hybrid coordinate is
            being used
   pzero  = logical flag that indicates if the pressure at the top of
            the model is zero (can only be non-zero if hybrid=true).

   pref   = reference pressure for the eta coordinate (101325. Pa)
              [real]
   tref   = reference temperature for the eta coordinate (288. degK)
              [real]
   gamma  = reference lapse rate for the eta coordinate (.0065 degK/m)
              [real]

</pre></a><!-- END DATA_TYPES -->
<!--------------------------------------------------------------------->
 <a name="ROUTINES" id="ROUTINES">
<hr>
<h2>PUBLIC ROUTINES</h2>
<!-- BEGIN ROUTINES -->
<pre>
<a name="vert_grid_init" id="vert_grid_init">

Vgrid = <b>vert_grid_init</b> (eta, peta, verbose) INPUT eta The sigma/eta values at model layer interfaces (half levels). The size of this array will determine the vertical resolution (i.e., number of levels = size(eta)-1). [real, dimension(:)] OPTIONAL INPUT peta A profile of reference pressure values at model layer interfaces (half levels) used to define the HYBRID vertical coordinate. If all values of "peta" are zero (the default), then the vertical coordinate used will be a pure sigma/eta coordinate. This variable along with "eta" determine the pressure at half model levels. The size of this array must be the same as eta. [real, dimension(:), default: peta=0.] verbose Integer flag that controls the amount of printed output. [integer, default: verbose = 0] RETURNS Vgrid Derived-type variable that contains all necessary vertical grid information needed by the dynamical core. [type(vert_grid_type)] ------------------------------------------------------------------------ </a><a name="compute_pres_depth"
id="compute_pres_depth">

<b>call compute_pres_depth</b> ( Vgrid, pssl, pdepth ) INPUT Vgrid Derived-type variable containing vertical grid information returned by a previous call to vert_grid_init. [type(vert_grid_type)] pssl The (surface) pressure at eta = 1. [real, dimension(:,:)] OUTPUT pdepth The pressure depth (i.e., weight) for all model layers. The first 2 dimensions of pdepth must be the same as the dimensions of pssl. The third dimension of pdepth must equal the number of full model levels. [real, dimension(:,:,nlev)] ------------------------------------------------------------------------ </a><a name="compute_pres_full"
id="compute_pres_full">

<b>call compute_pres_full</b> ( Vgrid, pssl, pfull [,phalf, dpde] ) INPUT Vgrid Derived-type variable containing vertical grid information returned by a previous call to vert_grid_init. [type(vert_grid_type)] pssl The (surface) pressure at eta = 1. [real, dimension(:,:)] OUTPUT pfull The pressure at full model levels. The first 2 dimensions of pfull must be the same as the dimensions of pssl. The third dimension of pfull must equal the number of full model levels. [real, dimension(:,:,nlev)] OPTIONAL INPUT phalf The pressure at half model levels (the interface between model layers). The first 2 dimensions of phalf must be the same as the dimensions of pssl. The third dimension of phalf must equal the number of half model levels. [real, dimension(:,:,nlev+1)] dpde The pressure depth (i.e., weight) for all model layers. The first 2 dimensions of dpde must be the same as the dimensions of pssl. The third dimension of dpde must equal the number of full model levels. [real, dimension(:,:,nlev)] ------------------------------------------------------------------------ </a><a name="compute_pres_half"
id="compute_pres_half">

<b>call compute_pres_half</b> ( Vgrid, pssl, phalf ) INPUT Vgrid Derived-type variable containing vertical grid information returned by a previous call to vert_grid_init. [type(vert_grid_type)] pssl The (surface) pressure at eta = 1. [real, dimension(:,:)] OUTPUT phalf The pressure at half model levels (the interface between model layers). The first 2 dimensions of phalf must be the same as the dimensions of pssl. The third dimension of phalf must equal the number of half model levels. [real, dimension(:,:,nlev+1)] ------------------------------------------------------------------------ </a><a name="compute_pres_weights"
id="compute_pres_weights">

<b>call compute_pres_weights</b> ( Vgrid, phalf, pfull, wta, wtb ) INPUT Vgrid Derived-type variable containing vertical grid information returned by a previous call to vert_grid_init. [type(vert_grid_type)] phalf The pressure at half model levels (the interface between model layers). [real, dimension(:,:,nlev+1)] pfull The pressure at full model levels. [real, dimension(:,:,nlev)] OUTPUT wta Weighting term used to average data from half levels to full model levels. wta is applied to data at the half level ABOVE. [real, dimension(:,:,nlev)] wtb Weighting term used to average data from half levels to full model levels. wtb is applied to data at the half level BELOW. [real, dimension(:,:,nlev)] NOTES For input and output arrays: phalf, pfull, wta, wtb. 1) The first 2 dimensions must all be the same. 2) The third dimension must agree with the vertical resolution of the model (as stated above). ------------------------------------------------------------------------ </a><a name="compute_pressures"
id="compute_pressures">

<b>call compute_pressures</b> ( Vgrid, pssl, phalf, pfull [, dpde, wta, wtb] ) INPUT Vgrid Derived-type variable containing vertical grid information returned by a previous call to vert_grid_init. [type(vert_grid_type)] pssl The (surface) pressure at eta = 1. [real, dimension(:,:)] OUTPUT phalf The pressure at half model levels (the interface between model layers). [real, dimension(:,:,nlev+1)] pfull The pressure at full model levels. [real, dimension(:,:,nlev)] OPTIONAL OUTPUT dpde The pressure depth (i.e., weight) for all model layers. [real, dimension(:,:,nlev)] wta Weighting term used to average data from half levels to full model levels. wta is applied to data at the half level ABOVE. [real, dimension(:,:,nlev)] wtb Weighting term used to average data from half levels to full model levels. wtb is applied to data at the half level BELOW. [real, dimension(:,:,nlev)] NOTES For output arrays: phalf, pfull, dpde, wta, wtb. 1) The first 2 dimensions must be the same as the dimensions of input argument pssl. 2) The third dimension must agree with the vertical resolution of the model (as stated above). ------------------------------------------------------------------------ </a><a name="compute_geop_height"
id="compute_geop_height">

<b>call compute_geop_height</b> ( Vgrid, fssl, vtemp, wta, wtb, zfull [, zhalf, mask] ) INPUT Vgrid Derived-type variable containing vertical grid information returned by a previous call to vert_grid_init. [type(vert_grid_type)] fssl The geopotential height (in m2/s2) at eta = 1. Note that for the eta coordinate fssl = 0. [real, dimension(:,:)] vtemp The virtual temperature (degK) at full model levels. [real, dimension(:,:,nlev)] wta Weighting term used to average data from half levels to full model levels. wta is applied to data at the half level ABOVE. [real, dimension(:,:,nlev)] wtb Weighting term used to average data from half levels to full model levels. wtb is applied to data at the half level BELOW. [real, dimension(:,:,nlev)] OUTPUT zfull The geopotential height (in m2/s2) at full model levels. [real, dimension(:,:,nlev)] OPTIONAL OUTPUT zhalf The geopotential height (in m2/s2) at half model levels (the interface between model layers). [real, dimension(:,:,nlev+1)] OPTIONAL INPUT mask A topography mask (0. or 1.) at full model levels for the eta/step-mountain vertical coordinate. [real, dimension(:,:,nlev)] NOTES For input and output arrays: vtemp, wta, wtb, zfull, zhalf, mask. 1) The first 2 dimensions must all be the same as the dimensions of input argument fssl. 2) The third dimension must agree with the vertical resolution of the model (as stated above). ------------------------------------------------------------------------ </a><a name="compute_height"
id="compute_height">

<b>call compute_height</b> ( Vgrid, fssl, temp, sphum, pfull, phalf, zfull, zhalf [, mask] ) INPUT Vgrid Derived-type variable containing vertical grid information returned by a previous call to vert_grid_init. [type(vert_grid_type)] fssl The geopotential height (in m2/s2) at eta = 1. Note that for the eta coordinate fssl = 0. [real, dimension(:,:)] temp The temperature (degK) at full model levels. [real, dimension(:,:,nlev)] sphum The specific humidity (Kg/Kg) at full model levels. [real, dimension(:,:,nlev)] pfull The pressure at full model levels. [real, dimension(:,:,nlev)] phalf The pressure at half model levels (the interface between model layers). [real, dimension(:,:,nlev+1)] OUTPUT zfull The geopotential height (in meters) at full model levels. [real, dimension(:,:,nlev)] zhalf The geopotential height (in meters) at half model levels (the interface between model layers). [real, dimension(:,:,nlev+1)] OPTIONAL INPUT mask A topography mask (0. or 1.) at full model levels for the eta/step-mountain vertical coordinate. [real, dimension(:,:,nlev)] NOTES For input and output arrays: temp, sphum, pfull, phalf, zfull, zhalf, mask. 1) The first 2 dimensions must all be the same as the dimensions of input argument fssl. 2) The third dimension must agree with the vertical resolution of the model (as stated above). ------------------------------------------------------------------------ </a><a name="compute_height_bottom"
id="compute_height_bottom">

<b>call compute_height_bottom</b> ( Vgrid, pssl, tbot, qbot, zbot, pbot [, kbot] ) INPUT Vgrid Derived-type variable containing vertical grid information returned by a previous call to vert_grid_init. [type(vert_grid_type)] pssl The (surface) pressure at eta = 1. [real, dimension(:,:)] tbot The temperature (degK) at the lowest full model level. [real, dimension(:,:)] qbot The specific humidity (Kg/Kg) at the lowest full model level. [real, dimension(:,:)] OUTPUT zbot Height (in meters) between the surface and the lowest full model level. [real, dimension(:,:)] pbot Pressure (in pascals) between the surface and the lowest full model level. [real, dimension(:,:)] OPTIONAL INPUT kbot Vertical index of the model closest to the surface. This argument does not need to be passed for the sigma coordinate system (at all grid points it will equal the number of levels as determined from Vgrid). [integer, dimension(:,:)] NOTES 1) The size of all 2D arrays must all be the same. 2) The number of vertical levels is determined from Vgrid. </a></pre></a><!-- END ROUTINES -->
<!--------------------------------------------------------------------->
 <a name="ERRORS" id="ERRORS">
<hr>
<h2>ERROR MESSAGES</h2>
<!-- BEGIN ERRORS -->
<pre>

<b>Fatal errors in compute_pres_depth</b>

     <b>incorrect dimension 3 for pdepth</b>
          The third dimension of output argument pdepth must equal
          the number of model levels.

     <b>pressure depth &lt;= 0.0</b>
          The pressure depth of a model layer has a mass &lt;= 0.
          This can only happen when running a hybrid vertical
          coordinate and the surface pressure becomes too small.
          This error may result as the model goes unstable.

<b>Fatal errors in compute_pres_full or compute_pressures</b>

     <b>incorrect dimension 3 for pfull</b>
          The third dimension of output argument pfull must equal
          the number of model levels.

<b>Fatal errors in compute_pres_half</b>

     <b>incorrect dimension 3 for phalf</b>
          The third dimension of output argument phalf must equal
          the number of model levels plus one.

<b>Fatal errors in compute_geop_height</b>

     <b>incorrect dimension 3 for zfull</b>
          The third dimension of output argument zfull must equal
          the number of model levels.

     <b>incorrect dimension 3 for zhalf</b>
          The third dimension of output argument zhalf must equal
          the number of model levels plus one.

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

    There are no known bugs.
   
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
