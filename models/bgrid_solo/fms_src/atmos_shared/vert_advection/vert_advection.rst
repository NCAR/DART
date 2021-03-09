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
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">
<html>
<head>
    <title>module vert_advection_mod</title>
   <link rel="stylesheet" href="http://www.gfdl.noaa.gov/~fms/style/doc.css" type="text/css">
   <meta http-equiv="Content-Type" content="text/html; charset=EUC-JP">
</head>
<body>
<font size=1 class="header">
<a href="#PUBLIC INTERFACE">PUBLIC INTERFACE </a>~
<a href="#PUBLIC ROUTINES">PUBLIC ROUTINES </a>~
<a href="#ERROR MESSAGES">ERROR MESSAGES </a>~
<a href="#REFERENCES">REFERENCES </a>~ 
<a href="#NOTES">NOTES</a>
</font><hr>


<h2>module vert_advection_mod</h2>
<a name="HEADER"></a>
<!-- BEGIN HEADER -->
<div>
     <b>Contact:</b> &nbsp;  Bruce Wyman <br>
     <b>Reviewers:</b>&nbsp; <br>
     <b>Change History:&nbsp; </b><a HREF="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/atmos/shared/vert_advection/vert_advection.f90">WebCVS Log</a> <br>
     <b>Last Modified:</b>&nbsp; $Date$
</div><br>

<!-- END HEADER -->
<!-------------------------------------------------------------------->
<a name="OVERVIEW"></a>
<hr>
<h4>OVERVIEW</h4>
<!-- BEGIN OVERVIEW -->
<p class="text">

     Computes the tendency due to vertical advection
     for an arbitrary quantity.

</p>
<!-- END OVERVIEW -->
<!-------------------------------------------------------------------->
<a name="DESCRIPTION"></a>
<!-- BEGIN DESCRIPTION -->
<div>
   <p>  The advective tendency may be computed in <i>advective form</i>
     (for use in spectral models) or <i>flux form</i>.  
     The spatial differencing may be <i>centered</i> (second or fourth order)
     or <i>finite volume</i> (van Leer) using a piecewise linear method.</p>

</div>
<!-- END DESCRIPTION -->
<!-------------------------------------------------------------------->
<a name="OTHER MODULES USED"></a>
<hr>
<h4>OTHER MODULES USED</h4>
<!-- BEGIN OTHER MODULES USED -->
<div><pre>
     fms_mod
</pre></div>
<!-- END OTHER MODULES USED -->
<!-------------------------------------------------------------------->
<a name="PUBLIC INTERFACE"></a>
<hr>
<h4>PUBLIC INTERFACE</h4>
<!-- BEGIN PUBLIC INTERFACE -->
<div>
<pre>

<b>use vert_advection_mod</b> [, only:  vert_advection,
                                 SECOND_CENTERED, FOURTH_CENTERED, VAN_LEER_LINEAR,
                                 FLUX_FORM, ADVECTIVE_FORM  ]

</pre>
<dl>
 <dt>  <a href="vert_advection.html#vert_advection">vert_advection</a>:
     <dd>   Computes the vertical advective tendency for an arbitrary quantity.
        There is no initialization routine necessary.

</dl></div><br>
<!-- END INTERFACE -->
<!-------------------------------------------------------------------->
<a name="PUBLIC ROUTINES"></a>
<hr>
<h4>PUBLIC ROUTINES</h4>
<!-- BEGIN PUBLIC ROUTINES -->
<ol type="a">
<li><a name="vert_advection"></a><h4>vert_advection</h4>
<pre>
<b>call vert_advection</b> ( dt, w, dz, r, rdt <span class="type">[, mask, scheme, form]</span> )

DESCRIPTION
   This routine computes the vertical advective tendency for
   an arbitrary quantity. The tendency can be computed using
   one of several different choices for spatial differencing
   and numerical form of the equations. These choices are
   controlled through optional arguments. 
   There is no initialization routine necessary.

INPUT
   dt   time step in seconds <span class="type">[real]</span>

   w    advecting velocity at the vertical boundaries of the grid boxes
        does not assume velocities at top and bottom are zero
        units = [units of dz / second]
        <span class="type">[real, dimension(:,:,:)]</span>
        <span class="type">[real, dimension(:,:)]</span>
        <span class="type">[real, dimension(:)]</span>

   dz   depth of model layers in arbitrary units (usually pressure)
        <span class="type">[real, dimension(:,:,:)]</span>
        <span class="type">[real, dimension(:,:)]</span>
        <span class="type">[real, dimension(:)]</span>

   r    advected quantity in arbitrary units
        <span class="type">[real, dimension(:,:,:)]</span>
        <span class="type">[real, dimension(:,:)]</span>
        <span class="type">[real, dimension(:)]</span>

OUTPUT
   rdt  advective tendency for quantity "r" weighted by depth of layer
        units = [units of r * units of dz / second]
        <span class="type">[real, dimension(:,:,:)]</span>
        <span class="type">[real, dimension(:,:)]</span>
        <span class="type">[real, dimension(:)]</span>

OPTIONAL INPUT
   mask    mask for below ground layers,
           where mask > 0 for layers above ground
           <span class="type">[real, dimension(:,:,:)]</span>
           <span class="type">[real, dimension(:,:)]</span>
           <span class="type">[real, dimension(:)]</span>

   scheme  spatial differencing scheme, use one of these values:
              SECOND_CENTERED = second-order centered
              FOURTH_CENTERED = fourth-order centered
              VAN_LEER_LINEAR = piecewise linear, finite volume (van Leer)
           <span class="type">[integer, default=VAN_LEER_LINEAR]</span>

   form    form of equations, use one of these values:
              FLUX_FORM      = solves for -d(w*r)/dt
              ADVECTIVE_FORM = solves for -w*d(r)/dt
           <span class="type">[integer, default=FLUX_FORM]</span>

NOTE
   The input/output arrays may be 1d, 2d, or 3d.
   The last dimension is always the vertical dimension.
   For the vertical dimension the following must be true:
   size(w,3) == size(dz,3)+1 == size(r,3)+1 == size(rdt,3)+1 == size(mask,3)+1
   All horizontal dimensions must have the same size (no check is done).

</pre>
</li>
</ol>
<!-- END PUBLIC ROUTINES -->
<!-------------------------------------------------------------------->
<a name="ERROR MESSAGES"></a>
<hr>
<h4>ERROR MESSAGES</h4>
<!-- BEGIN ERROR MESSAGES -->
<div>
<dl>
<dt><b>Errors in vert_advection_mod</b>
    <dd><span class="errmsg">vertical dimension of input arrays inconsistent</span></dd>
    <dd>  The following was not true: size(w,3) = size(r,3)+1.
    <br><br></dd>
    
    <dd><span class="errmsg">invalid value for optional argument scheme</span></dd>
    <dd>  The value of optional argument scheme must be one of the
        public parameters SECOND_CENTERED, FOURTH_CENTERED, or VAN_LEER_LINEAR.
    <br><br></dd>
    
    <dd><span class="errmsg">invalid value for optional argument form</span></dd>
    <dd>  The value of optional argument form must be one of
        the public parameters FLUX_FORM or ADVECTIVE_FORM.
    <br><br></dd>
</dl>
</div><br>
<!-- END ERROR MESSAGES -->
<!-------------------------------------------------------------------->
<a name="REFERENCES"></a>
<hr>
<h4>REFERENCES</h4>
<!-- BEGIN REFERENCES -->
<div>
<ol>

<li>     Lin, S.-J., W.C. Chao, Y.C. Sud, and G.K. Walker, 1994:
     A class of the van Leer-type transport schemes and its application
     to the moisture in a general circulation model. <i>Mon. Wea. Rev.</i>,
     <b>122</b>, 1575-1593.</li>

</ol>
</div><br>
<!-- END REFERENCES -->
<!-------------------------------------------------------------------->
<a name="KNOWN BUGS"></a>
<hr>
<h4>KNOWN BUGS</h4>
<!-- BEGIN KNOWN BUGS -->
<div>

     None.

</div><br>
<!-- END KNOWN BUGS -->
<!-------------------------------------------------------------------->
<a name="NOTES"></a>
<hr>
<h4>NOTES</h4>
<!-- BEGIN NOTES -->
<div>

     None.

</div><br>
<!-- END NOTES -->
<!-------------------------------------------------------------------->
<a name="FUTURE PLANS"></a>
<hr>
<h4>FUTURE PLANS</h4>
<!-- BEGIN FUTURE PLANS -->
<div>

      Add option for finite volume, piecewise parabolic method.

</div><br>
<!-- END FUTURE PLANS -->
<!-------------------------------------------------------------------->

<hr>
</body>
</html>
