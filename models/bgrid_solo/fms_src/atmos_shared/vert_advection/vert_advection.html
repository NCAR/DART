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
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>module vert_advection_mod</title>
<link rel="stylesheet" href=
"http://www.gfdl.noaa.gov/~fms/style/doc.css" type="text/css">
<meta http-equiv="Content-Type" content="text/html; charset=utf-8">
</head>
<body>
<font size="1" class="header"><a href="#PUBLIC%20INTERFACE">PUBLIC
INTERFACE</a> ~ <a href="#PUBLIC%20ROUTINES">PUBLIC ROUTINES</a> ~
<a href="#ERROR%20MESSAGES">ERROR MESSAGES</a> ~ <a href=
"#REFERENCES">REFERENCES</a> ~ <a href="#NOTES">NOTES</a></font>
<hr>
<h1>module vert_advection_mod</h1>
<!-- END HEADER -->
<!--================================================================-->
<a name="OVERVIEW" id="OVERVIEW"></a>
<hr>
<h2>OVERVIEW</h2>
<!-- BEGIN OVERVIEW -->
<p class="text">Computes the tendency due to vertical advection for
an arbitrary quantity.</p>
<!-- END OVERVIEW -->
<!--================================================================-->
<a name="DESCRIPTION" id="DESCRIPTION"></a> 
<!-- BEGIN DESCRIPTION -->
<div>
<p>The advective tendency may be computed in <i>advective form</i>
(for use in spectral models) or <i>flux form</i>. The spatial
differencing may be <i>centered</i> (second or fourth order) or
<i>finite volume</i> (van Leer) using a piecewise linear
method.</p>
</div>
<!-- END DESCRIPTION -->
<!--================================================================-->
<a name="OTHER MODULES USED"></a>
<hr>
<h2>OTHER MODULES USED</h2>
<!-- BEGIN OTHER MODULES USED -->
<div>
<pre>
     fms_mod
</pre></div>
<!-- END OTHER MODULES USED -->
<!--================================================================-->
<a name="PUBLIC INTERFACE"></a>
<hr>
<h2>PUBLIC INTERFACE</h2>
<!-- BEGIN PUBLIC INTERFACE -->
<div>
<pre>

<b>use vert_advection_mod</b> [, only:  vert_advection,
                                 SECOND_CENTERED, FOURTH_CENTERED, VAN_LEER_LINEAR,
                                 FLUX_FORM, ADVECTIVE_FORM  ]

</pre>
<dl>
<dt><a href=
"vert_advection.html#vert_advection">vert_advection</a>:</dt>
<dd>Computes the vertical advective tendency for an arbitrary
quantity. There is no initialization routine necessary.</dd>
</dl>
</div>
<br>
<!-- END INTERFACE -->
<!--================================================================-->
<a name="PUBLIC ROUTINES"></a>
<hr>
<h2>PUBLIC ROUTINES</h2>
<!-- BEGIN PUBLIC ROUTINES -->
<ol type="a">
<li><a name="vert_advection" id="vert_advection"></a>
<h2>vert_advection</h2>
<pre>
<b>call vert_advection</b> ( dt, w, dz, r, rdt <span class=
"type">[, mask, scheme, form]</span> )

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
           where mask &gt; 0 for layers above ground
           <span class="type">[real, dimension(:,:,:)]</span>
           <span class="type">[real, dimension(:,:)]</span>
           <span class="type">[real, dimension(:)]</span>

   scheme  spatial differencing scheme, use one of these values:
              SECOND_CENTERED = second-order centered
              FOURTH_CENTERED = fourth-order centered
              VAN_LEER_LINEAR = piecewise linear, finite volume (van Leer)
           <span class=
"type">[integer, default=VAN_LEER_LINEAR]</span>

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

</pre></li>
</ol>
<!-- END PUBLIC ROUTINES -->
<!--================================================================-->
<a name="ERROR MESSAGES"></a>
<hr>
<h2>ERROR MESSAGES</h2>
<!-- BEGIN ERROR MESSAGES -->
<div>
<dl>
<dt><b>Errors in vert_advection_mod</b></dt>
<dd><span class="errmsg">vertical dimension of input arrays
inconsistent</span></dd>
<dd>The following was not true: size(w,3) = size(r,3)+1.<br>
<br></dd>
<dd><span class="errmsg">invalid value for optional argument
scheme</span></dd>
<dd>The value of optional argument scheme must be one of the public
parameters SECOND_CENTERED, FOURTH_CENTERED, or
VAN_LEER_LINEAR.<br>
<br></dd>
<dd><span class="errmsg">invalid value for optional argument
form</span></dd>
<dd>The value of optional argument form must be one of the public
parameters FLUX_FORM or ADVECTIVE_FORM.<br>
<br></dd>
</dl>
</div>
<br>
<!-- END ERROR MESSAGES -->
<!--================================================================-->
<a name="REFERENCES" id="REFERENCES"></a>
<hr>
<h2>REFERENCES</h2>
<!-- BEGIN REFERENCES -->
<div>
<ol>
<li>Lin, S.-J., W.C. Chao, Y.C. Sud, and G.K. Walker, 1994: A class
of the van Leer-type transport schemes and its application to the
moisture in a general circulation model. <i>Mon. Wea. Rev.</i>,
<b>122</b>, 1575-1593.</li>
</ol>
</div>
<br>
<!-- END REFERENCES -->
<!--================================================================-->
<a name="KNOWN BUGS"></a>
<hr>
<h2>KNOWN BUGS</h2>
<!-- BEGIN KNOWN BUGS -->
<div>None.</div>
<br>
<!-- END KNOWN BUGS -->
<!--================================================================-->
<a name="NOTES" id="NOTES"></a>
<hr>
<h2>NOTES</h2>
<!-- BEGIN NOTES -->
<div>None.</div>
<br>
<!-- END NOTES -->
<!--================================================================-->
<a name="FUTURE PLANS"></a>
<hr>
<h2>FUTURE PLANS</h2>
<!-- BEGIN FUTURE PLANS -->
<div>Add option for finite volume, piecewise parabolic
method.</div>
<br>
<!-- END FUTURE PLANS -->
<!--================================================================-->
<hr>
</body>
</html>
