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
<title>module atmosphere_mod</title>
</head>
<body bgcolor="#AABBCC" text="#332211">
<div align="center"><font size="1"><a href="#INTERFACE">PUBLIC
INTERFACE</a> / <a href="#ROUTINES">ROUTINES</a> / <a href=
"#NAMELIST">NAMELIST</a> / <a href="#ERRORS">ERRORS</a></font><br>
<br></div>
<hr>
<h1>module atmosphere_mod</h1>
<a name="OVERVIEW" id="OVERVIEW">
<hr>
<h2>OVERVIEW</h2>
<!-- BEGIN OVERVIEW -->
<pre>

     Atmospheric driver for the dry B-grid dynamical core and Held-Suarez
     benchmark (aka simple_phyhsics).
</pre></a><!-- END OVERVIEW -->
<!-- ------------------------------------------------------------------>
 <a name="DESCRIPTION" id="DESCRIPTION"><!-- BEGIN DESCRIPTION -->
<pre>
     This module provides a standard interface to the B-grid dynamical
     core and B-grid interface to the simple physics.
     A similar interface for the spectral dynamical core exists that
     may be easily switched with this interface.

</pre></a><!-- END DESCRIPTION -->
<!-- ------------------------------------------------------------------>
 <a name="OTHER MODULES USED">
<hr>
<h2>OTHER MODULES USED</h2>
<!-- BEGIN OTHER MODULES USED -->
<pre>

   bgrid_horiz_mod
   bgrid_vert_mod
   bgrid_prog_var_mod
   bgrid_halo_mod
   bgrid_grid_change_mod
   bgrid_core_driver_mod
   hs_forcing_mod
   time_manager_mod
   fms_mod

</pre></a><!-- END MODULES_USED -->
<!-- ------------------------------------------------------------------>
 <a name="PUBLIC INTERFACE">
<hr>
<h2>PUBLIC INTERFACE</h2>
<!-- BEGIN PUBLIC INTERFACE -->
<pre>

  <b>use atmosphere_mod</b> [,only: atmosphere_init,       atmosphere_end,
                             atmosphere,
                             atmosphere_resolution, atmosphere_boundary,
                             get_atmosphere_axes  ]

  NOTES

     1)  Optional namelist interface <b>&amp;atmosphere_nml</b> may be
         read from file <b>input.nml</b>.
                                
</pre></a><!-- END PUBLIC INTERFACE -->
<!-- ------------------------------------------------------------------>
 <a name="PUBLIC ROUTINES">
<hr>
<h2>PUBLIC ROUTINES</h2>
<!-- BEGIN PUBLIC ROUTINES -->
<pre>

<b>call atmosphere_init</b> ( Time_init, Time, Time_step )

DESCRIPTION
   Initialization call for running the bgrid dynamical with the
   Held-Suarez GCM forcing.

INPUT
   Time_init   The initial (or base) time.  <font size="-1" color=
"#000099">[time_type]</font>

   Time        The current time.  <font size="-1" color=
"#000099">[time_type]</font>

   Time_step   The atmospheric model/physics time step.  <font size="-1"
color="#000099">[time_type]</font>

<!------------------------></pre>
<hr width="50%" align="left">
<pre>

<b>call atmosphere_end</b>

DESCRIPTION
   Termination call for the bgrid dynamical with Held-Suarez
   GCM forcing.  There are no arguments to this routine.

<!------------------------></pre>
<hr width="50%" align="left">
<pre>

<b>call atmosphere</b> ( Time )

DESCRIPTION
   Advances the B-grid prognostic variables one time step forward.
   The dynamical core, Held-Suarez forcing, diagnostics, and time
   differencing are all called.  This routine should only be called
   once per time step.

INPUT
   Time    The current time.  <font size="-1" color=
"#000099">[time_type]</font>

NOTE
   The prognostic variables are stored in a private derived-type 
   variabledefined in the header section of this module.

---------------------------------------------------------------------

<b>call get_atmosphere_axes</b> ( axes )

OUTPUT

   axes        axis identifiers for the atmospheric grids
               The size of axes at least 1 but not greater than 4.
               The axes returned are ordered (/ x, y, p_full, p_half /).
                 [integer, dimension(:)]

<!------------------------></pre>
<hr width="50%" align="left">
<pre>

<b>call atmosphere_resolution</b> ( nlon, nlat <font color=
"#007700">[, global]</font> )

DESCRIPTION
   Returns the resolution of compute domain for either the
   current processor or the global domain.

OUTPUT
   nlon   The number of longitude points in the compute domain.
             <font size="-1" color="#000099">[integer]</font>

   nlat   The number of latitude points in the compute domain.
             <font size="-1" color="#000099">[integer]</font>

OPTIONAL INPUT

   global  Flag that specifies whether the returned compute domain size is
           for the global grid (TRUE) or for the current processor (FALSE).
              <font size="-1" color=
"#000099">[logical, default: FALSE]</font>
           
<!------------------------></pre>
<hr width="50%" align="left">
<pre>

<b>call atmosphere_boundary</b> ( blon, blat <font color=
"#007700">[, global]</font> )

DESCRIPTION
   Returns the grid box edges of compute domain for either the
   current processor or the global domain.

OUTPUT
   blon    The west-to-east longitude edges of grid boxes (in radians).
              <font size="-1" color=
"#000099">[real, dimension(nlon+1)]</font>

   blat    The south-to-north latitude edges of grid boxes (in radians).
              <font size="-1" color=
"#000099">[real, dimension(nlat+1)]</font>

OPTIONAL INPUT
   global  Flag that specifies whether the returned grid box edges are
           for the global grid (TRUE) or for the current processor (FALSE).
              <font size="-1" color=
"#000099">[logical, default: FALSE]</font>
           
NOTE
   The size of the output arguments, blon and blat, must be +1 more than the
   output arguments for call atmosphere_resolution, nlon+1 and nlat+1, respectively.

</pre></a><!-- END PUBLIC ROUTINES -->
<!-- ------------------------------------------------------------------>
 <a name="NAMELIST" id="NAMELIST">
<hr>
<h2>NAMELIST</h2>
<!-- BEGIN NAMELIST -->
<pre>

<b>&amp;atmosphere_nml</b>

 physics_window  The number of "i" and "j" rows processed each time
                 the modular physics is called. To process the entire
                 domain use physics_window = 0,0.
                    <font size="-1" color=
"#000099">[integer, default: physics_window = 0,0]</font>

</pre></a><!-- END NAMELIST -->
<!-- ------------------------------------------------------------------>
 <a name="ERRORS" id="ERRORS">
<hr>
<h2>ERROR MESSAGES</h2>
<!-- BEGIN ERRORS -->
<pre>

<b>FATAL errors from get_atmosphere_axes in atmosphere_mod</b>

    <b>size of argument is incorrect</b>
        The size of the argument to get_atmosphere_axes must be
        between 1 and 4.


</pre></a><!-- END ERRORS -->
<!-- ------------------------------------------------------------------>
<hr>
</body>
</html>
