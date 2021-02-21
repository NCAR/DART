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
<title>module bgrid_prog_var_mod</title>
</head>
<body bgcolor="#AABBCC" text="#332211">
<div align="center"><font size="-2"><a href="#INTERFACE">PUBLIC
INTERFACE</a> / <a href="#DATA_TYPES">DATA</a> / <a href=
"#ROUTINES">ROUTINES</a> / <a href="#ERRORS">ERRORS</a> / <a href=
"#NOTES">NOTES</a></font><br>
<br></div>
<hr>
<h1>module bgrid_prog_var_mod</h1>
<a name="OVERVIEW" id="OVERVIEW">
<hr>
<h2>OVERVIEW</h2>
<!-- BEGIN OVERVIEW -->
<pre>

     Initializes a derived-type variable that contains prognostic
     variables for the B-grid dynamical core.

</pre></a><!-- END OVERVIEW -->
<!--------------------------------------------------------------------->
 <a name="DESCRIPTION" id="DESCRIPTION"><!-- BEGIN DESCRIPTION -->
<pre>
     This modules defines and initializes a derived-type variable
     (prog_var_type) that contains the following prognostic fields:
       surface pressure, temperature, momentum, and tracers.

     The main initialization interface allocates storage for
     a prog_var_type variable. Additional interfaces allocate
     storage for separate fields. All array storage includes the
     halo region when the allocation is done.

     There are also operators and interfaces for performing some
     simple operations between prog_var_type variables.
      
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
   bgrid_cold_start_mod
   fms_mod
   field_manager_mod
   tracer_manager_mod
 
</pre></a><!-- END MODULES_USED -->
<!--------------------------------------------------------------------->
 <a name="INTERFACE" id="INTERFACE">
<hr>
<h2>PUBLIC INTERFACE</h2>
<!-- BEGIN INTERFACE -->
<pre>

   <b>use bgrid_prog_var_mod</b> [, only: prog_var_type, 
                                   prog_var_init,
                                   var_init,
                                   prog_var_times_scalar,
                                   prog_var_equals_scalar,
                                   prog_var_time_diff,
                                   open_prog_var_file,
                                   read_prog_var,
                                   write_prog_var   ]

   <a href="#DATA_TYPES">prog_var_type</a>
        derived-type that contains horizontal and vertical
        grid sizes, and the prognostic variables for
        pressure, momentum, temperature, and tracers.
 
   <a href="#prog_var_init">prog_var_init</a>
        initializes prog_var_type variable (all data arrays are set to zero)
 
   <a href="#var_init">var_init</a>
        initializes storage for real arrays (sets to zero)
 
   <a href="#prog_var_times_scalar">prog_var_times_scalar</a>
        multiplies the prognostic variables within a prog_var_type
        variable by a constant real scalar
 
   <a href="#prog_var_time_diff">prog_var_time_diff</a>
        performs explicit time differencing
 
   <a href="#prog_var_equals_scalar">prog_var_equals_scalar</a>
        assigns a scalar to all prognostic variables in a prog_var_type
 
   <a href="#open_prog_var_file">open_prog_var_file</a>
        opens the restart file for bgrid prognostic variables
 
   <a href="#read_prog_var">read_prog_var</a>
        reads the restart file for bgrid prognostic variables
 
   <a href="#write_prog_var">write_prog_var</a>
        reads the restart file for bgrid prognostic variables
 
</pre></a><!-- END INTERFACE -->
<!--------------------------------------------------------------------->
 <a name="DATA_TYPES" id="DATA_TYPES">
<hr>
<h2>PUBLIC DATA</h2>
<!-- BEGIN DATA_TYPES -->
<pre>

<b>type prog_var_type</b>

--- integers (scalar) ---

     nlon = number of longitude points (first dimension)
            excludes halo points 
     nlat = number of latitude points (second dimension)
            excludes halo points
     nlev = number of vertical levels

     ntrace = number of tracers

     ilb = lower bound  (1st dimension) includes halo points
     iub = upper bound  (1st dimension) includes halo points
     jlb = lower bound  (2nd dimension) includes halo points
     jub = upper bound  (2nd dimension) includes halo points
     klb = lower bound  (3rd dimension)
     kub = upper bound  (3rd dimension)

--- prognostic fields (real arrays) ---

     ps   = surface pressure
              [real, dimension (ilb:iub, jlb:jub) ]
     pssl = surface pressure at eta=1 (sea level)
              [real, dimension (ilb:iub, jlb:jub) ]

     u    = zonal wind component
              [real, dimension (ilb:iub, jlb:jub, klb:kub) ]
     v    = meridional wind component
              [real, dimension (ilb:iub, jlb:jub, klb:kub) ]
     t    = temperature
              [real, dimension (ilb:iub, jlb:jub, klb:kub) ]
     r    = arbitrary number of tracers (includes specific humidity)
              [real, dimension (ilb:iub, jlb:jub, klb:kub, 1:ntrace) ]

</pre></a><!-- END DATA_TYPES -->
<!--------------------------------------------------------------------->
 <a name="ROUTINES" id="ROUTINES">
<hr>
<h2>PUBLIC ROUTINES</h2>
<!-- BEGIN ROUTINES -->
<pre>
<a name="prog_var_init" id="prog_var_init">

call prog_var_init ( Hgrid, nlev, ntrs, Vars )

  INPUT

     Hgrid   Derived-type variable containing horizontal grid constants.
               [type(horiz_grid_type), see horiz_grid_mod]

     nlev    The number of full model levels for the prognostic variables.
               [integer]

     ntrs    The total number of tracers.
               [integer]

  INPUT/OUTPUT

     Vars    Derived-type variable containing the model's prognostic
             fields (see above).
               [type(prog_var_type)]

---------------------------------------------------------------------
</a><a name="var_init" id="var_init">

The interface <b>var_init</b> can take several forms. var =&gt; <b>var_init</b> ( Hgrid ) var =&gt; <b>var_init</b> ( Hgrid, kdim ) var =&gt; <b>var_init</b> ( Hgrid, kdim, ntrace ) var =&gt; <b>var_init</b> ( ilb, iub, jlb, jub ) var =&gt; <b>var_init</b> ( ilb, iub, jlb, jub, kdim ) var =&gt; <b>var_init</b> ( ilb, iub, jlb, jub, kdim, ntrace ) INPUT Hgrid Derived-type variable containing horizontal grid constants. [type(horiz_grid_type), see horiz_grid_mod] ilb Lower bound/index for first dimension. [integer] iub Upper bound/index for first dimension. [integer] jlb Lower bound/index for second dimension. [integer] jub Upper bound/index for second dimension. [integer] kdim The size of the third dimension (or level dimension). [integer] ntrace The size of the fourth dimension (or tracer dimension). [integer] RETURNS The returned value is a pointer to memory. <b>Fields that are initialized this way must be declared at pointers.</b> Use the following syntax: real, pointer :: var(:,:) or var(:,:,:) or var(:,:,:,:) --------------------------------------------------------------------- </a><a name="prog_var_times_scalar"
id="prog_var_times_scalar">

<b>call prog_var_times_scalar</b> ( Var, scalar ) INPUT/OUTPUT Var prog_var_type which on output will have the prognostic variable components (u,v,t,r,ps,pssl) multiplied by scalar INPUT scalar a real scalar quantity --------------------------------------------------------------------- </a><a name="prog_var_equals_scalar"
id="prog_var_equals_scalar">

<b>call prog_var_equals_scalar</b> ( Var, scalar ) INPUT/OUTPUT Var prog_var_type which on output will have the prognostic variable components (u,v,t,r,ps,pssl) multiplied by scalar INPUT scalar a real scalar quantity --------------------------------------------------------------------- </a><a name="prog_var_time_diff"
id="prog_var_time_diff">

<b>call prog_var_time_diff</b> ( dt, Var_dt, Var, nt ) INPUT dt time step [real] INPUT/OUTPUT Var_dt input value is the tendency for prognostic variables, the output value is zero [prog_var_type] Var the prognostic variables, the input values are at time level n, and the output values are at time level n+1 [prog_var_type] OPTIONAL INPUT nt number of tracers to be advanced in time, by default all tracers will be advanced from time n to n+1 --------------------------------------------------------------------- </a><a name="open_prog_var_file"
id="open_prog_var_file">

<b>call open_prog_var_file</b> ( ix, jx, kx ) OUTPUT ix, jx, kx The 3-dimensional size of a prognostic field. [integer] ------------------------------------------------------------- </a><a name="read_prog_var"
id="read_prog_var">

<b>call read_prog_var</b> ( Hgrid, Var, eta, peta, fis, res ) INPUT/OUTPUT Hgrid Var OUTPUT eta peta fis res ------------------------------------------------------------- </a><a name="write_prog_var"
id="write_prog_var">

<b>call write_prog_var</b> ( Var, Hgrid, Vgrid, fis, res ) INPUT Var Hgrid Vgrid fis res </a></pre></a><!-- END ROUTINES -->
<!--------------------------------------------------------------------->
 <a name="ERRORS" id="ERRORS">
<hr>
<h2>ERROR MESSAGES</h2>
<!-- BEGIN ERRORS -->
<pre>

     None.

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

     None.

</pre></a><!-- END NOTES -->
<!--------------------------------------------------------------------->
 <a name="PLANS" id="PLANS">
<hr>
<h2>FUTURE PLANS</h2>
<!-- BEGIN PLANS -->
<pre>

     Need routines to release allocated memory.
     These may be called prog_var_end and var_end.

</pre></a><!-- END PLANS -->
<!--------------------------------------------------------------------->
<hr>
</body>
</html>
