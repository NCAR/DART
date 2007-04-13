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
<TITLE>module bgrid_prog_var_mod</TITLE>
<BODY BGCOLOR="#AABBCC" TEXT="#332211" >

<DIV ALIGN="CENTER"> <FONT SIZE=-2>
<A HREF="#INTERFACE">PUBLIC INTERFACE</A> / 
<A HREF="#DATA_TYPES">DATA</A> / 
<A HREF="#ROUTINES">ROUTINES</A> / 
<A HREF="#ERRORS">ERRORS</A> / 
<A HREF="#NOTES">NOTES</A> 
</FONT>
<BR><BR></DIV><HR>


<H2>module bgrid_prog_var_mod</H2>
<A NAME="HEADER">
<PRE>
     <B>Contact:</B>   B. Wyman
     <B>Reviewers:</B>
     <B>Change history:</B> <A HREF="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS
/atmos/bgrid/tools/bgrid_prog_var.f90">WebCVS Log for bgrid_prog_var.f90</A>

</PRE>
</A><!-- END HEADER -->
<!--------------------------------------------------------------------->
<A NAME="OVERVIEW">
<HR>
<H4>OVERVIEW</H4>
<!-- BEGIN OVERVIEW -->
<PRE>

     Initializes a derived-type variable that contains prognostic
     variables for the B-grid dynamical core.

</PRE>
</A><!-- END OVERVIEW -->
<!--------------------------------------------------------------------->
<A NAME="DESCRIPTION">
<!-- BEGIN DESCRIPTION -->
<PRE>
     This modules defines and initializes a derived-type variable
     (prog_var_type) that contains the following prognostic fields:
       surface pressure, temperature, momentum, and tracers.

     The main initialization interface allocates storage for
     a prog_var_type variable. Additional interfaces allocate
     storage for separate fields. All array storage includes the
     halo region when the allocation is done.

     There are also operators and interfaces for performing some
     simple operations between prog_var_type variables.
      
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
   bgrid_cold_start_mod
   fms_mod
   field_manager_mod
   tracer_manager_mod
 
</PRE>
</A><!-- END MODULES_USED -->
<!--------------------------------------------------------------------->
<A NAME="INTERFACE">
<HR>
<H4>PUBLIC INTERFACE</H4>
<!-- BEGIN INTERFACE -->
<PRE>

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
 
</PRE>
</A><!-- END INTERFACE -->
<!--------------------------------------------------------------------->
<A NAME="DATA_TYPES">
<HR>
<H4>PUBLIC DATA</H4>
<!-- BEGIN DATA_TYPES -->
<PRE>

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

</PRE>
</A><!-- END DATA_TYPES -->
<!--------------------------------------------------------------------->
<A NAME="ROUTINES">
<HR>
<H4>PUBLIC ROUTINES</H4>
<!-- BEGIN ROUTINES -->
<PRE>
<a name="prog_var_init">

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
<a name="var_init">

The interface <b>var_init</b> can take several forms.

var =&gt; <b>var_init</b> ( Hgrid )
var =&gt; <b>var_init</b> ( Hgrid, kdim )
var =&gt; <b>var_init</b> ( Hgrid, kdim, ntrace )

var =&gt; <b>var_init</b> ( ilb, iub, jlb, jub )
var =&gt; <b>var_init</b> ( ilb, iub, jlb, jub, kdim )
var =&gt; <b>var_init</b> ( ilb, iub, jlb, jub, kdim, ntrace )

  INPUT

     Hgrid   Derived-type variable containing horizontal grid constants.
               [type(horiz_grid_type), see horiz_grid_mod]

     ilb     Lower bound/index for first dimension.
               [integer]

     iub     Upper bound/index for first dimension.
               [integer]

     jlb     Lower bound/index for second dimension.
               [integer]

     jub     Upper bound/index for second dimension.
               [integer]

     kdim    The size of the third dimension (or level dimension).
               [integer]

     ntrace  The size of the fourth dimension (or tracer dimension).
               [integer]

  RETURNS

     The returned value is a pointer to memory.
     <b>Fields that are initialized this way must be declared at pointers.</b>

     Use the following syntax:

           real, pointer :: var(:,:) or var(:,:,:) or var(:,:,:,:)

---------------------------------------------------------------------
<a name="prog_var_times_scalar">

<b>call prog_var_times_scalar</b> ( Var, scalar )

  INPUT/OUTPUT

     Var      prog_var_type which on output will have the
              prognostic variable components (u,v,t,r,ps,pssl)
              multiplied by scalar

  INPUT

     scalar   a real scalar quantity

---------------------------------------------------------------------
<a name="prog_var_equals_scalar">

<b>call prog_var_equals_scalar</b> ( Var, scalar )

  INPUT/OUTPUT

     Var      prog_var_type which on output will have the
              prognostic variable components (u,v,t,r,ps,pssl)
              multiplied by scalar

  INPUT

     scalar   a real scalar quantity

---------------------------------------------------------------------
<a name="prog_var_time_diff">

<b>call prog_var_time_diff</b> ( dt, Var_dt, Var, nt )

  INPUT

     dt      time step [real]

  INPUT/OUTPUT

     Var_dt  input value is the tendency for prognostic variables,
             the output value is zero [prog_var_type]

     Var     the prognostic variables, the input values are at time
             level n, and the output values are at time level n+1
                [prog_var_type]

  OPTIONAL INPUT

     nt      number of tracers to be advanced in time, by default
             all tracers will be advanced from time n to n+1

---------------------------------------------------------------------
<a name="open_prog_var_file">

<b>call open_prog_var_file</b> ( ix, jx, kx )

OUTPUT

  ix, jx, kx   The 3-dimensional size of a prognostic field.
                 [integer]

-------------------------------------------------------------
<a name="read_prog_var">

<b>call read_prog_var</b> ( Hgrid, Var, eta, peta, fis, res )

INPUT/OUTPUT
   Hgrid
   Var

OUTPUT
   eta
   peta
   fis
   res

-------------------------------------------------------------
<a name="write_prog_var">

<b>call write_prog_var</b> ( Var, Hgrid, Vgrid, fis, res )

INPUT
   Var
   Hgrid
   Vgrid
   fis
   res

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

     None.

</PRE>
</A><!-- END NOTES -->
<!--------------------------------------------------------------------->
<A NAME="PLANS">
<HR>
<H4>FUTURE PLANS</H4>
<!-- BEGIN PLANS -->
<PRE>

     Need routines to release allocated memory.
     These may be called prog_var_end and var_end.

</PRE>
</A><!-- END PLANS -->
<!--------------------------------------------------------------------->

<HR>
</BODY>
</HTML>
