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
<TITLE>module bgrid_core_mod</TITLE>
<BODY BGCOLOR="#AABBCC" TEXT="#332211" >

<DIV ALIGN="CENTER"> <FONT SIZE=-2>
<A HREF="#INTERFACE">PUBLIC INTERFACE</A> / 
<A HREF="#DATA_TYPES">DATA</A> / 
<A HREF="#ROUTINES">ROUTINES</A> / 
<A HREF="#ERRORS">ERRORS</A> / 
<A HREF="#REFERENCES">REFERENCES</A> / 
<A HREF="#NOTES">NOTES</A> 
</FONT>
<BR><BR></DIV><HR>


<H2>module bgrid_core_mod</H2>
<A NAME="HEADER">
<PRE>
     <B>Contact:</B>   B. Wyman
     <B>Reviewers:</B>
     <B>Change history:</B> <A HREF="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS
/atmos/bgrid/model/bgrid_core.f90">WebCVS Log for bgrid_core.f90</A>

</PRE>
</A><!-- END HEADER -->
<!--------------------------------------------------------------------->
<A NAME="OVERVIEW">
<HR>
<H4>OVERVIEW</H4>
<!-- BEGIN OVERVIEW -->
<PRE>

     The basic driver module for the B-grid atmospheric dynamical core.
     This version can be run on distributed memory machines.

</PRE>
</A><!-- END OVERVIEW -->
<!--------------------------------------------------------------------->
<A NAME="DESCRIPTION">
<!-- BEGIN DESCRIPTION -->
<PRE>
     The B-grid dynamical core updates the tendencies of the 
     prognostic variables, the variables themselves are not modified.
     The prognostic variables are surface pressure, the zonal and 
     meridional momentum components, temperature, and optional
     tracers. Specific humidity is considered a tracer.

     This version of the B-grid dynamics use a two time level scheme.
     The gravity waves are integrated using a forward-backward scheme
     with a relatively short time step and the advection terms are
     integrated using a modified Euler backward scheme and longer
     time step. Horizontal diffusion is done on the advective time step.

     Derived-type variables for the horizontal and vertical constants
     must be initialized before initializing the B-grid dynamical core.
     The initialization routine of the B-grid dynamical core returns
     a derived-type variable that is needed when updating the tendencies
     of the prognostic dynamical variables. These prognostic variables
     are also contained as a derived-type variable.

     Here is a detailed write-up of the finite differencing used in the
     <A HREF="bgrid_core.tech.ps">B-grid atmospheric dynamical core</A>.

</PRE>
</A><!-- END DESCRIPTION -->
<!--------------------------------------------------------------------->
<A NAME="MODULES_USED">
<HR>
<H4>OTHER MODULES USED</H4>
<!-- BEGIN MODULES_USED -->
<PRE>

   bgrid_prog_var_mod
   bgrid_horiz_mod
   bgrid_vert_mod
   bgrid_masks_mod
   bgrid_advection_mod
   bgrid_horiz_diff_mod
   bgrid_horiz_adjust_mod
   bgrid_vert_adjust_mod
   bgrid_polar_filter_mod
   bgrid_boundary_mod
   bgrid_sponge_mod
   fms_mod
   constants_mod

</PRE>
</A><!-- END MODULES_USED -->
<!--------------------------------------------------------------------->
<A NAME="INTERFACE">
<HR>
<H4>PUBLIC INTERFACE</H4>
<!-- BEGIN INTERFACE -->
<PRE>

   <b>use bgrid_core_mod</b> [ ,only: bgrid_dynam_type, 
                               bgrid_core_init,
                               update_bgrid_core,
                               bgrid_core_end    ]

   <a href="#DATA_TYPES">bgrid_dynam_type</a>
        A derived-type variable that contains grid constants, time step information,
        and other pre-computed terms that are needed by the dynamical core.
        This includes other derived-type variables and pointers to the physical
        grids of the model, but does not include the prognostic variables.

   <a href="#bgrid_core_init">bgrid_core_init</a>
        Initializes the derived-type variable, bgrid_dynam_type, and initializes
        other modules used. There is no namelist for bgrid_core_mod, all runtime
        options are controlled with optional arguments to this subroutine.
        This subroutine must be called before update_bgrid_core.

   <a href="#update_bgrid_core">update_bgrid_core</a>
        Called at every (atmospheric) time step to update the tendencies of the
        prognostic variables.  No time differencing is done. 

   <a href="#bgrid_core_end">bgrid_core_end</a>
        Called at the end of the model run to terminate the module.
        Currently this call does nothing.

</PRE>
</A><!-- END INTERFACE -->
<!--------------------------------------------------------------------->
<A NAME="DATA_TYPES">
<HR>
<H4>PUBLIC DATA</H4>
<!-- BEGIN DATA_TYPES -->
<PRE>

<b>type bgrid_dynam_type</b>

-- derived data types contained in bgrid_dynam_type --
   (see the appropriate module for details)

    Hgrid = horizontal grid constants, initialized by horiz_grid_mod
              [horiz_grid_type]
    Vgrid = vertical grid constants, initialized by vert_grid_mod
              [vert_grid_type]
    Masks = eta coordinate topography masks and indices,
            initialized by grid_masks_mod  [grid_mask_type]
    Pfilt = polar filter constants, initialized by polar_filter_mod
              [pfilt_control_type]
    Advec = advection constants, initialized by advection_mod
              [advec_control_type]
    Hdiff = horizontal diffusion constants, initialized by horiz_diff_mod
              [hdiff_control_type]

-- real, dimension(:,:) --

    fis  = geopotential height of the surface (m2/s2)
    fisl = geopotential height at eta=1. (for eta coord = 0.0,
    res  = reciprocal of eta at the surface

-- scalars --

    nt_adj = no. of adjustment time steps per advection step [integer]
    nt_adv = no. of advection  time steps per atmospheric step [integer]
    dt_adj = adjustment time step in seconds [real]

    dry_model = if true then water vapor will not be considered in
                the equation of state [logical]
    verbose   = verbose flag [integer]

    wcorr     = coefficient for the grid separation smoothing operator [real]
    fopt      = filtering option [integer]


 NOTE

   Hgrid, Vgrid, fis, and res are pointers to the arguments passed
   via the initialization call (bgrid_core_init).


</PRE>
</A><!-- END DATA_TYPES -->
<!--------------------------------------------------------------------->
<A NAME="ROUTINES">
<HR>
<H4>PUBLIC ROUTINES</H4>
<!-- BEGIN ROUTINES -->
<PRE>
<a name="bgrid_core_init">

Dynam = <b>bgrid_core_init</b> ( Hgrid, Vgrid, fis, res, dt, ntadj, ntadv
                  <FONT COLOR="#007700">[, damp_order_vel,  damp_order_tmp,  damp_order_trs,
                     damp_coeff_vel,  damp_coeff_tmp,  damp_coeff_trs,
                     damp_scheme, damp_slope_coeff_vel, damp_slope_coeff_tmp,
                     num_horiz_fill, num_vert_fill,
                     advec_order_vel, advec_order_tmp, advec_order_trs,
                     advec_coeff_vel, advec_coeff_tmp, advec_coeff_trs,
                     grid_sep_coeff, filter_option, filter_weight,
                     ref_lat_filter, num_sponge_levels,
                     sponge_coeff_vel, sponge_coeff_tmp, sponge_coeff_trs,
                     dry_model, verbose ]</FONT>  )

INPUT

   Hgrid     Derived-type variable containing horizontal grid constants (see horiz_grid_mod).
                [type(horiz_grid_type)]

   Vgrid     Derived-type variable containing vertical grid constants (see vert_grid_mod).
                [type(vert_grid_type)]

   fis       Geopotential height of the surface (m2/s2).
             Should have horizontal indexing consistent with the B-grid core.
                [real, dimension(:,:)]

   res       Reciprocal of eta at the surface. Used at a switch for sigma vs. eta.
             For the sigma coordinate, all res = 0.
             Should have horizontal indexing consistent with the B-grid core.
                [real, dimension(:,:)]

   dtadv     Time step in seconds for each call to update_bgrid_core.
             This should be the atmospheric time step.
                [real]

   ntadj     The number of adjustment time steps for each advective time step.
                [integer]

   ntadv     The number of advection time steps for each call to update bgrid_core.
                [integer]


RETURNS

   Dynam      Derived-type variable that contains quantities needed
              by the dynamical core.   [type(bgrid_dynam_type)]


OPTIONAL INPUT  (Note: These argument can be modified through a
                        namelist interface in the bgrid_core_driver_mod)

   damp_order_vel       The horizontal damping order for momentum,
   damp_order_tmp       temperature, and default order for all 
   damp_order_trs       prognostic tracers. The damping order must be 
                        an even number; damp_order = 0 turns off damping.
                           <FONT SIZE=-1 COLOR="#000099">[integer, default: damp_order = 4]</FONT>

   damp_coeff_vel       The horizontal damping coefficients for
   damp_coeff_tmp       momentum, temperature, and default value for
   damp_coeff_trs       all prognostic tracers. The coefficients are
                        expressed as non-dimensional values for the 
                        second-order diffusion operator (range = 0,1). 
                           <FONT SIZE=-1 COLOR="#000099">[real, default: damp_coeff = 0.50]</FONT>

   damp_scheme          Determines how horizontal damping coefficients
                        vary with latitude.
                           = 1, constant
                           = 2, varies as inverse of diagonal grid distance
                           = 3, varies as inverse of x-axis grid distance
                        Note: damp_scheme = 1 is recommended, 
                        damp_scheme = 2,3 is experimental.
                            <FONT SIZE=-1 COLOR="#000099">[integer, default: damp_scheme = 1]</FONT>



   damp_slope_corr_vel  The topography slope correction applied to horizontal
   damp_slope_corr_tmp  damping of momentum and temperature (including all
                        prognostic tracers).  The coefficients (with range = 0,1)
                        are expressed as arrays of size 4.  The first 3 values are
                        coefficients for the lowest 3 model layers, the last value
                        represents the remaining uppermost layers.  A NON-ZERO
                        value turns the correction ON.  Typical values might be
                        (/ .25, .50, .75, .95 /).
                           <FONT SIZE=-1 COLOR="#000099">[real, dimension(4), default: damp_slope_corr = 0.,0.,0.,0.]</FONT>

   advec_order_vel      The advection order for momentum, temperature,
   advec_order_tmp      and default order for all prognostic tracers.
   advec_order_trs      The advection order must be an even number.
                           <FONT SIZE=-1 COLOR="#000099">[integer, default: advec_order = 2]</FONT>

   advec_coeff_vel      Coefficients for modified Euler-backward advection
   advec_coeff_tmp      scheme for momentum, temperature, and all
   advec_coeff_trs      prognostic tracers.
                        NOTE: advec_coeff=0 is the Euler-forward scheme which
                        is unstable, advec_coeff=1 is the Euler-backward scheme
                        which is highly dissipative.
                           <FONT SIZE=-1 COLOR="#000099">[real, default: advec_coeff = 0.7]</FONT>

   num_horiz_fill       The number of successive horizontal and vertical passes
   num_vert_fill        applied in the tracer borrowing/filling scheme.  This
                        conservative scheme is used to fill negative tracer values.
                        It is applied in either the vertical and horizontal directions.
                        Each successive pass should remove more negative values,
                        however an optimum number of passes is probably between 1-3.
                        This is applied after advection to all prognostic tracers.
                           <FONT SIZE=-1 COLOR="#000099">[integer, default: num_fill = 1]</FONT>

   grid_sep_coeff       Coefficient to suppress grid-separation problem
                        associated with the B-grid. Currently, this option has been
                        disabled within the model, so that this coefficient does nothing.
                           <FONT SIZE=-1 COLOR="#000099">[real, default: grid_sep_coeff = 0.00]</FONT>

   filter_option        Determines how polar filtering is performed.
                        filter_option = 0,  NO filtering
                                      = 1,  not implemented
                                      = 2,  filter horiz OMG/DIV,
                                            advec mass tendencies,
                                            and momentum
                           <FONT SIZE=-1 COLOR="#000099">[integer, default: filter_option = 2]</FONT>

   filter_weight        Weight applied to the polar filter that will
                        increase (or decrease) the strength of the standard
                        polar filter response function.
                        SS(new) = SS(std)**filter_weight,
                        where SS(std) is the Arakawa and Lamb response function.
                           <FONT SIZE=-1 COLOR="#000099">[integer, default: filter_weight = 1 ]</FONT>

   ref_lat_filter       The reference latitude at which polar filtering
                        (in each hemisphere) will begin to be applied.
                        Setting this argument &gt;= 90. will turn off
                        polar filtering.
                           <FONT SIZE=-1 COLOR="#000099">[real, default: ref_lat_filter = 60.]</FONT>

   num_sponge_levels    Number of uppermost model level where a band-pass
                        filter is applied to damp undesirable waves.
                        Currently num_sponge_levels &gt; 1 is not allowed.
                        If num_sponge_levels = 0, no damping is done.
                           <FONT SIZE=-1 COLOR="#000099">[integer, default: num_sponge_levels = 0 ]</FONT>

   sponge_coeff_vel     Damping coefficients for the sponge layer(s) in
   sponge_coeff_tmp     the uppermost model levels. Coefficients have been
   sponge_coeff_trs     normalized and must be in the range [0,1].
                        If num_sponge_levels = 0, the value of the coefficients
                        is ignored.  There is no option to specify coefficients
                        that vary with level, although currently
                        num_sponge_levels &gt; 1 is not allowed.
                           <FONT SIZE=-1 COLOR="#000099">[real, default: sponge_coeff = 0.]</FONT>

   dry_model            Flag that determines whether or not water vapor effects
                        are included in the hydrostatic equation.
                           <FONT SIZE=-1 COLOR="#000099">[logical, default: dry_model = true]</FONT>

   verbose              Flag that control additional printed output.
                        Currently, this option is not being used.
                          <FONT SIZE=-1 COLOR="#000099">[integer, default: verbose = 0]</FONT>


NOTES

    fis and res should be dimension by the size of the global grid,
    the number of longitude points by number of latitude points.
    <b>The declaration of fis and res must have the target attribute
    and the storage must be static.</b>


---------------------------------------------------------------
<a name="update_bgrid_core">

<b>call update_bgrid_core</b> ( Var, Var_dt, Dynam, omega )

  INPUT

     Var     A derived-type variable that contains the B-grid's
             prognostic variables.
                [prog_var_type, see prog_var_mod]

  INPUT/OUTPUT

     Var_dt  A derived-type variable that contains the TENDENCIES
             for the B-grid's prognostic variables.
                [prog_var_type, see prog_var_mod]

     Dynam   The derived-type variable returned by a previous call
             to bgrid_core_init (see above).
                [type(bgrid_dynam_type)]

  OUTPUT

     omega   The omega diagnostic (from the thermodynamic equation) in
             pascals per second. The array should have horizontal dimensions
             that are consistent with the B-grid dynamical core.
                [real, dimension(ilb:,jlb:,:)]

---------------------------------------------------------------
<a name="bgrid_core_end">

<b>call bgrid_core_end</b> ( Dynam )

  INPUT

     Dynam    Derived-type variable returned by a previous call to bgrid_core_init.
                [type(bgrid_dynam_type)]

</PRE>
</A><!-- END ROUTINES -->
<!-- ------------------------------------------------------------------>
<A NAME="ERRORS">
<HR>
<H4>ERROR MESSAGES</H4>
<!-- BEGIN ERRORS -->
<PRE>

<b>FATAL ERRORS in bgrid_core_init</b>

    <b>input argument ntadj must be &gt;= 1</b>
        If you were using the namelist interface &amp;bgrid_core_driver_nml
        then check the namelist variable corresponding to this variable.

    <b>input argument ntadv must be &gt;= 1</b>
        If you were using the namelist interface &amp;bgrid_core_driver_nml
        then check the namelist variable corresponding to this variable.

    <b>input argument dtadv must be &gt; 0.</b>
        The model time step usually is set at the highest program level.
        Check the namelist for the main program.

</PRE>
</A><!-- END ERRORS -->
<!-- ------------------------------------------------------------------>
<A NAME="REFERENCES">
<HR>
<H4>REFERENCES</H4>
<!-- BEGIN REFERENCES -->
<PRE>

     Here is a detailed write-up of the finite differencing used in the
     <A HREF="bgrid_core.tech.ps">B-grid atmospheric dynamical core</A>.

</PRE>
</A><!-- END REFERENCES -->
<!-- ------------------------------------------------------------------>
<A NAME="BUGS">
<HR>
<H4>KNOWN BUGS</H4>
<!-- BEGIN BUGS -->
<PRE>

  There are no known bugs.

</PRE>
</A><!-- END BUGS -->
<!-- ------------------------------------------------------------------>
<A NAME="NOTES">
<HR>
<H4>NOTES</H4>
<!-- BEGIN NOTES -->
<PRE>

 internal options:

    alpha_implicit   determines how the coriolis and press grad force
                     terms are solved
                         = 0.5  trapezoidal implicit
                         = 1.0        fully implicit
                     [real, default: alpha_implicit = 0.5]

</PRE>
</A><!-- END NOTES -->
<!-- ------------------------------------------------------------------>
<A NAME="PLANS">
<HR>
<H4>FUTURE PLANS</H4>
<!-- BEGIN PLANS -->
<PRE>

  1) Addition polar filter options 
     (e.g., filter prognostic variables or tendencies)

</PRE>
</A><!-- END PLANS -->
<!-- ------------------------------------------------------------------>

<HR>
</BODY>
</HTML>
