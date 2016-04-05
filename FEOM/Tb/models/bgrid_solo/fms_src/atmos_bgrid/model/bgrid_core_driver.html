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
<TITLE>module bgrid_core_driver_mod</TITLE>
<BODY BGCOLOR="#AABBCC" TEXT="#332211" >

<DIV ALIGN="CENTER"> <FONT SIZE=-2>
<A HREF="#INTERFACE">PUBLIC INTERFACE</A> / 
<A HREF="#DATA_TYPES">DATA</A> / 
<A HREF="#ROUTINES">ROUTINES</A> / 
<A HREF="#NAMELIST">NAMELIST</A> / 
<A HREF="#ERRORS">ERRORS</A> / 
<A HREF="#NOTES">NOTES</A> 
</FONT>
<BR><BR></DIV><HR>

<H2>module bgrid_core_driver_mod</H2>
<A NAME="HEADER">
<PRE>
     <B>Contact:</B>   B. Wyman
     <B>Reviewers:</B>
     <B>Change history:</B> <A HREF="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS
/atmos/bgrid/model/bgrid_core_driver.f90">WebCVS Log for bgrid_core_driver.f90</A>

</PRE>
</A><!-- END HEADER -->
<!--------------------------------------------------------------------->
<A NAME="OVERVIEW">
<HR>
<H4>OVERVIEW</H4>
<!-- BEGIN OVERVIEW -->
<PRE>

     Provides high-level interfaces to the B-grid dynamical core
     that allow easy initialization, integration, and diagnostics.

</PRE>
</A><!-- END OVERVIEW -->
<!--------------------------------------------------------------------->
<A NAME="DESCRIPTION">
<!-- BEGIN DESCRIPTION -->
<PRE>
     There is a namelist interface for initializing the optional
     arguments to subroutine bgrid_core_init.

</PRE>
</A><!-- END DESCRIPTION -->
<!--------------------------------------------------------------------->
<A NAME="MODULES_USED">
<HR>
<H4>OTHER MODULES USED</H4>
<!-- BEGIN MODULES_USED -->
<PRE>

   bgrid_core_mod
   bgrid_horiz_mod
   bgrid_vert_mod
   bgrid_prog_var_mod
   bgrid_halo_mod
   bgrid_diagnostics_mod
   bgrid_integrals_mod
   bgrid_conserve_energy_mod
   time_manager_mod
   fms_mod

</PRE>
</A><!-- END MODULES_USED -->
<!--------------------------------------------------------------------->
<A NAME="INTERFACE">
<HR>
<H4>PUBLIC INTERFACE</H4>
<!-- BEGIN INTERFACE -->
<PRE>

   <b>use bgrid_core_driver_mod</b> [ ,only: bgrid_dynam_type, 
                                      bgrid_core_driver_init,
                                      bgrid_core_driver,
                                      bgrid_core_driver_end,
                                      bgrid_core_time_diff,
                                      get_bottom_data,
                                      put_bottom_data   ]

   <a href="#DATA_TYPES">bgrid_dynam_type</a>
        A defined-type variable that contains constants needed
        by the B-grid dynamical core (see bgrid_core_mod)

   <a href="#bgrid_core_driver_init">bgrid_core_driver_init</a>
        Initializes the B-grid dynamical core.
        This interface returns values for the "bgrid_dynam_type" and
        prog_var_type" derived-type variables. Also internally
        initialized are other B-grid derived-type variables for
        the horizontal and vertical grid constants.

   <a href="#bgrid_core_driver">bgrid_core_driver</a>
        A wrapper for integrating the dynamical core one (atmospheric) time step.
        Note that only the tendencies of prognostic variables are updated.

   <a href="#bgrid_core_time_diff">bgrid_core_time_diff</a>
        Performs the time differencing of the prognostics variables 
        and outputs diagnostics for the B-grid dynamical core.

   <a href="#bgrid_core_driver_end">bgrid_core_driver_end</a>
        A wrapper for terminating the dynamical core.

   <a href="#get_bottom_data">get_bottom_data</a>, <a href="#put_bottom_data">put_bottom_data</a>
        Routines for getting and putting data at the model level
        closest to the ground.


   NOTES

      * A namelist interface called <a href="#NAMELIST">bgrid_core_driver_nml</a> is read
        from file <b>input.nml</b>.

</PRE>
</A><!-- END INTERFACE -->
<!--------------------------------------------------------------------->
<A NAME="DATA_TYPES">
<HR>
<H4>PUBLIC DATA</H4>
<!-- BEGIN DATA_TYPES -->
<PRE>

<b>type (bgrid_dynam_type)</b>

     See <A HREF="bgrid_core.html#DATA_TYPES">bgrid_core_mod</A> for details.

</PRE>
</A><!-- END DATA_TYPES -->
<!--------------------------------------------------------------------->
<A NAME="ROUTINES">
<HR>
<H4>PUBLIC ROUTINES</H4>
<!-- BEGIN ROUTINES -->
<PRE>
<a name="bgrid_core_driver_init">

<b>call bgrid_core_driver_init</b> ( Time_init, Time, Time_step, 
                               Var, Var_dt, Dynam, phys_axes )

DESCRIPTION
   Returns initialized/allocated values for the "bgrid_dynam_type"
   and "prog_var_type" derived-type variables. Also internally
   initialized are other B-grid derived-type variables for the
   horizontal and vertical grid constants.

INPUT
   Time_init  The initial (or base) time.  <FONT SIZE=-1 COLOR="#000099">[time_type]</FONT>

   Time       The current time.  <FONT SIZE=-1 COLOR="#000099">[time_type]</FONT>

   Time_step  The atmospheric model/physics time step.  <FONT SIZE=-1 COLOR="#000099">[time_type]</FONT>

INPUT/OUTPUT
   Var        A derived-type variable that contains the prognostic
              variables for the B-grid dynamical core.
              The returned values will have been initialized
              by prog_var_mod (most likely read from a restart file).
                 <FONT SIZE=-1 COLOR="#000099">[type(prog_var_type)]</FONT>

   Var_dt     A derived-type variable that contains the prognostic
              variable time tendencies. The returned value is zero.
                 <FONT SIZE=-1 COLOR="#000099">[type(prog_var_type)]</FONT>

   Dynam      A derived-type variable that contains almost everything
              needed by the dynamical core.
                 <FONT SIZE=-1 COLOR="#000099">[type(bgrid_dynam_type)]</FONT>

OUTPUT
   phys_axes  Axis identifiers as returned by the diagnostics manager
              and needed for subsequent calls to the diagnostics manager.
                 <FONT SIZE=-1 COLOR="#000099">[integer, dimension(4)]</FONT>
              
<!------------------------>
<HR WIDTH="50%" ALIGN=LEFT>
<a name="bgrid_core_driver">

<b>call bgrid_core_driver</b> ( Time_diag, Var, Var_dt, Dynam, omega )

DESCRIPTION
   Updates the prognostic variable tendencies with the dynamical
   core tendencies for the current atmospheric time step.
   Also calls diagnostics routines for outputting the dynamical
   core tendencies.

INPUT
   Time_diag  The diagnostics time, usually the current time + time step.
                 <FONT SIZE=-1 COLOR="#000099">[type(time_type)]</FONT>

   Var        A derived-type variable that contains the B-grid's
              prognostic variables.
                 <FONT SIZE=-1 COLOR="#000099">[type(prog_var_type)]</FONT>

INPUT/OUTPUT
   Var_dt     A derived-type variable that contains the TENDENCIES
              for the B-grid's prognostic variables.
                 <FONT SIZE=-1 COLOR="#000099">[type(prog_var_type)]</FONT>

   Dynam      The derived-type variable returned by a previous call
              to bgrid_core_driver_init (see above).
                 <FONT SIZE=-1 COLOR="#000099">[type(bgrid_dynam_type)]</FONT>

OUTPUT
   omega      The omega diagnostic (from the thermodynamic equation) in
              pascals per second. The array should have horizontal dimensions that
              are consistent with the data domain of the B-grid dynamical core.
                  <FONT SIZE=-1 COLOR="#000099">[real, dimension(ilb:,jlb:,:)]</FONT>

<!------------------------>
<HR WIDTH="50%" ALIGN=LEFT>
<a name="bgrid_core_time_diff">

<b>call bgrid_core_time_diff</b> ( omega, Time_diag, Dynam, Var, Var_dt )

DESCRIPTION
     Performs the time differencing of the prognostics variables 
     and outputs diagnostics for the B-grid dynamical core.

INPUT
   omega      The pressure vertical velocity in Pascals/second.
              This is only needed for diagnostic purposes.
                 <FONT SIZE=-1 COLOR="#000099">[real, dimension(:,:,:)]</FONT>

   Time_diag  The diagnostics time, usually the current time + time step.
                 <FONT SIZE=-1 COLOR="#000099">[type(time_type)]</FONT>

   Dynam      The derived-type variable returned by a previous call
              to bgrid_core_driver_init (see above).
                 <FONT SIZE=-1 COLOR="#000099">[type(bgrid_dynam_type)]</FONT>

INPUT/OUTPUT
   Var        The prognostic variables. The input quantities are at the
              current and on output they are at the next time step.
                 <FONT SIZE=-1 COLOR="#000099">[type(prog_var_type)]</FONT>

   Var_dt     The time tendencies for the prognostic variables.
              The output tendencies will have been set to zero.
                 <FONT SIZE=-1 COLOR="#000099">[type(prog_var_type)]</FONT>

<!------------------------>
<HR WIDTH="50%" ALIGN=LEFT>
<a name="bgrid_core_driver_end">

<b>call bgrid_core_driver_end</b> (Dynam)

DESCRIPTION
   Termination routine for the B-grid dynamical core.

INPUT
   Dynam   The derived-type variable returned by a previous call
           to bgrid_core_driver_init (see above).
              <FONT SIZE=-1 COLOR="#000099">[type(bgrid_dynam_type)]</FONT>

<!------------------------>
<HR WIDTH="50%" ALIGN=LEFT>
<a name="get_bottom_data">

<b>call get_bottom_data</b> ( a, b, a_bot, b_bot, <FONT COLOR="#007700">[, k_bot]</FONT> )

DESCRIPTION
   Given a pair of 3-dimensional model fields this interface returns
   the 2-dimensional fields at the model level closest to the ground.
   If optional argument "kbot" is NOT present the returned field
   will be the 2-d field at k = size(a,3).

INPUT
   a, b            Three-dimension fields on the model grid.
                   The last dimension varies from the top of the atmosphere
                   towards the surface.
                      <FONT SIZE=-1 COLOR="#000099">[real, dimension(:,:,:)]</FONT>

OUTPUT
   a_bot, b_bot    Data located at the model level closest to the ground.
                   Must have the same size as the first two dimensions of a and b.
                      <FONT SIZE=-1 COLOR="#000099">[real, dimension(:,:)]</FONT>

OPTIONAL INPUT
   k_bot           The vertical index for the model level closest to
                   the ground. Must have the same size as a_bot and b_bot.
                      <FONT SIZE=-1 COLOR="#000099">[integer, dimension(:,:)]</FONT>

<!------------------------>
<HR WIDTH="50%" ALIGN=LEFT>
<a name="put_bottom_data">

<b>call put_bottom_data</b> ( a_bot, b_bot, a, b <FONT COLOR="#007700">[, k_bot]</FONT> )

DESCRIPTION
   Puts 2-dimensional data given at the lowest model level
   into their 3-dimensional model fields.

INPUT
   a_bot, b_bot   Data located at the model level closest to the ground.
                  This data will be inserted into arrays a and b.
                     <FONT SIZE=-1 COLOR="#000099">[real, dimension(:,:)]</FONT>

INPUT/OUTPUT
   a, b           Three-dimension fields on the model grid.
                     <FONT SIZE=-1 COLOR="#000099">[real, dimension(:,:,:)]</FONT>

OPTIONAL INPUT
   k_bot          The vertical index for the model level closest to
                  the ground. Must have the same size as a_bot and b_bot.
                     <FONT SIZE=-1 COLOR="#000099">[integer, dimension(:,:)]</FONT>

</PRE>
</A><!-- END ROUTINES -->
<!--------------------------------------------------------------------->
<A NAME="NAMELIST">
<HR>
<H4>NAMELIST</H4>
<!-- BEGIN NAMELIST -->
<PRE>

 <b>&bgrid_core_driver_nml</b>

    damp_scheme          Determines how horizontal damping coefficients
                         vary with latitude.
                            = 1, constant
                            = 2, varies as inverse of diagonal grid distance
                            = 3, varies as inverse of x-axis grid distance
                         Note: damp_scheme = 1 is recommended, 
                         damp_scheme = 2,3 is experimental.
                             <FONT SIZE=-1 COLOR="#000099">[integer, default: damp_scheme = 1]</FONT>
   
    damp_order_wind      The horizontal damping order for momentum,
    damp_order_temp      temperature, and default order for all
    damp_order_tracer    prognostic tracers. The damping order must be
                         an even number; damp_order = 0 turns off damping.
                            <FONT SIZE=-1 COLOR="#000099">[integer, default: damp_order = 4]</FONT>

    damp_coeff_wind      The horizontal damping coefficients for
    damp_coeff_temp      momentum, temperature, and default value for
    damp_coeff_tracer    all prognostic tracers. The coefficients are
                         expressed as non-dimensional values for the
                         second-order diffusion operator (range = 0,1).
                            <FONT SIZE=-1 COLOR="#000099">[real, default: damp_coeff = 0.50]</FONT>
 
    slope_corr_wind      The topography slope correction applied to horizontal
    slope_corr_temp      damping of momentum and temperature (including all
                         prognostic tracers).  The coefficients (with range = 0,1)
                         are expressed as arrays of size 4.  The first 3 values are
                         coefficients for the lowest 3 model layers, the last value
                         represents the remaining uppermost layers.  A NON-ZERO
                         value turns the correction ON.  Typical values might be
                         (/ .25, .50, .75, .95 /).
                           <FONT SIZE=-1 COLOR="#000099">[real, dimension(4), default: slope_corr = 0.,0.,0.,0.]</FONT>

    advec_order_wind     The advection order for momentum, temperature,
    advec_order_temp     and default order for all prognostic tracers.
    advec_order_tracer   The advection order must be an even number.
                           <FONT SIZE=-1 COLOR="#000099">[integer, default: advec_order = 2]</FONT>

    advec_coeff_wind     Coefficients for modified Euler-backward advection
    advec_coeff_temp     scheme for momentum, temperature, and all
    advec_coeff_tracer   prognostic tracers.
                         NOTE: advec_coeff=0 is the Euler-forward scheme which
                         is unstable, advec_coeff=1 is the Euler-backward scheme
                         which is highly dissipative.
                           <FONT SIZE=-1 COLOR="#000099">[real, default: advec_coeff = 0.7]</FONT>

    num_fill_pass        The number of successive passes applied in the tracer
                         borrowing/filling scheme.  This conservative scheme is
                         used to fill negative tracer values. It is applied in
                         both the vertical and horizontal directions.
                         Each successive pass should remove more negative values,
                         however an optimum number of passes is probably between 1-3.
                         This is applied after advection to all prognostic tracers.
                           <FONT SIZE=-1 COLOR="#000099">[integer, default: num_fill_pass = 1]</FONT>

    grid_sep_coeff      Coefficient to suppress grid-separation problem 
                        associated with the B-grid. Currently, this option has been
                        disabled within the model, so that this coefficient does nothing.
                           <FONT SIZE=-1 COLOR="#000099">[real, default: grid_sep_coeff = 0.00]</FONT>

    filter_option       Determines how polar filtering is performed.
                        filter_option = 0,  NO filtering
                                      = 1,  not implemented
                                      = 2,  filter horiz OMG/DIV,
                                            advec mass tendencies,
                                            and momentum
                           <FONT SIZE=-1 COLOR="#000099">[integer, default: filter_option = 2]</FONT>

    filter_weight       Weight applied to the polar filter that will
                        increase (or decrease) the strength of the standard
                        polar filter response function.
                        SS(new) = SS(std)**filter_weight, 
                        where SS(std) is the Arakawa and Lamb response function.
                           <FONT SIZE=-1 COLOR="#000099">[integer, default: filter_weight = 1 ]</FONT>

    ref_lat_filter      The reference latitude at which polar filtering
                        (in each hemisphere) will begin to be applied.
                        Setting this argument &gt;= 90. will turn off
                        polar filtering.
                          <FONT SIZE=-1 COLOR="#000099">[real, default: ref_lat_filter = 60.]</FONT>

    num_sponge_levels   Number of uppermost model level where a band-pass
                        filter is applied to damp undesirable waves.
                        Currently num_sponge_levels &gt; 1 is not allowed.
                        If num_sponge_levels = 0, no damping is done.
                          <FONT SIZE=-1 COLOR="#000099">[integer, default: num_sponge_levels = 0 ]</FONT>

    sponge_coeff_wind     Damping coefficients for the sponge layer(s) in
    sponge_coeff_temp     the uppermost model levels. Coefficients have been 
    sponge_coeff_tracer   normalized and must be in the range [0,1].
                          If num_sponge_levels = 0, the value of the coefficients
                          is ignored.  There is no option to specify coefficients
                          that vary with level, although currently 
                          num_sponge_levels &gt; 1 is not allowed.
                             <FONT SIZE=-1 COLOR="#000099">[real, default: sponge_coeff = 0.]</FONT>

    halo                The number of halo rows along all (NEWS) boundaries.
                        There is currently no namelist option that allows unequal
                        halo boundary.  NOTE: Additional halo rows are not 
                        necessary when using higher order horizontal damping or
                        advection, and may in fact result in poorer cpu performance.
                           <FONT SIZE=-1 COLOR="#000099">[integer, default, halo = 1]</FONT>

    num_adjust_dt       The number of adjustment time steps for each advection
                        time step, where num_adjust_dt &gt;= 1. 
                          <FONT SIZE=-1 COLOR="#000099">[integer, default: num_adjust_dt = 3]</FONT>

    num_advec_dt        The number of advection/dynamics time steps for each
                        atmospheric/physics time step, where num_advec_dt &gt;= 1.
                          <FONT SIZE=-1 COLOR="#000099">[integer, default: num_advec_dt = 1]</FONT>

    decomp              The domain decomposition, where decomp(1) = x-axis
                        decomposition, decomp(2) = y-axis decomposition.
                        * If decomp(1)*decomp(2) does not equal the number
                          of processors the model will fail.
                        * If decomp(1)=decomp(2)=0 then default rules apply.
                        * By default, one-dimensional decomposition (in Y) is used.
                          When there is fewer than 2 points per processor, then 2-D
                          decomposition is used.
                            <FONT SIZE=-1 COLOR="#000099">[integer, dimension(2), default: decomp = 0,0]</FONT>

    do_conserve_energy  If TRUE the temperature tendency will be updated to
                        guarantee that the dynamical core conserves total energy.
                        The correction is applied to a uniform global value.
                          <FONT SIZE=-1 COLOR="#000099">[logical, default: do_conserve_energy=.false.]</FONT>

    verbose             Flag that control additional printed output.
                        Currently, this option is not being used.
                          <FONT SIZE=-1 COLOR="#000099">[integer, default: verbose = 0]</FONT>

 NOTES


</PRE>
</A><!-- END NAMELIST -->
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

     None.

</PRE>
</A><!-- END PLANS -->
<!--------------------------------------------------------------------->

<HR>
</BODY>
</HTML>
