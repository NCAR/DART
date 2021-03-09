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
<TITLE>module bgrid_integrals</TITLE>
<BODY BGCOLOR="#AABBCC" TEXT="#332211" >

<DIV ALIGN="CENTER"> <FONT SIZE=-2>
<A HREF="#INTERFACE">PUBLIC INTERFACE</A> / 
<A HREF="#ROUTINES">ROUTINES</A> / 
<A HREF="#NAMELIST">NAMELIST</A> / 
<A HREF="#ERRORS">ERRORS</A> / 
<A HREF="#NOTES">NOTES</A> 
</FONT>
<BR><BR></DIV><HR>


<H2>module bgrid_integrals</H2>
<A NAME="HEADER">
<PRE>
     <B>Contact:</B>   Bruce Wyman
     <B>Reviewers:</B>
     <B>Change history:</B> <A HREF="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS
/atmos/bgrid/tools/bgrid_integrals.f90">WebCVS Log for bgrid_integrals.f90</A>

</PRE>
</A><!-- END HEADER -->
<!--------------------------------------------------------------------->
<A NAME="OVERVIEW">
<HR>
<H4>OVERVIEW</H4>
<!-- BEGIN OVERVIEW -->
<PRE>

     Computes and outputs global integrals of various quantities
     for the B-grid dynamical core.
</PRE>
</A><!-- END OVERVIEW -->
<!--------------------------------------------------------------------->
<A NAME="DESCRIPTION">
<!-- BEGIN DESCRIPTION -->
<PRE>

     Global integrals of the following quantities (except wind speed)
     are computed and output as ascii to either standard output or
     a user supplied file name:

        time                 (n)
        surface pressure     (ps)
        temperature          (tavg)
        minimum temperature  (tmin)
        maximum wind speed   (vmax)
        kinetic energy       (ke)
        total energy         (te)
        enstrophy            (ens)
        tracers 1-4 (trs --->)

     An internal variable allows for backward compatibility with
     the previous version. Zonal and eddy kinetic energy averages
     can replace tmin and te.
     
        zonal kinetic energy (zke)
        eddy kinetic energy  (eke)

     Notes:

        1) All quantities are instantaneous, EXCEPT tmin and vmax
           which are determined over the output interval.
        2) Global integrals are area and pressure weighted, then
           normalized by the global average pressure (see notes below).
        3) A header record is output before the data, the header record
           names are in parentheses.
        4) Zonal and eddy kinetic energy require computing a zonal mean value.
           The zonal mean will not be computed correctly when there is
           domain decomposition along the x-axis.

</PRE>
</A><!-- END DESCRIPTION -->
<!--------------------------------------------------------------------->
<A NAME="MODULES_USED">
<HR>
<H4>OTHER MODULES USED</H4>
<!-- BEGIN MODULES_USED -->
<PRE>

    bgrid_change_grid_mod
    bgrid_horiz_mod
    bgrid_vert_mod
    bgrid_masks_mod
    bgrid_prog_var_mod
    time_manager_mod
    fms_mod
    constants_mod
    mpp_mod
    mpp_io_mod
    mpp_domains_mod
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

   <b>use bgrid_integrals_mod</b> [, only: bgrid_integrals_init,
                                       bgrid_integrals,
                                       bgrid_integrals_end ]

   <a href="#bgrid_integrals_init">bgrid_integrals_init</a>
        The initialization routine that must be called once before
        calling bgrid_integrals.

   <a href="#bgrid_integrals">bgrid_integrals</a>
         Routine that computes and writes integrals to the desired output file.
         It is called every time step.  An internal alarm determines whether
         the integrals should be computed and written.

   <a href="#bgrid_integrals_end">bgrid_integrals_end</a>
         Called outside the timeloop at the end of an integration.
         Closes all open IO units.

   Notes:

     1) Namelist <a href="#NAMELIST">&bgrid_integrals_nml</a> controls the output frequency of 
        B-grid global integrals.

     2) When the model is stopped, there is no saving of restart data
        such as the alarm information. The alarm will be reset when the
        model is restarted.

</PRE>
</A><!-- END INTERFACE -->
<!--------------------------------------------------------------------->
<A NAME="ROUTINES">
<HR>
<H4>PUBLIC ROUTINES</H4>
<!-- BEGIN ROUTINES -->
<PRE>
<a name="bgrid_integrals_init">

<b>call bgrid_integrals_init</b> ( Time_init, Time )

 INPUT

   Time_init    Reference (base) time.     [time_type]

   Time         The current time, must always be greater or equal to Time_init.
                The base time (Time_init) will be subtracted from this time
                when labeling the integral output.   [time_type]

 NOTES

   Time_init is saved internally by this module as the base time.

------------------------------------------------------------------------
<a name="bgrid_integrals">

<b>call bgrid_integrals</b> ( Time, Hgrid, Vgrid, Var, Masks )

 INPUT

   Time         The current time, must always be greater or equal to Time_init.
                   [time_type]

   Hgrid        derived-type variable containing horizontal grid constants
                for the B-grid dynamical core   [horiz_grid_type]

   Vgrid        derived-type variable containing vertical grid constants
                for the B-grid dynamical core     [vert_grid_type]

   Var          derived-type variable containing the prognostic variables
                for the B-grid dynamical core        [prog_var_type]

   Masks        derived-type variable containing the grid box masks for
                the eta coordinate   [grid_mask_type]


 NOTES

   1) If it is not time to output integrals, bgrid_integrals will do nothing.
   2) Integrals are computed only when they will be written to standard 
      output or the requested file name.

------------------------------------------------------------------------
<a name="bgrid_integrals_end">

<b>call bgrid_integrals_end</b>

</PRE>
</A><!-- END ROUTINES -->
<!--------------------------------------------------------------------->
<A NAME="NAMELIST">
<HR>
<H4>NAMELIST</H4>
<!-- BEGIN NAMELIST -->
<PRE>

<b>&bgrid_integrals_nml</b>

 file_name  = optional file name for output (max length of 32 chars);
              if no name is specified (the default) then
              standard output will be used
                 [character, default: filename  = ' ']

 time_units = specifies the time units used for time,
              the following values are valid character strings
                 time_units = 'seconds'
                            = 'minutes'
                            = 'hours'   (default)
                            = 'days'

 output_interval = time interval in units of "time_units" for
                   global b-grid integral diagnostics;
                   * if an interval of zero is specified then no
                     diagnostics will be generated
                   * a negative value tries to use a value from a
                     restart file
                       [real, default: output_interval = -1.0]

 chksum_file_name = File name for output integrals in hexadecimal format.
                    If chksum_file_name is specified, integrals will 
                    be computed from the global data field for exact
                    reproducibility (both in the standard output file and
                    in chksum_file_name). Reproducibility is only an issue
                    when running on multiple processors.
                    If chksum_file_name is not specified, integrals may not
                    reproduce on multiple processors.
                    The standard output file is always produced regardless 
                    of chksum_file_name.

 trsout           = The names of the tracers written to the standard
                    output file. 
                      [character, dimension(4), 
                       default: 'sphum  ','liq_wat','ice_wat','cld_amt']

 tracer_file_name = The output format can accommodate 4 tracers. If output
                    for  more than 4 tracers is needed, "tracer_file_name"
                    specifies the file for additional tracer integrals.

</PRE>
</A><!-- END NAMELIST -->
<!--------------------------------------------------------------------->
<A NAME="ERRORS">
<HR>
<H4>ERROR MESSAGES</H4>
<!-- BEGIN ERRORS -->
<PRE>

<b>FATAL Errors in bgrid_integrals_mod</b>

    <b>must call bgrid_integrals_init</b>
        The initialization routine must be called before calling
        bgrid_integrals.

    <b>wsum=0</b>
        A global weighted average has a weight sum of zero.

<b>NOTES in bgrid_integrals_mod</b>

    <b>checksum integrals of zonal and eddy KE will not be exact
    with x-axis decomposition</b>
        Do not expect the output integrals called before calling
        bgrid_integrals.

    <b>end of the output period did not coincide with the end of the model run</b>
        Quantities that are buffered between output intervals (tmin,vmax)
        may not reproduce between runs.

</PRE>
</A><!-- END ERRORS -->
<!--------------------------------------------------------------------->
<A NAME="BUGS">
<HR>
<H4>KNOWN BUGS</H4>
<!-- BEGIN BUGS -->
<PRE>

     No known bugs.

</PRE>
</A><!-- END BUGS -->
<!--------------------------------------------------------------------->
<A NAME="NOTES">
<HR>
<H4>NOTES</H4>
<!-- BEGIN NOTES -->
<PRE>

 FORMULAS USED FOR COMPUTING INTEGRAL QUANTITIES

   1) pressure

                sum { dp dA }
        p_avg = ------------- 
                  sum {dA}

   2) temperature and tracers

                 sum { T dp dA }
        t_avg = -----------------
                sum {dA}  p_avg

   3)  kinetic energy

                 0.5 * sum { V**2 dp dA }
        ke_avg = ------------------------
                   sum {dA}  p_avg

        where V =  u  or  v   for total kinetic energy
                = [u] or [v]  for zonal mean kinetic energy
                =  u* or  v*  for eddy kinetic energy

                  [ ]  = zonal mean
                  ( )* = deviation from the zonal mean

   4) enstrophy

                   sum { R**2 dp dA }
        ens_avg = --------------------
                   sum {dA}  p_avg

        R = relative vorticity

</PRE>
</A><!-- END NOTES -->
<!--------------------------------------------------------------------->
<A NAME="PLANS">
<HR>
<H4>FUTURE PLANS</H4>
<!-- BEGIN PLANS -->
<PRE>

    1) Netcdf capabilities?
    2) Switch to a very general diagnostics handler?

</PRE>
</A><!-- END PLANS -->
<!--------------------------------------------------------------------->

<HR>
</BODY>
</HTML>
