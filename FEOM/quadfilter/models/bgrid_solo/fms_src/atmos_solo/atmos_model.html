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
<TITLE>program atmos_model</TITLE>
<BODY BGCOLOR="#AABBCC" TEXT="#332211" >

<DIV ALIGN="CENTER"> <FONT SIZE=-2>
<A HREF="#PUBLIC INTERFACE">PUBLIC INTERFACE</A> / 
<A HREF="#NAMELIST">NAMELIST</A> / 
<A HREF="#ERROR MESSAGES">ERRORS</A> / 
<A HREF="#NOTES">NOTES</A> 
</FONT>
<BR><BR></DIV><HR>

<H2>program atmos_model</H2>
<A NAME="HEADER">
<PRE>
     <B>Contact:</B>   B. Wyman
     <B>Reviewers:</B>
     <B>Change history:</B> <A HREF="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS
/atmos/driver/solo/atmos_model.f90">WebCVS Log for solo/atmos_model.f90</A>

</PRE>
</A><!-- END HEADER -->
<!--------------------------------------------------------------------->
<A NAME="OVERVIEW">
<HR>
<H4>OVERVIEW</H4>
<!-- BEGIN OVERVIEW -->
<PRE>

     A main program for running a stand-alone atmospheric model.

</PRE>
</A><!-- END OVERVIEW -->
<!--------------------------------------------------------------------->
<A NAME="DESCRIPTION">
<!-- BEGIN DESCRIPTION -->
<PRE>
     This version is suitable for running a dry dynamical core (Held-Suarez
     benchmark, aka simple_physics) and shallow water models.

</A><!-- END DESCRIPTION -->
<!--------------------------------------------------------------------->
<A NAME="OTHER MODULES USED">
<HR>
<H4>OTHER MODULES USED</H4>
<!-- BEGIN OTHER MODULES USED -->
<PRE>

       atmosphere_mod
     time_manager_mod
     diag_manager_mod
              fms_mod
    field_manager_mod
   tracer_manager_mod

</PRE>
</A><!-- END OTHER MODULES USED -->
<!--------------------------------------------------------------------->
<A NAME="PUBLIC INTERFACE">
<HR>
<H4>PUBLIC INTERFACE</H4>
<!-- BEGIN PUBLIC INTERFACE -->
<PRE>

   This is a main program. There are no callable interfaces.

   A namelist interface called <b>&main_nml</b> must reside
   in file <b>input.nml</b>. See the details below.

</PRE>
</A><!-- END PUBLIC INTERFACE -->
<!--------------------------------------------------------------------->
<A NAME="NAMELIST">
<HR>
<H4>NAMELIST</H4>
<!-- BEGIN NAMELIST -->
<PRE>

<b>&main_nml</b>

  current_time  The time (day,hour,minute,second) that the current
                integration starts with.
                  [integer, dimension(4), default: current_time=0]

  override      Flag that determines whether the namelist variable
                current_time should override the time in the
                restart file INPUT/atmos_model.res. If the restart file
                does not exist then override has not effect, the
                value of current_date will be used.
                  [logical, default: override=false]

  days          The number of days that the current integration will
                be run for.   [integer, default: days=0]

  hours         The number of hours that the current integration will
                be run for.   [integer, default: hours=0]

  minutes       The number of minutes that the current integration will
                be run for.   [integer, default: minutes=0]

  seconds       The number of seconds that the current integration will
                be run for.   [integer, default: seconds=0]


  dt_atmos      Time step in seconds for the atmospheric model.
                Must be specified.
                   [integer, default: dt_atmos=0]

  Notes:

    1) If no value is set for current_time (or default value specified)
       then the value from restart file "INPUT/atmos_model.res" will
       be used. If neither a namelist value or restart file value exist
       the program will fail.

    2) The actual run length will be the sum of days, hours,
       minutes, and seconds.  A run length of zero is not a valid option.

</PRE>
</A><!-- END NAMELIST -->
<!--------------------------------------------------------------------->
<A NAME="ERROR MESSAGES">
<HR>
<H4>ERROR MESSAGES</H4>
<!-- BEGIN ERROR MESSAGES -->
<PRE>

<b>FATAL ERRORS in program atmos_model</b>

    <b>dt_atmos has not been specified</b>
        A value must be specified for variable "dt_atmos" in
        namelist &main_nml. See the namelist documentation for details.

    <b>invalid base date - must have year = month = 0</b>
        There is no calendar associated with this model.
        The base date retrieved from the diagnostics manager assumes
        that a year and month exist, this is not allowed.
        Set the base date year and month to zero in diag_table.

    <b>initial time is greater than current time</b>
        If a restart file is present, then the namelist value for either
        current_time or base time (from diag_table) was incorrectly set.

    <b>run length must be multiple of atmosphere time step</b>
        There must be an even number of atmospheric time steps for the
        requested run length.

<b>WARNINGS in program atmos_model</b>

    <b>final time does not match expected ending time</b>
        This error should probably not occur because of checks done at
        initialization time.

</PRE>
</A><!-- END ERROR MESSAGES -->
<!--------------------------------------------------------------------->
<A NAME="KNOWN BUGS">
<HR>
<H4>KNOWN BUGS</H4>
<!-- BEGIN KNOWN BUGS -->
<PRE>

     None.

</PRE>
</A><!-- END KNOWN BUGS -->
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
<A NAME="FUTURE PLANS">
<HR>
<H4>FUTURE PLANS</H4>
<!-- BEGIN FUTURE PLANS -->
<PRE>

     None.

</PRE>
</A><!-- END FUTURE PLANS -->
<!--------------------------------------------------------------------->

<HR>
</BODY>
</HTML>
