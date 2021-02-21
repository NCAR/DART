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
<title>program atmos_model</title>
</head>
<body bgcolor="#AABBCC" text="#332211">
<div align="center"><font size="-2"><a href=
"#PUBLIC%20INTERFACE">PUBLIC INTERFACE</a> / <a href=
"#NAMELIST">NAMELIST</a> / <a href="#ERROR%20MESSAGES">ERRORS</a> /
<a href="#NOTES">NOTES</a></font><br>
<br></div>
<hr>
<h1>program atmos_model</h1>
<a name="OVERVIEW" id="OVERVIEW">
<hr>
<h2>OVERVIEW</h2>
<!-- BEGIN OVERVIEW -->
<pre>

     A main program for running a stand-alone atmospheric model.

</pre></a><!-- END OVERVIEW -->
<!--------------------------------------------------------------------->
 <a name="DESCRIPTION" id="DESCRIPTION"><!-- BEGIN DESCRIPTION -->
<pre>
     This version is suitable for running a dry dynamical core (Held-Suarez
     benchmark, aka simple_physics) and shallow water models.

</pre></a><!-- END DESCRIPTION -->
<!--------------------------------------------------------------------->
 <a name="OTHER MODULES USED">
<hr>
<h2>OTHER MODULES USED</h2>
<!-- BEGIN OTHER MODULES USED -->
<pre>

       atmosphere_mod
     time_manager_mod
     diag_manager_mod
              fms_mod
    field_manager_mod
   tracer_manager_mod

</pre></a><!-- END OTHER MODULES USED -->
<!--------------------------------------------------------------------->
 <a name="PUBLIC INTERFACE">
<hr>
<h2>PUBLIC INTERFACE</h2>
<!-- BEGIN PUBLIC INTERFACE -->
<pre>

   This is a main program. There are no callable interfaces.

   A namelist interface called <b>&amp;main_nml</b> must reside
   in file <b>input.nml</b>. See the details below.

</pre></a><!-- END PUBLIC INTERFACE -->
<!--------------------------------------------------------------------->
 <a name="NAMELIST" id="NAMELIST">
<hr>
<h2>NAMELIST</h2>
<!-- BEGIN NAMELIST -->
<pre>

<b>&amp;main_nml</b>

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

</pre></a><!-- END NAMELIST -->
<!--------------------------------------------------------------------->
 <a name="ERROR MESSAGES">
<hr>
<h2>ERROR MESSAGES</h2>
<!-- BEGIN ERROR MESSAGES -->
<pre>

<b>FATAL ERRORS in program atmos_model</b>

    <b>dt_atmos has not been specified</b>
        A value must be specified for variable "dt_atmos" in
        namelist &amp;main_nml. See the namelist documentation for details.

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

</pre></a><!-- END ERROR MESSAGES -->
<!--------------------------------------------------------------------->
 <a name="KNOWN BUGS">
<hr>
<h2>KNOWN BUGS</h2>
<!-- BEGIN KNOWN BUGS -->
<pre>

     None.

</pre></a><!-- END KNOWN BUGS -->
<!--------------------------------------------------------------------->
 <a name="NOTES" id="NOTES">
<hr>
<h2>NOTES</h2>
<!-- BEGIN NOTES -->
<pre>

     None.

</pre></a><!-- END NOTES -->
<!--------------------------------------------------------------------->
 <a name="FUTURE PLANS">
<hr>
<h2>FUTURE PLANS</h2>
<!-- BEGIN FUTURE PLANS -->
<pre>

     None.

</pre></a><!-- END FUTURE PLANS -->
<!--------------------------------------------------------------------->
<hr>
</body>
</html>
