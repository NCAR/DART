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
<html>
<head>
<META http-equiv="Content-Type" content="text/html; charset=EUC-JP">
<title>Module time_manager_mod</title>
<link type="text/css" href="http://www.gfdl.noaa.gov/~fms/style/doc.css" rel="stylesheet">
<STYLE TYPE="text/css">
          .fixed {
            font-size:medium;
            font-family:monospace;
            border-style:none;
            border-width:0.1em;
            padding:0.1em;
            color:#663366;
          }
        </STYLE>
</head>
<body>
<a name="TOP"></a><font class="header" size="1"><a href="#PUBLIC INTERFACE">PUBLIC INTERFACE </a>~
          <a href="#PUBLIC DATA">PUBLIC DATA </a>~
          <a href="#PUBLIC ROUTINES">PUBLIC ROUTINES </a>~
          <a href="#NAMELIST">NAMELIST </a>~
          <a href="#DIAGNOSTIC FIELDS">DIAGNOSTIC FIELDS </a>~
          <a href="#ERROR MESSAGES">ERROR MESSAGES </a>~
          <a href="#REFERENCES">REFERENCES </a>~ 
          <a href="#NOTES">NOTES</a></font>
<hr>
<h2>module time_manager_mod</h2>
<a name="HEADER"></a>
<!-- BEGIN HEADER -->
<div>
<b>Contact:</b>&nbsp;<a href="mailto:fms@gfdl.noaa.gov">   fms </a>
<br>
<b>Reviewers:</b>&nbsp;<br>
<b>Change History: </b>&nbsp;<a href="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/">WebCVS Log</a>
<br>
<b>Last Modified:</b>&nbsp;2002/04/16 21:00:57<br>
<br>
</div>
<!-- END HEADER -->
<a name="OVERVIEW"></a>
<hr>
<h4>OVERVIEW</h4>
<!-- BEGIN OVERVIEW -->
<p class="text">   A software package that provides a set of simple interfaces for
   modelers to perform computations related to time and dates. </p>
<!-- END OVERVIEW -->
<a name="DESCRIPTION"></a>
<!-- BEGIN DESCRIPTION -->
<div>   The module defines a type that can be used to represent discrete
   times (accurate to one second) and to map these times into dates
   using a variety of calendars. A time is mapped to a date by
   representing the time with respect to an arbitrary base date (refer
   to &lt;B&gt;NOTES&lt;/B&gt; section for the <a href="#base date">base date</a>   setting).
   <br>
<br>
   The time_manager provides a single defined type, time_type, which is
   used to store time and date quantities. A time_type is a positive
   definite quantity that represents an interval of time. It can be
   most easily thought of as representing the number of seconds in some
   time interval. A time interval can be mapped to a date under a given
   calendar definition by using it to represent the time that has passed
   since some base date. A number of interfaces are provided to operate
   on time_type variables and their associated calendars. Time intervals
   can be as large as n days where n is the largest number represented by
   the default integer type on a compiler. This is typically considerably
   greater than 10 million years (assuming 32 bit integer representation)
   which is likely to be adequate for most applications. The description
   of the interfaces is separated into two sections. The first deals with
   operations on time intervals while the second deals with operations
   that convert time intervals to dates for a given calendar. </div>
<br>
<!-- END DESCRIPTION -->
<a name="OTHER MODULES USED"></a>
<hr>
<h4>OTHER MODULES USED</h4>
<!-- BEGIN OTHER MODULES USED -->
<div>
<pre>fms_mod</pre>
</div>
<!-- END OTHER MODULES USED -->
<!-- BEGIN PUBLIC INTERFACE -->
<a name="PUBLIC INTERFACE"></a>
<hr>
<h4>PUBLIC INTERFACE</h4>
<div>
<pre>use time_manager_mod [, only:  set_time,<br>                               get_time,<br>                               increment_time,<br>                               decrement_time,<br>                               time_gt,<br>                               time_ge,<br>                               time_lt,<br>                               time_le,<br>                               time_eq,<br>                               time_ne,<br>                               time_plus,<br>                               time_minus,<br>                               time_scalar_mult,<br>                               scalar_time_mult,<br>                               time_divide,<br>                               time_real_divide,<br>                               time_scalar_divide,<br>                               interval_alarm,<br>                               repeat_alarm,<br>                               set_calendar_type,<br>                               get_calendar_type,<br>                               get_date,<br>                               set_date,<br>                               increment_date,<br>                               decrement_date,<br>                               days_in_month,<br>                               leap_year,<br>                               length_of_year,<br>                               days_in_year,<br>                               month_name,<br>                               time_manager_init,<br>                               print_time,<br>                               print_date ]</pre>
<dl>
<dt>
<a href="#set_time">set_time</a>:</dt>
<dd>   Given some number of seconds and days, returns the
   corresponding time_type. </dd>
<dt>
<a href="#get_time">get_time</a>:</dt>
<dd>   Given a time interval, returns the corresponding seconds and days. </dd>
<dt>
<a href="#increment_time">increment_time</a>:</dt>
<dd>   Given a time and an increment of days and seconds, returns
   a time that adds this increment to an input time. </dd>
<dt>
<a href="#decrement_time">decrement_time</a>:</dt>
<dd>   Given a time and a decrement of days and seconds, returns
   a time that subtracts this decrement from an input time. </dd>
<dt>
<a href="#time_gt">time_gt</a>:</dt>
<dd>   Returns true if time1 &gt; time2. </dd>
<dt>
<a href="#time_ge">time_ge</a>:</dt>
<dd>   Returns true if time1 &gt;= time2. </dd>
<dt>
<a href="#time_lt">time_lt</a>:</dt>
<dd>   Returns true if time1 &lt; time2. </dd>
<dt>
<a href="#time_le">time_le</a>:</dt>
<dd>   Returns true if time1 &lt;= time2. </dd>
<dt>
<a href="#time_eq">time_eq</a>:</dt>
<dd>   Returns true if time1 == time2. </dd>
<dt>
<a href="#time_ne">time_ne</a>:</dt>
<dd>   Returns true if time1 /= time2. </dd>
<dt>
<a href="#time_plus">time_plus</a>:</dt>
<dd>   Returns sum of two time_types. </dd>
<dt>
<a href="#time_minus">time_minus</a>:</dt>
<dd>   Returns difference of two time_types. </dd>
<dt>
<a href="#time_scalar_mult">time_scalar_mult</a>:</dt>
<dd>   Returns time multiplied by integer factor n. </dd>
<dt>
<a href="#scalar_time_mult">scalar_time_mult</a>:</dt>
<dd>   Returns time multiplied by integer factor n. </dd>
<dt>
<a href="#time_divide">time_divide</a>:</dt>
<dd>   Returns the largest integer, n, for which time1 &gt;= time2 * n. </dd>
<dt>
<a href="#time_real_divide">time_real_divide</a>:</dt>
<dd>   Returns the double precision quotient of two times. </dd>
<dt>
<a href="#time_scalar_divide">time_scalar_divide</a>:</dt>
<dd>   Returns the largest time, t, for which n * t &lt;= time. </dd>
<dt>
<a href="#interval_alarm">interval_alarm</a>:</dt>
<dd>   Given a time, and a time interval, this function returns true
   if this is the closest time step to the alarm time. </dd>
<dt>
<a href="#repeat_alarm">repeat_alarm</a>:</dt>
<dd>   Repeat_alarm supports an alarm that goes off with
   alarm_frequency and lasts for alarm_length. </dd>
<dt>
<a href="#set_calendar_type">set_calendar_type</a>:</dt>
<dd>   Sets the default calendar type for mapping time intervals to dates. </dd>
<dt>
<a href="#get_calendar_type">get_calendar_type</a>:</dt>
<dd>   Returns the value of the default calendar type for mapping
   from time to date. </dd>
<dt>
<a href="#get_date">get_date</a>:</dt>
<dd>   Given a time_interval, returns the corresponding date under
   the selected calendar. </dd>
<dt>
<a href="#set_date">set_date</a>:</dt>
<dd>   Given an input date in year, month, days, etc., creates a
   time_type that represents this time interval from the
   internally defined base date. </dd>
<dt>
<a href="#increment_date">increment_date</a>:</dt>
<dd>   Increments the date represented by a time interval and the
   default calendar type by a number of seconds, etc. </dd>
<dt>
<a href="#decrement_date">decrement_date</a>:</dt>
<dd>   Decrements the date represented by a time interval and the
   default calendar type by a number of seconds, etc. </dd>
<dt>
<a href="#days_in_month">days_in_month</a>:</dt>
<dd>   Given a time interval, gives the number of days in the
   month corresponding to the default calendar. </dd>
<dt>
<a href="#leap_year">leap_year</a>:</dt>
<dd>   Returns true if the year corresponding to the date for the
   default calendar is a leap year. Returns false for
   THIRTY_DAY_MONTHS and NO_LEAP. </dd>
<dt>
<a href="#length_of_year">length_of_year</a>:</dt>
<dd>   Returns the mean length of the year in the default calendar setting. </dd>
<dt>
<a href="#days_in_year">days_in_year</a>:</dt>
<dd>   Returns the number of days in the calendar year corresponding to
   the date represented by time for the default calendar. </dd>
<dt>
<a href="#month_name">month_name</a>:</dt>
<dd>   Returns a character string containing the name of the
   month corresponding to month number n. </dd>
<dt>
<a href="#time_manager_init">time_manager_init</a>:</dt>
<dd>   Write the version information to the log file. </dd>
<dt>
<a href="#print_time">print_time</a>:</dt>
<dd>   Prints the given time_type argument as a time (using days and seconds). </dd>
<dt>
<a href="#print_date">print_date</a>:</dt>
<dd>   prints the time to standard output (or optional unit) as a date. </dd>
</dl>
</div>
<br>
<!-- END PUBLIC INTERFACE -->
<a name="PUBLIC DATA"></a>
<hr>
<h4>PUBLIC DATA</h4>
<!-- BEGIN PUBLIC DATA -->
<div>
<table align="center" cellspacing="2" CELLPADDING="2" BORDER="2">
<tr>
<th> Name  </th><th> Type  </th><th> Value  </th><th> Units  </th><th> Description  </th>
</tr>
<tr>
<td> time_type  </td><td> derived type  </td><td> ---  </td><td> ---  </td><td>    Derived-type data variable used to store time and date quantities. It
   contains two PRIVATE variables: seconds and days.   </td>
</tr>
</table>
<br>
</div>
<!-- END PUBLIC DATA -->
<a name="PUBLIC ROUTINES"></a>
<hr>
<h4>PUBLIC ROUTINES</h4>
<!-- BEGIN PUBLIC ROUTINES -->
<ol type="a">
<li>
<a name="set_time"></a>
<h4>set_time</h4>
<pre>&lt;B&gt; <b>set_time</b> &lt;/B&gt;(seconds, days)</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   Given some number of seconds and days, returns the
   corresponding time_type. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>seconds&nbsp;&nbsp;&nbsp;</tt></td><td>   A number of seconds (can be greater than 86400),  must be positive. <br>&nbsp;&nbsp;&nbsp;<span class="type">[integer, dimension(scalar)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>days&nbsp;&nbsp;&nbsp;</tt></td><td>   A number of days, must be positive. <br>&nbsp;&nbsp;&nbsp;<span class="type">[integer, dimension(scalar)]</span></td>
</tr>
</table>
</dd>
<br>
<dt>
<b>OUTPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>&nbsp;&nbsp;&nbsp;</tt></td><td>   A time interval corresponding to this number of days and seconds. <br>&nbsp;&nbsp;&nbsp;<span class="type">[, dimension]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="get_time"></a>
<h4>get_time</h4>
<pre>
<b>call get_time </b>&lt;/B&gt;(time, seconds, days)</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   Given a time interval, returns the corresponding seconds and days. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>time&nbsp;&nbsp;&nbsp;</tt></td><td>   A time interval. <br>&nbsp;&nbsp;&nbsp;<span class="type">[time_type]</span></td>
</tr>
</table>
</dd>
<br>
<dt>
<b>OUTPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>seconds&nbsp;&nbsp;&nbsp;</tt></td><td>   A number of seconds (&lt; 86400). <br>&nbsp;&nbsp;&nbsp;<span class="type">[integer, dimension(scalar)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>days&nbsp;&nbsp;&nbsp;</tt></td><td>   A number of days, must be positive. <br>&nbsp;&nbsp;&nbsp;<span class="type">[integer, dimension(scalar)]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="increment_time"></a>
<h4>increment_time</h4>
<pre> 
<b>increment_time</b> (time, seconds, days)</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   Given a time and an increment of days and seconds, returns
   a time that adds this increment to an input time.
   Increments a time by seconds and days; increments cannot be negative. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>time&nbsp;&nbsp;&nbsp;</tt></td><td>   A time interval. <br>&nbsp;&nbsp;&nbsp;<span class="type">[time_type, dimension]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>seconds&nbsp;&nbsp;&nbsp;</tt></td><td>   Increment of seconds (can be greater than 86400);  must be positive. <br>&nbsp;&nbsp;&nbsp;<span class="type">[integer, dimension(scalar)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>days&nbsp;&nbsp;&nbsp;</tt></td><td>   Increment of days;  must be positive. <br>&nbsp;&nbsp;&nbsp;<span class="type">[integer, dimension(scalar)]</span></td>
</tr>
</table>
</dd>
<br>
<dt>
<b>OUTPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>&nbsp;&nbsp;&nbsp;</tt></td><td>   A time that adds this increment to the input time. <br>&nbsp;&nbsp;&nbsp;<span class="type">[, dimension]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="decrement_time"></a>
<h4>decrement_time</h4>
<pre> 
<b>decrement_time</b> (time, seconds, days)</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   Decrements a time by seconds and days; decrements cannot be negative. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>time&nbsp;&nbsp;&nbsp;</tt></td><td>   A time interval. <br>&nbsp;&nbsp;&nbsp;<span class="type">[time_type, dimension]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>seconds&nbsp;&nbsp;&nbsp;</tt></td><td>   Decrement of seconds (can be greater than 86400);  must be positive. <br>&nbsp;&nbsp;&nbsp;<span class="type">[integer, dimension(scalar)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>days&nbsp;&nbsp;&nbsp;</tt></td><td>   Decrement of days;  must be positive. <br>&nbsp;&nbsp;&nbsp;<span class="type">[integer, dimension(scalar)]</span></td>
</tr>
</table>
</dd>
<br>
<dt>
<b>OUTPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>&nbsp;&nbsp;&nbsp;</tt></td><td>   A time that subtracts this decrement from an input time. If
   the result is negative, it is considered a fatal error. <br>&nbsp;&nbsp;&nbsp;<span class="type">[, dimension]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="time_gt"></a>
<h4>time_gt</h4>
<pre>&lt;B&gt; <b>time_gt</b> &lt;/B&gt;(time1, time2)</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   Returns true if time1 &gt; time2. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>time1&nbsp;&nbsp;&nbsp;</tt></td><td>   A time interval. <br>&nbsp;&nbsp;&nbsp;<span class="type">[time_type, dimension]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>time2&nbsp;&nbsp;&nbsp;</tt></td><td>   A time interval. <br>&nbsp;&nbsp;&nbsp;<span class="type">[time_type, dimension]</span></td>
</tr>
</table>
</dd>
<br>
<dt>
<b>OUTPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>&nbsp;&nbsp;&nbsp;</tt></td><td>   Returns true if time1 &gt; time2 <br>&nbsp;&nbsp;&nbsp;<span class="type">[logical, dimension]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="time_ge"></a>
<h4>time_ge</h4>
<pre>&lt;B&gt; <b>time_ge</b> &lt;/B&gt;(time1, time2)</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   Returns true if time1 &gt;= time2. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>time1&nbsp;&nbsp;&nbsp;</tt></td><td>   A time interval. <br>&nbsp;&nbsp;&nbsp;<span class="type">[time_type, dimension]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>time2&nbsp;&nbsp;&nbsp;</tt></td><td>   A time interval. <br>&nbsp;&nbsp;&nbsp;<span class="type">[time_type, dimension]</span></td>
</tr>
</table>
</dd>
<br>
<dt>
<b>OUTPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>&nbsp;&nbsp;&nbsp;</tt></td><td>   Returns true if time1 &gt;= time2 <br>&nbsp;&nbsp;&nbsp;<span class="type">[logical, dimension]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="time_lt"></a>
<h4>time_lt</h4>
<pre>&lt;B&gt; <b>time_lt</b> &lt;/B&gt;(time1, time2)</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   Returns true if time1 &lt; time2. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>time1&nbsp;&nbsp;&nbsp;</tt></td><td>   A time interval. <br>&nbsp;&nbsp;&nbsp;<span class="type">[time_type, dimension]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>time2&nbsp;&nbsp;&nbsp;</tt></td><td>   A time interval. <br>&nbsp;&nbsp;&nbsp;<span class="type">[time_type, dimension]</span></td>
</tr>
</table>
</dd>
<br>
<dt>
<b>OUTPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>&nbsp;&nbsp;&nbsp;</tt></td><td>   Returns true if time1 &lt; time2 <br>&nbsp;&nbsp;&nbsp;<span class="type">[logical, dimension]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="time_le"></a>
<h4>time_le</h4>
<pre>&lt;B&gt; <b>time_le</b> &lt;/B&gt;(time1, time2)</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   Returns true if time1 &lt;= time2. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>time1&nbsp;&nbsp;&nbsp;</tt></td><td>   A time interval. <br>&nbsp;&nbsp;&nbsp;<span class="type">[time_type, dimension]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>time2&nbsp;&nbsp;&nbsp;</tt></td><td>   A time interval. <br>&nbsp;&nbsp;&nbsp;<span class="type">[time_type, dimension]</span></td>
</tr>
</table>
</dd>
<br>
<dt>
<b>OUTPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>&nbsp;&nbsp;&nbsp;</tt></td><td>   Returns true if time1 &lt;= time2 <br>&nbsp;&nbsp;&nbsp;<span class="type">[logical, dimension]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="time_eq"></a>
<h4>time_eq</h4>
<pre>&lt;B&gt; <b>time_eq</b> &lt;/B&gt;(time1, time2)</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   Returns true if time1 == time2. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>time1&nbsp;&nbsp;&nbsp;</tt></td><td>   A time interval. <br>&nbsp;&nbsp;&nbsp;<span class="type">[time_type, dimension]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>time2&nbsp;&nbsp;&nbsp;</tt></td><td>   A time interval. <br>&nbsp;&nbsp;&nbsp;<span class="type">[time_type, dimension]</span></td>
</tr>
</table>
</dd>
<br>
<dt>
<b>OUTPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>&nbsp;&nbsp;&nbsp;</tt></td><td>   Returns true if time1 == time2 <br>&nbsp;&nbsp;&nbsp;<span class="type">[logical, dimension]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="time_ne"></a>
<h4>time_ne</h4>
<pre>&lt;B&gt; <b>time_ne</b> &lt;/B&gt;(time1, time2)</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   Returns true if time1 /= time2. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>time1&nbsp;&nbsp;&nbsp;</tt></td><td>   A time interval. <br>&nbsp;&nbsp;&nbsp;<span class="type">[time_type, dimension]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>time2&nbsp;&nbsp;&nbsp;</tt></td><td>   A time interval. <br>&nbsp;&nbsp;&nbsp;<span class="type">[time_type, dimension]</span></td>
</tr>
</table>
</dd>
<br>
<dt>
<b>OUTPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>&nbsp;&nbsp;&nbsp;</tt></td><td>   Returns true if time1 /= time2 <br>&nbsp;&nbsp;&nbsp;<span class="type">[logical, dimension]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="time_plus"></a>
<h4>time_plus</h4>
<pre>&lt;B&gt; <b>time_plus</b> &lt;/B&gt;(time1, time2)</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   Returns sum of two time_types. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>time1&nbsp;&nbsp;&nbsp;</tt></td><td>   A time interval. <br>&nbsp;&nbsp;&nbsp;<span class="type">[time_type, dimension]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>time2&nbsp;&nbsp;&nbsp;</tt></td><td>   A time interval. <br>&nbsp;&nbsp;&nbsp;<span class="type">[time_type, dimension]</span></td>
</tr>
</table>
</dd>
<br>
<dt>
<b>OUTPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>&nbsp;&nbsp;&nbsp;</tt></td><td>   Returns sum of two time_types. <br>&nbsp;&nbsp;&nbsp;<span class="type">[time_type, dimension]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="time_minus"></a>
<h4>time_minus</h4>
<pre>&lt;B&gt; <b>time_minus</b> &lt;/B&gt;(time1, time2)</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   Returns difference of two time_types. WARNING: a time type is positive 
   so by definition time1 - time2  is the same as time2 - time1. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>time1&nbsp;&nbsp;&nbsp;</tt></td><td>   A time interval. <br>&nbsp;&nbsp;&nbsp;<span class="type">[time_type, dimension]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>time2&nbsp;&nbsp;&nbsp;</tt></td><td>   A time interval. <br>&nbsp;&nbsp;&nbsp;<span class="type">[time_type, dimension]</span></td>
</tr>
</table>
</dd>
<br>
<dt>
<b>OUTPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>&nbsp;&nbsp;&nbsp;</tt></td><td>   Returns difference of two time_types. <br>&nbsp;&nbsp;&nbsp;<span class="type">[time_type, dimension]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="time_scalar_mult"></a>
<h4>time_scalar_mult</h4>
<pre>&lt;B&gt; <b>time_scalar_mult</b> &lt;/B&gt;(time, n)</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   Returns time multiplied by integer factor n. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>time&nbsp;&nbsp;&nbsp;</tt></td><td>   A time interval. <br>&nbsp;&nbsp;&nbsp;<span class="type">[time_type, dimension]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>n&nbsp;&nbsp;&nbsp;</tt></td><td>   A time interval. <br>&nbsp;&nbsp;&nbsp;<span class="type">[integer, dimension]</span></td>
</tr>
</table>
</dd>
<br>
<dt>
<b>OUTPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>&nbsp;&nbsp;&nbsp;</tt></td><td>   Returns time multiplied by integer factor n. <br>&nbsp;&nbsp;&nbsp;<span class="type">[time_type, dimension]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="scalar_time_mult"></a>
<h4>scalar_time_mult</h4>
<pre>&lt;B&gt; <b>scalar_time_mult</b> &lt;/B&gt;(n, time)</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   Returns time multiplied by integer factor n. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>time&nbsp;&nbsp;&nbsp;</tt></td><td>A time interval.<br>&nbsp;&nbsp;&nbsp;<span class="type">[time_type, dimension]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>n&nbsp;&nbsp;&nbsp;</tt></td><td>   An integer. <br>&nbsp;&nbsp;&nbsp;<span class="type">[integer, dimension]</span></td>
</tr>
</table>
</dd>
<br>
<dt>
<b>OUTPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>&nbsp;&nbsp;&nbsp;</tt></td><td>   Returns time multiplied by integer factor n. <br>&nbsp;&nbsp;&nbsp;<span class="type">[time_type, dimension]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="time_divide"></a>
<h4>time_divide</h4>
<pre>&lt;B&gt; <b>time_divide</b> &lt;/B&gt;(time1, time2)</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   Returns the largest integer, n, for which time1 &gt;= time2 * n. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>time1&nbsp;&nbsp;&nbsp;</tt></td><td>   A time interval. <br>&nbsp;&nbsp;&nbsp;<span class="type">[time_type, dimension]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>time2&nbsp;&nbsp;&nbsp;</tt></td><td>   A time interval. <br>&nbsp;&nbsp;&nbsp;<span class="type">[time_type, dimension]</span></td>
</tr>
</table>
</dd>
<br>
<dt>
<b>OUTPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>&nbsp;&nbsp;&nbsp;</tt></td><td>   Returns the largest integer, n, for which time1 &gt;= time2 * n. <br>&nbsp;&nbsp;&nbsp;<span class="type">[integer, dimension]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="time_real_divide"></a>
<h4>time_real_divide</h4>
<pre>&lt;B&gt; <b>time_real_divide</b> &lt;/B&gt;(time1, time2)</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   Returns the double precision quotient of two times. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>time1&nbsp;&nbsp;&nbsp;</tt></td><td>   A time interval. <br>&nbsp;&nbsp;&nbsp;<span class="type">[time_type, dimension]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>time2&nbsp;&nbsp;&nbsp;</tt></td><td>   A time interval. <br>&nbsp;&nbsp;&nbsp;<span class="type">[time_type, dimension]</span></td>
</tr>
</table>
</dd>
<br>
<dt>
<b>OUTPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>&nbsp;&nbsp;&nbsp;</tt></td><td>   Returns the double precision quotient of two times <br>&nbsp;&nbsp;&nbsp;<span class="type">[integer, dimensiondouble precision]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="time_scalar_divide"></a>
<h4>time_scalar_divide</h4>
<pre>&lt;B&gt; <b>time_scalar_divide</b> &lt;/B&gt;(time, n)</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   Returns the largest time, t, for which n * t &lt;= time. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>time&nbsp;&nbsp;&nbsp;</tt></td><td>   A time interval. <br>&nbsp;&nbsp;&nbsp;<span class="type">[time_type, dimension]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>n&nbsp;&nbsp;&nbsp;</tt></td><td>   An integer factor. <br>&nbsp;&nbsp;&nbsp;<span class="type">[integer, dimension]</span></td>
</tr>
</table>
</dd>
<br>
<dt>
<b>OUTPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>&nbsp;&nbsp;&nbsp;</tt></td><td>   Returns the largest time, t, for which n * t &lt;= time. <br>&nbsp;&nbsp;&nbsp;<span class="type">[integer, dimensiondouble precision]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="interval_alarm"></a>
<h4>interval_alarm</h4>
<pre> 
<b>interval_alarm</b> (time, time_interval, alarm, alarm_interval)</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   This is a specialized operation that is frequently performed in models.
   Given a time, and a time interval, this function is true if this is the
   closest time step to the alarm time. The actual computation is:
   <br>
<br>
   if((alarm_time - time) &lt;= (time_interval / 2))
   <br>
<br>
   If the function is true, the alarm time is incremented by the
   alarm_interval; WARNING, this is a featured side effect. Otherwise, the
   function is false and there are no other effects. CAUTION: if the
   alarm_interval is smaller than the time_interval, the alarm may fail to
   return true ever again.  Watch
   for problems if the new alarm time is less than time + time_interval </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>time&nbsp;&nbsp;&nbsp;</tt></td><td>   Current time. <br>&nbsp;&nbsp;&nbsp;<span class="type">[time_type]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>time_interval&nbsp;&nbsp;&nbsp;</tt></td><td>   A time interval. <br>&nbsp;&nbsp;&nbsp;<span class="type">[time_type]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>alarm_interval&nbsp;&nbsp;&nbsp;</tt></td><td>   A time interval. <br>&nbsp;&nbsp;&nbsp;<span class="type">[time_type]</span></td>
</tr>
</table>
</dd>
<br>
<dt>
<b>INPUT/OUTPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>alarm&nbsp;&nbsp;&nbsp;</tt></td><td>   An alarm time, which is incremented by the alarm_interval
   if the function is true. <br>&nbsp;&nbsp;&nbsp;<span class="type">[time_type]</span></td>
</tr>
</table>
</dd>
<br>
<dt>
<b>OUTPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>interval_alarm&nbsp;&nbsp;&nbsp;</tt></td><td>   Returns either True or false. <br>&nbsp;&nbsp;&nbsp;<span class="type">[logical]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="repeat_alarm"></a>
<h4>repeat_alarm</h4>
<pre> 
<b>repeat_alarm</b> 
</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   Repeat_alarm supports an alarm that goes off with alarm_frequency and
   lasts for alarm_length.  If the nearest occurence of an alarm time
   is less than half an alarm_length from the input time, repeat_alarm
   is true.  For instance, if the alarm_frequency is 1 day, and the 
   alarm_length is 2 hours, then repeat_alarm is true from time 2300 on 
   day n to time 0100 on day n + 1 for all n. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>time&nbsp;&nbsp;&nbsp;</tt></td><td>   Current time. <br>&nbsp;&nbsp;&nbsp;<span class="type">[time_type]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>alarm_frequency&nbsp;&nbsp;&nbsp;</tt></td><td>   A time interval for alarm_frequency. <br>&nbsp;&nbsp;&nbsp;<span class="type">[time_type]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>alarm_length&nbsp;&nbsp;&nbsp;</tt></td><td>   A time interval for alarm_length. <br>&nbsp;&nbsp;&nbsp;<span class="type">[time_type]</span></td>
</tr>
</table>
</dd>
<br>
<dt>
<b>OUTPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>repeat_alarm&nbsp;&nbsp;&nbsp;</tt></td><td>   Returns either True or false. <br>&nbsp;&nbsp;&nbsp;<span class="type">[logical]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="set_calendar_type"></a>
<h4>set_calendar_type</h4>
<pre>
<b>call set_calendar_type </b>(type)</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   A constant number for setting the calendar type. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>type&nbsp;&nbsp;&nbsp;</tt></td><td>   A constant number for setting the calendar type. <br>&nbsp;&nbsp;&nbsp;<span class="type">[integer, dimension]</span></td>
</tr>
</table>
</dd>
<br>
<dt>
<b>OUTPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>calendar_type&nbsp;&nbsp;&nbsp;</tt></td><td>   A constant number for default calendar type. <br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
</table>
</dd>
<br>
<dt>
<b>NOTE</b>
</dt>
<dd>   At present, four integer constants are defined for setting
   the calendar type: THIRTY_DAY_MONTHS, JULIAN, NO_LEAP, and
   GREGORIAN. However, GREGORIAN CALENDAR is not completely
   implemented. Selection of this type will result in illegal
   type error. </dd>
<br>
<br>
</dl>
</li>
<li>
<a name="get_calendar_type"></a>
<h4>get_calendar_type</h4>
<pre> 
<b>get_calendar_type</b> ()</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   There are no arguments in this function. It returns the value of
   the default calendar type for mapping from time to date. </dd>
<br>
<br>
</dl>
</li>
<li>
<a name="get_date"></a>
<h4>get_date</h4>
<pre>
<b>call get_date </b>(time, year, month, day, hour, minute, second)</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   Given a time_interval, returns the corresponding date under
   the selected calendar. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>time&nbsp;&nbsp;&nbsp;</tt></td><td>   A time interval.<br>&nbsp;&nbsp;&nbsp;<span class="type">[time_type]</span></td>
</tr>
</table>
</dd>
<br>
<dt>
<b>OUTPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>day&nbsp;&nbsp;&nbsp;</tt></td><td>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>month&nbsp;&nbsp;&nbsp;</tt></td><td>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>year&nbsp;&nbsp;&nbsp;</tt></td><td>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>second&nbsp;&nbsp;&nbsp;</tt></td><td>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>minute&nbsp;&nbsp;&nbsp;</tt></td><td>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>hour&nbsp;&nbsp;&nbsp;</tt></td><td>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
</table>
</dd>
<br>
<dt>
<b>NOTE</b>
</dt>
<dd>   For all but the thirty_day_months calendar, increments to months
   and years must be made separately from other units because of the
   non-associative nature of the addition. All the input increments
   must be positive. </dd>
<br>
<br>
</dl>
</li>
<li>
<a name="set_date"></a>
<h4>set_date</h4>
<pre> 
<b>set_date</b> (year, month, day, hours, minutes, seconds)</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   Given a date, computes the corresponding time given the selected
   date time mapping algorithm. Note that it is possible to specify
   any number of illegal dates; these should be checked for and generate
   errors as appropriate. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>time&nbsp;&nbsp;&nbsp;</tt></td><td>   A time interval.<br>&nbsp;&nbsp;&nbsp;<span class="type">[time_type]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>day&nbsp;&nbsp;&nbsp;</tt></td><td>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>month&nbsp;&nbsp;&nbsp;</tt></td><td>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>year&nbsp;&nbsp;&nbsp;</tt></td><td>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>second&nbsp;&nbsp;&nbsp;</tt></td><td>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>minute&nbsp;&nbsp;&nbsp;</tt></td><td>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>hour&nbsp;&nbsp;&nbsp;</tt></td><td>
<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
</table>
</dd>
<br>
<dt>
<b>OUTPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>set_date&nbsp;&nbsp;&nbsp;</tt></td><td>   A time interval.<br>&nbsp;&nbsp;&nbsp;<span class="type">[time_type]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="increment_date"></a>
<h4>increment_date</h4>
<pre> 
<b>increment_date</b> (time, years, months, days, hours, minutes, seconds)</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   Given a time and some date increment, computes a new time.  Depending
   on the mapping algorithm from date to time, it may be possible to specify
   undefined increments (i.e. if one increments by 68 days and 3 months in
   a Julian calendar, it matters which order these operations are done and
   we don't want to deal with stuff like that, make it an error). </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>time&nbsp;&nbsp;&nbsp;</tt></td><td>   A time interval.<br>&nbsp;&nbsp;&nbsp;<span class="type">[time_type]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>day&nbsp;&nbsp;&nbsp;</tt></td><td>An increment of days.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>month&nbsp;&nbsp;&nbsp;</tt></td><td>An increment of months.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>year&nbsp;&nbsp;&nbsp;</tt></td><td>An increment of years.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>second&nbsp;&nbsp;&nbsp;</tt></td><td>An increment of seconds.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>minute&nbsp;&nbsp;&nbsp;</tt></td><td>An increment of minutes.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>hour&nbsp;&nbsp;&nbsp;</tt></td><td>An increment of hours.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
</table>
</dd>
<br>
<dt>
<b>OUTPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>increment_date&nbsp;&nbsp;&nbsp;</tt></td><td>   A new time based on the input 
   time interval and the default calendar type. <br>&nbsp;&nbsp;&nbsp;<span class="type">[time_type]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="decrement_date"></a>
<h4>decrement_date</h4>
<pre> 
<b>decrement_date</b> (time, years, months, days, hours, minutes, seconds)</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   Given a time and some date decrement, computes a new time.  Depending
   on the mapping algorithm from date to time, it may be possible to specify
   undefined decrements (i.e. if one decrements by 68 days and 3 months in
   a Julian calendar, it matters which order these operations are done and
   we don't want to deal with stuff like that, make it an error). </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>time&nbsp;&nbsp;&nbsp;</tt></td><td>   A time interval.<br>&nbsp;&nbsp;&nbsp;<span class="type">[time_type]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>day&nbsp;&nbsp;&nbsp;</tt></td><td>A decrement of days.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>month&nbsp;&nbsp;&nbsp;</tt></td><td>A deincrement of months.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>year&nbsp;&nbsp;&nbsp;</tt></td><td>A deincrement of years.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>second&nbsp;&nbsp;&nbsp;</tt></td><td>A deincrement of seconds.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>minute&nbsp;&nbsp;&nbsp;</tt></td><td>A deincrement of minutes.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>hour&nbsp;&nbsp;&nbsp;</tt></td><td>A deincrement of hours.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
</table>
</dd>
<br>
<dt>
<b>OUTPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>decrement_date&nbsp;&nbsp;&nbsp;</tt></td><td>   A new time based on the input 
   time interval and the default calendar type. <br>&nbsp;&nbsp;&nbsp;<span class="type">[time_type]</span></td>
</tr>
</table>
</dd>
<br>
<dt>
<b>NOTE</b>
</dt>
<dd>   For all but the thirty_day_months calendar, decrements to months
   and years must be made separately from other units because of the
   non-associative nature of addition. All the input decrements must
   be positive. If the result is a negative time (i.e. date before the
   base date) it is considered a fatal error. </dd>
<br>
<br>
</dl>
</li>
<li>
<a name="days_in_month"></a>
<h4>days_in_month</h4>
<pre>&lt;B&gt; <b>days_in_month</b> (time)</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   Given a time, computes the corresponding date given the selected
   date time mapping algorithm. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>time&nbsp;&nbsp;&nbsp;</tt></td><td>A time interval.<br>&nbsp;&nbsp;&nbsp;<span class="type">[time_type, dimension]</span></td>
</tr>
</table>
</dd>
<br>
<dt>
<b>OUTPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>days_in_month&nbsp;&nbsp;&nbsp;</tt></td><td>   The number of days in the month given the selected time
   mapping algorithm. <br>&nbsp;&nbsp;&nbsp;<span class="type">[integer, dimension]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="leap_year"></a>
<h4>leap_year</h4>
<pre> 
<b>leap_year</b> (time)</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   Is this date in a leap year for default calendar? Returns true if the 
   year corresponding to the date for the default calendar is a leap year. 
   Returns false for THIRTY_DAY_MONTHS and NO_LEAP. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>time&nbsp;&nbsp;&nbsp;</tt></td><td>A time interval.<br>&nbsp;&nbsp;&nbsp;<span class="type">[time_type, dimension]</span></td>
</tr>
</table>
</dd>
<br>
<dt>
<b>OUTPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>leap_year&nbsp;&nbsp;&nbsp;</tt></td><td>   True if the year corresponding to the date for the default
   calendar is a leap year. False for THIRTY_DAY_MONTHS and
   NO_LEAP and otherwise. <br>&nbsp;&nbsp;&nbsp;<span class="type">[calendar_type, dimension]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="length_of_year"></a>
<h4>length_of_year</h4>
<pre> 
<b>length_of_year</b> ()</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   There are no arguments in this function. It returns the mean
   length of the year in the default calendar setting. </dd>
<br>
<br>
</dl>
</li>
<li>
<a name="days_in_year"></a>
<h4>days_in_year</h4>
<pre> 
<b>days_in_year</b> ()</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   Returns the number of days in the calendar year corresponding to
   the date represented by time for the default calendar. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>time&nbsp;&nbsp;&nbsp;</tt></td><td>A time interval.<br>&nbsp;&nbsp;&nbsp;<span class="type">[time_type]</span></td>
</tr>
</table>
</dd>
<br>
<dt>
<b>OUTPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>&nbsp;&nbsp;&nbsp;</tt></td><td>   The number of days in this year for the default calendar type. </td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="month_name"></a>
<h4>month_name</h4>
<pre> 
<b>month_name</b> (n)</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   Returns a character string containing the name of the
   month corresponding to month number n. Definition is the
   same for all calendar types. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>n&nbsp;&nbsp;&nbsp;</tt></td><td>Month number.<br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
</table>
</dd>
<br>
<dt>
<b>OUTPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>month_name&nbsp;&nbsp;&nbsp;</tt></td><td>   The character string associated with a month. For now all
   calendars have 12 months and will return standard names. <br>&nbsp;&nbsp;&nbsp;<span class="type">[character]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="time_manager_init"></a>
<h4>time_manager_init</h4>
<pre> 
<b>time_manager_init</b> ()</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   Initialization routine. This routine does not have to be called, all 
   it does is write the version information to the log file. </dd>
<br>
<br>
</dl>
</li>
<li>
<a name="print_time"></a>
<h4>print_time</h4>
<pre> 
<b>print_time</b> (time,str,unit)</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   Prints the given time_type argument either as a time (using days and
   seconds). NOTE: there is no check for PE number. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>time&nbsp;&nbsp;&nbsp;</tt></td><td>   Time that will be printed. <br>&nbsp;&nbsp;&nbsp;<span class="type">[time_type]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>str&nbsp;&nbsp;&nbsp;</tt></td><td>   Character string that precedes the printed time or date. <br>&nbsp;&nbsp;&nbsp;<span class="type">[character (len=*)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>unit&nbsp;&nbsp;&nbsp;</tt></td><td>   Unit number for printed output. The default unit is stdout. <br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
<li>
<a name="print_date"></a>
<h4>print_date</h4>
<pre> 
<b>print_date</b> (time,str,unit)</pre>
<dl>
<dt>
<b>DESCRIPTION</b>
</dt>
<dd>   Prints the given time_type argument as a date (using year,month,day,
   hour,minutes and seconds). NOTE: there is no check for PE number. </dd>
<br>
<br>
<dt>
<b>INPUT</b>
</dt>
<dd>
<table border="0">
<tr>
<td valign="top" align="left"><tt>time&nbsp;&nbsp;&nbsp;</tt></td><td>   Time that will be printed. <br>&nbsp;&nbsp;&nbsp;<span class="type">[time_type]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>str&nbsp;&nbsp;&nbsp;</tt></td><td>   Character string that precedes the printed time or date. <br>&nbsp;&nbsp;&nbsp;<span class="type">[character (len=*)]</span></td>
</tr>
<tr>
<td valign="top" align="left"><tt>unit&nbsp;&nbsp;&nbsp;</tt></td><td>   Unit number for printed output. The default unit is stdout. <br>&nbsp;&nbsp;&nbsp;<span class="type">[integer]</span></td>
</tr>
</table>
</dd>
<br>
</dl>
</li>
</ol>
<!-- END PUBLIC ROUTINES -->
<a name="NAMELIST"></a>
<!-- BEGIN NAMELIST -->
<!-- END NAMELIST --><a name="DIAGNOSTIC FIELDS"></a>
<!-- BEGIN DIAGNOSTIC FIELDS -->
<!-- END DIAGNOSTIC FIELDS --><a name="DATA SETS"></a>
<!-- BEGIN DATA SETS -->
<hr>
<h4>DATA SETS</h4>
<div>None.<br>
<br>
</div>
<!-- END DATA SETS -->
<a name="ERROR MESSAGES"></a>
<!-- BEGIN ERROR MESSAGES -->
<hr>
<h4>ERROR MESSAGES</h4>
<div>None.<br>
<br>
</div>
<!-- END ERROR MESSAGES -->
<a name="REFERENCES"></a>
<hr>
<h4>REFERENCES</h4>
<!-- BEGIN REFERENCES -->
<div>
        None.
      </div>
<br>
<!-- END REFERENCES -->
<a name="COMPILER SPECIFICS"></a>
<hr>
<h4>COMPILER SPECIFICS</h4>
<!-- BEGIN COMPILER SPECIFICS -->
<div>
        None.
      </div>
<br>
<!-- END COMPILER SPECIFICS -->
<a name="PRECOMPILER OPTIONS"></a>
<hr>
<h4>PRECOMPILER OPTIONS</h4>
<!-- BEGIN PRECOMPILER OPTIONS -->
<div>
        None.
      </div>
<br>
<!-- END PRECOMPILER OPTIONS -->
<a name="LOADER OPTIONS"></a>
<hr>
<h4>LOADER OPTIONS</h4>
<!-- BEGIN LOADER -->
<div>None.<br>
<br>
</div>
<!-- END LOADER OPTIONS -->
<a name="TEST PROGRAM"></a>
<hr>
<h4>TEST PROGRAM</h4>
<!-- BEGIN TEST PROGRAM -->
<div>
<dl>
<dt>time_main2</dt>
<dd> 
<pre>        use time_manager_mod
        implicit none
        type(time_type) :: dt, init_date, astro_base_date, time, final_date
        type(time_type) :: next_rad_time, mid_date
        type(time_type) :: repeat_alarm_freq, repeat_alarm_length
        integer :: num_steps, i, days, months, years, seconds, minutes, hours
        integer :: months2, length
        real :: astro_days
   
Set calendar type
    call set_calendar_type(THIRTY_DAY_MONTHS)
        call set_calendar_type(JULIAN)
    call set_calendar_type(NO_LEAP)
   
 Set timestep
        dt = set_time(1100, 0)
   
 Set initial date
        init_date = set_date(1992, 1, 1)
   
 Set date for astronomy delta calculation
        astro_base_date = set_date(1970, 1, 1, 12, 0, 0)
   
 Copy initial time to model current time
        time = init_date
   
 Determine how many steps to do to run one year
        final_date = increment_date(init_date, years = 1)
        num_steps = (final_date - init_date) / dt
        write(*, *) 'Number of steps is' , num_steps
   
 Want to compute radiation at initial step, then every two hours
        next_rad_time = time + set_time(7200, 0)
   
 Test repeat alarm
        repeat_alarm_freq = set_time(0, 1)
        repeat_alarm_length = set_time(7200, 0)
   
 Loop through a year
        do i = 1, num_steps
   
 Increment time
        time = time + dt
   
 Test repeat alarm
        if(repeat_alarm(time, repeat_alarm_freq, repeat_alarm_length)) &amp;
        write(*, *) 'REPEAT ALARM IS TRUE'
   
 Should radiation be computed? Three possible tests.
 First test assumes exact interval; just ask if times are equal
     if(time == next_rad_time) then
 Second test computes rad on last time step that is &lt;= radiation time
     if((next_rad_time - time) &lt; dt .and. time &lt; next_rad) then
 Third test computes rad on time step closest to radiation time
         if(interval_alarm(time, dt, next_rad_time, set_time(7200, 0))) then
           call get_date(time, years, months, days, hours, minutes, seconds)
           write(*, *) days, month_name(months), years, hours, minutes, seconds
   
 Need to compute real number of days between current time and astro_base
           call get_time(time - astro_base_date, seconds, days)
           astro_days = days + seconds / 86400.
       write(*, *) 'astro offset ', astro_days
        end if
   
 Can compute daily, monthly, yearly, hourly, etc. diagnostics as for rad
   
 Example: do diagnostics on last time step of this month
        call get_date(time + dt, years, months2, days, hours, minutes, seconds)
        call get_date(time, years, months, days, hours, minutes, seconds)
        if(months /= months2) then
           write(*, *) 'last timestep of month'
           write(*, *) days, months, years, hours, minutes, seconds
        endif
   
 Example: mid-month diagnostics; inefficient to make things clear
        length = days_in_month(time)
        call get_date(time, years, months, days, hours, minutes, seconds)
        mid_date = set_date(years, months, 1) + set_time(0, length) / 2
   
        if(time &lt; mid_date .and. (mid_date - time) &lt; dt) then
           write(*, *) 'mid-month time'
           write(*, *) days, months, years, hours, minutes, seconds
        endif
   
        end do</pre>   end program time_main2 </dd>
</dl>
</div>
<br>
<!-- END TEST PROGRAM -->
<a name="KNOWN BUGS"></a>
<hr>
<h4>KNOWN BUGS</h4>
<!-- BEGIN KNOWN BUGS -->
<div>
        None.
      </div>
<br>
<!-- END KNOWN BUGS -->
<a name="NOTES"></a>
<hr>
<h4>NOTES</h4>
<!-- BEGIN NOTES -->
<div>   The Gregorian calendar type is not completely implemented, and currently
   no effort is put on it since it doesn't differ from Julian in use between
   1901 and 2099. <br>
<br>   The &lt;a name="base date"&gt;base date&lt;/a&gt; is implicitly defined so users don't 
   need to be concerned with it. For the curious, the base date is defined as 
   0 seconds, 0 minutes, 0 hours, day 1, month 1, year 1 for the Julian and 
   thirty_day_months calendars, and 1 January, 1900, 0 seconds, 0 minutes, 
   0 hour for the Gregorian calendar. <br>
<br>   Please note that a time is a positive definite quantity. <br>
<br>   See the <a href="TEST PROGRAM">Test Program </a>   for a simple program 
   that shows some of the capabilities of the time manager. </div>
<br>
<!-- END NOTES -->
<a name="FUTURE PLANS"></a>
<hr>
<h4>FUTURE PLANS</h4>
<!-- BEGIN FUTURE PLANS -->
<div>
        None.
      </div>
<br>
<!-- END FUTURE PLANS -->
<hr>
<div align="right">
<font size="-1"><a href="#TOP">top</a></font>
</div>
</body>
</html>
