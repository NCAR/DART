<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>program obs_seq_coverage</title>
<link rel="stylesheet" type="text/css" href=
"../../../docs/html/doc.css">
<link href="../../../docs/images/dart.ico" rel="shortcut icon">
</head>
<body>
<h1>program <em class="program">obs_seq_coverage</em></h1>
<table border="0" summary="" cellpadding="5">
<tr>
<td valign="middle"><img src="../../../docs/images/Dartboard7.png"
alt="DART project logo" height="70"></td>
<td>Jump to <a href="../../../docs/index.html">DART Documentation
Main Index</a></td>
</tr>
</table>
<a href="#Namelist">NAMELIST</a> / <a href="#Modules">MODULES</a> /
<a href="#FilesUsed">FILES</a> / <a href="#Usage">USAGE</a> /
<a href="#References">REFERENCES</a> / <a href="#Errors">ERRORS</a>
/ <a href="#FuturePlans">PLANS</a> / <a href="#Legalese">TERMS OF
USE</a>
<h2>Overview</h2>
<p><em class="program">obs_seq_coverage</em> queries a set of
observation sequence files to determine which observation locations
report frequently enough to be useful for a verification study. The
big picture is to be able to pare down a large set of observations
into a compact observation sequence file to run through <a href=
"../filter/filter.html">filter</a> with all of the intended
observation types flagged as <em class="option">evaluate_only</em>.
DART's forward operators then get applied and all the forecasts are
preserved in a standard <em class="file">obs_seq.final</em> file -
perhaps more appropriately called <em class=
"file">obs_seq.forecast</em>! Paring down the input observation
sequence file cuts down on the unnecessary application of the
forward operator to create observation copies that will not be used
anyway ...</p>
<img src="../../../docs/images/forecasting_diagram.png" alt=
"forecast evaluation schematic" width="90%">
<p><em class="program">obs_seq_coverage</em> results in two output
files:</p>
<ul>
<li><em class="file">obsdef_mask.txt</em> contains the list of
observation definitions (but not the observations themselves) that
are desired. The observation definitions include the locations and
times for each of the desired observation types. This file is read
by <a href=
"../../../assimilation_code/programs/obs_selection/obs_selection.html">
obs_selection</a> and combined with the raw observation sequence
files to create the observation sequence file appropriate for use
in a forecast.<br></li>
<li><em class="file">obsdef_mask.nc</em> contains information
needed to be able to plot the times and locations of the
observations in a manner to help explore the design of the
verification locations/network. <em class=
"file">obsdef_mask.nc</em> is <em class="strong">required</em> by
<a href=
"../../../assimilation_code/programs/obs_seq_verify/obs_seq_verify.html">
obs_seq_verify</a>, the program that reorders the observations into
a structure that makes it easy to calculate statistics like ROC,
etc.</li>
</ul>
<p>The following section explains the strategy and requirements for
determining what observations will be used to verify a forecast.
Since it is 'standard practice' to make several forecasts to build
statistical strength, it is important to use the SAME set of
observation locations for all the forecasts that will be verified
together. To make the discussion easier, let's define the
<em class="option">verification network</em> as the set of
locations and times for a particular observation type.<br>
<br>
The entire discussion about finding locations that are repeatedly
observed through time boils down to the simple statement that if
the observation is within about 500cm of a previous observation,
they are treated as co-located observations. For some very high
resolution applications, this may be insufficient, but there it is.
For observations at pressure levels, see the <a href="#Usage">Word
about vertical levels</a>.<br>
<br>
The only complicated part of determining the verification network
is the temporal component. The initial time (usually an <em class=
"option">analysis time</em> from a previous assimilation), the
<em class="option">verification interval</em>, and the <em class=
"option">forecast length</em> completely specify the temporal
aspect of a forecast. The following example has a verification
interval of 6 hours and a forecast length of 24 hours. We adopt the
convention of also including the initial conditions (a "nowcast")
in the "forecast", so there are 5 times of interest - which we will
call <em class="option">verification times</em> and are represented
by <img src="../../../docs/images/verification_time_icon.png" alt=
"verification icon" width="3%">. The candidate observation sequence
files are scanned to select all the observations that are
<strong>closest</strong> to the verification times. The difference
in time between the "nowcast" and the "forecast" is the <em class=
"option">forecast lead</em>.<br>
<br>
<img src="../../../docs/images/simple_forecast.png" alt=
"simple forecast" width="60%"><br>
<br>
So - that is simple enough if there is only one forecast, but this
is rarely the case. Let's say we have a second forecast. Ideally,
we'd like to verify at exactly the same locations and forecast
leads - otherwise we're not really comparing the same things. If
the second verification network happens to be at locations that are
easy to predict, we're comparing apples and oranges. The <em class=
"italic">fair</em> way to proceed is to determine the verification
network that is the same for all forecasts. This generally results
in a pretty small set of observations - a problem we will deal with
later.<br>
<br>
The diagram below illustrates the logic behind determining the list
of verification times for a pretty common scenario: a 24-hour
forecast with a forecast lead of 6 hours, repeated the next day.
The <em class="option">first_analysis</em> is at VT1 - let's call
it 00Z day 1. We need to have observations available at:<br>
VT1 (00Z day1), VT2 (06Z day1), VT3 (12Z day1), VT4 (18Z day1), and
VT5 (24Z day1 / 00Z day2). The <em class=
"option">last_analysis</em> starts at VT5 00Z day 2 and must verify
at<br>
VT5 (00Z day2), VT6 (06Z day2), VT7 (12Z day2), VT8 (18Z day2), and
VT9 (24Z day2 / 00Z day3).<br>
<br>
<a name="coverage" id="coverage"></a> <img src=
"../../../docs/images/obs_seq_coverage_diagram.png" alt=
"coverage timetable" width="100%"><br>
<br>
Note that, if you wanted to, you could launch forecasts at VT2,
VT3, and VT4 without adding extra constraints on the verification
network. <em class="program">obs_seq_coverage</em> simply provides
these possible forecasts "for free", there is no assumption about
<strong>needing</strong> them. We will use the variable <em class=
"option">verification_times</em> to describe the complete set of
times for all possible forecasts. In our example above, there are 5
possible forecasts, each forecast consisting of 5 verification
times (the analysis time and the 4 forecast lead times). As such,
there are 9 unique verification times.<br>
<br>
Note that no attempt is made at checking the QC value of the
candidate observations. One of the common problems is that the
region definition does not mesh particularly well with the model
domain and the DART forward operator fails because it would have to
extrapolate (which is not allowed). Without checking the QC value,
this can mean there are a lot of 'false positives'; observations
that seemingly could be used to validate, but are actually just
outside the model domain. I'm working on that ....<br>
<br>
The <a href="#Usage">USAGE</a> section has more on the actual use
of <em class="program">obs_seq_coverage</em>.</p>
<!--==================================================================-->
<!--=================== DESCRIPTION OF A NAMELIST  ===================-->
<!--==================================================================-->
<a name="Namelist" id="Namelist"></a>
<hr>
<h2>NAMELIST</h2>
<p>This namelist is read from the file <em class=
"file">input.nml</em>. Namelists start with an ampersand '&amp;'
and terminate with a slash '/'. Character strings that contain a
'/' must be enclosed in quotes to prevent them from prematurely
terminating the namelist.</p>
<div class="namelist">
<pre>
&amp;obs_seq_coverage_nml
   obs_sequences     = ''
   obs_sequence_list = ''
   obs_of_interest   = ''
   textfile_out      = 'obsdef_mask.txt'
   netcdf_out        = 'obsdef_mask.nc'
   calendar          = 'Gregorian'
   first_analysis    =  2003, 1, 1, 0, 0, 0
   last_analysis     =  2003, 1, 2, 0, 0, 0
   forecast_length_days          = 1
   forecast_length_seconds       = 0
   verification_interval_seconds = 21600
   temporal_coverage_percent     = 100.0
   lonlim1                       =  -888888.0
   lonlim2                       =  -888888.0
   latlim1                       =  -888888.0
   latlim2                       =  -888888.0
   verbose                       = .false.
   debug                         = .false.
  /
</pre></div>
<br>
<br>
<p>Note that -888888.0 is not a useful number. To use the defaults
delete these lines from the namelist, or set them to 0.0, 360.0 and
-90.0, 90.0.</p>
<p>The date-time integer arrays in this namelist have the form
(YYYY, MM, DD, HR, MIN, SEC).</p>
<p>The allowable ranges for the region boundaries are: latitude
[-90.,90], longitude [0.,Inf.]</p>
<p>You can specify <strong>either</strong> <em class=
"option">obs_sequences</em> <strong>or</strong> <em class=
"option">obs_sequence_list</em> -- not both. One of them has to be
an empty string ... i.e. <em class="option">''</em>.</p>
<div>
<table border="0" cellpadding="10" width="100%" summary=
'namelist description'>
<thead align="left">
<tr>
<th>Item</th>
<th>Type</th>
<th>Description</th>
</tr>
</thead>
<tbody valign="top">
<tr>
<td>obs_sequences</td>
<td>character(len=256)</td>
<td>Name of the observation sequence file(s).<br>
This may be a relative or absolute filename. If the filename
contains a '/', the filename is considered to be comprised of
everything to the right, and a directory structure to the left. The
directory structure is then queried to see if it can be incremented
to handle a sequence of observation files. The default behavior of
<em class="program">obs_seq_coverage</em> is to look for additional
files to include until the files are exhausted or an <em class=
"file">obs_seq.final</em> file is found that contains observations
beyond the timeframe of interest.<br>
e.g. 'obsdir_001/obs_seq.final' will cause <em class=
"program">obs_seq_coverage</em> to look for
'obsdir_002/obs_seq.final', and so on.<br>
If this is set, <em class="option">obs_sequence_list</em> must be
set to ' '.</td>
</tr>
<tr>
<td>obs_sequence_list</td>
<td>character(len=256)</td>
<td>Name of an ascii text file which contains a list of one or more
observation sequence files, one per line. If this is specified,
<em class="option">obs_sequences</em> must be set to ' '. Can
be created by any method, including sending the output of the 'ls'
command to a file, a text editor, or another program.</td>
</tr>
<tr>
<td>obs_of_interest</td>
<td>character(len=32), dimension(:)</td>
<td>These are the observation types that will be verified. It is an
array of character strings that must match the standard DART
observation types. Simply add as many or as few observation types
as you need. Could be 'METAR_U_10_METER_WIND',
'METAR_V_10_METER_WIND',..., for example.</td>
</tr>
<tr>
<td>textfile_out</td>
<td>character(len=256)</td>
<td>The name of the file that will contain the observation
definitions of the verfication observations. Only the metadata from
the observations (location, time, obs_type) are preserved in this
file. They are in no particular order. <a href=
"../../../assimilation_code/programs/obs_selection/obs_selection.html">
obs_selection</a> will use this file as a 'mask' to extract the
real observations from the candidate observation sequence
files.</td>
</tr>
<tr>
<td>netcdf_out</td>
<td>character(len=256)</td>
<td>The name of the file that will contain the observation
definitions of the unique locations that match <strong>any</strong>
of the verification times. This file is used in conjunction with
<a href=
"../../../assimilation_code/programs/obs_seq_verify/obs_seq_verify.html">
obs_seq_verify</a> to reorder the <em class=
"file">obs_seq.forecast</em> into a structure that will facilitate
calculating the statistics and scores of the forecasts.</td>
</tr>
<tr>
<td>calendar</td>
<td>character(len=129)</td>
<td>The type of the calendar used to interpret the dates.</td>
</tr>
<tr>
<td>first_analysis</td>
<td>integer, dimension(6)</td>
<td>The start time of the first forecast. Also known as the
analysis time of the first forecast. The six integers are: year,
month, day, hour, hour, minute, second -- in that order.</td>
</tr>
<tr>
<td>last_analysis</td>
<td>integer, dimension(6)</td>
<td>The start time of the last forecast. The six integers are:
year, month, day, hour, hour, minute, second -- in that order. This
needs to be a perfect multiple of the <em class=
"option">verification_interval_seconds</em> from the start of
<em class="option">first_analysis</em>.</td>
</tr>
<tr>
<td>forecast_length_days<br>
forecast_length_seconds</td>
<td>integer</td>
<td>both values are used to determine the <strong>total</strong>
length of any single forecast.</td>
</tr>
<tr>
<td>verification_interval_seconds</td>
<td>integer</td>
<td>The number of seconds between each verification.
<ul style="list-style: none;">
<li> 1 h == 3600s</li>
<li> 2 h == 7120s</li>
<li> 3 h == 10800s</li>
<li> 6 h == 21600s</li>
<li>12 h == 43200s</li>
</ul>
</td>
</tr>
<tr>
<td>temporal_coverage_percent</td>
<td>real</td>
<td>While it is possible to specify that you do not need an
observation at <strong>every</strong> time, it makes the most
sense. This is not actually <strong>required</strong> to be 100%
but 100% results in the most robust comparison.</td>
</tr>
<tr>
<td>lonlim1</td>
<td>real</td>
<td>Westernmost longitude of desired region.</td>
</tr>
<tr>
<td>lonlim2</td>
<td>real</td>
<td>Easternmost longitude of desired region. <em class="new">If
this value is <b>less than</b> the westernmost value, it defines a
region that spans the prime meridian.</em> It is perfectly
acceptable to specify lonlim1 = 330 , lonlim2 = 50 to identify a
region like "Africa".</td>
</tr>
<tr>
<td>latlim1</td>
<td>real</td>
<td>Southernmost latitude of desired region.</td>
</tr>
<tr>
<td>latlim2</td>
<td>real</td>
<td>Northernmost latitude of desired region.</td>
</tr>
<tr>
<td>verbose</td>
<td>logical</td>
<td>Print extra run-time information.</td>
</tr>
<tr>
<td>debug</td>
<td>logical</td>
<td>Enable debugging messages. May generate a lot of output.</td>
</tr>
</tbody>
</table>
</div>
<br>
<br>
<p>For example:</p>
<pre>
&amp;obs_seq_coverage_nml
   obs_sequences     = ''
   obs_sequence_list = 'obs_coverage_list.txt'
   obs_of_interest   = 'METAR_U_10_METER_WIND',
                       'METAR_V_10_METER_WIND'
   textfile_out      = 'obsdef_mask.txt'
   netcdf_out        = 'obsdef_mask.nc'
   calendar          = 'Gregorian'
   first_analysis    =  2003, 1, 1, 0, 0, 0
   last_analysis     =  2003, 1, 2, 0, 0, 0
   forecast_length_days          = 1
   forecast_length_seconds       = 0
   verification_interval_seconds = 21600
   temporal_coverage_percent     = 100.0
   lonlim1    =    0.0
   lonlim2    =  360.0
   latlim1    =  -90.0
   latlim2    =   90.0
   verbose    = .false.
   /
</pre>
<br>
<!--==================================================================-->
 <a name="Modules" id="Modules"></a>
<hr>
<h2>OTHER MODULES USED</h2>
<pre>
assim_model_mod
types_mod
location_mod
model_mod
null_mpi_utilities_mod
obs_def_mod
obs_kind_mod
obs_sequence_mod
random_seq_mod
time_manager_mod
utilities_mod
</pre>
<!--==================================================================-->
<!-- Describe the Files Used by this module.                          -->
<!--==================================================================-->
<a name="FilesUsed" id="FilesUsed"></a>
<hr>
<h2>FILES</h2>
<ul>
<li><em class="file">input.nml</em> is used for <em class=
"option">obs_seq_coverage_nml</em><br></li>
<li>A text file containing the metadata for the observations to be
used for forecast evaluation is created. This file is subsequently
required by <a href=
"../../../assimilation_code/programs/obs_selection/obs_selection.html">
obs_selection</a> to subset the set of input observation sequence
files into a single observation sequence file (<em class=
"file">obs_seq.evaluate</em>) for the forecast step.<br>
(<em class="file">obsdef_mask.txt</em> is the default
name)<br></li>
<li>A netCDF file containing the metadata for a much larger set of
observations that may be used is created. This file is subsequently
required by <a href=
"../../../assimilation_code/programs/obs_seq_coverage/obs_seq_coverage.html">
obs_seq_coverage</a> to define the desired times and locations for
the verification.<br>
(<em class="file">obsdef_mask.nc</em> is the default name)</li>
</ul>
<!--==================================================================-->
<!-- Discuss  typical usage of obs_seq_coverage.                      -->
<!--==================================================================-->
<a name="Usage" id="Usage"></a>
<hr>
<h2>USAGE</h2>
<p><em class="program">obs_seq_coverage</em> is built in
.../DART/models/<em class="italic">your_model</em>/work, in the
same way as the other DART components.<br>
<br>
There is no requirement on the reporting time/frequence of the
candidate voxels. Once the verification times have been defined,
the observation <strong>closest in time</strong> to the
verification time is selected, the others are ignored. Only
observations within half the verification interval are eligible to
be considered "close".<br>
<br>
<strong>A word about vertical levels.</strong> If the desired
observation type has UNDEFINED or SURFACE for the vertical
coordinate system, there is no concern about trying to match the
vertical. If the desired observation types use PRESSURE; the
following 14 levels are used as the standard levels: 1000, 925,
850, 700, 500, 400, 300, 250, 200, 150, 100, 70, 50, 10 (all hPa).
<strong>No other vertical coordinate system is
supported.</strong></p>
<h3 class="indent1">Example: a single 48-hour forecast that is
evaluated every 6 hours.</h3>
<p><img src="../../../docs/images/verification_48hrX6hr.png" alt=
"Example 1" width="75%"><br>
<br>
In this example, we are generating an <em class=
"file">obsdef_mask.txt</em> file for a single forecast. All the
required input observation sequence filenames will be contained in
a file referenced by the <em class="option">obs_sequence_list</em>
variable. We'll also restrict the observations to a specific
rectangular (in Lat/Lon) region at a particular level. It is
convenient to turn on the verbose option the first time to get a
feel for the logic. Here are the namelist settings if you want to
verify the METAR_U_10_METER_WIND and METAR_V_10_METER_WIND
observations over the entire globe every 6 hours for 2 days
starting 18Z 8 Jun 2008:</p>
<div class="routine">
<pre>
&amp;obs_seq_coverage_nml
   obs_sequences      = ''
   obs_sequence_list  = <em class="input">'obs_file_list.txt'</em>
   obs_of_interest    = <em class=
"input">'METAR_U_10_METER_WIND'</em>,
                        <em class=
"input">'METAR_V_10_METER_WIND'</em>
   textfile_out       = 'obsdef_mask.txt'
   netcdf_out         = 'obsdef_mask.nc'
   calendar           = 'Gregorian'
   first_analysis     = <em class=
"input"> 2008, 6, 8, 18, 0, 0 </em>
   last_analysis      = <em class=
"input"> 2008, 6, 8, 18, 0, 0 </em>
   forecast_length_days          = <em class="input">2</em>
   forecast_length_seconds       = <em class="input">0</em>
   verification_interval_seconds = <em class="input">21600</em>
   temporal_coverage_percent     = 100.0
   lonlim1            = <em class="input">   0.0</em>
   lonlim2            = <em class="input"> 360.0</em>
   latlim1            = <em class="input"> -90.0</em>
   latlim2            = <em class="input">  90.0</em>
   verbose            = <em class="input">.true.</em>
   /
</pre></div>
<p>The first step is to create a file containing the list of
observation sequence files you want to use. This can be done with
the unix command 'ls' with the -1 option (that's a number one) to
put one file per line, particularly if the files are organized in a
nice fashion. If your observation sequence are organized like
this:</p>
<pre>
/Exp1/Dir20080101/obs_seq.final
/Exp1/Dir20080102/obs_seq.final
/Exp1/Dir20080103/obs_seq.final
...
/Exp1/Dir20081231/obs_seq.final
</pre>
<p>then</p>
<div class="unix">ls -1 /Exp1/Dir*/obs_seq.final &gt;
obs_file_list.txt</div>
<p>creates the desired file. Then, simply run <em class=
"program">obs_seq_coverage</em> - you may want to save the run-time
output to a file. It is convenient to turn on the verbose option
the first time. Here is a portion of the run-time output:</p>
<div class="unix">
<pre>
[thoar@mirage2 work]$ ./obs_seq_coverage | &amp; tee my.log
 Starting program obs_seq_coverage
 Initializing the utilities module.
 Trying to log to unit           10
 Trying to open file dart_log.out
 
 --------------------------------------
 Starting ... at YYYY MM DD HH MM SS = 
                 2011  2 22 13 15  2
 Program obs_seq_coverage
 --------------------------------------
 
 set_nml_output Echo NML values to log file only
 Trying to open namelist log dart_log.nml
 location_mod: Ignoring vertical when computing distances; horizontal only
 ------------------------------------------------------
 
 
 -------------- ASSIMILATE_THESE_OBS_TYPES --------------
 RADIOSONDE_TEMPERATURE
 RADIOSONDE_U_WIND_COMPONENT
 RADIOSONDE_V_WIND_COMPONENT
 SAT_U_WIND_COMPONENT
 SAT_V_WIND_COMPONENT
 -------------- EVALUATE_THESE_OBS_TYPES --------------
 RADIOSONDE_SPECIFIC_HUMIDITY
 ------------------------------------------------------
 
 METAR_U_10_METER_WIND is type           36
 METAR_V_10_METER_WIND is type           37
 
 There are            9  verification times per forecast.
 There are            1  supported forecasts.
 There are            9  total times we need observations.
 
 At least           9  observations times are required at:
 verification #            1 at 2008 Jun 08 18:00:00
 verification #            2 at 2008 Jun 09 00:00:00
 verification #            3 at 2008 Jun 09 06:00:00
 verification #            4 at 2008 Jun 09 12:00:00
 verification #            5 at 2008 Jun 09 18:00:00
 verification #            6 at 2008 Jun 10 00:00:00
 verification #            7 at 2008 Jun 10 06:00:00
 verification #            8 at 2008 Jun 10 12:00:00
 verification #            9 at 2008 Jun 10 18:00:00
 
 obs_seq_coverage  opening obs_seq.final.2008060818
 QC index           1  NCEP QC index
 QC index           2  DART quality control
 
First observation time day=148812, sec=64380
First observation date 2008 Jun 08 17:53:00
 Processing obs        10000  of        84691
 Processing obs        20000  of        84691
 Processing obs        30000  of        84691
 Processing obs        40000  of        84691
 Processing obs        50000  of        84691
 Processing obs        60000  of        84691
 Processing obs        70000  of        84691
 Processing obs        80000  of        84691
 obs_seq_coverage  doneDONEdoneDONE does not exist. Finishing up.
 
 There were          442  voxels matching the input criterion.
...
</pre></div>
<h4 class="indent1">Discussion</h4>
<p>Note that the values of <em class=
"code">ASSIMILATE_THESE_OBS_TYPES</em> and <em class=
"code">EVALUATE_THESE_OBS_TYPES</em> are completely irrelevant -
since we're not actually doing an assimilation. The
<strong>BIG</strong> difference between the two output files is
that <em class="file">obsdef_mask.txt</em> contains the metadata
for just the matching observations while <em class=
"file">obsdef_mask.nc</em> contains the metadata for all candidate
locations as well as a lot of information about the desired
verification times. It is possible to explore <em class=
"file">obsdef_mask.nc</em> to review the selection criteria to
include observations/"voxels" that do not perfectly match the
original selection criteria.<br>
<br>
Now that you have the <em class="file">obsdef_mask.nc</em>, you can
explore it with <a href=
"http://www.unidata.ucar.edu/software/netcdf/old_docs/docs_4_1/netcdf/ncdump.html">
ncdump</a>.</p>
<div class="unix">
<pre>
netcdf obsdef_mask {
dimensions:
        voxel = UNLIMITED ; // (512 currently)
        time = 9 ;
        analysisT = 1 ;
        forecast_lead = 9 ;
        nlevels = 14 ;
        linelen = 256 ;
        nlines = 446 ;
        stringlength = 32 ;
        location = 3 ;
variables:
        int voxel(voxel) ;
                voxel:long_name = "desired voxel flag" ;
                voxel:description = "1 == good voxel" ;
        double time(time) ;
                time:long_name = "verification time" ;
                time:units = "days since 1601-1-1" ;
                time:calendar = "GREGORIAN" ;
        double analysisT(analysisT) ;
                analysisT:long_name = "analysis (start) time of each forecast" ;
                analysisT:units = "days since 1601-1-1" ;
                analysisT:calendar = "GREGORIAN" ;
        int forecast_lead(forecast_lead) ;
                forecast_lead:long_name = "current forecast length" ;
                forecast_lead:units = "seconds" ;
        double verification_times(analysisT, forecast_lead) ;
                verification_times:long_name = "verification times during each forecast run" ;
                verification_times:units = "days since 1601-1-1" ;
                verification_times:calendar = "GREGORIAN" ;
                verification_times:rows = "each forecast" ;
                verification_times:cols = "each verification time" ;
        float mandatory_level(nlevels) ;
                mandatory_level:long_name = "mandatory pressure levels" ;
                mandatory_level:units = "Pa" ;
        char namelist(nlines, linelen) ;
                namelist:long_name = "input.nml contents" ;
        char obs_type(voxel, stringlength) ;
                obs_type:long_name = "observation type string at this voxel" ;
        double location(voxel, location) ;
                location:description = "location coordinates" ;
                location:location_type = "loc3Dsphere" ;
                location:long_name = "threed sphere locations: lon, lat, vertical" ;
                location:storage_order = "Lon Lat Vertical" ;
                location:units = "degrees degrees which_vert" ;
        int which_vert(voxel) ;
                which_vert:long_name = "vertical coordinate system code" ;
                which_vert:VERTISUNDEF = -2 ;
                which_vert:VERTISSURFACE = -1 ;
                which_vert:VERTISLEVEL = 1 ;
                which_vert:VERTISPRESSURE = 2 ;
                which_vert:VERTISHEIGHT = 3 ;
                which_vert:VERTISSCALEHEIGHT = 4 ;
        int ntimes(voxel) ;
                ntimes:long_name = "number of observation times at this voxel" ;
        double first_time(voxel) ;
                first_time:long_name = "first valid observation time at this voxel" ;
                first_time:units = "days since 1601-1-1" ;
                first_time:calendar = "GREGORIAN" ;
        double last_time(voxel) ;
                last_time:long_name = "last valid observation time at this voxel" ;
                last_time:units = "days since 1601-1-1" ;
                last_time:calendar = "GREGORIAN" ;
        double ReportTime(voxel, time) ;
                ReportTime:long_name = "time of observation" ;
                ReportTime:units = "days since 1601-1-1" ;
                ReportTime:calendar = "GREGORIAN" ;
                ReportTime:missing_value = 0. ;
                ReportTime:_FillValue = 0. ;

// global attributes:
                :creation_date = "YYYY MM DD HH MM SS = 2011 03 01 09 28 40" ;
                :obs_seq_coverage_source = "$URL$" ;
                :obs_seq_coverage_revision = "$Revision$" ;
                :obs_seq_coverage_revdate = "$Date$" ;
                :min_steps_required = 9 ;
                :forecast_length_days = 2 ;
                :forecast_length_seconds = 0 ;
                :verification_interval_seconds = 21600 ;
                :obs_of_interest_001 = "METAR_U_10_METER_WIND" ;
                :obs_of_interest_002 = "METAR_V_10_METER_WIND" ;
                :obs_seq_file_001 = "obs_seq.final.2008060818" ;
data:

 time = 148812.75, 148813, 148813.25, 148813.5, 148813.75, 148814, 148814.25, 
    148814.5, 148814.75 ;

 forecast_lead = 0, 21600, 43200, 64800, 86400, 108000, 129600, 151200, 172800 ;
}
</pre></div>
<p>The first thing to note is that there are more voxels (512) than
reported during the run-time output (442). Typically, there will be
many more voxels in the netCDF file than will meet the selection
criteria - but this is just an example. Some of the voxels in the
netCDF file do not meet the selection criteria - meaning they do
not have observations at all 9 required times. Furthermore, there
are 512 locations for ALL of the desired observation types. In
keeping with the DART philosophy of scalar observations, each
observation type gets a separate voxel. There are
<strong>not</strong> 512 METAR_U_10_METER_WIND observations and 512
METAR_V_10_METER_WIND observations. There are N
METAR_U_10_METER_WIND observations and M METAR_V_10_METER_WIND
observations where N+M = 512. And only 442 of them have
observations at all the times required for the verification. Dump
the <em class="option">obs_type</em> variable to see what voxel has
what observation type.<br>
<br>
The <em class="option">voxel</em> variable is fundamentally a flag
that indicates if the station has all of the desired verification
times. Combine that information with the <em class=
"option">obs_type</em> and <em class="option">location</em> to
determine where your verifications of any particular observation
type will take place.<br>
<br>
Now that you have the <em class="file">obsdef_mask.txt</em>, you
can run <a href=
"../../../assimilation_code/programs/obs_selection/obs_selection.html">
obs_selection</a> to subset the observation sequence files into one
compact file to use in your ensemble forecast.</p>
<!--==================================================================-->
<!-- Cite references, if need be.                                     -->
<!--==================================================================-->
<a name="References" id="References"></a>
<hr>
<h2>REFERENCES</h2>
<ul>
<li>none - but this seems like a good place to start: <a href=
"http://www.cawcr.gov.au/projects/verification/">The Centre for
Australian Weather and Climate Research - Forecast Verification
Issues, Methods and FAQ</a></li>
</ul>
<!--==================================================================-->
<!-- Describe all the error conditions and codes.                     -->
<!--==================================================================-->
<a name="Errors" id="Errors"></a>
<hr>
<h2>ERROR CODES and CONDITIONS</h2>
<div class="errors">
<table border="1" cellspacing="1" cellpadding="10" width="100%">
<tr>
<th>Routine</th>
<th>Message</th>
<th>Comment</th>
</tr>
<tr><!-- routine -->
<td valign="top">obs_seq_coverage</td>
<!-- message -->
<td valign="top">'namelist: temporal_coverage_percent (xxxx) must
be == 100.0 for now.)'</td>
<!-- comment -->
<td valign="top">it is required that ALL verification times be
present for all forecasts</td>
</tr>
<tr><!-- routine -->
<td valign="top">obs_seq_coverage</td>
<!-- message -->
<td valign="top">'specify "obs_sequences" or
"obs_sequence_list"'</td>
<!-- comment -->
<td valign="top">one of these namelist variables MUST be an empty
string</td>
</tr>
<tr><!-- routine -->
<td valign="top">obs_seq_coverage</td>
<!-- message -->
<td valign="top">'xxxxxx ' is not a known observation type.'</td>
<!-- comment -->
<td valign="top">one of the <em class="option">obs_of_interest</em>
namelist entries specifies an observation type that is not
supported. Perhaps you need to rerun <em class=
"program">preprocess</em> with support for the observation, or
perhaps it is spelled incorrectly. All DART observation types are
strictly uppercase.</td>
</tr>
<tr><!-- routine -->
<td valign="top">obs_seq_coverage</td>
<!-- message -->
<td valign="top">'need at least 1 qc and 1 observation copy'</td>
<!-- comment -->
<td valign="top">an observation sequence does not have all the
metadata necessary. Cannot use "<em class=
"file">obs_seq.in</em>"-class sequences.</td>
</tr>
<tr><!-- routine -->
<td valign="top">obs_seq_coverage</td>
<!-- message -->
<td valign="top">'num_copies ##### does not match #####'</td>
<!-- comment -->
<td valign="top">ALL observation sequences must contain the same
'copy' information. At some point it may be possible to mix
"<em class="file">obs_seq.out</em>"-class sequences with
"<em class="file">obs_seq.final</em>"-class sequences, but this
seems like it can wait.</td>
</tr>
<tr><!-- routine -->
<td valign="top">obs_seq_coverage</td>
<!-- message -->
<td valign="top">'No location had at least ### reporting
times.'</td>
<!-- comment -->
<td valign="top">The input selection criteria did not result in any
locations that had observations at all of the required verification
times.</td>
</tr>
<tr><!-- routine -->
<td valign="top">set_required_times</td>
<!-- message -->
<td valign="top">'namelist: forecast length is not a multiple of
the verification interval'</td>
<!-- comment -->
<td valign="top">The namelist settings for <em class=
"option">forecast_length_[days,seconds]</em> and <em class=
"option">verification_interval_seconds</em> do not make sense.
Refer to <a href="#coverage">the forecast time diagram</a>.</td>
</tr>
<tr><!-- routine -->
<td valign="top">set_required_times</td>
<!-- message -->
<td valign="top">'namelist: last analysis time is not a multiple of
the verification interval'</td>
<!-- comment -->
<td valign="top">The namelist settings for <em class=
"option">first_analysis</em> and <em class=
"option">last_analysis</em> are not separated by a multiple of
<em class="option">verification_interval_seconds</em>. Refer to
<a href="#coverage">the forecast time diagram</a>.</td>
</tr>
</table>
</div>
<h2>KNOWN BUGS</h2>
<p>none</p>
<!--==================================================================-->
<!-- Describe Future Plans.                                           -->
<!--==================================================================-->
<a name="FuturePlans" id="FuturePlans"></a>
<hr>
<h2>FUTURE PLANS</h2>
<p><em class="removed">Relax the restriction requiring 100.0%
temporal coverage.</em><br>
Sensibly require that we only require against observations that
DART can compute. i.e. the prior forward operator must complete
successfully - almost all cases of the operator failing are
extrapolation issues; the observation is outside the domain.<br>
<br>
Note that no attempt is made at checking the QC value of the
candidate observations. One of the common problems is that the
region definition does not mesh particularly well with the model
domain and the DART forward operator fails because it would have to
extrapolate (which is not allowed). Without checking the QC value,
this can mean there are a lot of 'false positives'; observations
that seemingly could be used to validate, but are actually just
outside the model domain. I'm working on that ....</p>
<!--==================================================================-->
<!-- Legalese & Metadata                                              -->
<!--==================================================================-->
<a name="Legalese" id="Legalese"></a>
<hr>
<h2>Terms of Use</h2>
<p>DART software - Copyright UCAR. This open source software is
provided by UCAR, "as is", without charge, subject to all terms of
use at <a href=
"http://www.image.ucar.edu/DAReS/DART/DART_download">http://www.image.ucar.edu/DAReS/DART/DART_download</a></p>
<!--==================================================================-->
</body>
</html>
