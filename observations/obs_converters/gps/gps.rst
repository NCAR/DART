<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN"
          "http://www.w3.org/TR/html4/strict.dtd">
<HTML>
<HEAD>
<TITLE>GPS Observations</TITLE>
<link rel="stylesheet" type="text/css" href="../../../docs/html/doc.css" />
<link href="../../../docs/images/dart.ico" rel="shortcut icon" />
</HEAD>
<BODY>
<A NAME="TOP"></A>

<H1>GPS Observations</H1>

<table border=0 summary="" cellpadding=5>
<tr>
    <td valign=middle>
    <img src="../../../docs/images/Dartboard7.png" alt="DART project logo" height=70 />
    </td>
    <td>Jump to <a href="../../../docs/index.html">DART Documentation Main Index</a></td>
</tr>
</table>

<A HREF="#DataSources">DATA SOURCES</A> /
<A HREF="#Programs">PROGRAMS</A> / 
<A HREF="#Namelist">NAMELIST</A> /
<A HREF="#Workflow">WORKFLOW</A> / 
<A HREF="#Modules">MODULES</A> /
<A HREF="#Errors">ERRORS</A> /
<A HREF="#FuturePlans">FUTURE PLANS</A> /
<A HREF="#Legalese">TERMS OF USE</A>

<H2>Overview</H2>

<P>
The <a href="http://www.cosmic.ucar.edu">COSMIC</a> project provides data from 
a series of satellites.  There are two forms of the data that are used by 
DART: GPS Radio Occultation data and Electron Density.  The programs in this 
directory extract the data from the distribution files and put them into 
DART observation sequence (obs_seq) file format.
</P>
<H3>Radio Occultation</H3>
<P>
The COSMIC satellites measure the phase delay caused by deviation of the
straight-line path of the GPS satellite signal as it passes through the 
Earth's atmosphere when the GPS and COSMIC
satellites rise and set relative to each other.
This deviation results from changes in the angle of refraction of light
as it passes through regions of varying density of atmosphere.
These changes are a result of variations in the temperature,
pressure, and moisture content.
Vertical profiles of temperature and moisture can be derived as 
the signal passes through more and more atmosphere 
until it is obscured by the earth's horizon.
There are thousands of observations each day
distributed around the globe, including in areas which
previously were poorly observed.
These data are converted with the <em class=code>convert_cosmic_gps_cdf.f90</em>
program and create DART observations of GPSRO_REFRACTIVITY.
</P>

<H3>Electron Density</H3>
<P>
The COSMIC satellites also provide ionospheric profiles of electron density.
The accuracy is generally about 10<sup>-4</sup> &sim; 10<sup>-5</sup> cm<sup>-3</sup>.
These data are converted with the <em class=code>convert_cosmic_ionosphere.f90</em>
program and create DART observations tagged as COSMIC_ELECTRON_DENSITY.
</P>

<!--==================================================================-->

<A NAME="DataSources"></A>
<div class="top">[<a href="#">top</a>]</div><hr />
<H2>DATA SOURCES</H2>

<P>
Data from the
<a href="http://www.cosmic.ucar.edu">COSMIC</a> Program
are available by signing up on the
<a href="http://cosmic-io.cosmic.ucar.edu/cdaac">data access</a>
web page.  We prefer delivery in
<a href="http://www.unidata.ucar.edu/software/netcdf">netCDF</a>
file format.
<P>

<H3>Radio Occultation</H3>
<P>
The files we use as input to these conversion programs are
the Level 2 data, Atmospheric Profiles (filenames include
the string 'atmPrf').
<br />
<br />
Each vertical profile is stored in a separate netCDF file,
and there are between 1000-3000 profiles/day, so converting a day's
worth of observations used to involve downloading many individual files.
There are now daily tar files available which makes it simpler
to download the raw data all in a single file and then untar it to
get the individual profiles.
<br />
<br />
The scripts in the <em class=file>shell_scripts</em> directory can
now download profiles from any of the available
satellites that return GPS RO data to the CDAAC web site.  See the
<em class=file>gpsro_to_obsseq.csh</em> or
<em class=file>convert_many_gpsro.csh</em> script for where to 
specify the satellites to be included.
</P>

<H3>Electron Density</H3>
<P>
The files we have used as input to these conversion programs are from 
the COSMIC 2013 Mission and have a data type of 'ionPrf'.
<br />
<br />
The file naming convention and file format are described by COSMIC 
<a href=' http://cdaac-www.cosmic.ucar.edu/cdaac/cgi_bin/fileFormats.cgi?type=ionPrf'>here</a> 
and there can be more than 1000 profiles/day. Like the GPS radio
occultation data, the profiles are now available in a single daily
tar file which can be downloaded then be unpacked into the individual files.
COSMIC has instructions on ways to download the data at<br />
<a href='http://cdaac-www.cosmic.ucar.edu/cdaac/tar/rest.html'
>http://cdaac-www.cosmic.ucar.edu/cdaac/tar/rest.html</a>
</P>

<!--==================================================================-->

<A NAME="Programs"></A>
<div class="top">[<a href="#">top</a>]</div><hr />
<H2>PROGRAMS</H2>

<H3>convert_cosmic_gps_cdf</H3>
<P>
The data are distributed in
<a href="http://www.unidata.ucar.edu/software/netcdf">netCDF</a>
file format.  DART requires all observations to be in a proprietary
format often called DART "obs_seq" format.
The files in this directory (a combination
of C shell scripts and a Fortran source executable)
do this data conversion.
<br />
<br />
The shell_scripts directory contains several example scripts, including one 
which downloads the raw data files a day at a time (<em class=program>download_script.csh</em>), 
and one which executes the conversion program (<em class=program>convert_script.csh</em>).  
These scripts make 6 hour files by default, but have options for other times.  
Each profile is stored in a separate netcdf file and there are usually
between 1000-3000 files/day, so the download process can be
lengthy.  You probably want to download as a separate preprocess step
and do not use the script options to automatically delete the input files.
Keep the files around until you are sure you are satisified with the
output files and then delete them by hand.
<br />
<br />
The conversion executable <em class=program>convert_cosmic_gps_cdf</em>, 
reads the namelist
<A HREF="#Namelist"> <em class=code>&amp;convert_cosmic_gps_nml</em> </A>
from the file <em class=file>input.nml</em>.
<br />
<br />
The namelist lets you select from one of two different forward
operators.  The 'local' forward operator computes the expected
observation value at a single point: the requested height at the
tangent point of the ray between satellites.  The 'non-local' operator
computes values along the ray-path and does an integration to get
the expected value.  The length of the integration segments and
height at which to end the integration are given in the namelist.
In some experiments the difference between the two types of
operators was negligible.  This choice is made at the time of
the conversion, and the type of operator is stored in the observation,
so at runtime the corresponding forward operator will be used to
compute the expected observation value.
<br />
<br />
The namelist also lets you specify at what heights you want
observations to be extracted.  The raw data is very dense in
the vertical; using all values would not results in a set of
independent observations.  The current source code no longer
does an intermediate interpolation; the original profiles
appear to be smooth enough that this is not needed.  The requested
vertical output heights are interpolated directly from the full profile.
</P>

<H3>convert_cosmic_ionosphere</H3>
<P>
Each profile is interpolated to a set of desired levels that are specified at run time.
During the conversion process, each profile is checked for negative values of 
electron density above the minimum desired level. If negative values are found,
the entire profile is discarded. If an observation sequence file already exists,
the converter will simply add the new observations to it. Multiple profiles may
be converted in a single execution, so it is easy to consolidate all the profiles
for a single day into a single observation sequence file, for example.
<br />
<br />
<em class=program>convert_cosmic_ionosphere</em>
reads the namelist
<A HREF="#Namelist"> <em class=code>&amp;convert_cosmic_ionosphere_nml</em> </A>
from the file <em class=file>input.nml</em>.
<br />
<br />
The original observation times are preserved in the conversion process. If it is
desired to subset the observation sequence file such that observations too far away
from desired assimilation times are rejected, a separate post-processing step
using the <a href='obs_sequence_tool.html'>obs_sequence_tool</a> is required.
A script will be necessary to take a start date, an end date, an assimilation
time step, and a desired time 'window' - and strip out the unwanted observations from
a series of observation sequence files.
<br />
<br />
There are multiple ways of specifying the observation error variance at run 
time. They are implemented in a routine named 
<em class=code>electron_density_error()</em> and are selected by the namelist
variable <em class=code>observation_error_method</em>.
<P>
<TABLE border=0 cellpadding=1 width=100% summary='namelist description'>
<tr><td>'constant'<td>
    <td>a scalar value for all observations</td>
</tr>
<tr><td>'scaled'<td>
    <td>the electron density is multiplied by a scalar value</td>
</tr>
<tr><td>'lookup'<td>
    <td>a lookup table is read</td>
</tr>
<tr><td>'scaled_lookup'<td>
    <td>the lookup table value is multiplied by a scalar value and the electron density value</td>
</tr>
</TABLE>

<blockquote>
I-Te Lee: " ... the original idea for error of ionospheric observation is 1%. Thus, I 
put the code as "oerr = 0.01_r8 * obsval". Liu et. al and Yue et al investigated 
the Abel inversion error of COSMIC ionosphere profile, both of them figure out 
the large error would appear at the lower altitude and push model toward wrong 
direction at the lower ionosphere while assimilating these profiles.
On the other hand, the Abel inversion error depends on the ionospheric electron 
density structure, which is a function of local time, altitude and geomagnetic latitude.
To simplify the procedure to define observation error of profiles, Xinan Yue help 
me to estimate an error matrix and saved in the file which named 'f3coerr.nc'. 
...
The number in the matrix is error percentage (%), which calculated by OSSE.  
Here are two reference papers.  In the end, the observation error consists of 
instrumentation error (10%) and Abel error."

<ul><li>X.&nbsp;Yue, W.S.&nbsp;Schreiner, J.&nbsp;Lei, S.V.&nbsp;Sokolovskiy, C.&nbsp;Rocken, D.C.&nbsp;Hunt, and Y.-H.&nbsp;Kuo (2010),<br />
<a href="https://www.ann-geophys.net/28/217/2010/"
>Error analysis of Abel retrieved electron density profiles from radio occultation measurements.</a><br />
<i>Annales Geophysicae: Atmospheres, Hydrospheres and Space Sciences</i>. <b>28</b> No. 1, pp 217-222,
doi:10.5194/angeo-28-217-2010
</li>

<li> J.Y.&nbsp;Liu, C.Y.&nbsp;Lin, C.H.&nbsp;Lin, H.F.&nbsp;Tsai, S.C.&nbsp;Solomon, Y.Y.&nbsp;Sun, I.T.&nbsp;Lee, W.S.&nbsp;Schreiner, and Y.H.&nbsp;Kuo (2010), <br>
<a href="http://dx.doi.org/10.1029/2009JA015079"
>Artificial plasma cave in the low-latitude ionosphere results from the radio occultation inversion of the FORMOSAT-3/COSMIC}</a>,
<i>Journal of Geophysical Research: Space Physics</i>. <b>115</b> No. A7, pp 2156-2202,
doi:10.1029/2009JA015079
</li>
</ul>
</blockquote>

<P>
It is possible to create observation sequence files for perfect model 
experiments that have realistic observation sampling patterns and observation
error variances that <strong>do not have any actual electron densities</strong>. 
The COSMIC data files are read, but the electron density information is not 
written.  Keep in mind that some methods of specifying the observation error 
variance require knowledge of the observation value. If the observation 
value is bad or the entire profile is bad, no observation locations are 
created for the profile.
<P>

<!--==================================================================-->
<!--=================== DESCRIPTION OF A NAMELIST ====================-->
<!--==================================================================-->

<A NAME="Namelist"></A>
<div class="top">[<a href="#">top</a>]</div><hr />
<H2>NAMELIST</H2>
<P>
This namelist is read from the file <em class=file>input.nml</em>.
Namelists start with an ampersand
'&amp;' and terminate with a slash '/'.
Character strings that contain a '/' must be
enclosed in quotes to prevent them from 
prematurely terminating the namelist.
</P>

<div class=namelist>
<pre>
&amp;convert_cosmic_gps_nml
   obs_levels             = -1.0
   use_original_kuo_error = .false.
   local_operator         = .true.
   ray_ds                 = 5000.0
   ray_htop               = 15000.0
   gpsro_netcdf_file      = 'cosmic_gps_input.nc'
   gpsro_netcdf_filelist  = ''
   gpsro_out_file         = 'obs_seq.gpsro'
 /
</pre>
</div>

<br />
<br />

<div>
<TABLE border=0 cellpadding=10 width=100% summary='namelist description'>
<THEAD align=left>
<TR><TH> Item </TH>
    <TH> Type </TH>
    <TH> Description </TH> </TR>
</THEAD>

<TBODY valign=top>
  
<TR><TD> obs_levels </TD>
    <TD> integer(200) </TD>
    <TD>A series of heights, in kilometers, where observations from this profile
should be interpolated.  (Note that the other distances and heights in the
namelist are specified in meters.)  The values should be listed in increasing
height order.
</TD></TR>
  
<TR><TD> use_original_kuo_error </TD>
    <TD> logical </TD>
   <TD>If .true. use the observation error variances for a refractivity observation
that come from a Kuo paper and were implied to be used for the CONUS domain.
If .false. use observation error variances similar to what is used in GSI.
</TD></TR>
  
<TR><TD> local_operator </TD>
    <TD> logical </TD>
   <TD>If .true. compute the observation using a method which assumes all
effects occur at the tangent point.  If .false. integrate along the tangent
line and do ray-path reconstruction.
</TD></TR>
  
<TR><TD> ray_ds </TD>
    <TD> real(r8) </TD>
    <TD>For the non-local operator only, the delta stepsize, in meters, to use
for the along-path integration in each direction out from the tangent point.
</TD></TR>
  
<TR><TD> ray_htop </TD>
    <TD> real(r8) </TD>
    <TD>For the non-local operator only, stop the integration when one of the
endpoints of the next integration step goes above this height.  Specify in meters.
</TD></TR>
  
<TR><TD> gpsro_netcdf_file </TD>
    <TD> character(len=128) </TD>
    <TD>The input filename when converting a single profile.  Only one of the 
file or filelist items can have a valid value, so to use the single filename
set the list name 'gpsro_netcdf_filelist' to the empty string ('&nbsp;').
</TD></TR>

<TR><TD> gpsro_netcdf_filelist </TD>
    <TD> character(len=128) </TD> 
    <TD>To convert a series of profiles in a single execution create a text file
which contains each input file, in ascii, one filename per line.  Set this
item to the name of that file, and set 'gpsro_netcdf_file' to the empty
string ('&nbsp;').  
</TD></TR>

<TR><TD> gpsro_out_file </TD>
    <TD> character(len=128) </TD> 
    <TD>The output file to be created.  To be compatible with earlier
versions of this program, if this file already exists it will be read in and
the new data will be appended to that file.  
</TD></TR>

</TBODY> 
</TABLE> 

<P>A more useful example follows:</P>
<pre>
&amp;convert_cosmic_gps_nml
  gpsro_netcdf_file      = ''
  gpsro_netcdf_filelist  = 'flist'
  gpsro_out_file         = 'obs_seq.gpsro'
  local_operator         = .true.
  use_original_kuo_error = .false.
  ray_ds                 = 5000.0
  ray_htop               = 13000.1
  obs_levels =        0.2,  0.4,  0.6,  0.8,
                1.0,  1.2,  1.4,  1.6,  1.8,
                2.0,  2.2,  2.4,  2.6,  2.8,
                3.0,  3.2,  3.4,  3.6,  3.8,
                4.0,  4.2,  4.4,  4.6,  4.8,
                5.0,  5.2,  5.4,  5.6,  5.8,
                6.0,  6.2,  6.4,  6.6,  6.8,
                7.0,  7.2,  7.4,  7.6,  7.8,
                8.0,  8.2,  8.4,  8.6,  8.8,
                9.0,  9.2,  9.4,  9.6,  9.8,
               10.0, 10.2, 10.4, 10.6, 10.8,
               11.0, 11.2, 11.4, 11.6, 11.8,
               12.0, 12.2, 12.4, 12.6, 12.8,
               13.0, 13.2, 13.4, 13.6, 13.8,
               14.0, 14.2, 14.4, 14.6, 14.8,
               15.0, 15.2, 15.4, 15.6, 15.8,
               16.0, 16.2, 16.4, 16.6, 16.8,
               17.0, 17.2, 17.4, 17.6, 17.8,
               18.0, 19.0, 20.0, 21.0, 22.0,
               23.0, 24.0, 25.0, 26.0, 27.0,
               28.0, 29.0, 30.0, 31.0, 32.0,
               33.0, 34.0, 35.0, 36.0, 37.0,
               38.0, 39.0, 40.0, 41.0, 42.0,
               43.0, 44.0, 45.0, 46.0, 47.0,
               48.0, 49.0, 50.0, 51.0, 52.0,
               53.0, 54.0, 55.0, 56.0, 57.0,
               58.0, 59.0, 60.0,
 /
</pre>
</div>


<div class=namelist>
<pre>
&amp;convert_cosmic_ionosphere_nml
  input_file               = ''
  input_file_list          = 'input_file_list.txt'
  output_file              = 'obs_seq.out'
  observation_error_file   = 'none'
  observation_error_method = 'scaled_lookup'
  locations_only           = .false.
  obs_error_factor         = 1.0
  verbose                  = 0
  obs_levels               = -1.0
 /
</pre>
</div>

<br />
<br />

<div>
<TABLE border=0 cellpadding=10 width=100% summary='namelist description'>
<THEAD align=left>
<TR><TH> Item </TH>
    <TH> Type </TH>
    <TH> Description </TH> </TR>
</THEAD>

<TBODY valign=top>
  
<TR><TD> input_file </TD>
    <TD> character(len=256) </TD>
    <TD>The input filename when converting a single profile.  Only one of 
        the <em class=code>input_file</em> or 
        <em class=code>input_file_list</em> items can have a valid value,
        so to use a single filename set 
        <em class=code>input_file_list&nbsp;=&nbsp;''</em>
</TD></TR>

<TR><TD> input_file_list </TD>
    <TD> character(len=256) </TD> 
    <TD>To convert a series of profiles in a single execution create a text 
        file which contains one filename per line.
        Set this item to the name of that file, and set 
        <em class=code>input_file&nbsp;=&nbsp;''</em>
</TD></TR>

<TR><TD> output_file </TD>
    <TD> character(len=256) </TD> 
    <TD>The output file to be created. If this file already exists the new 
        data will be added to that file. DART observation sequences are 
        linked lists. When the list is traversed, the observations are in 
        ascending time order. The order they appear in the file is 
        completely irrelevant. 
</TD></TR>

<TR><TD> observation_error_file </TD>
    <TD> character(len=256) </TD> 
    <TD>This specifies a lookup table. The table created by I-Te Lee and 
        Xinan Yue is called <em class=file>f3coerr.nc</em>.
</TD></TR>

<TR><TD> observation_error_method </TD>
    <TD> character(len=128) </TD> 
    <TD> There are multiple ways of specifying the observation error
         variance. This character string allows you to select the method.
         The selection is not case-sensitive.  Allowable values are: 
         'constant', 'scaled', 'lookup', or 'scaled_lookup'.
         Anything else will result in an error. Look in the 
         <em class=code>electron_density_error()</em> routine for specifics.
</TD></TR>

<TR><TD> locations_only </TD>
    <TD> logical </TD> 
    <TD> If <em class=code>locations_only&nbsp;=&nbsp;.true.</em> then the actual
         observation values are not written to the output observation 
         sequence file. This is useful for designing an OSSE that has a
         realistic observation sampling pattern. Keep in mind that some
         methods of specifying the observation error variance require
         knowledge of the observation value. If the observation value
         is bad or the entire profile is bad, this profile is rejected -
         even if <em class=code>locations_only&nbsp;=&nbsp;.true.</em> 
</TD></TR>

<TR><TD> obs_error_factor </TD>
    <TD> real(r8) </TD> 
    <TD> This is the scalar that is used in several of the methods
         specifying the observation error variance.
</TD></TR>

<TR><TD> verbose </TD>
    <TD> integer </TD> 
    <TD> controls the amount of run-time output echoed to the screen.
         0 is nearly silent, higher values write out more. The filenames 
         of the profiles that are skipped are ALWAYS printed.
</TD></TR>

<TR><TD> obs_levels </TD>
    <TD> integer(200) </TD>
    <TD>A series of heights, in kilometers, where observations from this profile
        should be interpolated.  (Note that the other distances and heights 
        in the namelist are specified in meters.)  The values must be listed
        in increasing height order.
</TD></TR>
  
</TBODY> 
</TABLE> 

<P>A more useful example follows:</P>
<pre>
&amp;convert_cosmic_ionosphere_nml
   input_file               = ''
   input_file_list          = 'file_list.txt'
   output_file              = 'obs_seq.out'
   observation_error_file   = 'f3coeff.dat'
   observation_error_method = 'scaled'
   locations_only           = .false.
   obs_error_factor         = 0.01
   verbose                  = 1
   obs_levels = 160.0, 170.0, 180.0, 190.0, 200.0,
                210.0, 220.0, 230.0, 240.0, 250.0,
                260.0, 270.0, 280.0, 290.0, 300.0,
                310.0, 320.0, 330.0, 340.0, 350.0,
                360.0, 370.0, 380.0, 390.0, 400.0,
                410.0, 420.0, 430.0, 440.0, 450.0
  /
</pre>
</div>

<P><!-- useless divider for 'top' --></P>

<!--==================================================================-->

<A NAME="Workflow"></A>
<div class="top">[<a href="#">top</a>]</div><hr />
<H2>WORKFLOW FOR BATCH CONVERSIONS</H2>

<!-- THIS NEEDS MORE WORK.  IT'S BASICALLY THE CONTENTS OF
     THE README DROPPED IN HERE WITH A </pre> PREPENDED.  
     BUT IT'S BETTER THAN NOTHING.  -->

<P>
If you are converting only a day or two of observations
you can download the files by hand and call the converter 
programs from the command line.  However if you are going
convert many days/months/years of data you need an automated
script, possibly submitted to a batch queue on a large machine.
The following instructions describe shell scripts we provide
as a guide in the <em class=file>shell_scripts</em> directory.
You will have to adapt them for your own system
unless you are running on an NCAR superscomputer.
</P>

<P>
<b>Making DART Observations from Radio Occultation atmPrf Profiles:</b>
<br />

<pre>
Description of the scripts provided to process the COSMIC and 
CHAMP GPS radio occultation data.

Summary of workflow:  
1) cd to the ../work directory and run ./quickbuild.csh to compile everything.  
2) Edit ./gpsro_to_obsseq.csh once to set the directory where the DART
    code is installed, and your CDAAC web site user name and password.
3) Edit ./convert_many_gpsro.csh to set the days of data to download/convert/remove.
4) Run ./convert_many_gpsro.csh either on the command line or submit to a batch system.


More details:

1) quickbuild.csh:

Make sure your $DART/mkmf/mkmf.template is one that matches the
platform and compiler for your system.  It should be the same as
how you have it set to build the other DART executables.

Run quickbuild.csh and it should compile all the executables needed
to do the GPS conversion into DART obs_sequence files.


2) gpsro_to_obsseq.csh:

Edit gpsro_to_obsseq.csh once to set the DART_DIR to where you have
downloaded the DART distribution.  (There are a few additional options
in this script, but the distribution version should be good for most users.)
If you are downloading data from the CDAAC web site, set your
web site user name and password.  After this you should be able to 
ignore this script.


3) convert_many_gpsro.csh:

A wrapper script that calls the converter script a day at a time.
Set the days of data you want to download/convert/remove.  See the
comments at the top of this script for the various options to set.  
Rerun this script for all data you need.  This script depends on
the advance_time executable, which should automatically be built
in the ../work directory, but you may have to copy or link to a
version from this dir.  you also need a minimal input.nml here:

&amp;utilities_nml
 /

is all the contents it needs.


It can be risky to use the automatic delete/cleanup option - if there are
any errors in the script or conversion (file system full, bad file format,
etc) and the script doesn't exit, it can delete the input files before 
the conversion has succeeded.  But if you have file quota concerns
this allows you to keep the total disk usage lower.
</pre>
</P>

<br />
<br />

<P>
<b>Making DART Observations from Ionospheric ionPrf Profiles:</b>
<br />

<pre>
0) run quickbuild.csh as described above

1) iono_to_obsseq.csh

set the start and stop days.  downloads from the CDAAC and
untars into 100s of files per day.  runs the converter to
create a single obs_seq.ion.YYYYMMDD file per day.

2) split_obs_seq.csh

split the daily files into X minute/hour files - set the
window times at the top of the file before running.
</pre>
</P>

<br />
<br />

<P>
<b>Notes on already converted observations on the NCAR supercomputers</b>
<br />
<br />

<b>GPS Radio Occultation Data:</b>
<br />
<br />

<pre>
See /glade/p/image/Observations/GPS

These are DART observation sequence files that contain
radio-occultation measurements from the COSMIC
(and other) satellites.  

Uses temperature/moisture bending of the signals as they
pass through the atmosphere between GPS source satellites
and low-earth-orbit receiving satellites to compute the 
delay in the arrival of data. the files also contain the
bending angle data, but we are not using that currently.


the subdirectories include:

local -- original processed files, single obs at nadir
local-cosmic2013 -- reprocessed by CDAAC in 2013
local-test2013 -- 2013 data, denser in vertical, diff errors
local-complete2013 - all satellites available for that time, 
 new errors (from lydia c), 2013 cosmic reprocessed data
nonlocal -- original processed files, ray-path integrated
rawdata -- netcdf data files downloaded from the CDAAC

local: the ob is at a single location (the tangent point
of the ray and earth) and the entire effect is assumed 
to be impacting the state at that point.

non-local: computes the ob value by doing a line integral
along the ray path to accumulate the total effect.

(in our experiments we have compared both and did not see 
a large difference between the two methods, and so have 
mistly used the local version because it's faster to run.)


some directories contain only the gps obs and must be
merged (with the obs_sequence_tool) with the rest of
the conventional obs before assimilation.

some directories contain both the gps-only files and
the obs merged with NCEP and ACARS data.


if a directory exists but is empty, the files are
likely archived on the HPSS.  see the README files
in the next level directory down for more info on
where they might be.

nsc
jan 2016

</pre>

<br />
<br />

<b>Ionosphere Data:</b>
<br />
<br />

<pre>
See /glade/p/image/Observation/ionosphere

These are COSMIC 'ionPrf' ionospheric profile observations.

They are downloaded from the CDAAC website as daily tar files
and unpacked into the 'raw' directory.  They distribute these
observations with one profile per netcdf file.  Each profile has 
data at ~500-1000 different levels.

Our converter has a fixed number of levels in the namelist
and we interpolate between the two closest levels to get the
data for that level.  If you give the converter a list of
input netcdf files it will convert all of them into a 
single output file.

The 'daily' directory is a collection of all the profiles for
that day.

The 'convert' directory has the executables and scripting
for breaking up the daily files into 10 minute files which
are put in the '10min' directory.  Change the 'split_obs_seq.csh'
script to change the width of this window, or the names of
the output files.

The 'verify.csh' script prints out any missing files, which
happens if there are no profiles in the given window.

Our convention is to make a 0 length file for missing intervals
and we expect the filter run script to look at the file size
and loop if there is a file but with no contents.  This will
allow us to distinguish between a time where we haven't converted
the observations and a time where there are no observations.
In that case the script should add time to the next model
advance request and loop to the next interval.

</pre>
</P>

<br />
<br />

<P><!-- useless divider for 'top' --></P>

<!--==================================================================-->

<A NAME="Modules"></A>
<div class="top">[<a href="#">top</a>]</div><hr />
<H2>MODULES USED</H2>
<P>
<em class=program>convert_cosmic_gps_cdf</em> and
<em class=program>convert_cosmic_ionosphere</em> use the same set of modules.
</P>
<PRE>
assimilation_code/location/threed_sphere/location_mod.f90
assimilation_code/modules/assimilation/adaptive_inflate_mod.f90
assimilation_code/modules/assimilation/assim_model_mod.f90
assimilation_code/modules/io/dart_time_io_mod.f90
assimilation_code/modules/io/direct_netcdf_mod.f90
assimilation_code/modules/io/io_filenames_mod.f90
assimilation_code/modules/io/state_structure_mod.f90
assimilation_code/modules/io/state_vector_io_mod.f90
assimilation_code/modules/observations/obs_kind_mod.f90
assimilation_code/modules/observations/obs_sequence_mod.f90
assimilation_code/modules/utilities/distributed_state_mod.f90
assimilation_code/modules/utilities/ensemble_manager_mod.f90
assimilation_code/modules/utilities/netcdf_utilities_mod.f90
assimilation_code/modules/utilities/null_mpi_utilities_mod.f90
assimilation_code/modules/utilities/null_win_mod.f90
assimilation_code/modules/utilities/options_mod.f90
assimilation_code/modules/utilities/random_seq_mod.f90
assimilation_code/modules/utilities/sort_mod.f90
assimilation_code/modules/utilities/time_manager_mod.f90
assimilation_code/modules/utilities/types_mod.f90
assimilation_code/modules/utilities/utilities_mod.f90
models/template/model_mod.f90
models/utilities/default_model_mod.f90
observations/forward_operators/obs_def_mod.f90
observations/forward_operators/obs_def_utilities_mod.f90
observations/obs_converters/utilities/obs_utilities_mod.f90
</PRE>

<!--==================================================================-->

<A NAME="KnownBugs"></A>
<div class="top">[<a href="#">top</a>]</div><hr />
<H2>KNOWN BUGS</H2>
<P>
Some COSMIC files seem to have internal times which differ from the
times encoded in the filenames by as much as 2-3 minutes.  If it is
important to get all the observations within a particular time window
files with filenames from a few minutes before and after the window 
should be converted.
Times really outside the window can be excluded in a separate step using the
<a href="../../../assimilation_code/programs/obs_sequence_tool/obs_sequence_tool.html">obs_sequence_tool</a>.
</P>

<!--==================================================================-->
<!-- Describe Errors.                                                 -->
<!--==================================================================-->

<A NAME="Errors"></A>
<div class="top">[<a href="#">top</a>]</div><hr />
<H2>ERRORS</H2>
<P>
The converters have a parameter declaring the maximum number of desired 
levels as 200. If more than 200 levels are entered as input
(to <em class=code>obs_levels</em>), a rather uninformative run-time 
error is generated:
</P>
<pre>
 ERROR FROM:
  routine: check_namelist_read
  message:  INVALID NAMELIST ENTRY:  / in namelist convert_cosmic_ionosphere_nml
</pre>
<P>
Your error may be different if <em class=code>obs_levels</em> is not the
last namelist item before the slash '/'
</P>

<!--==================================================================-->
<!-- Describe Future Plans.                                           -->
<!--==================================================================-->

<A NAME="FuturePlans"></A>
<div class="top">[<a href="#">top</a>]</div><hr />
<H2>FUTURE PLANS</H2>
<P>
none
</P>

<!--==================================================================-->
<!-- Legalese & Metadata                                              -->
<!--==================================================================-->

<A NAME="Legalese"></A>
<div class="top">[<a href="#">top</a>]</div><hr />
<H2>Terms of Use</H2>

<P>
DART software - Copyright UCAR. This open source software is provided
by UCAR, "as is", without charge, subject to all terms of use at
<a href="http://www.image.ucar.edu/DAReS/DART/DART_download">
http://www.image.ucar.edu/DAReS/DART/DART_download</a>
</P>

<!--==================================================================-->

</BODY>
</HTML>
