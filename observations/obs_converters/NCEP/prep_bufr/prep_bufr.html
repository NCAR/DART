<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>program prepbufr</title>
<link rel="stylesheet" type="text/css" href=
"../../../../docs/html/doc.css">
<link href="../../../../docs/images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>PROGRAM <em class="program">prepbufr</em></h1>
<table border="0" summary="" cellpadding="5">
<tr>
<td valign="middle"><img src=
"../../../../docs/images/Dartboard7.png" alt="DART project logo"
height="70"></td>
<td>Jump to <a href="../../../../docs/index.html">DART
Documentation Main Index</a></td>
</tr>
</table>
<a href="#DataSources">DATA SOURCES</a> / <a href=
"#Programs">PROGRAMS</a> / <a href="#Modules">MODULES</a> /
<a href="#Namelist">NAMELIST</a> / <a href="#Errors">ERRORS</a> /
<a href="#FuturePlans">FUTURE PLANS</a> / <a href="#Legalese">TERMS
OF USE</a>
<h2>Overview</h2>
<p>Translating NCEP PREPBUFR files into DART obs_seq.out files
(input file to filter) is a 2 stage process. The first stage uses
NCEP software to translate the PREPBUFR file into an intermediate
text file. This is described in this document. The second step is
to translate the intermediate files into obs_seq.out files, which
is done by create_real_obs, as described in <a href=
"../ascii_to_obs/create_real_obs.html">create_real_obs</a> .</p>
<!--==================================================================-->
<a name="Instructions" id="Instructions"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>INSTRUCTIONS</h2>
<p>The prep_bufr package is free-standing and has not been
completely assimilated into the DART architecture. It also requires
adaptation of the sources codes and scripts to the computing
environment where it will be run. It is not so robust that it can
be controlled just with input parameters. It may not have the same
levels of error detection and warning that the rest of DART has, so
the user should very careful about checking the end product for
correctness.</p>
<h3 class="indent1">Overview of what needs to be built and run</h3>
<p>More detailed instructions follow, but this section describes a
quick overview of what programs you will be building and
running.</p>
<h4 class="indent1">Building</h4>
<p>Running the install.sh script will build the library and main
executable. You will probably have to edit this script to set which
fortran compiler is available on your system.</p>
<p>If you have raw unblocked PREPBUFR files you will need to
convert them to blocked format (what prepbufr expects as input).
The blk/ublk section of the build script compiles the <em class=
"code">cword.x</em> converter program.</p>
<p>If you are running on an Intel (little-endian) based machine you
will need the <em class="code">grabbufr</em> byte swapping program
that is also built by this script.</p>
<h4 class="indent1">One-shot execution</h4>
<p>If you are converting a single obs file, or are walking through
the process by hand for the first time, you can follow the more
detailed build instructions below, and then run the prep_bufr.x
program by hand. This involves the following steps:</p>
<ul>
<li>building the executables.</li>
<li>running the blocker if needed (generally not if you have
downloaded the blocked format PREPBUFR files).</li>
<li>running the binary format converter if you are on an Intel
(little-endian) machine.</li>
<li>linking the input file to a fixed input filename</li>
<li>running prepbufr.x to convert the file</li>
<li>copying the fixed output filename to the desired output
filename</li>
</ul>
<h4 class="indent1">Production mode</h4>
<p>If you have multiple days (or months) of observations that you
are intending to convert, there is a script in the work
subdirectory which is set up to run the converter on a sequence of
raw data files, and concatenate the output files together into one
output file per day. Edit the work/prepbufr.csh script and set the
necessary values in the 'USER SET PARAMETERS' section near the top.
This script can either be run from the command line, or it can be
submitted to a batch queue for a long series of conversion
runs.</p>
<h3 class="indent1">Installation of the NCEP PREPBUFR decoding
program</h3>
<p>This package is currently organized into files under the
DART/observations/NCEP/prep_bufr directory:</p>
<pre>
src           Source code of the NCEP PREPBUFR decoder
lib           NCEP BUFR library source
install.sh    A script to install the NCEP PREPBUFR decoder and the NCEP BUFR library.
exe           Executables of the decoder and converter.
data          Where the NCEP PREPBUFR files (prepqm****) could be loaded into
              from the NCAR Mass Store (the script assumes this is the default location).
work          Where we run the script to do the decoding.
convert_bufr  Source code (grabbufr) to convert the binary big-endian PREPBUFR files to 
              little-endian files, and a script to compile the program.
blk_ublk      Source code (cwordsh) to convert between blocked and unblocked format.
docs          Some background information about NCEP PREPBUFR observations.
</pre>
<h4 class="indent1">The decoding program: src/prepbufr.f</h4>
<p>The program prepbufr.f is used to decode the NCEP reanalysis
PREPBUFR data into intermediate text files. This program was
originally developed by NCEP. It has been modified to output
surface pressure, dry temperature, specific humidity, and wind
components (U/V) of conventional radiosonde, aircraft reports, and
satellite cloud motion derived wind. There are additional
observation types on the PREPBUFR files, but using them they would
require significant modifications of prepbufr and require detailed
knowledge of the NCEP PREPBUFR files. The NCEP quality control
indexes for these observations based on NCEP forecasts are also
output and used in DART observation sequence files. The NCEP
PREPBUFR decoding program is written in Fortran 77 and has been
successfully compiled on Linux computers using pgi90, SGI®
computers with f77, IBM® SP® systems with xlf, and Intel® based
Mac® with gfortran.</p>
<p>If your operating system uses modules you may need to remove the
default compiler and add the one desired for this package. For
example</p>
<ul>
<li>which pgf90 (to see if pgf90 is available.)</li>
<li>module rm intel64 netcdf64 mpich64</li>
<li>module add pgi32</li>
</ul>
<p>To compile the BUFR libraries and the decoding program, set the
CPLAT variable in the install.sh script to match the compilers
available on your system. CPLAT = linux is the default. Execute the
install.sh script to complete the compilations for the main
decoding program, the NCEP BUFR library, and the conversion
utilities.</p>
<p>The executables (i.e., prepbufr.x, prepbufr_03Z.x) are placed in
the ../exe directory.</p>
<p>Platforms tested:</p>
<ul>
<li>Linux clusters with Intel, PGI, Pathscale, GNU Fortran,</li>
<li>Mac OS X with Intel, GNU Fortran,</li>
<li>SGI Altix with Intel</li>
<li>Cray with Intel, Cray Fortran.</li>
</ul>
<h4 class="indent1">The byte-swapping program
convert_bufr/grabbufr.f</h4>
<p>For platforms with little-endian binary file format (e.g. Intel,
AMD®, and non-MIPS SGI processors) the program grabbufr.f is used
to convert the big-endian format NCEP PREPBUFR data into
little-endian format. The grabbufr.f code is written in Fortran 90,
and has been compiled can be compiled with the pgf90 compiler on a
Linux system, with gfortran on an Intel based Mac, and the ifort
compiler on other Linux machines. More detailed instructions for
building it can be found in convert_bufr/README, but the base
install script should build this by default. In case of problems,
cd into the convert_bufr subdirectory, edit convert_bufr.csh to set
your compiler, and run it to compile the converter code
(grabbufr).</p>
<p>This program reads the whole PREPBUFR file into memory, and
needs to know the size of the file (in bytes). Unfortunately, the
system call STAT() returns this size as one number in an array, and
the index into that array differs depending on the system and
sometimes the word size (32 vs 64) of the compiler. To test that
the program is using the right offset into this array, you can
compile and run the stat_test.f program. It takes a single filename
argument and prints out information about that file. One of the
numbers will be the file size in bytes. Compare this to the size
you see with the 'ls -l' command for that same file. If the numbers
do not agree, find the right index and edit the grabbufr.f source
file. Look for the INDEXVAL line near the first section of
executable code.</p>
<p>If grabbufr.f does not compile because the getarg() or iargc()
subroutines are not found or not available, then either use the
arg_test.f program to debug how to get command line arguments into
a fortran program on your system, or simply go into the grabbufr.f
source and comment out the section which tries to parse command
line arguments and comment in the hardcoded input and output
filenames. Now to run this program you must either rename the data
files to these predetermined filenames, or you can use links to
temporarily give the files the names needed.</p>
<h4 class="indent1">The blocking program blk_ublk/cword.x</h4>
<p>The prepbufr.x program expects to read a blocked input file,
which is generally what is available for download. However, if you
have an unblocked file that you need to convert, there is a
conversion program. The install.sh script will try to build this by
default, but in case of problems you can build it separately.
Change directories into the blk_ublk subdirectory and read the
README_cwordsh file for more help. The cwordsh shell-script wrapper
shows how to run the executable cwordsh.x executable.</p>
<p>Note that if you can get the blocked file formats to begin with,
this program is not needed.</p>
<h3 class="indent1">Getting the NCEP Reanalysis PREPBUFR format
data from NCAR HPSS.</h3>
<p>The NCEP PREPBUFR files (prepqmYYMMDDHH) can be found within the
NCEP reanalysis dataset, ds090.0, on NCAR Mass Store System
(HPSS).</p>
<p>To find the files:</p>
<ul>
<li>go to the <a href="http://rda.ucar.edu/datasets/ds090.0/"
target="_blank">NCAR/NCEP reanalysis archive.</a></li>
<li>Click on the "Inventories" tab.</li>
<li>Select the year you are interested in.</li>
<li>Search for files with the string "prepqm" in the name.</li>
<li>Depending on the year the format of the filenames change, but
they should contain the year, usually as 2 digits, the month, and
then either the start/stop day for weekly files, or the letters A
and B for semi-monthly files.</li>
</ul>
<p>Depending on the year you select, the prepqm files can be
weekly, monthly, or semi-monthly. Each tar file has a unique
dataset number of the form "A#####". For example, for January of
2003, the 4 HPSS TAR files are: A21899, A21900, A21901, A21902.
After September 2003, these files include AIRCRAFT data (airplane
readings taken at cruising elevation) but not ACARS data (airplane
readings taken during takeoff and landing). There are different
datasets which include ACARS data but their use is restricted and
you must contact the RDA group to get access.</p>
<p>If you are running on a machine with direct access to the NCAR
HPSS, then change directories into the prep_bufr/data subdirectory
and run:<br>
<br>
<em class="input">&gt; hsi get /DSS/A##### rawfile</em><br>
<br>
where ##### is the data set number you want.</p>
<p>These files may be readable tar files, or they may require
running the <em class="file">cosconvert</em> program first. See if
the <em class="file">tar</em> command can read them:<br>
<br>
<em class="input">&gt; tar -tvf rawfile</em><br>
<br>
If you get a good table of contents then simply rename the file and
untar it:<br>
<br>
<em class="input">&gt; mv rawfile data.tar</em><br>
<em class="input">&gt; tar -xvf data.tar</em><br>
<br>
However, if you get an error from the tar command you will need to
run the <em class="file">cosconvert</em> program to convert the
file into a readable tar file. On the NCAR machine <em class=
"machine">yellowstone</em>, run:<br>
<br>
<em class="input">&gt; /glade/u/home/rdadata/bin/cosconvert -b
rawfile data.tar</em><br>
<br>
On other platforms, download the appropriate version from: <a href=
"http://rda.ucar.edu/libraries/io/cos_blocking/utils/" target=
"_blank">http://rda.ucar.edu/libraries/io/cos_blocking/utils/</a> .
Build and run the converter and then you should have a tar file you
can unpack.</p>
<p>The output of tar should yield individual 6-hourly NCEP PREPBUFR
data files for the observations in the +/- 3-hour time windows of
00Z, 06Z, 12Z, and 18Z of each day. Note that DART obs_seq files
are organized such that a 24 hour file with 4 observation times
would contain observations from 3:01Z to 3:00Z of the next day,
centered on 6Z, 12Z, 18Z and "24Z". In addition, there are some
observations at 3:00Z on the PREPBUFR file labelled with 06Z. Then,
in order to make a full day intermediate file incorporating all the
required obs from the "next" day, you'll need the PREPBUFR files
through 6Z of the day after the last day of interest. For example,
to generate the observation sequence for Jan 1, 2003, the decoded
NCEP PREPBUFR text files for Jan 1 and 2, 2003 are needed, and
hence the PREPBUFR files</p>
<ul>
<li>prepqm03010106</li>
<li>prepqm03010112</li>
<li>prepqm03010118</li>
<li>prepqm03010200</li>
<li>prepqm03010206</li>
</ul>
<p>are needed.</p>
<h3 class="indent1">Running the NCEP PREPBUFR decoding program</h3>
<p>In prep_bufr/work/prepbufr.csh set the appropriate values of the
year, month, first day, and last day of the period you desire, and
the variable "convert" to control conversion from big- to
little-endian. Confirm that the raw PREPBUFR files are in ../data,
or that prepbufr.csh has been changed to find them. Execute
prepbufr.csh in the work directory. It has code for running in the
LSF batch environment, but not PBS.</p>
<p>Currently, this script generates decoded PREPBUFR text data each
24 hours which contains the observations within the time window of
-3:01 hours to +3:00Z within each six-hour synoptic time. These
daily output text files are named as temp_obs.yyyymmdd. These text
PREPBUFR data files can then be read by
DART/observations/NCEP/ascii_to_obs/work/<a href=
"../ascii_to_obs/create_real_obs.html">create_real_obs</a> to
generate the DART daily observation sequence files.</p>
<p>There is an alternate section in the script which creates a
decoded PREPBUFR text data file each 6 hours (so they are 1-for-1
with the original PREPBUFR files). Edit the script prepbufr.csh and
look for the commented out code which outputs 4 individual files
per day. Note that if you chose this option, you will have to make
corresponding changes in the create_obs_seq.csh script in step
2.</p>
<!--==================================================================-->
<a name="OtherModulesUsed" id="OtherModulesUsed"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>OTHER MODULES USED</h2>
<p>This is a piece of code that is intended to be 'close' to the
original, as such, we have not modified it to use the DART build
mechanism. This code does not use any DART modules.</p>
<!--==================================================================-->
<!--=================== DESCRIPTION OF A NAMELIST  ===================-->
<!--==================================================================-->
<a name="Namelist" id="Namelist"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>NAMELIST</h2>
<p>This namelist is read from the file <em class=
"file">input.nml</em>. Namelists start with an ampersand '&amp;'
and terminate with a slash '/'. Character strings that contain a
'/' must be enclosed in quotes to prevent them from prematurely
terminating the namelist.</p>
<div class="namelist">
<pre>
&amp;prep_bufr_nml
   obs_window       = 1.5,
   obs_window_upa   = 1.5,
   obs_window_air   = 1.5,
   obs_window_sfc   = 0.8,
   obs_window_cw    = 1.5,
   land_temp_error  = 2.5,
   land_wind_error  = 3.5,
   land_moist_error = 0.2,
   otype_use        = missing,
   qctype_use       = missing,
/
</pre></div>
<br>
<br>
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
<td>obs_window</td>
<td>real</td>
<td>Window of time to include observations. If &gt; 0, overrides
all the other more specific window sizes. Set to -1.0 to use
different time windows for different obs types. The window is +/-
this number of hours, so the total window size is twice this
value.</td>
</tr>
<tr>
<td>obs_window_upa</td>
<td>real</td>
<td>Window of time to include sonde observations (+/- hours) if
obs_window is &lt; 0, otherwise ignored.</td>
</tr>
<tr>
<td>obs_window_air</td>
<td>real</td>
<td>Window of time to include aircraft observations (+/- hours) if
obs_window is &lt; 0, otherwise ignored.</td>
</tr>
<tr>
<td>obs_window_sfc</td>
<td>real</td>
<td>Window of time to include surface observations (+/- hours) if
obs_window is &lt; 0, otherwise ignored.</td>
</tr>
<tr>
<td>obs_window_cw</td>
<td>real</td>
<td>Window of time to include cloud wind observations (+/- hours)
if obs_window is &lt; 0, otherwise ignored.</td>
</tr>
<tr>
<td>otype_use</td>
<td>real(300)</td>
<td>Report Types to extract from bufr file. If unspecified, all
types will be converted.</td>
</tr>
<tr>
<td>qctype_use</td>
<td>integer(300)</td>
<td>QC types to include from the bufr file. If unspecified, all QC
values will be accepted.</td>
</tr>
<tr>
<td>land_temp_error</td>
<td>real</td>
<td>observation error for land surface temperature observations
when none is in the input file.</td>
</tr>
<tr>
<td>land_wind_error</td>
<td>real</td>
<td>observation error for land surface wind observations when none
is in the input file.</td>
</tr>
<tr>
<td>land_moisture_error</td>
<td>real</td>
<td>observation error for land surface moisture observations when
none is in the input file.</td>
</tr>
</tbody>
</table>
</div>
<br>
<!--==================================================================-->
<!-- Describe the Files Used by this module.                          -->
<!--==================================================================-->
 <a name="FilesUsed" id="FilesUsed"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>FILES</h2>
<ul>
<li>input file(s); NCEP PREPBUFR observation files named using
ObsBase with the "yymmddhh" date tag on the end. Input to grabbufr
if big- to little-endian is to be done. Input to prepbufr if
not.</li>
<li>intermediate (binary) prepqm.little; output from grabbufr,
input to prepbufr.</li>
<li>intermediate (text) file(s) "temp_obs.yyyymmddhh"; output from
prepbufr, input to create_real_obs</li>
</ul>
<!--==================================================================-->
<!-- Cite references, if need be.                                     -->
<!--==================================================================-->
<a name="References" id="References"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>REFERENCES</h2>
DART/observations/NCEP/prep_bufr/docs/* (NCEP text files describing
the PREPBUFR files) 
<!--==================================================================-->
<!-- Describe all the error conditions and codes.                     -->
<!--==================================================================-->
 <a name="Errors" id="Errors"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>ERROR CODES and CONDITIONS</h2>
<p>Various, see the source code, doc directory, and README files
for more help.</p>
<h2>KNOWN BUGS</h2>
<p>none at this time</p>
<!--==================================================================-->
<!-- Describe Future Plans.                                           -->
<!--==================================================================-->
<a name="FuturePlans" id="FuturePlans"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>FUTURE PLANS</h2>
<p>This converter could be combined with the DART library code to
go directly from PREPBUFR files to DART obs_sequence files.</p>
<p>The converter should make an output file with a fully qualified
date in the name, so constructing files that start at 21Z and
continue to 3Z the next day are not as difficult to process. There
are currently two separate versions of the converter only because
one deals with the date wrap. This shouldn't be needed.</p>
<!--==================================================================-->
<!-- Legalese & Metadata                                              -->
<!--==================================================================-->
<a name="Legalese" id="Legalese"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>Terms of Use</h2>
<p>DART software - Copyright UCAR. This open source software is
provided by UCAR, "as is", without charge, subject to all terms of
use at <a href=
"http://www.image.ucar.edu/DAReS/DART/DART_download">http://www.image.ucar.edu/DAReS/DART/DART_download</a></p>
<!--==================================================================-->
</body>
</html>
