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
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">
<html>
<head>
   <title>Table diag_table_tk</title>
   <meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1">
   <link rel="stylesheet" href="http://www.gfdl.noaa.gov/~fms/style/doc.css" type="text/css">
</head>
<body>
<div class="header"> <font size=1>
<a href="#INSTALLATION">INSTALLATION </a>~
<a href="#USAGE">USAGE</a> ~ 
<a href="#BUGS AND FUTURE PLANS">BUGS AND FUTURE PLANS</a>
</font></div>
<hr>
<h2>diag_table_tk</h2>
<a NAME="HEADER"></a>
<!-- BEGIN HEADER -->
<div>
<b>Contact:</b>&nbsp; Matt Harrison <br>
<b>Reviewers:</b>&nbsp;         <br>
<b>Change History:</b> <a href="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/shared/diag_manager/diag_table_tk">WebCVS Log</a> <br>
<b>Last Modified:</b>&nbsp; $Date$ <br>
<b>Language:</b>&nbsp; Perl/Tk
<br><br></div>
<!-- END HEADER -->
<!-------------------------------------------------------------------->
<a NAME="OVERVIEW"></a>
<hr>
<h4> OVERVIEW</h4>
<!-- BEGIN OVERVIEW -->
<div>
<p>The script diag_table_tk is a GUI written in Perl/Tk for building diagnostics
tables, which are used by the <a href="models/bgrid_solo/fms_src/shared/diag_manager/diag_manager.html">Diagnostics Manager</a> for run-time specification of diagnostics.</p>
</div>
<!-- END OVERVIEW -->
<!-------------------------------------------------------------------->
<a NAME="DESCRIPTION"></a>
<!-- BEGIN DESCRIPTION -->
<div><p>The diagnostics table allows users to specify sampling rates and the
choice of fields at run time. The table consists of comma-separated ASCII
values and may be hand-edited. The preferred method of building a table
is to use the provided GUI interface diag_table_tk. A default diag table
is provided with each runscript.</p>
<p>The table is separated into three sections.</p>
<ol>
<li>
<b>Global section:</b> The first two lines of the table contain the experiment
title and base date. The base date is the reference time used for the time
units. The base date must be greater than or equal to the model start date.
The date consists of six space-separated integers: year, month, day, hour,
minute, and second.
<br><br></li>
<li>
<b>File section:</b> File lines contain 6 fields - file name, output frequency,
output frequency units, file format (currently only support NetCDF), time
units and long name for time axis. The format is as follows:

<pre>
"file_name", output_freq, "output_freq_units", format, "time_units", "time_long_name"


output_freq:  
         &#62; 0  output frequency in "output_units"
         = 0  output frequency every time step
         =-1  output frequency at end of run

output_freq_units = units used for output frequency
        (years, months, days, minutes, hours, seconds)

format:   1 NetCDF

time_units   = units used to label the time axis
         (days, minutes, hours, seconds)
</pre>
</li>
<li>
<b>Field section:</b> Field lines contain 8 fields - module name, field
name, output field name, file name, time sampling (for averaging, currently
only support all timesteps), time average, other operations (currently
not implemented) and pack value (1,2,4 or 8). The format is as follows:
<pre>
"module_name", "field_name", "output_name", "file_name" "time_sampling", 
time_avg, "other_opts", packing

module_name :  e.g. "atmos_mod", "land_mod"

time_avg = .true. or .false.

packing  = 1  double precision
         = 2  float
         = 4  packed 16-bit integers
         = 8  packed 1-byte (not tested?)
</pre>
</li>
</ol>
</div>
<!-- END DESCRIPTION -->
<!-------------------------------------------------------------------->
<a NAME="INSTALLATION"></a>
<hr>
<h4> INSTALLATION</h4>
<!-- BEGIN INSTALLATION -->
<div>
diag_table_tk requires the following perl modules:
<pre>
use English;
use Tk;
use Cwd;
require Tk::FileSelect;
require Tk::Text;
use Tk::widgets qw/Dialog ErrorDialog ROText/;
use Tk::FileDialog;
use Tk::Balloon;
use File::Find;
</pre>
<p>Most of these are built by default in perl 5.004 and above; however, you
may need to install the perl Tk modules.
<p>Obtain Tk and Tk-FileDialog from:
<br><a href="http://www.cpan.org">http://www.cpan.org</a>
<p>Obtain Tk and Tk-FileDialog in RPM (red hat package manager) format
from:
<br><a href="http://rpmfind.net/linux/rpm2html/search.php?query=perl-Tk">
http://rpmfind.net/linux/rpm2html/search.php?query=perl-Tk</a>
<br><a href="http://rpmfind.net/linux/rpm2html/search.php?query=perl-Tk-FileDialog">
http://rpmfind.net/linux/rpm2html/search.php?query=perl-Tk-FileDialog</a>
<br><br></div>
<!-- END INSTALLATION -->
<!-------------------------------------------------------------------->
<a NAME="USAGE"></a>
<hr>
<h4>USAGE</h4>
<!-- BRGIN USAGE -->
<ol>
<li>
<b>Load and edit previously saved diag tables.</b>
<br>Choose "Load Table" from the "File" menu.
<br><br></li>
<li>
<b>Quick parsing of f90 source code for fields which may be registered.
Fields are grouped by module name.</b>
<br>
To obtain a list of available diagnostic fields, choose "Modify Output
Field Entry" from the main menu. Enter the path to the directory containing
your source code, and click "Search". After the search is complete, you
can look in the "Field" menu for the list of available fields.
<br><br></li>
<li>
<b>Easy table editing, including ability to delete or edit selected lines.</b>
<br>
To edit the text of an entry, click the "Show Table" button, Select
the entry you wish to edit by clicking on the "Entry List" button, then
click "Edit Entry". A new window will open in which you can make changes.
Click "Save Changes" when you are finished.
<br><br></li>
<li>
<b>Error checks to help ensure that your diag table will work properly.</b>
<br>
Ensures proper spacing and formatting.
<br><br></li>
<li>
<b>Online Help is available.</b>
<br>
Choose "Help" from the menubar.
<br><br>
</li>
</ol>
<!-- END USAGE -->
<!-------------------------------------------------------------------->
<a NAME="BUGS AND FUTURE PLANS"></a>
<hr>
<h4>BUGS AND FUTURE PLANS</h4>
<!-- BEGIN BUGS AND FUTURE PLANS -->
<div>The "cancel" button doesn't seem to work.
<p>The Show Table window should be opened by default when you click "modify
table a" or "modify table b". Visual feedback is good.
<p>It should warn you if you make changes and quit without saving.
<br><br></div>
<!-- END BUGS AND FUTURE PLANS -->
<hr>
</body>
</html>
