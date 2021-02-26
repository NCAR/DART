<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>program gen_sampling_err_table</title>
<link rel="stylesheet" type="text/css" href=
"../../../docs/html/doc.css">
<link href="../../../docs/images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>PROGRAM <em class="program">gen_sampling_err_table</em></h1>
<table border="0" summary="" cellpadding="5">
<tr>
<td valign="middle"><img src="../../../docs/images/Dartboard7.png"
alt="DART project logo" height="70"></td>
<td>Jump to <a href="../../../docs/index.html">DART Documentation
Main Index</a></td>
</tr>
</table>
<a href="#Namelist">NAMELIST</a> / <a href="#Modules">MODULES</a> /
<a href="#FilesUsed">FILES</a> / <a href=
"#References">REFERENCES</a> / <a href="#Errors">ERRORS</a> /
<a href="#FuturePlans">PLANS</a> / <a href="#Legalese">TERMS OF
USE</a> <a name="Overview" id="Overview"></a>
<h2>Overview</h2>
<p>Utility program which computes a table of values needed to apply
Sampling Error Correction (SEC) during assimilation. These values
are used to correct covariances based on small sample size
statistics. See <a href="#References">reference</a> below.</p>
<p>The name of the SEC table is always <em class=
"file">sampling_error_correction_table.nc</em>. This is a NetCDF
format file. If this file already exists in the current directory
any tables for new ensemble sizes will be appended to the existing
file. If the file does not exist a new file will be created by this
tool. The resulting file should be copied into the current working
directory when <em class="program">filter</em> is run.</p>
<p>A file with 40 common ensemble sizes is distributed with the
system. Any new ensemble sizes can be generated on demand. Be aware
that the computation can be time consuming. The job may need to be
submitted to a batch system if many new ensemble sizes are being
generated, or start the job on a laptop and leave it to run
overnight.</p>
<p>The file contains a "sparse array" of ensemble sizes. Only sizes
which have an existing table are stored in the file so large
ensemble sizes do not require a large jump in the size of the
output file.</p>
<p>This program uses the random number generator to compute the
correction factors. The generator is seeded with the ensemble size
so repeated runs of the program will generate the same values for
the tables.</p>
<!--==================================================================-->
<!--=================== DESCRIPTION OF A NAMELIST ====================-->
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
&amp;gen_sampling_error_table_nml
   ens_sizes = -1
   debug = .false.
   /
</pre></div>

<h3>Description of each namelist entry</h3>

<dl>
<dt>ens_sizes</dt>
<dd>
<details><i>type: </i>integer(200)</details>
<p>
   List of ensemble sizes to compute Sampling Error Correction
   tables for. These do not need to be in any particular order.
   Duplicates will be removed and any sizes which already have tables
   computed in the output file will be skipped. The file which comes
   with the system already has tables computed for these ensemble
   sizes:
</p>
<pre>
<code>
ens_sizes = 5, 6, 7, 8, 9, 10, 12, 14, 15, 16, 18, 20,
            22, 24, 28, 30, 32, 36, 40, 44, 48, 49, 50,
            52, 56, 60, 64, 70, 72, 80, 84, 88, 90, 96,
            100, 120, 140, 160, 180, 200
</code>
</pre>
</dd>
<dt>debug</dt>
<dd>
<details><i>type: </i>logical</details>
<p>
   If true print out debugging info.
</p>
</dd>
</dl>

<!--==================================================================-->
<a name="Examples" id="Examples"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>EXAMPLES</h2>
<p>To add tables for ensemble sizes 128 and 256 run the program
with this namelist:</p>
<div>
<pre>
&amp;gen_sampling_error_table_nml
   ens_sizes = 128, 256,
   debug = .false.
   /
</pre></div>
<!--==================================================================-->
<a name="Modules" id="Modules"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>MODULES USED</h2>
<pre>
types_mod
utilities_mod
random_seq_mod
netcdf
</pre>
<!--==================================================================-->
<!-- Describe the Files Used by this module.                          -->
<!--==================================================================-->
<a name="FilesUsed" id="FilesUsed"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>FILES</h2>
<ul>
<li>output file is always <em class=
"file">sampling_error_corrrection_table.nc</em> If one exists new
ensemble sizes will be appended. If it doesn't exist a new file
will be created. This is a NetCDF format file.</li>
</ul>
<!--==================================================================-->
<!-- Cite references, if need be.                                     -->
<!--==================================================================-->
<a name="References" id="References"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>REFERENCES</h2>
<ul>
<li>Ref: Anderson, J., 2012: Localization and Sampling Error
Correction in Ensemble Kalman Filter Data Assimilation. Mon. Wea.
Rev., 140, 2359-2371, doi: 10.1175/MWR-D-11-00013.1.</li>
</ul>
<!--==================================================================-->
<!-- Describe all the error conditions and codes.                     -->
<!--==================================================================-->
<a name="Errors" id="Errors"></a>
<div class="top">[<a href="#">top</a>]</div>
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
<td valign="top">gen_sampling_err_tool</td>
<!-- message -->
<td valign="top"></td>
<!-- comment -->
<td valign="top"></td>
</tr>
</table>
</div>
<h2>KNOWN BUGS</h2>
<p>This table was regenerated using a different random number
generator than the Lanai release of DART used. If you are doing a
regression test with Lanai and have Sampling Error Correction
enabled you will not get bitwise reproducible results. <a href=
"mailto:dart@ucar.edu">Contact the DART team</a> if you are
interested in a table that includes the exact tables used in the
Lanai version.</p>
<!--==================================================================-->
<!-- Describe Future Plans.                                           -->
<!--==================================================================-->
<a name="FuturePlans" id="FuturePlans"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>FUTURE PLANS</h2>
<p>none</p>
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
