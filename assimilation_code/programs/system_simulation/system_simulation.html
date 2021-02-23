<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>system simulation programs</title>
<link rel="stylesheet" type="text/css" href=
"../../../docs/html/doc.css">
<link href="../../../docs/images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>system simulation programs</h1>
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
USE</a>
<h2>Overview</h2>
<p>A collection of standalone programs for simulating various
properties of ensembles.</p>
<ul>
<li><em class="file">gen_sampling_err_table.f90</em></li>
<li><em class="file">full_error.f90</em></li>
<li><em class="file">obs_sampling_err.f90</em></li>
<li><em class="file">sampling_error.f90</em></li>
<li><em class="file">system_simulation.f90</em></li>
<li><em class="file">test_sampling_err_table.f90</em></li>
<li><em class="file">correl_error.f90</em></li>
</ul>
<p><strong>The program of most interest here is <em class=
"file">gen_sampling_err_table.f90</em> which generates the lookup
table needed when using sampling error correction in <em class=
"program">filter</em>.</strong> Talk to Jeff Anderson about the
other programs in this directory.</p>
<p>To enable the sampling error correction algorithm in <em class=
"program">filter</em>, set the namelist item <a href=
"../../modules/assimilation/assim_tools_mod.html#Namelist">&amp;assim_tools_nml : sampling_error_correction</a>
to <em class="input">.true.</em>, and copy the netCDF file
<a class="file">system_simulation/sampling_error_correction_table.nc
into the run directory.<br>
<br>
The supported set of precomputed ensemble sizes can be found by
exploring the <em class="code">ens_sizes</em> variable in</a>
<a class="file">sampling_error_correction_table.nc. To add support
for another ensemble size, build the executables in the</a>
<a href="../system_simulation/work">work</a> directory, (usually by
running <em class="program">quickbuild.csh</em>) set the <em class=
"code">ens_sizes</em> (it takes a list, but keep it short) namelist
item in <em class="file">work/input.nml</em>, and run <em class=
"program">gen_sampling_err_table</em>. It generates a LARGE number
of samples <em>per ensemble size</em> for statistical rigor. Larger
ensemble sizes take longer to generate, and compiler optimizations
matter - perhaps significantly. For example, the numbers below come
from calculating one ensemble size at a time on my desktop machine
with gfortran and basic optimization:</p>
<table width="50%">
<tr>
<th align="left">ensemble size</th>
<th align="left">run-time (seconds)</th>
</tr>
<tr>
<td>10</td>
<td>57</td>
</tr>
<tr>
<td>50</td>
<td>273</td>
</tr>
<tr>
<td>100</td>
<td>548</td>
</tr>
</table>
<p>The basic structure of <a class=
"file">sampling_error_correction_table.nc is shown below.</a></p>
<div>
<pre>
0[1095] desktop:system_simulation/work % <em class=
"input">ncdump -v ens_sizes *nc</em>
netcdf sampling_error_correction_table {
dimensions:
        bins = 200 ;
        ens_sizes = UNLIMITED ; // (40 currently)
variables:
        int count(ens_sizes, bins) ;
                count:description = "number of samples in each bin" ;
        double true_corr_mean(ens_sizes, bins) ;
        double alpha(ens_sizes, bins) ;
                alpha:description = "sampling error correction factors" ;
        int ens_sizes(ens_sizes) ;
                ens_sizes:description = "ensemble size used for calculation" ;

// global attributes:
                :num_samples = 100000000 ;
                :title = "Sampling Error Corrections for fixed ensemble sizes." ;
                :reference = "Anderson, J., 2012: Localization and Sampling Error 
                              Correction in Ensemble Kalman Filter Data Assimilation.
                              Mon. Wea. Rev., 140, 2359-2371, doi: 10.1175/MWR-D-11-00013.1." ;
                :version = "" ;
data:

<em>These ensemble sizes are already supported!</em>
 ens_sizes = 5,  6,  7,  8,  9, 10, 12, 14, 15, 16, 18, 20, 22, 24, 28, 30, 32, 36, 40, 44,
            48, 49, 50, 52, 56, 60, 64, 70, 72, 80, 84, 88, 90, 96, 100, 120, 140, 160, 180, 200
}
</pre></div>
<!--==================================================================-->
<!--============== DESCRIPTION OF A NAMELIST ========================-->
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
   ens_sizes = 5,  6,  7,  8,  9, 10, 12, 14, 15, 16, 18, 20, 22, 24, 28, 30, 32, 36, 40, 44,
              48, 49, 50, 52, 56, 60, 64, 70, 72, 80, 84, 88, 90, 96, 100, 120, 140, 160, 180, 200
   debug = .false.
   /
</pre></div>
<br>
<br>
<table border="0" cellpadding="10" width="100%" summary=
'namelist description'>
<thead align="left">
<tr>
<th>Item</th>
<th>Type</th>
<th>Description</th>
</tr>
</thead>
<tr>
<td valign="top">ens_sizes</td>
<td valign="top">integer(200)</td>
<td>An array of ensemble sizes to compute. Any new size gets
appended to the variables in the netCDF file. Any order is fine,
the array does not have to be monotonic. The numbers listed in the
example exist in the file distributed with DART. <em class=
"italic">Do not get carried away by generating a lot of new
ensemble sizes in one execution.</em> The table of run-time above
should give you some indication of how long it takes to create a
new entry.</td>
</tr>
<tr>
<td valign="top">debug</td>
<td valign="top">logical</td>
<td>A switch to add some run-time output. Generally not
needed.</td>
</tr>
</table>
<!--==================================================================-->
<a name="Modules" id="Modules"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>MODULES USED</h2>
<pre>
types_mod
utilities_mod
random_seq_mod
</pre>
<!--==================================================================-->
<!-- Describe the Files Used by this module.                          -->
<!--==================================================================-->
<a name="FilesUsed" id="FilesUsed"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>FILES for <em class="program">gen_sampling_err_table</em></h2>
<ul>
<li><em class="file">input.nml</em> for the run-time input</li>
<li><em class="file">sampling_error_correction_table.nc</em> is
both read and written. Any new ensemble sizes are simply appended
to the file.</li>
<li><em class="file">dart_log.out</em> has the run-time
output.</li>
</ul>
<h2>FILES for <em class="program">full_error</em></h2>
<ul>
<li><em class="file">input.nml</em> for the run-time input</li>
<li><em class="file">final_full.N</em> are created - N is the
ensemble size.</li>
<li><em class="file">dart_log.out</em> has the run-time
output.</li>
</ul>
<!--==================================================================-->
<!-- Cite references, if need be.                                     -->
<!--==================================================================-->
<a name="References" id="References"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>REFERENCES</h2>
<ul>
<li><b>Anderson, J. L.</b>, 2012: Localization and Sampling Error
Correction in Ensemble Kalman Filter Data Assimilation. <i>Mon.
Wea. Rev.</i>, <b>140</b>, 2359-2371 <a href=
"http://dx.doi.org/doi:10.1175/MWR-D-11-00013.1">doi:
10.1175/MWR-D-11-00013.1</a></li>
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
<td valign="top"><em class=
"program">gen_sampling_err_table</em></td>
<!-- message -->
<td valign="top">duplicate ensemble size found</td>
<!-- comment -->
<td valign="top">no need to recompute an alpha for an ensemble size
already supported.</td>
</tr>
<tr><!-- routine -->
<td valign="top"><em class=
"program">gen_sampling_err_table</em></td>
<!-- message -->
<td valign="top">existing file used a different bin size</td>
<!-- comment -->
<td valign="top">The code has been modified to use a different
number of bins than the existing netCDF file. If that's what you
intend, you need to make a new file.</td>
</tr>
<tr><!-- routine -->
<td valign="top"><em class=
"program">gen_sampling_err_table</em></td>
<!-- message -->
<td valign="top">existing file uses <em class="it">N</em> samples,
the program has <em class="it">Y</em> samples.</td>
<!-- comment -->
<td valign="top">The code has been modified to use a different
number of replicates used to estimate the <em class=
"code">alphas</em>. If that's what you intend, you need to make a
new file.</td>
</tr>
<tr><!-- routine -->
<td valign="top"><em class="program">full_error</em></td>
<!-- message -->
<td valign="top">cannot handle task counts &gt; 99999</td>
<!-- comment -->
<td valign="top">Ensemble size must be less than 100,000.</td>
</tr>
<tr><!-- routine -->
<td valign="top"><em class="program">full_error</em></td>
<!-- message -->
<td valign="top">empty bin</td>
<!-- comment -->
<td valign="top">The sample size must be large enough for all bins
to have counts</td>
</tr>
</table>
</div>
<h2>KNOWN BUGS</h2>
<p>none</p>
<!--==================================================================-->
<!-- Descibe Future Plans.                                            -->
<!--==================================================================-->
<a name="FuturePlans" id="FuturePlans"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>FUTURE PLANS</h2>
<p>none at this time.</p>
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
