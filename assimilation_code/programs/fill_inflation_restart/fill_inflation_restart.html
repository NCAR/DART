<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>program fill_inflation_restart</title>
<link rel="stylesheet" type="text/css" href=
"../../../docs/html/doc.css">
<link href="../../../docs/images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>PROGRAM <em class="program">fill_inflation_restart</em></h1>
<table border="0" summary="" cellpadding="5">
<tr>
<td valign="middle"><img src="../../../docs/images/Dartboard7.png"
alt="DART project logo" height="70"></td>
<td>Jump to <a href="../../../docs/index.html">DART Documentation
Main Index</a></td>
</tr>
</table>
<a href="#Namelist">NAMELIST</a> / 
<!--A HREF="#Modules">MODULES</A-->
<a href="#FilesUsed">FILES</a> / <a href=
"#References">REFERENCES</a> / <a href="#Errors">ERRORS</a> /
<a href="#FuturePlans">PLANS</a> / <a href="#Legalese">TERMS OF
USE</a>
<h2>Overview</h2>
<p>Utility program to create inflation restart files with constant
values.</p>
<p>These files can be used as input for the first step of a
multi-step assimilation when adaptive inflation is being used. This
allows the namelist items <em class=
"code">inf_initial_from_restart</em> and <em class=
"code">inf_sd_initial_from_restart</em> in the <em class=
"code">&amp;filter_nml</em> namelist to be <em class=
"code">.TRUE.</em> for all steps of the assimilation including the
very first one. (These items control whether inflation values are
read from an input file or read from constants in the
namelist.)</p>
<p>Adaptive inflation restart files are written at the end of a
<em class="program">filter</em> run and are needed as input for the
next timestep. This program creates files that can be used for the
initial run of filter when no inflation restart files have been
created by filter but are required to be read as input.</p>
<p>This program reads the inflation values to use from the
<em class="code">&amp;fill_inflation_restart_nml</em> namelist for
setting the prior inflation mean and standard deviation, and/or the
posterior inflation mean and standard deviation. It does not use
the inflation values in the <em class="code">&amp;filter</em>
namelist.</p>
<p>This program uses the information from the model_mod code to
determine the number of items in the state vector. It must be
compiled with the right model's model_mod, and if the items in the
state vector are selectable by namelist options, the namelist when
running this program must match exactly the namelist used during
the assimilation run.</p>
<p>It creates files with names consistent with the input names
expected by filter:</p>
<pre>
<em class="file">input_priorinf_mean.nc</em>
<em class="file">input_priorinf_sd.nc</em>
<em class="file">input_postinf_mean.nc</em>
<em class="file">input_postinf_sd.nc</em>
</pre>
<p>An older (and deprecated) alternative to running <em class=
"program">fill_inflation_restart</em> is to create inflation netcdf
files by using one of the NCO utilities like "<em class=
"program">ncap2</em>" on a copy of another restart file to set the
initial inflation mean, and another for the initial inflation
standard deviation. Inflation mean and sd values look exactly like
restart values, arranged by variable type like T, U, V, etc.</p>
<p>Depending on your version of the NCO utilities, you can use
<em class="program">ncap2</em> to set the T,U and V inf values
using one of two syntaxes:</p>
<div class="unix">
<pre>
  ncap2 -s 'T=1.0;U=1.0;V=1.0' wrfinput_d01 input_priorinf_mean.nc
  ncap2 -s 'T=0.6;U=0.6;V=0.6' wrfinput_d01 input_priorinf_sd.nc
  -or-
  ncap2 -s 'T(:,:,:)=1.0;U(:,:,:)=1.0;V(:,:,:)=1.0' wrfinput_d01 input_priorinf_mean.nc
  ncap2 -s 'T(:,:,:)=0.6;U(:,:,:)=0.6;V(:,:,:)=0.6' wrfinput_d01 input_priorinf_sd.nc
</pre></div>
<p>Some versions of the NCO utilities change the full 3D arrays
into a single scalar. If that's your result (check your output with
<tt>ncdump -h</tt>) use the alternate syntax or a more recent
version of the NCO tools.</p>
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
&amp;fill_inflation_restart_nml

   write_prior_inf    = .FALSE.
   prior_inf_mean     = -88888.8888
   prior_inf_sd       = -88888.8888

   write_post_inf     = .FALSE.
   post_inf_mean      = -88888.8888
   post_inf_sd        = -88888.8888

   single_file        = .FALSE.
   input_state_files  = ''
   verbose            = .FALSE.
/
</pre></div>
<p>The namelist controls which files are created and what values
are written to the restart files.</p>
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
<td>write_prior_inf</td>
<td>logical</td>
<td>Setting this to .TRUE. writes both the prior inflation mean and
standard deviation files: <em class=
"file">input_priorinf_mean.nc</em>, <em class=
"file">input_priorinf_sd.nc</em>.</td>
</tr>
<tr>
<td>prior_inf_mean</td>
<td>real(r8)</td>
<td>Prior inflation mean value.</td>
</tr>
<tr>
<td>prior_inf_sd</td>
<td>real(r8)</td>
<td>Prior inflation standard deviation value.</td>
</tr>
<tr>
<td>write_post_inf</td>
<td>logical</td>
<td>Setting this to .TRUE. writes both the posterior inflation mean
and standard deviation files <em class=
"file">input_postinf_mean.nc</em>, <em class=
"file">input_postinf_sd.nc</em>.</td>
</tr>
<tr>
<td>post_inf_mean</td>
<td>real(r8)</td>
<td>Posterior inflation mean value.</td>
</tr>
<tr>
<td>post_inf_sd</td>
<td>real(r8)</td>
<td>Posterior inflation standard deviation value.</td>
</tr>
<tr>
<td>single_file</td>
<td>logical</td>
<td>Currently not supported, but would be used in the case where
you have a single restart file that contains all of the ensemble
members. Must be .false.</td>
</tr>
<tr>
<td>input_state_files</td>
<td>character(:)</td>
<td>List one per domain, to be used as a template for the output
inflation files.</td>
</tr>
<tr>
<td>verbose</td>
<td>logical</td>
<td>Setting this to .TRUE. gives more output, and is generally used
for debugging</td>
</tr>
</tbody>
</table>
</div>
<br>
<br>
<p>Here is an example of a typical namelist for <em class=
"program">fill_inflation_restart</em> :</p>
<div class="namelist">
<pre>
&amp;fill_inflation_restart_nml

   write_prior_inf    = .TRUE.
   prior_inf_mean     = 1.01
   prior_inf_sd       = 0.6

   write_post_inf     = .FALSE.
   post_inf_mean      = 1.0
   post_inf_sd        = 0.6

   single_file        = .FALSE.
   input_state_files  = ''
   verbose            = .FALSE.
/
</pre></div>
<!--==================================================================-->
<!--A NAME="Modules"></A>
<div class="top">[<a href="#">top</a>]</div><hr />
<H2>MODULES USED</H2>
<PRE>
types_mod
utilities_mod
ensemble_manager_mod
assim_model_mod
model_mod
mpi_utilities_mod
</PRE-->
<!--==================================================================-->
<p><!-- stupid break to put top at end of page --></p>
<!--==================================================================-->
<!-- Describe the Files Used by this module.                          -->
<!--==================================================================-->
<a name="FilesUsed" id="FilesUsed"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>FILES</h2>
<p>Creates:</p>
<pre>
<em class="file">input_priorinf_mean.nc</em>
<em class="file">input_priorinf_sd.nc</em>
<em class="file">input_postinf_mean.nc</em>
<em class="file">input_postinf_sd.nc</em>
</pre>
<p>based on the template file from the specific model this code is
compiled for.</p>
<!--==================================================================-->
<!-- Cite references, if need be.                                     -->
<!--==================================================================-->
<a name="References" id="References"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>REFERENCES</h2>
<ul>
<li>none</li>
</ul>
<!--==================================================================-->
<!-- Describe all the error conditions and codes.                     -->
<!--==================================================================-->
<a name="Errors" id="Errors"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>ERROR CODES and CONDITIONS</h2>
<p>Only works for models which have individual restart files and
not the 'single_file' format, where all the ensemble members are
contained in one file.</p>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>KNOWN BUGS</h2>
<p>none</p>
<!--==================================================================-->
<!-- Describe Future Plans.                                           -->
<!--==================================================================-->
<a name="FuturePlans" id="FuturePlans"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>FUTURE PLANS</h2>
<p>If requested we can implement the 'single_file' version of
fill_inflation_restart.</p>
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
