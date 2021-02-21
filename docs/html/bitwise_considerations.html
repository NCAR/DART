<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>bitwise considerations</title>
<link rel="stylesheet" type="text/css" href="../html/doc.css">
<link href="../images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>Bitwise Considerations</h1>
<table border="0" summary="" cellpadding="5">
<tr>
<td valign="middle"><img src="../images/Dartboard7.png" alt=
"DART project logo" height="70"></td>
<td>Jump to <a href="../index.html">DART Documentation Main
Index</a></td>
</tr>
</table>
By bitwise we mean bit for bit identical results in the output
obs_sequence file and the restarts (netcdf-to-netcdf or DART
format-to-dart format) when comparing one version of the code to
another. For testing the code to be bitwise with Lanai there are
several things to change/set in Manhattan and your mkmf:
<ol>
<li>assim_tools_mod.f90:<br>
<code>lanai_bitwise = .true.</code> This is hard coded to false
because it is very slow to be bitwise with Lanai.
<code>lanai_bitwise = .true.</code> causes the vertical conversion
of observations to be done in get_close_obs inside the sequential
obs do loop. See <a href="vertical_conversion.html">vertical
conversion</a> for details about this change.</li>
<li>filter_nml:<br>
<code>output_forward_op_errors = .true.</code> This will cause the
forward operator code to calculate all forward operators, even if
some ensemble members fail.</li>
<li>mkmf:<br>
ifort fp-model precise<br>
In general use the debugging FLAGS in the mkmfs provided with
DART.</li>
<li>sampling_error_correction:<br>
Use the sampling_error_correction_table.Lanai.nc from the
assimilation_code/programs/system_simulation/work directory. These
values are identical to the Lanai release.</li>
<li>posterior inflation:<br>
Try to avoid testing cases which use posterior inflation. The
posterior inflation has additional code in the Manhattan version
compared to Lanai. If you need to test a case that has posterior,
copy assim_tools/assim_tools_mod.f90 from the 'classic' release to
your Lanai build, and run your test against that version.</li>
</ol>
<h3>Important</h3>
<p>The <em>CAM</em> and <em>bgrid_solo</em> model_mods have been
altered so the state is in a different order inside filter. Thus
DART format restarts will <b>not</b> be bitwise with Lanai DART
format restarts, but netcdf files will be (after running
dart_to_cam).</p>
<p><!-- make sure the 'top' is aligned correctly --></p>
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
