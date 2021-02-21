<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>program obs_impact_tool</title>
<link rel="stylesheet" type="text/css" href=
"../../../docs/html/doc.css">
<link href="../../../docs/images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>PROGRAM <em class="program">obs_impact_tool</em></h1>
<table border="0" summary="" cellpadding="5">
<tr>
<td valign="middle"><img src="../../../docs/images/Dartboard7.png"
alt="DART project logo" height="70"></td>
<td>Jump to <a href="../../../docs/index.html">DART Documentation
Main Index</a></td>
</tr>
</table>
<a href="#Namelist">NAMELIST</a> / <a href="#Examples">EXAMPLES</a>
/ <a href="#Modules">MODULES</a> / <a href="#FilesUsed">FILES</a> /
<a href="#References">REFERENCES</a> / <a href="#Errors">ERRORS</a>
/ <a href="#FuturePlans">PLANS</a> / <a href="#Legalese">TERMS OF
USE</a> <a name="Overview" id="Overview"></a>
<h2>Overview</h2>
<p>The standard DART algorithms compute increments for an
observation and then compute corresponding increments for each
model state variable due to that observation. To do this, DART
computes a sample regression coefficient using the prior ensemble
distributions of a state variable and the observation. The
increments for each member of the observation are multiplied by
this regression coefficient and then added to the corresponding
prior ensemble member for the state variable. However, in many
cases, it is appropriate to reduce the impact of an observation on
a state variable; this is called localization. The standard DART
algorithms allow users to specify a localization that is a function
of the horizontal (and optionally vertical) distance between the
observation and the state variable. The localization is a value
between 0 and 1 and multiplies the regression coefficient when
updating state ensemble members.</p>
<p>Sometimes, it may be desirable to do an additional localization
that is a function of the type of observation and the state vector
quantity. This program allows users to construct a table that is
read by filter at run-time to localize the impact of sets of
observation types on sets of state vectorquantities. Users can
create named sets of observation types and sets of state vector
quantities and specify a localization for the impact of the
specified observation types on the state vector quantities.</p>
<p>An example would be to create a subset of observations of tracer
concentration for a variety of tracers, and a subset of dynamic
state variable quantities like temperatures and wind components. It
has been common to set this localization value to 0 so that tracer
observations have no impact on dynamic state quantities, however,
the tool allows values between 0 and 1 to be specified.</p>
<p>This tool allows related collections of observation types and
state vector quantities to be named and then express the
relationship of the named groups to each other in a concise way. It
can also define relationships by exceptions.</p>
<p>All the listed observation types and state vector quantities
must be known by the system. If they are not, look at the
&amp;preprocess_nml :: input_items namelist which specifies which
obs_def_xxx_mod.f90 files are included, which is where observation
types are defined. Quantities are defined in <em class=
"file">assimilation_code/modules/observations/DEFAULT_obs_kinds_mod.F90</em>.
(Note you must add new quantities in 2 places if you do alter this
file.)</p>
<p>Format of the input file can be any combination of these types
of sections:</p>
<div>
<pre>


# hash mark starts a comment.

# the GROUP keyword starts a group and must be followed
# by a name.  All types or quantities listed before the END
# line becomes members of this group.

# GROUPs cannot contain nested groups.

GROUP groupname1
 QTY_xxx  QTY_xxx  QTY_xxx
 QTY_xxx                          # comments can be here
END GROUP

GROUP groupname2
 QTY_xxx  
 QTY_xxx  
 QTY_xxx
 QTY_xxx
END GROUP

# GROUPs can also be defined by specifying ALL, ALLQTYS,
# or ALLTYPES and then EXCEPT and listing the types or
# quantities which should be removed from this group.
# ALL EXCEPT must be the first line in a group, and all
# subsequent items are removed from the list.
# The items listed after EXCEPT can include the names
# of other groups.

GROUP groupnameM
ALL EXCEPT QTY_xxx QTY_xxx
QTY_xxx
END GROUP

GROUP groupnameN
ALL EXCEPT groupnameY
END GROUP


# once any groups have been defined, a single instance
# of the IMPACT table is specified by listing a TYPE,
# QTY, or group in column 1, then a QTY or GROUP
# in column 2 (the second name cannot be a specific type).
# column 3 must be 0.0 or 1.0.  subsequent entries
# that overlap previous entries have precedence
# (last entry wins).

IMPACT
 QTY_xxx     QTY_xxx      0.0
 QTY_xxx     groupname1   0.0
 groupname1  QTY_xxx      0.0
 groupname1  groupname1   0.0
END IMPACT

</pre></div>
<p>Namelist interface <a href="#Namelist"><em class=
"code">&amp;obs_impact_tool_nml</em></a> must be read from file
<em class="file">input.nml</em>.</p>
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
&amp;obs_impact_tool_nml
  input_filename          = 'cross_correlations.txt'
  output_filename         = 'control_impact_runtime.txt'
  debug                   = .false.
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
<td>input_filename</td>
<td>character(len=512)</td>
<td>Name of an ascii text file which describes how the interaction
of observations to state vector values and observations to other
observations should be controlled. See the <a href=
"#Overview">Overview</a> section for details about the format of
the input file entries.</td>
</tr>
<tr>
<td>output_filename</td>
<td>character(len=512)</td>
<td>Name of an ascii text file which created by this tool. It can
be read at filter run time to control the impact of observations on
state vector items and other observation values. The format of this
file is set by this tool and should not be modified by hand. Rerun
this tool to recreate the file.</td>
</tr>
<tr>
<td>debug</td>
<td>logical</td>
<td>If true print out debugging info.</td>
<td></td>
</tr>
</tbody>
</table>
</div>
<br>
<br>
<p><!-- --></p>
<!--==================================================================-->
<a name="Examples" id="Examples"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>EXAMPLES</h2>
<p>To prevent chemistry species from impacting the meterological
variables in the model state, and vice versa:</p>
<div>
<pre>
GROUP chem
 QTY_CO QTY_NO QTY_C2H4
END GROUP

GROUP met
 ALLQTYS EXCEPT chem
END GROUP

IMPACT
 chem   met    0.0
 met    chem   0.0
END IMPACT

</pre></div>
<p><!-- --></p>
<!--==================================================================-->
<a name="Modules" id="Modules"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>MODULES USED</h2>
<pre>
types_mod
utilities_mod
parse_args_mod
</pre>
<p><!-- --></p>
<!--==================================================================-->
<!-- Describe the Files Used by this module.                          -->
<!--==================================================================-->
<a name="FilesUsed" id="FilesUsed"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>FILES</h2>
<ul>
<li>two text files, one input and one output.</li>
<li>obs_impact_tool.nml</li>
</ul>
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
<p><!-- --></p>
<!--==================================================================-->
<!-- Describe all the error conditions and codes.                     -->
<!--==================================================================-->
<a name="Errors" id="Errors"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>ERROR CODES and CONDITIONS</h2>
<div class="errors">
<table border="1" cellspacing="1" cellpadding="10" width="100%"
summary="error msgs">
<tr>
<th>Routine</th>
<th>Message</th>
<th>Comment</th>
</tr>
<tr><!-- routine -->
<td valign="top">obs_impact_tool</td>
<!-- message -->
<td valign="top">Only use single process</td>
<!-- comment -->
<td valign="top">Only a single mpi process can be used with this
program</td>
</tr>
<tr><!-- routine -->
<td valign="top">obs_impact_tool</td>
<!-- message -->
<td valign="top">cannot nest groups</td>
<!-- comment -->
<td valign="top">Groups cannot contain other groups. You can
exclude a group from another group.</td>
</tr>
<tr><!-- routine -->
<td valign="top">obs_impact_tool</td>
<!-- message -->
<td valign="top">Impact must be 0.0 or 1.0</td>
<!-- comment -->
<td valign="top">Currently the impact can be either full or
nothing. Contact the DART developers if you want to experiment with
other values.</td>
</tr>
</table>
</div>
<h2>KNOWN BUGS</h2>
<p>none</p>
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
