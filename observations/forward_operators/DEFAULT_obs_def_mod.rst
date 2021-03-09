<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN"
          "http://www.w3.org/TR/html4/strict.dtd">
<HTML>
<HEAD>
<TITLE>module DEFAULT_obs_def_mod</TITLE>
<link rel="stylesheet" type="text/css" href="../../docs/html/doc.css" />
<link href="../../docs/images/dart.ico" rel="shortcut icon" />
</HEAD>
<BODY>
<A NAME="TOP"></A>

<H1>MODULE DEFAULT_obs_def_mod</H1>

<table border=0 summary="" cellpadding=5>
<tr>
    <td valign=middle>
    <img src="../../docs/images/Dartboard7.png" alt="DART project logo" height=70 />
    </td>
    <td>Jump to <a href="../../docs/index.html">DART Documentation Main Index</a></td>
</tr>
</table>

<A HREF="#Legalese">TERMS OF USE</A>

<H2>Overview</H2>

<P>
DEFAULT_obs_def.F90 is a template used by the program <em class="unix">preprocess</em> to create
<em class="file">obs_def_mod.f90</em>.
<br>
<br>
To read more detailed instructions on how to add new
observation types, see the documentation for
<a href="obs_def_mod.html">obs_def_mod</a>.
<em class="file">obs_def_*_mod.f90</em> files are specified as
input to the <em class="unix">preprocess</em> program by namelist,
and a new <em class="file">obs_def_mod.f90</em> file is generated which contains
all the selected observation types.


<br>
<br>
Information from zero or more special obs_def modules, such as
<em class="file">obs_def_1d_state_mod.f90</em> or
<em class="file">obs_def_reanalyis_bufr_mod.f90</em>,
(also documented in this directory) are incorporated into the
DEFAULT_obs_def_mod.F90 template by <em class="unix">preprocess</em>. If no special
obs_def files are included in the preprocess namelist, a minimal
<em class="file">obs_def_mod.f90</em> is created which can only
support identity forward observation operators.  Any identity observations on
the obs_seq.out file will be assimilated, regardless of the obs types specified in
assimilate_these_obs_types.
<br>
<br>
The documentation below describes the special
formatting that is included in the <em class="file">DEFAULT_obs_def_mod.F90</em> in order to guide the
<em class="unix">preprocess</em> program.
<p>

Up to seven sections of code are inserted into <em class="file">DEFAULT_obs_def_mod.F90</em> from each of
the special <em class="file">obs_def_*_mod.f90</em> files. The insertion point for each
section is denoted by a special comment line that must be included <em>verbatim</em> in
<em class="file">DEFAULT_obs_def_mod.F90</em>.
These special comment lines and their significance are:
</P>
<ol>
<li>
<code>! DART PREPROCESS MODULE CODE INSERTED HERE </code> <BR>
Some special observation definition modules (see for instance
<em class="file">obs_def_1d_state_mod.f90</em>)
contain code for evaluating forward observation
operators, reading or writing special information about an observation
definition to an obs sequence file, or for interactive definition of an
observation. The entire module code section is inserted here,
so the resulting output file will be completely self-contained.
Fortran 90 allows multiple modules to be defined in a single source
file, and subsequent module code can use previously defined modules,
so this statement must preceed the rest of the other comment lines.<br><br>
</li>
<li>
<code>! DART PREPROCESS USE FOR OBS_QTY_MOD INSERTED HERE </code> <BR>
The quantities available to DART are defined by passing quantity files
from <em class="file">DART/assimilation_code/modules/observations</em> to
<em class="unix">preprocess</em>.
Unique integer values for each quantity are assigned by <em class="unix">preprocess</em> and
the use statements for these entries are inserted here.
<br><br>
</li>
<li>
<code>! DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE INSERTED HERE </code> <BR>
Some special observation definition modules (see for instance
<em class="file">obs_def_1d_state_mod.f90</em>)
contain code for evaluating forward observation
operators, reading or writing special information about an observation
definition to an obs sequence file, or for interactive definition of an
observation. The use statements for these routines from the special
observation definition modules are inserted here.<br><br>
</li>
<li>
<code>! DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF INSERTED HERE </code><BR>
Special observation definition modules must contain case statement code saying
what to do to evaluate a forward observation operator for each observation type
that they define.
This code is inserted here.<br><br>
</li>
<li>
<code>! DART PREPROCESS READ_OBS_DEF INSERTED HERE  </code><BR>
Special observation definition modules must contain case statement code saying
what to do to read any additional information required for each observation type
that they define from an observation sequence file.
This code is inserted here.<br><br>
</li>
<li>
<code>! DART PREPROCESS WRITE_OBS_DEF INSERTED HERE </code><BR>
Special observation definition modules must contain case statement code saying
what to do to write any additional information required for each observation
type that they define to an observation sequence file.
This code is inserted here.<br><br>
</li>
<li>
<code>! DART PREPROCESS INTERACTIVE_OBS_DEF INSERTED HERE </code><BR>
Special observation definition modules must contain case statement code saying
what to do to interactively create any additional information required for each
observation type that they define.
This code is inserted here.<br><br>
</li>
</ol>

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
