<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN"
          "http://www.w3.org/TR/html4/strict.dtd">
<HTML>
<HEAD>
<TITLE>program preprocess</TITLE>
<link rel="stylesheet" type="text/css" href="../../../docs/html/doc.css" />
<link href="../../../docs/images/dart.ico" rel="shortcut icon" />
</HEAD>
<BODY>
<A NAME="TOP"></A>

<H1>PROGRAM preprocess</H1>

<table border=0 summary="" cellpadding=5>
<tr>
    <td valign=middle>
    <img src="../../../docs/images/Dartboard7.png" alt="DART project logo" height=70 />
    </td>
    <td>Jump to <a href="../../../docs/index.html">DART Documentation Main Index</a></td>
</tr>
</table>

<A HREF="#Namelist">NAMELIST</A> /
<A HREF="#Modules">MODULES</A> /
<A HREF="#FilesUsed">FILES</A> /
<A HREF="#References">REFERENCES</A> /
<A HREF="#Errors">ERRORS</A> /
<A HREF="#FuturePlans">PLANS</A> /
<A HREF="#Legalese">TERMS OF USE</A>

<H2>Overview</H2>

<P>
Preprocess is a DART-supplied preprocessor program.  Preprocess is used to
insert observation specific code into DART at compile time.

<p>
In DART, forward operators are not specific to any one model. To achieve this
separation between models and forward operators DART makes a distinction between
an observation <em>type</em> and a physical <em>quantity</em>.
For example, a radiosonde used to measure windspeed would be a <em>type</em> of
observation. Zonal wind and meridional wind are <em>quantities</em> used
to calculate windspeed. Specifying many observation types allows DART to be able
to evaluate some observations and assimilate others even if the instruments
measure the same quantity.

<p>
Preprocess takes user supplied observation and quantity files and combines them
with template files to produce code for DART. Use the namelist option
'obs_type_files' to specify the input observation files and  the namelist option
'quantity_files' to specify the input quantity files.

<p>
<ul>
<li>If no quantity files are given, a default list of quantities is used.
<li>If no obs_type_files are given, only identity observations can be used in
  the filter (i.e. the state variable values are directly observed;
forward operator is an identity)
</ul>

<p>
The template files <em class="file">DEFAULT_obs_def_mod.F90</em> and
<em class="file">DEFAULT_obs_kind_mod.F90</em> contain specially formatted
comment lines. These comment lines are used as markers to insert observation
specific information. Prepreocess relies these comment lines being used <em>verbatim</em>.

<p>
There is no need to to alter <em class="file">DEFAULT_obs_def_mod.F90</em> or
<em class="file">DEFAULT_obs_kind_mod.F90</em>. Detailed instructions for adding
new observation types can be found in <a href="obs_def_mod.html">obs_def_mod.html</a>. New
quantities should be added to a quantity file, for example a new atmosphere quantity
should be added to <em class="file">atmosphere_quantities_mod.f90</em>.
<p>
Every line in a quantity file between the start and end markers must be a
comment or a quantity definition (QTY_string). Multiple name-value pairs can be
specified for a quantity but are not required. For example, temperature may be defined:
<code> ! QTY_TEMPERATURE units="K" minval=0.0</code>. Comments are allowed between
quantity definitions or on the same line as the definition.  The code snippet below
shows acceptable formats for quantity definitions
<p>
<code>
! BEGIN DART PREPROCESS QUANTITY DEFINITIONS <br>
! <br>
! Formats accepted:<br>
!<br>
! QTY_string<br>
! QTY_string name=value<br>
! QTY_string name=value name2=value2<br>
!<br>
! QTY_string ! comments<br>
!<br>
! ! comment<br>
!<br>
! END DART PREPROCESS QUANTITY DEFINITIONS<br>

</code>

<p>
The output files produced by preprocess are named
<em class="file">assimilation_code/modules/observations/obs_kind_mod.f90</em> and
<em class="file">observations/forward_operators/obs_def_mod.f90</em>, but can
be renamed by namelist control if needed. Be aware that if you change the name
of these output files, you will need to change the path_names files for
DART executables.
<br>
<br>

<!--==================================================================-->
<!--=================== DESCRIPTION OF A NAMELIST  ===================-->
<!--==================================================================-->

<A NAME="Namelist"></A>
<div class="top">[<a href="#">top</a>]</div><hr />
<H2>NAMELIST</H2>
<P>
When you run preprocess, the namelist is read from the file <em class=file>input.nml</em>
in the directory where preprocess is run.
<p>
Namelists start with an ampersand
'&amp;' and terminate with a slash '/'.
Character strings that contain a '/' must be
enclosed in quotes to prevent them from
prematurely terminating the namelist.
</P>

<div class=namelist>
<pre>
&amp;preprocess_nml
  overwrite_output        = .true.,
  input_obs_def_mod_file  = '../../../observations/forward_operators/DEFAULT_obs_def_mod.F90',
  output_obs_def_mod_file = '../../../observations/forward_operators/obs_def_mod.f90',
  input_obs_qty_mod_file  = '../../../assimilation_code/modules/observations/DEFAULT_obs_kind_mod.F90',
  output_obs_qty_mod_file = '../../../assimilation_code/modules/observations/obs_kind_mod.f90',
  quantity_files          = '../../../assimilation_code/modules/observations/atmosphere_quantities_mod.f90',
  obs_type_files          = '../../../observations/forward_operators/obs_def_reanalysis_bufr_mod.f90',
                            '../../../observations/forward_operators/obs_def_rel_humidity_mod.f90',
                            '../../../observations/forward_operators/obs_def_altimeter_mod.f90'
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

<TR><TD>input_obs_def_mod_file </TD>
    <TD>character(len=256)<BR>
    <TD> Path name of the template obs def module to be preprocessed.
      The default is <em class="file">../../../observations/forward_operators/DEFAULT_obs_def_mod.F90</em>.
This file must have the appropriate commented lines indicating where the different
parts of the input special obs definition modules are to be inserted.
</TD></TR>

<TR><TD>output_obs_def_mod_file </TD>
    <TD>character(len=256)<BR>
    <TD> Path name of output obs def module to be created by preprocess.
      The default is <em class="file">../../../observations/forward_operators/obs_def_mod.f90</em>.
</TD></TR>

<TR><TD>input_obs_qty_mod_file </TD>
    <TD>character(len=256)<BR>
    <TD> Path name of input obs quantity file to be preprocessed. The default path name is
<em class="file">../../../assimilation_code/modules/observations/DEFAULT_obs_kind_mod.F90</em>.
This file must have the appropriate commented lines indicating where the
different quantity modules are to be inserted.
</TD></TR>

<TR><TD>output_obs_qty_mod_file </TD>
    <TD>character(len=256)<BR>
    <TD> Path name of output obs quantity module to be created by preprocess.
      The default is <em class="file">../../../assimilation_code/modules/observations/obs_kind_mod.f90</em>.
</TD></TR>

<TR><TD>obs_type_files </TD>
    <TD>character(len=256)(:)<BR>
    <TD> A list of files containing observation definitions for the type of
      observations you want to use with DART. The maximum number
      of files is limited to MAX_OBS_TYPE_FILES = 1000. The DART obs_def files are
      in <em class="file">observations/forward_operators/obs_def_*.mod.f90</em>.
</TD></TR>

<TR><TD>overwrite_output </TD>
    <TD>logical<BR>
    <TD> By defualt, preprocess will overwrite the existing obs_kind_mod.f90 and
       obs_def_mod.f90 files. Set overwrite_output = .false. if you want to
       preprocess to not overwrite existing files.
</TD></TR>


</TBODY>
</TABLE>
</div>

<br />
<br />

<!--==================================================================-->

<A NAME="Modules"></A>
<div class="top">[<a href="#">top</a>]</div><hr />
<H2>MODULES USED</H2>
<PRE>
parse_arges_mod
types_mod
utilities_mod
</PRE>

<P>
Namelist interface <em class=code>&amp;preprocess_nml</em>
must be read from file <em class=file>input.nml</em>.
</P>

<!--==================================================================-->
<!-- Describe the Files Used by this module.                          -->
<!--==================================================================-->

<A NAME="FilesUsed"></A>
<div class="top">[<a href="#">top</a>]</div><hr />
<H2>FILES</H2>
<UL><LI>input_obs_def_mod_file, specified by namelist; usually <em class="file">DEFAULT_obs_def_mod.F90</em>.
    <LI>output_obs_def_mod_file, specified by namelist; usually <em class="file">obs_def_mod.f90</em>.
    <LI>input_obs_qty_mod_file, specified by namelist; usually <em class="file">DEFAULT_obs_kind_mod.F90</em>.
    <LI>output_obs_qty_mod_file, specified by namelist; usually <em class="file">obs_kind_mod.f90</em>.
    <LI>obs_type_files, specified by namelist; usually files like <em class="file">obs_def_reanalysis_bufr_mod.f90</em>.
    <LI>quantity_files, specified by namelist; usually files like <em class="file">atmosphere_quantities_mod.f90</em>.
    <LI>namelistfile
</UL>

<!--==================================================================-->
<!-- Cite references, if need be.                                     -->
<!--==================================================================-->

<A NAME="References"></A>
<div class="top">[<a href="#">top</a>]</div><hr />
<H2>REFERENCES</H2>
<ul>
<li> none </li>
</ul>

<!--==================================================================-->
<!-- Describe all the error conditions and codes.                     -->
<!--==================================================================-->

<A NAME="Errors"></A>
<div class="top">[<a href="#">top</a>]</div><hr />
<H2>ERROR CODES and CONDITIONS</H2>
<div class=errors>
<TABLE border=1 cellspacing=1 cellpadding=10 width=100%>
<TR><TH>Routine</TH><TH>Message</TH><TH>Comment</TH></TR>


<TR><!-- routine --><TD VALIGN=top>preprocess</TD>
    <!-- message --><TD VALIGN=top>file ____ does not exist (and must)</TD>
    <!-- comment --><TD VALIGN=top>The input obs_type_files and qty_files must exist.</TD>
</TR>


<TR><!-- routine --><TD VALIGN=top>preprocess</TD>
    <!-- message --><TD VALIGN=top>file _____ does NOT contain ! BEGIN DART PREPROCESS QUANTITY LIST</TD>
    <!-- comment --><TD VALIGN=top>Each special obs_def input file must contain this comment string.</TD>
</TR>

<TR><!-- routine --><TD VALIGN=top>preprocess</TD>
    <!-- message --><TD VALIGN=top>file _____ does NOT contain " END DART PREPROCESS KIND LIST</TD>
    <!-- comment --><TD VALIGN=top>Each special obs_def input file must contain this comment string.</TD>
</TR>

<TR><!-- routine --><TD VALIGN=top>preprocess</TD>
    <!-- message --><TD VALIGN=top>Input DEFAULT obs_kind file ended unexpectedly.</TD>
    <!-- comment --><TD VALIGN=top>Did not find strings indicating where to insert special obs_def sections in the input obs_kind module.</TD>
</TR>

<TR><!-- routine --><TD VALIGN=top>preprocess</TD>
    <!-- message --><TD VALIGN=top>Input DEFAULT obs_def file ended unexpectedly.</TD>
    <!-- comment --><TD VALIGN=top>Did not find strings indicating where to insert special obs_def sections in the input obs_def module.</TD>
</TR>

<TR><!-- routine --><TD VALIGN=top>preprocess</TD>
    <!-- message --><TD VALIGN=top>file _____ does NOT contain ! BEGIN DART PREPROCESS.</TD>
    <!-- comment --><TD VALIGN=top>Input special obs_def file must contain this comment string.</TD>
</TR>

<TR><!-- routine --><TD VALIGN=top>preprocess</TD>
    <!-- message --><TD VALIGN=top>file _____ does NOT contain ! END DART PREPROCESS.</TD>
    <!-- comment --><TD VALIGN=top>Input special obs_def file must contain this comment string.</TD>
</TR>

<TR><!-- routine --><TD VALIGN=top>preprocess</TD>
    <!-- message --><TD VALIGN=top>'Incompatible duplicate entry detected'</TD>
    <!-- comment --><TD VALIGN=top>A quantity has been defined more than once but with conflicting metadata in each definition.</TD>
</TR>




</TABLE>
</div>

<H2>KNOWN BUGS</H2>
<P>
none
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
