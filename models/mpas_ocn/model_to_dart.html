<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN"
"http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<HTML>
<HEAD>
<TITLE>program model_to_dart (MPAS OCN)</TITLE>
<link rel="stylesheet" type="text/css" href="../../docs/html/doc.css" />
<link href="../../docs/images/dart.ico" rel="shortcut icon" />
</HEAD>
<BODY>
<A NAME="TOP"></A>

<H1>PROGRAM <em class=program>model_to_dart</em> for MPAS OCN</H1>

<table border=0 summary="" cellpadding=5>
<tr>
    <td valign=middle>
    <img src="../../docs/images/Dartboard7.png" alt="DART project logo" height=70 />
    </td>
    <td>Jump to <a href="../../docs/index.html">DART Documentation Main Index</a></td>
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
   <em class=program>model_to_dart</em> is the program that reads an MPAS OCN
   analysis file (nominally named <em class=file>mpas_restart.nc</em>) and 
   creates a DART state vector file 
   (e.g. <em class=file>perfect_ics, filter_ics, ... </em>). The MPAS analysis
   files have a <strong>Time</strong> UNLIMITED Dimension, which indicates
   there may (at some point) be more than one timestep in the file. The DART
   routines are currently designed to use the LAST timestep. If the Time
   dimension of length 3, we use the third timestep. A warning message is
   issued and indicates exactly the time being used.
   <br />
   <br />
   <em class=file>input.nml</em><em class=code>&amp;mpas_vars_nml</em>
   defines the list of MPAS variables used to build the DART state vector.
   This namelist is more fully described in the 
   <a href="model_mod.html">MPAS model_mod.html</a> documentation. 
   For example:
</P>
   
<pre>&amp;mpas_vars_nml
   mpas_state_variables = 'temperature',  'QTY_TEMPERATURE',
                          'salinity',     'QTY_SALINITY',
                          'rho',          'QTY_DENSITY',
                          'u',            'QTY_EDGE_NORMAL_SPEED',
                          'h',            'QTY_SEA_SURFACE_HEIGHT'
                          'tracer1',      'QTY_TRACER_CONCENTRATION'
   /
</pre>  
 
<P>
   Conditions required for successful execution of <em class=program>model_to_dart</em> are:
</P>

<UL>
   <LI>a valid <em class=file>input.nml</em> namelist file for DART which contains</LI>
   <LI>a MPAS OCN analysis file (nominally named <em class=file>mpas_analysis.nc</em>).</LI>
</UL>

<P>
Since this program is called repeatedly for every ensemble member,
we have found it convenient to link the MPAS OCN analysis files
to a static input filename (e.g. <em class=file>mpas_analysis.nc</em>).
The default DART filename is <em class=file>dart_ics</em> -
this may be moved or linked as necessary. 
</P>

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
&amp;model_to_dart_nml
   model_to_dart_output_file = 'dart_ics'
   /
</pre>
</div>

<br />

<div class=namelist>
<pre>
&amp;model_nml
   model_analysis_filename  = 'mpas_analysis.nc'
   /

(partial namelist)
</pre>
</div>

<br />

<div class=namelist>
<pre>
&amp;mpas_vars_nml
   mpas_state_variables = '',
   mpas_state_bounds = '',
   /
</pre>
</div>

<br />
<br />


<P>The model_to_dart namelist includes:</P>

<div>
<TABLE border=0 cellpadding=10 width=100% summary='namelist description'>
<THEAD align=left>
<TR><TH> Item </TH>
    <TH> Type </TH>
    <TH> Description </TH> </TR>
</THEAD>

<TBODY valign=top>
<TR><TD> model_to_dart_output_file </TD>
    <TD> character(len=128) </TD>
    <TD>The name of the DART file containing the model state
derived from the MPAS analysis file.
</TD></TR>
</TABLE>
</div>

<br />

<P>Two more namelists need to be mentioned.
   The <a href="model_mod.html#Namelist">model_nml</a> namelist specifies the
   MPAS analysis file to be used as the source.
   The <a href="model_mod.html#mpas_vars_nml">mpas_vars_nml</a> namelist specifies the
   MPAS variables that will comprise the DART state vector.
</P>

<P>For example:</P>
<pre>
&amp;mpas_vars_nml
   mpas_state_variables = 'temperature',  'QTY_TEMPERATURE',
                          'salinity',     'QTY_SALINITY',
                          'rho',          'QTY_DENSITY',
                          'u',            'QTY_EDGE_NORMAL_SPEED',
                          'h',            'QTY_SEA_SURFACE_HEIGHT'
                          'tracer1',      'QTY_TRACER_CONCENTRATION'
   /
</pre>

<br />

<!--==================================================================-->

<A NAME="Modules"></A>
<div class="top">[<a href="#">top</a>]</div><hr />
<H2>MODULES USED</H2>
<PRE>
assim_model_mod.f90
types_mod.f90
location_mod.f90
model_to_dart.f90
model_mod.f90
null_mpi_utilities_mod.f90
obs_kind_mod.f90
random_seq_mod.f90
time_manager_mod.f90
utilities_mod.f90
</PRE>

<!--==================================================================-->
<!-- Describe the Files Used by this module.                          -->
<!--==================================================================-->

<A NAME="FilesUsed"></A>
<div class="top">[<a href="#">top</a>]</div><hr />
<H2>FILES Read</H2>
<UL><LI>MPAS analysis file; <em class=file>mpas_analysis.nc</em></LI>
    <LI>DART namelist file; <a href="work/input.nml">input.nml</a></LI>
</UL>

<H2>FILES Written</H2>
<UL><LI>DART initial conditions/restart file; e.g. <em class=file>dart_ics</em></LI>
</UL>

<!--==================================================================-->
<!-- Cite references, if need be.                                     -->
<!--==================================================================-->

<A NAME="References"></A>
<div class="top">[<a href="#">top</a>]</div><hr />
<H2>REFERENCES</H2>
<P>
none
</P>

<!--==================================================================-->
<!-- Describe all the error conditions and codes.                     -->
<!--==================================================================-->

<A NAME="Errors"></A>
<div class="top">[<a href="#">top</a>]</div><hr />
<H2>ERROR CODES and CONDITIONS</H2>
<P>
none - all error messages come from modules that have their own documentation.
</P>

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
None.
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
