<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN"
          "http://www.w3.org/TR/html4/strict.dtd">
<HTML>
<HEAD>
<TITLE>program fill_inflation_restart</TITLE>
<link rel="stylesheet" type="text/css" href="../../../docs/html/doc.css" />
<link href="../../../docs/images/dart.ico" rel="shortcut icon" />
</HEAD>
<BODY>
<A NAME="TOP"></A>

<H1>PROGRAM <em class=program>fill_inflation_restart</em></H1>

<table border=0 summary="" cellpadding=5>
<tr>
    <td valign=middle>
    <img src="../../../docs/images/Dartboard7.png" alt="DART project logo" height=70 />
    </td>
    <td>Jump to <a href="../../../docs/index.html">DART Documentation Main Index</a></td>
</tr>
</table>

<A HREF="#Namelist">NAMELIST</A> /
<!--A HREF="#Modules">MODULES</A-->
<A HREF="#FilesUsed">FILES</A> /
<A HREF="#References">REFERENCES</A> /
<A HREF="#Errors">ERRORS</A> /
<A HREF="#FuturePlans">PLANS</A> /
<A HREF="#Legalese">TERMS OF USE</A>

<H2>Overview</H2>

<P>
Utility program to create inflation restart files with constant values.
</P><P>
These files can be used as input for the first step of a 
multi-step assimilation when adaptive inflation is being used. 
This allows the namelist items
<em class=code>inf_initial_from_restart</em> and
<em class=code>inf_sd_initial_from_restart</em> in the 
<em class=code>&amp;filter_nml</em>
namelist to be <em class=code>.TRUE.</em> for all steps of the 
assimilation including the very first one.   
(These items control whether inflation values
are read from an input file or read from constants in the namelist.)
</P><P>
Adaptive inflation restart files are written
at the end of a <em class=program>filter</em> run
and are needed as input for the next timestep.  
This program creates files that can be used
for the initial run of filter when no inflation restart files have
been created by filter but are required to be read as input.
</P><P>
This program reads the inflation values to use from the 
<em class=code>&amp;fill_inflation_restart_nml</em>
namelist for setting the prior inflation mean and standard deviation,
and/or the posterior inflation mean and standard deviation.
It does not use the inflation values in the 
<em class=code>&amp;filter</em> namelist.
</P>

<P>
This program uses the information from the model_mod code to determine
the number of items in the state vector.  It must be compiled with the
right model's model_mod, and if the items in the state vector are selectable
by namelist options, the namelist when running this program must match
exactly the namelist used during the assimilation run.
</P>

<P>                   
It creates files with names consistent with the input names expected by filter:
</P>
<pre>
<em class=file>input_priorinf_mean.nc</em>
<em class=file>input_priorinf_sd.nc</em>
<em class=file>input_postinf_mean.nc</em>
<em class=file>input_postinf_sd.nc</em>
</pre>                 


<P>
An older (and deprecated) alternative to running <em class=program>fill_inflation_restart</em>
is to create inflation netcdf files by using one of the NCO utilities like
"<em class=program>ncap2</em>"
on a copy of another restart file to set the initial inflation mean, and another
for the initial inflation standard deviation.  Inflation mean and sd values
look exactly like restart values, arranged by variable type like T, U, V, etc.
</P>
<P>
Depending on your version of the NCO utilities, you can 
use <em class=program>ncap2</em> to set the T,U and V inf values using one of two syntaxes:
</P>
<div class=unix>
<pre>
  ncap2 -s 'T=1.0;U=1.0;V=1.0' wrfinput_d01 input_priorinf_mean.nc
  ncap2 -s 'T=0.6;U=0.6;V=0.6' wrfinput_d01 input_priorinf_sd.nc
  -or-
  ncap2 -s 'T(:,:,:)=1.0;U(:,:,:)=1.0;V(:,:,:)=1.0' wrfinput_d01 input_priorinf_mean.nc
  ncap2 -s 'T(:,:,:)=0.6;U(:,:,:)=0.6;V(:,:,:)=0.6' wrfinput_d01 input_priorinf_sd.nc
</pre>
</div>
<P>
Some versions of the NCO utilities change the full 3D arrays into a 
single scalar.  If that's your result (check your output with <tt>ncdump -h</tt>)
use the alternate syntax or a more recent version of the NCO tools.
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
</pre>
</div>

<P>
The namelist controls which files are created and
what values are written to the restart files.
</P>

<div>
<TABLE border=0 cellpadding=10 width=100% summary='namelist description'>
<THEAD align=left>
<TR><TH> Item </TH>
    <TH> Type </TH>
    <TH> Description </TH> </TR>
</THEAD>

<TBODY valign=top>

<TR><TD>write_prior_inf</TD>
    <TD>logical</TD>
    <TD>Setting this to .TRUE. writes both the prior inflation mean and
    standard deviation files: <em class=file>input_priorinf_mean.nc</em>, 
    <em class=file>input_priorinf_sd.nc</em>.
</TD></TR>

<TR><TD>prior_inf_mean</TD>
    <TD>real(r8)</TD>
    <TD> Prior inflation mean value.
</TD></TR>

<TR><TD>prior_inf_sd</TD>
    <TD>real(r8)</TD>
    <TD>Prior inflation standard deviation value.
</TD></TR>

<TR><TD>write_post_inf</TD>
    <TD>logical</TD>
    <TD>Setting this to .TRUE. writes both the posterior inflation mean and
    standard deviation files <em class=file>input_postinf_mean.nc</em>, 
    <em class=file>input_postinf_sd.nc</em>.
</TD></TR>

<TR><TD>post_inf_mean</TD>
    <TD>real(r8)</TD>
    <TD>Posterior inflation mean value.
</TD></TR>

<TR><TD>post_inf_sd</TD>
    <TD>real(r8)</TD>
    <TD>Posterior inflation standard deviation value.
</TD></TR>

<TR><TD>single_file</TD>
    <TD>logical</TD>
    <TD>Currently not supported, but would be used in the
    case where you have a single restart file that contains
    all of the ensemble members.  Must be .false.
</TD></TR>

<TR><TD>input_state_files</TD>
    <TD>character(:)</TD>
    <TD>List one per domain, to be used as a template for the output 
        inflation files.
</TD></TR>

<TR><TD>verbose</TD>
    <TD>logical</TD>
    <TD>Setting this to .TRUE. gives more output, and is 
    generally used for debugging
</TD></TR>

</TBODY> 
</TABLE>
</div>

<br />
<br />

<P>
Here is an example of a typical namelist for <em class=program>fill_inflation_restart</em> :
</P>

<div class=namelist>
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
</pre>
</div>

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
<P><!-- stupid break to put top at end of page --></P>

<!--==================================================================-->
<!-- Describe the Files Used by this module.                          -->
<!--==================================================================-->

<A NAME="FilesUsed"></A>
<div class="top">[<a href="#">top</a>]</div><hr />
<H2>FILES</H2>
<P>
Creates:
</P>
<pre>
<em class=file>input_priorinf_mean.nc</em>
<em class=file>input_priorinf_sd.nc</em>
<em class=file>input_postinf_mean.nc</em>
<em class=file>input_postinf_sd.nc</em>
</pre>
<P>
based on the template file from the specific model
this code is compiled for.
</P>

<!--==================================================================-->
<!-- Cite references, if need be.                                     -->
<!--==================================================================-->

<A NAME="References"></A>
<div class="top">[<a href="#">top</a>]</div><hr />
<H2>REFERENCES</H2>
<ul>
<li>none</li>
</ul>

<!--==================================================================-->
<!-- Describe all the error conditions and codes.                     -->
<!--==================================================================-->

<A NAME="Errors"></A>
<div class="top">[<a href="#">top</a>]</div><hr />
<H2>ERROR CODES and CONDITIONS</H2>
<P>
Only works for models which have individual restart files and not the
'single_file' format, where all the ensemble members are contained in
one file.
</P>

<div class="top">[<a href="#">top</a>]</div><hr />
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
If requested we can implement the 'single_file' version
of fill_inflation_restart.
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
