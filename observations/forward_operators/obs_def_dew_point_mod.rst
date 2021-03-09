<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN"
          "http://www.w3.org/TR/html4/strict.dtd">
<HTML>
<HEAD>
<TITLE>module obs_def_dew_point_mod</TITLE>
<link rel="stylesheet" type="text/css" href="../../docs/html/doc.css" />
<link href="../../docs/images/dart.ico" rel="shortcut icon" />
</HEAD>
<BODY>
<A NAME="TOP"></A>

<H1>MODULE obs_def_dew_point_mod</H1>

<table border=0 summary="" cellpadding=5>
<tr>
    <td valign=middle>
    <img src="../../docs/images/Dartboard7.png" alt="DART project logo" height=70 />
    </td>
    <td>Jump to <a href="../../docs/index.html">DART Documentation Main Index</a></td>
</tr>
</table>

<A HREF="#Interface">INTERFACES</A> /
<A HREF="#FilesUsed">FILES</A> /
<A HREF="#References">REFERENCES</A> /
<A HREF="#Errors">ERRORS</A> /
<A HREF="#FuturePlans">PLANS</A> /
<A HREF="#Legalese">TERMS OF USE</A>

<H2>Overview</H2>

<P>
Provides a subroutine to calculate the dew point temperature from model 
temperature, specific humidity, and pressure.
<br>
<br>
Revision 2801 implements a more robust method 
(based on Bolton's Approximation) for calculating dew point.
This has been further revised to avoid a numerical instability that could lead
to failed forward operators for dewpoints almost exactly 0 C.
</P>

<!--==================================================================-->

<A NAME="OtherModulesUsed"></A>
<div class="top">[<a href="#">top</a>]</div><hr />
<H2>OTHER MODULES USED</H2>
<PRE>
types_mod
utilities_mod
location_mod (most likely threed_sphere)
assim_model_mod
obs_kind_mod
</PRE>

<!--==================================================================-->
<!-- Declare all public entities ...                                  -->
<!-- duplicate public routines template as many times as necessary    -->
<!-- make sure you replace all yyyroutine?? strings                   -->
<!--==================================================================-->
<!--Note to authors. The first row of the table is different.         -->
<!--==================================================================-->

<A NAME="Interface"></A>
<div class="top">[<a href="#">top</a>]</div><hr />
<H2>PUBLIC INTERFACES</H2>

<TABLE>
<TR><TD><em class=call>use obs_def_dew_point_mod, only : </em></TD>
                   <TD><A HREF="#get_expected_dew_point">get_expected_dew_point</A></TD></TR>
</TABLE>

<P>
   A note about documentation style. 
   Optional arguments are enclosed in brackets 
   <em class=optionalcode>[like this]</em>.
</P>

<!--============= DESCRIPTION OF A SUBROUTINE =======================-->

<A NAME="get_expected_dew_point"></A>
<br>
<div class=routine>
<em class=call> call get_expected_dew_point(state_vector, location, key, td, istatus) </em>
<pre>
real(r8),            intent(in)  :: <em class=code>state_vector</em>
type(location_type), intent(in)  :: <em class=code>location</em>
integer,             intent(in)  :: <em class=code>key</em>
real(r8),            intent(out) :: <em class=code>td</em>
integer,             intent(out) :: <em class=code>istatus</em>
</pre>
</div>

<div class=indent1>
<!-- Description -->

<P>
Calculates the dew point temperature (Kelvin).
</P>

<TABLE width=100% border=0 summary="" cellpadding=3>
<TR><TD valign=top><em class=code>state_vector&nbsp;&nbsp;&nbsp;</em></TD>
    <TD>A one dimensional representation of the model state vector</TD></TR>
<TR><TD valign=top><em class=code>location</em></TD>
    <TD>Location for this obs</TD></TR>
<TR><TD valign=top><em class=code>key</em></TD>
    <TD>Controls whether upper levels (key = 1) or 2-meter (key = 2) 
         is required.</TD></TR>
<TR><TD valign=top><em class=code>td</em></TD>
    <TD>The returned dew point temperature value</TD></TR>
<TR><TD valign=top><em class=code>istatus</em></TD>
    <TD>Returned integer describing problems with applying
                forward operator</TD></TR>
</TABLE>

</div>
<br>

<!--==================================================================-->
<!-- Describe the Files Used by this module.                          -->
<!--==================================================================-->

<A NAME="FilesUsed"></A>
<div class="top">[<a href="#">top</a>]</div><hr />
<H2>FILES</H2>
<UL><LI>NONE</LI>
</UL>

<!--==================================================================-->
<!-- Cite references, if need be.                                     -->
<!--==================================================================-->

<A NAME="References"></A>
<div class="top">[<a href="#">top</a>]</div><hr />
<H2>REFERENCES</H2>
<ol>
<li> Bolton, David, 1980: The Computation of Equivalent Potential Temperature. Monthly Weather Review, 108, 1046-1053. </li>
</ol>

<!--==================================================================-->
<!-- Describe all the error conditions and codes.                     -->
<!--==================================================================-->

<A NAME="Errors"></A>
<div class="top">[<a href="#">top</a>]</div><hr />
<H2>ERROR CODES and CONDITIONS</H2>
<div class=errors>
<TABLE border=1 cellspacing=1 cellpadding=10 width=100%>
<TR><TH>Routine</TH><TH>Message</TH><TH>Comment</TH></TR>

<TR><!-- routine --><TD VALIGN=top>get_expected_dew_point</TD>
    <!-- message --><TD VALIGN=top>'key has to be 1 (upper levels) or 2 (2-meter), got ',key</TD>
    <!-- comment --><TD VALIGN=top>The input value of key is not allowed.</TD>
</TR>
</TABLE>
</div>

<H2>KNOWN BUGS</H2>
<P>
none at this time
</P>

<!--==================================================================-->
<!-- Describe Future Plans.                                           -->
<!--==================================================================-->

<A NAME="FuturePlans"></A>
<div class="top">[<a href="#">top</a>]</div><hr />
<H2>FUTURE PLANS</H2>
<P>
none at this time
</P>

<!--==================================================================-->
<!-- PrivateComponents                                                -->
<!--==================================================================-->


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
