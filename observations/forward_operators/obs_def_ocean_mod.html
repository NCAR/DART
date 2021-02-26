<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>module obs_def_ocean_mod</title>
<link rel="stylesheet" type="text/css" href=
"../../docs/html/doc.css">
<link href="../../docs/images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>MODULE obs_def_ocean_mod</h1>
<table border="0" summary="" cellpadding="5">
<tr>
<td valign="middle"><img src="../../docs/images/Dartboard7.png"
alt="DART project logo" height="70"></td>
<td>Jump to <a href="../../docs/index.html">DART Documentation Main
Index</a></td>
</tr>
</table>
<a href="#Interface">INTERFACES</a> / <a href=
"#Namelist">NAMELIST</a> / <a href="#FilesUsed">FILES</a> /
<a href="#References">REFERENCES</a> / <a href="#Errors">ERRORS</a>
/ <a href="#FuturePlans">PLANS</a> / <a href=
"#PrivateComponents">PRIVATE COMPONENTS</a> / <a href=
"#Legalese">TERMS OF USE</a>
<h2>Overview</h2>

<p>DART includes a flexible, powerful, and slightly complicated
mechanism for incorporating new types of observations. The
<em class="file">obs_def_ocean_mod</em> module being described here
is used by the program <em class="unix">preprocess</em> to insert
appropriate definitions of ocean observations into the <em class=
"file">DEFAULT_obs_def_mod.f90</em> template and generate the
source files <em class="file">obs_def_mod.f90</em> and <em class=
"file">obs_kind_mod.f90</em> that are used by <em class=
"unix">filter</em> and other DART programs.<br>
<br>
There are no code segments in this module, only definitions of
observation types that map specific observation types to generic
observation quantities. DART contains logic that supports a limited
inheritance of attributes. If you need to interpolate observations
of 'FLOAT_TEMPERATURE', DART will check to see if a specific
routine is provided for that type, if none exists, the
interpolation routine for the generic 'QTY_TEMPERATURE' will be
used; that way one interpolation routine may support many
observation types.<br>
<br>
The mandatory header line is followed by lines that have the
observation type name (an all caps Fortran 90 identifier) and their
associated generic quantity identifier from the obs_kind module. If
there is no special processing needed for an observation type, and
no additional data needed beyond the standard contents of an
observation, then a third word on the line, the <em class=
"unix">COMMON_CODE</em> will instruct the preprocess program to
automatically generate all stubs and code needed for this type. For
observation types needing any special code or additional data, this
word should not be specified and the user must supply the code
manually. One of the future extensions of this module will be to
support acoustic tomographic observations, which will necessitate
specific support routines.</p>
<h3 class="indent1">ocean variable types and their corresponding
quantities</h3>
<pre>
<code>
! BEGIN DART PREPROCESS KIND LIST
!SALINITY,                      QTY_SALINITY,              COMMON_CODE
!TEMPERATURE,                   QTY_TEMPERATURE,           COMMON_CODE
!U_CURRENT_COMPONENT,           QTY_U_CURRENT_COMPONENT,   COMMON_CODE
!V_CURRENT_COMPONENT,           QTY_V_CURRENT_COMPONENT,   COMMON_CODE
!SEA_SURFACE_HEIGHT,            QTY_SEA_SURFACE_HEIGHT,    COMMON_CODE
!ARGO_U_CURRENT_COMPONENT,      QTY_U_CURRENT_COMPONENT,   COMMON_CODE
!ARGO_V_CURRENT_COMPONENT,      QTY_V_CURRENT_COMPONENT,   COMMON_CODE
!ARGO_SALINITY,                 QTY_SALINITY,              COMMON_CODE
!ARGO_TEMPERATURE,              QTY_TEMPERATURE,           COMMON_CODE
!ADCP_U_CURRENT_COMPONENT,      QTY_U_CURRENT_COMPONENT,   COMMON_CODE
!ADCP_V_CURRENT_COMPONENT,      QTY_V_CURRENT_COMPONENT,   COMMON_CODE
!ADCP_SALINITY,                 QTY_SALINITY,              COMMON_CODE
!ADCP_TEMPERATURE,              QTY_TEMPERATURE,           COMMON_CODE
!FLOAT_SALINITY,                QTY_SALINITY,              COMMON_CODE
!FLOAT_TEMPERATURE,             QTY_TEMPERATURE,           COMMON_CODE
!DRIFTER_U_CURRENT_COMPONENT,   QTY_U_CURRENT_COMPONENT,   COMMON_CODE
!DRIFTER_V_CURRENT_COMPONENT,   QTY_V_CURRENT_COMPONENT,   COMMON_CODE
!DRIFTER_SALINITY,              QTY_SALINITY,              COMMON_CODE
!DRIFTER_TEMPERATURE,           QTY_TEMPERATURE,           COMMON_CODE
!GLIDER_U_CURRENT_COMPONENT,    QTY_U_CURRENT_COMPONENT,   COMMON_CODE
!GLIDER_V_CURRENT_COMPONENT,    QTY_V_CURRENT_COMPONENT,   COMMON_CODE
!GLIDER_SALINITY,               QTY_SALINITY,              COMMON_CODE
!GLIDER_TEMPERATURE,            QTY_TEMPERATURE,           COMMON_CODE
!MOORING_U_CURRENT_COMPONENT,   QTY_U_CURRENT_COMPONENT,   COMMON_CODE
!MOORING_V_CURRENT_COMPONENT,   QTY_V_CURRENT_COMPONENT,   COMMON_CODE
!MOORING_SALINITY,              QTY_SALINITY,              COMMON_CODE
!MOORING_TEMPERATURE,           QTY_TEMPERATURE,           COMMON_CODE
!SATELLITE_MICROWAVE_SST,       QTY_TEMPERATURE,           COMMON_CODE
!SATELLITE_INFRARED_SST,        QTY_TEMPERATURE,           COMMON_CODE
!SATELLITE_SSH,                 QTY_SEA_SURFACE_HEIGHT,    COMMON_CODE
!SATELLITE_SSS,                 QTY_SALINITY,              COMMON_CODE
! END DART PREPROCESS KIND LIST
</code>
</pre>

<p class="indent1">New observation types may be added to this list
with no loss of generality. Supporting the observations and
actually <strong>assimilating</strong> them are somewhat different
and is controlled by the <em class="file">input.nml</em><em class=
"unix">&amp;obs_kind_nml</em> <a href=
"../../assimilation_code/modules/observations/obs_kind_mod.html#Namelist">
assimilate_these_obs_types</a> variable. This provides the
flexibility to have an observation sequence file containing many
different observation types and being able to selectively choose
what types will be assimilated.</p>
<!--==================================================================-->
<a name="OtherModulesUsed" id="OtherModulesUsed"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>OTHER MODULES USED</h2>
<p class="indent1">none</p>
<!--==================================================================-->
<!-- Declare all public entities ...                                  -->
<!-- duplicate public routines template as many times as necessary    -->
<!-- make sure you replace all yyyroutine?? strings                   -->
<!--==================================================================-->
<!--Note to authors. The first row of the table is different.         -->
<!--==================================================================-->
<a name="Interface" id="Interface"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>PUBLIC INTERFACES</h2>
<p class="indent1">none</p>
<!--==================================================================-->
<!-- Declare all public entities ...                                  -->
<!-- duplicate public routines template as many times as necessary    -->
<!-- make sure you replace all yyyroutine?? strings                   -->
<!--==================================================================-->
<a name="PublicEntities" id="PublicEntities"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>PUBLIC COMPONENTS</h2>
<p class="indent1">none</p>
<!--==================================================================-->
<!-- Describe the Files Used by this module.                          -->
<!--==================================================================-->
<a name="FilesUsed" id="FilesUsed"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>FILES</h2>
<p class="indent1">none</p>
<!--==================================================================-->
<!-- Cite references, if need be.                                     -->
<!--==================================================================-->
<a name="References" id="References"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>REFERENCES</h2>
<p class="indent1">none</p>
<!--==================================================================-->
<!-- Describe all the error conditions and codes.                     -->
<!--==================================================================-->
<a name="Errors" id="Errors"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>ERROR CODES and CONDITIONS</h2>
<p class="indent1">none</p>
<!--==================================================================-->
<!-- Describe the bugs.                                               -->
<!--==================================================================-->
<h2>KNOWN BUGS</h2>
<p>none at this time</p>
<!--==================================================================-->
<!-- Describe Future Plans.                                           -->
<!--==================================================================-->
<a name="FuturePlans" id="FuturePlans"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>FUTURE PLANS</h2>
<ul class="indent1">
<li>support acoustic tomography observations - similar to GPS
soundings in the atmosphere. Take a look at <a href=
"obs_def_gps_mod.f90">obs_def_gps_mod.f90</a></li>
</ul>
<!--==================================================================-->
<!-- PrivateComponents                                                -->
<!--==================================================================-->
<a name="PrivateComponents" id="PrivateComponents"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr>
<h2>PRIVATE COMPONENTS</h2>
<p>N/A</p>
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
