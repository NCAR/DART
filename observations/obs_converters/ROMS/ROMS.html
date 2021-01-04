<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN"
          "http://www.w3.org/TR/html4/strict.dtd">
<HTML>
<HEAD>
<TITLE>ROMS verification observations to DART observation sequences</TITLE>
<link rel="stylesheet" type="text/css" href="../../../docs/html/doc.css" />
<link href="../../../docs/images/dart.ico" rel="shortcut icon" />
</HEAD>
<BODY>
<A NAME="TOP"></A>

<H1>ROMS observations to DART observation sequences</H1>

<table border=0 summary="" cellpadding=5>
<tr>
    <td valign=middle>
    <img src="../../../docs/images/Dartboard7.png" alt="DART project logo" height=70 />
    </td>
    <td>Jump to <a href="../../../docs/index.html">DART Documentation Main Index</a></td>
</tr>
</table>

<A HREF="#DataSources">DATA SOURCES</A> /
<A HREF="#Programs">PROGRAMS</A> / 
<A HREF="#Namelist">NAMELIST</A> /
<A HREF="#Modules">MODULES</A> /
<A HREF="#Errors">ERRORS</A> /
<A HREF="#FuturePlans">FUTURE PLANS</A> /
<A HREF="#Legalese">TERMS OF USE</A>

<H2>Overview</H2>

<P>
The relationship between ROMS and DART is slightly different
than most other models.  ROMS has the ability to apply its own forward operator 
as the model is advancing (a capability needed for variational assimilation) 
which produces something the ROMS community calls '<em class=italic>verification</em>' 
observations.  The observation file that is input to ROMS is specified by the 
<em class=file>s4dvar.in</em>:<em class=code>OBSname</em> variable.
The verification obs are written out to a netcdf file whose name is specified by the
<em class=file>s4dvar.in</em>:<em class=code>MODname</em> variable.
Since each ROMS model is advancing independently, a set of verification observation
files are created during a DART/ROMS assimilation cycle. This set of files
can be converted using <em class=program>convert_roms_obs</em> to produce a 
DART observation sequence file that has precomputed forward operators (FOs).  
<em class=program>convert_roms_obs</em> can also convert 
<em class=file>s4dvar.in</em>:<em class=code>OBSname</em>,<em class=code>MODname</em> 
files to a DART observation sequence file that does not have the precomputed FOs.
</P>

<P>
The ROMS verification observation files <strong>must</strong> contain
the <strong>obs_provenance as a global attribute</strong> and the following variables:
</P>
<ul>
   <li><em class=green>obs_lat, obs_lon, obs_depth</em></li>
   <li><em class=green>obs_value</em></li>
   <li><em class=green>obs_error</em></li>
   <li><em class=green>obs_time</em></li>
   <li><em class=green>NLmodel_value</em></li>
   <li><em class=green>obs_scale</em></li>
   <li><em class=green>obs_provenance</em></li>
</ul>
<P>Note that the 
   <em class=green>obs_provenance:flag_values</em>, and 
   <em class=green>obs_provenance:flag_meanings</em> 
   attributes are totally ignored - those relationships are specified by the
   global attribute <strong>obs_provenance</strong>.
</P>

<P>Locations only specified by <em class=green>obs_Xgrid, obs_Ygrid, obs_depth</em> 
are <strong>not</strong> supported.
</P>

<P>
The conversion of a (set of) ROMS verification observations requires metadata to
coordinate the relationship of the ROMS observation provenance to a DART observation 
TYPE.  ROMS provides significant flexibility when specifying the observation 
provenance and it is simply impractical for DART to try to support all of them. 
An example of the current practice is described in the 
<a href="#Programs">PROGRAMS</a> section below.
<br/>
<br/>
<strong>Important:</strong> <em class=program>filter</em> and 
<em class=program>perfect_model_obs</em> must also be informed which DART 
observation types use precomputed forward operators.  This is done by 
setting the <em class=file>input.nml</em><em class=program>&amp;obs_kind_nml</em> 
namelist. An example is shown at the end of the 
<a href="#obs_kind_section">PROGRAMS</a> section below.
</P>


<!--==================== DESCRIPTION OF A NAMELIST =====================-->

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
&amp;convert_roms_obs_nml
   ens_size               = 1
   roms_mod_obs_files     = ''
   roms_mod_obs_filelist  = 'filelist.txt'
   dart_output_obs_file   = 'obs_seq.out'
   append_to_existing     = .false.
   use_precomputed_values = .true.
   add_random_noise       = .false.
   pert_amplitude         = 0.01
   verbose                = 0
   type_translations      = 'NULL'
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

<TR><TD> ens_size </TD>
    <TD> integer </TD>
    <TD> Number of ensemble members which are expected to
be found when creating the expected obs values.  This must
match the number of ROMS "mod" files listed in either the
'roms_mod_obs_files' or 'roms_mod_obs_filelist' namelist
items.  It is an error if they are not the same length.
</TD></TR>

<TR><TD> roms_mods_obs_files </TD>
    <TD> character(len=256), dimension(100) </TD>
    <TD> List of filenames, one per ensemble member, that contain the
observation values for each ensemble member.  These are output from
the ROMS program.  If listing the files explicitly in this list,
'roms_mod_obs_filelist' must be '&nbsp;' (null).
</TD></TR>

<TR><TD> roms_mods_obs_filelist </TD>
    <TD> character(len=256) </TD>
    <TD> The name of an ASCII file which contains, one per line,
a list of filenames, one per ensemble member, that contain the
expected obs values for each ensemble member.  The
filenames should NOT be quoted. These are output from
the ROMS program.  If using a filelist, then
'roms_mod_obs_files' must be '&nbsp;' (null).
</TD></TR>

<TR><TD> dart_output_obs_file </TD>
    <TD> character(len=256) </TD>
    <TD> The name of the DART obs_seq file to create.  If a file
already exists with this name, it is either appended to or overwritten
depending on the 'append_to_existing' setting below.
</TD></TR>

<TR><TD> append_to_existing </TD>
    <TD> logical </TD>
    <TD> If an existing 'dart_output_obs_file' is found,
this namelist item controls how it is handled.  If .true.
the new observations are appended to the existing file.
If .false. the new observations overwrite the existing file.
</TD></TR>

<TR><TD> use_precomputed_values </TD>
    <TD> logical </TD>
    <TD> flag to indicate that the output DART observation 
         sequence file should include the verification observation 
         values from all of the ROMS observation files. 
         If <em class=mono>.true.</em> this will result in the 
         DART file having the precomputed FOs to be used in the 
         DART assimilation. 
         If <em class=mono>.false.</em> this will result in 
         DART files having the instrument values only.
</TD></TR>

<TR><TD> add_random_noise </TD>
    <TD> logical </TD>
    <TD> Almost always should be .false. . The exception is
the first cycle of an assimilation if all the ROMS input
files are identical (no ensemble currently exists).  
To create differences in the forward
operator values (since they are computed by ROMS), we can
add gaussian noise here to give them perturbed values.
This should be set as well as the 
"perturb_from_single_instance = .true."
namelist in the <em class=program>&amp;filter_nml</em> namelist.
After the first cycle, both these should be set back
to .false. .
</TD></TR>

<TR><TD> pert_amplitude </TD>
    <TD> real(r8) </TD>
    <TD> Ignored unless 'add_random_noise' is .true. .
Controls the range of random values added to the
expected obs values.  Sets the width of a gaussian.
</TD></TR>

<TR><TD> verbose </TD>
    <TD> integer </TD>
    <TD> If greater than 0, prints more information during
the conversion.  
</TD></TR>

<TR><TD> type_translations </TD>
    <TD> character(256), dimension(2, 100) </TD>
    <TD> A set of strings which control the mapping of ROMS
observation types to DART observation types.  These should
be specified in pairs.  The first column should be a string
that occurs in the global attribute '<em class=mono>obs_provenance</em>'.
Note that the <em class=mono>obs_provenance:flag_values</em> and 
<em class=mono>obs_provenance:flag_meanings</em> attributes are ignored.
The second column should be a DART specific obs type that
is found in <em class=file>DART/assimilation_code/modules/observations/obs_kind_mod.f90</em>, 
which is created by the DART <em class=program>preprocess</em> program.
</TD></TR>

</TBODY>
</TABLE>
</div>

<br />
<br />


<!--==================================================================-->

<A NAME="DataSources"></A>
<HR />
<H2>DATA SOURCES</H2>

<P>
The origin of the input observation files used by ROMS are completely unknown to me.
</P>

<!--==================================================================-->

<A NAME="Programs"></A>
<HR />
<H2>PROGRAMS</H2>

<ul><li><a href="../convert_roms_obs.html">convert_roms_obs</a>
    <li><a href="../../../assimilation_code/programs/obs_seq_to_netcdf/obs_seq_to_netcdf.html">obs_seq_to_netcdf</a>
    <li><a href="../../../assimilation_code/programs/obs_sequence_tool/assimilation_code/programs/obs_sequence_tool/obs_sequence_tool.html">obs_sequence_tool</a>
    <li><a href="../../../assimilation_code/programs/preprocess/preprocess.html">preprocess</a>
    <li><a href="../../../assimilation_code/programs/advance_time/advance_time.html">advance_time</a>
</ul>

<P>
Only <em class=program>convert_roms_obs</em> will be discussed here. 

<P> The <strong>global attribute</strong> <em class=mono>obs_provenance</em> 
is used to relate the observation provenance to DART observation TYPES. 
The ROMS 'MODname' netCDF file(s) must have both the <em class=mono>obs_provenance</em>
variable and a <em class=mono>obs_provenance</em> <strong>global attribute</strong>.
The <strong>exact</strong> strings must be repeated in the DART 
<em class=program>convert_roms_obs_nml</em>:<em class=code>type_translations</em> 
variable to be able to convert from the integer value of the obs_provenance to
th DART type in the following example:
</P>
<em class=unix>ncdump -h roms_mod_obs.nc</em> (the output has been pruned for clarity)
</P>

<pre>
netcdf roms_mod_obs {
dimensions:
        record = 2 ;
        survey = 5376 ;
        state_var = 8 ;
        datum = 2407217 ;
variables:
        <em>{snip}</em>
        int obs_provenance(datum) ;
                <del>obs_provenance:long_name = "observation origin" ;
                obs_provenance:flag_values = 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 ;</del>
        double obs_time(datum) ;
                obs_time:long_name = "time of observation" ;
                obs_time:units = "days since 1900-01-01 00:00:00 GMT" ;
                obs_time:calendar = "gregorian" ;
        double obs_lon(datum) ;
                obs_lon:long_name = "observation longitude" ;
                obs_lon:units = "degrees_east" ;
        double obs_lat(datum) ;
                obs_lat:long_name = "observation latitude" ;
                obs_lat:units = "degrees_north" ;
        double obs_depth(datum) ;
                obs_depth:long_name = "ROMS internal depth of observation variable" ;
                obs_depth:units = "meters or fractional z-levels" ;
                obs_depth:negative_value = "downwards" ;
                obs_depth:missing_value = 1.e+37 ;
        double obs_error(datum) ;
                obs_error:long_name = "observation error covariance" ;
        double obs_value(datum) ;
                obs_value:long_name = "observation value" ;
        double obs_scale(datum) ;
                obs_scale:long_name = "observation screening/normalization scale" ;
                obs_scale:_FillValue = 0. ;
        double NLmodel_value(datum) ;
                NLmodel_value:long_name = "nonlinear model at observation locations" ;
                NLmodel_value:_FillValue = 1.e+37 ;
        <em>{snip}</em>
     :obs_provenance = "\n",
             "1: gridded AVISO sea level anomaly (zeta)\n",
             "2: gridded Aquarius SSS (salinity)\n",
             "3: XBT from Met Office (temperature)\n",
             "4: CTD from Met Office (temperature)\n",
             "5: CTD from Met Office (salinity)\n",
             "6: ARGO floats (temperature)\n",
             "7: ARGO floats (salinity)\n",
             "8: glider UCSD (temperature)\n",
             "9: glider UCSD (salinity)\n",
             "10: blended satellite SST (temperature)" ;
        <em>{snip}</em>
</pre>

<P>Note the integer values that start the obs_provenance strings are used to 
interpret the integer contents of the obs_provenance variable. They need not 
be consecutive, nor in any particular order, but they must not appear more 
than once.
<br />
<br />
The following is the relevent section of the DART <em class=file>input.nml</em>:
</P>
<pre>
&amp;convert_roms_obs_nml
   ens_size               = 32
   roms_mod_obs_filelist  = 'precomputed_files.txt'
   dart_output_obs_file   = 'obs_seq.out'
   append_to_existing     = .false.
   use_precomputed_values = .true.
   add_random_noise       = .false.
   verbose                = 1
   type_translations = "gridded AVISO sea level anomaly (zeta)", "SATELLITE_SSH",
                       "gridded Aquarius SSS (salinity)",        "SATELLITE_SSS",
                       "XBT from Met Office (temperature)",      "XBT_TEMPERATURE",
                       "CTD from Met Office (temperature)",      "CTD_TEMPERATURE",
                       "CTD from Met Office (salinity)",         "CTD_SALINITY",
                       "ARGO floats (temperature)",              "ARGO_TEMPERATURE",
                       "ARGO floats (salinity)",                 "ARGO_SALINITY",
                       "glider UCSD (temperature)",              "GLIDER_TEMPERATURE",
                       "glider UCSD (salinity)",                 "GLIDER_SALINITY",
                       "blended satellite SST (temperature)",    "SATELLITE_BLENDED_SST"
  /
</pre>

<P>A complete list of DART observation TYPES is available in <a href="../../forward_operators/obs_def_ocean_mod.f90">obs_def_ocean_mod.f90</a>
</P>

<A NAME="obs_kind_section"></A>
<P>
Any or all of the DART observation types that appear in the second column of 
<em class=mono>type_translations</em> must also be designated as 
observations that have precomputed forward operators. This is done by 
setting the <em class=file>input.nml</em><em class=program>&amp;obs_kind_nml</em> 
namelist as follows:
</P>

<pre>
&amp;obs_kind_nml
   assimilate_these_obs_types =          'SATELLITE_SSH',
                                         'SATELLITE_SSS',
                                         'XBT_TEMPERATURE',
                                         'CTD_TEMPERATURE',
                                         'CTD_SALINITY',
                                         'ARGO_TEMPERATURE',
                                         'ARGO_SALINITY',
                                         'GLIDER_TEMPERATURE',
                                         'GLIDER_SALINITY',
                                         'SATELLITE_BLENDED_SST'
   use_precomputed_FOs_these_obs_types = 'SATELLITE_SSH',
                                         'SATELLITE_SSS',
                                         'XBT_TEMPERATURE',
                                         'CTD_TEMPERATURE',
                                         'CTD_SALINITY',
                                         'ARGO_TEMPERATURE',
                                         'ARGO_SALINITY',
                                         'GLIDER_TEMPERATURE',
                                         'GLIDER_SALINITY',
                                         'SATELLITE_BLENDED_SST'
  /
</pre>

<!--==================================================================-->
<!-- Describe the bugs.                                               -->
<!--==================================================================-->

<A NAME="KnownBugs"></A>
<HR />
<H2>KNOWN BUGS</H2>
<P>
none
</P>

<!--==================================================================-->
<!-- Describe Future Plans.                                           -->
<!--==================================================================-->

<A NAME="FuturePlans"></A>
<HR />
<H2>FUTURE PLANS</H2>
<P>
<OL><li>either remove the need for <em class=green>ens_size</em> or use
        it to check the consistency of the number of files to convert. </li>
    <li>remove the <em class=green>roms_obs_files()</em> mechanism.</li>
    <li>remove the <em class=green>roms_obs_file</em> mechanism -or- IF
        the <em class=green>ens_size</em> is specified, just
        use the first file and ignore the <em class=green>roms_obs_file</em>.</li>
    <li>extend the <em class=green>type_translations</em> table to be able to
        support a superset of obs_provenance-to-DART_TYPE relationships.</li>
    <li>support the ability to use the IJK locations. 
        </li>
</OL>

<!--==================================================================-->
<!-- Legalese & Metadata                                              -->
<!--==================================================================-->

<A NAME="Legalese"></A>
<HR />
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
