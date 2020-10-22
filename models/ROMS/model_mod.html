<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN"
          "http://www.w3.org/TR/html4/strict.dtd">
<HTML>
<HEAD>
<TITLE>module model_mod (ROMS)</TITLE>
<link rel="stylesheet" type="text/css" href="../../docs/html/doc.css" />
<link href="../../docs/images/dart.ico" rel="shortcut icon" />
</HEAD>
<BODY>
<A NAME="TOP"></A>

<H1>MODULE model_mod (ROMS)</H1>

<table border=0 summary="" cellpadding=5>
<tr>
    <td valign=middle>
    <img src="../../docs/images/Dartboard7.png" alt="DART project logo" height=70 />
    </td>
    <td>Jump to <a href="../../docs/index.html">DART Documentation Main Index</a></td>
</tr>
</table>

<A HREF="#Namelist">NAMELIST</A> /
<A HREF="#Interface">INTERFACES</A> /
<A HREF="#FilesUsed">FILES</A> /
<A HREF="#References">REFERENCES</A> /
<A HREF="#Errors">ERRORS</A> /
<A HREF="#FuturePlans">PLANS</A> /
<A HREF="#Legalese">TERMS OF USE</A>

<H2>Overview</H2>

<P>
This is the DART interface to the 
<a href="https://www.myroms.org">Regional Ocean Modeling System</a> - <strong>ROMS</strong>.
This document describes the relationship between ROMS and DART and provides an overview of
how to perform ensemble data assimilation with ROMS to provide ocean states that are 
consistent with the information provided by various ocean observations.
<br />
<br />
Running ROMS is complicated. It is <strong>strongly</strong> recommended that you become 
very familiar with running ROMS before you attempt a ROMS-DART assimilation experiment.
Running DART is complicated. It is <strong>strongly</strong> recommended that you become
very familiar with running DART before you attempt a ROMS-DART assimilation experiment.
Running ROMS-DART takes expertise in both areas.
<br />
<br />
We recommend working through the 
<a href="../../docs/tutorial/index.html">DART tutorial</a> to learn the concepts
of ensemble data assimilation and the capabilities of DART.
<br />
<br />
The ROMS code is not distributed with DART, it can be obtained from the ROMS website 
<a href="https://www.myroms.org">https://www.myroms.org</a>. There you will also find 
instructions on how to compile and run ROMS. DART can use the 'verification observations'
from ROMS (basically the estimate of the observation at the location and time 
computed as the model advances) so it would be worthwhile to become familiar with 
that capability of ROMS. DART calls these 'precomputed forward operators'.  
<br />
<br />
DART can also use observations from the
<a href="https://www.nodc.noaa.gov/OC5/indprod.html">World Ocean Database</a> - WOD. 
The conversion from the WOD formats to the DART observation sequence format is 
accomplished by the converters in the
<em class=file>observations/obs_converters/WOD</em> directory. They are described
by <a href="../../observations/obs_converters/WOD/WOD.html">WOD.html</a>.
The DART forward operators require interpolation from the ROMS 
terrain-following and horizontally curvilinear orthogonal coordinates to the 
observation location.
<em>As of revision 12153</em>, this capability is still under development. 
</P>

<P>
<strong>A note about file names:</strong> During the course of an experiment, many files
are created. To make them unique, the <em class=code>ocean_time</em> is converted from
<blockquote>"seconds&nbsp;since&nbsp;1900-01-01&nbsp;00:00:00"</blockquote> to the equivalent
number of DAYS. An <em class=it>integer</em> number of days. The intent is to tag the 
filename to reflect the valid time of the model state. This could be used as the DSTART 
for the next cycle, so it makes sense to me.
<br />
<br />
The confusion comes when applied to the observation files.
The input observation files for the ROMS 4DVAR system typically have a DSTART that 
designates the start of the forecast cycle and the file must contain observation from DSTART
to the end of the forecast. Makes sense. The model runs to the end of the forecast,
harvesting the verification observations along the way.
<br />
<br />
So then DART converts all those verification observations and tags that file ... with the same
time tag as all the other output files ... which reflects the <em class=code>ocean_time</em>
(converted to days). The input observation file to ROMS will have a different DSTART time
in the filename than the corresponding verification files. Ugh. You are free to come up
with a better plan. These are just examples, after all. Hopefully good examples.
</P>

<P>
The procedure to perform an assimilation experiment is outlined in the following steps:
</P>

<ol><li>Compile ROMS (as per the ROMS instructions).</li>
    <li>Compile all the DART executables (in the normal fashion).</li>
    <li>Stage a directory with all the files required to advance an ensemble of ROMS models and DART.</li>
    <li>Modify the run-time controls in <em class=file>ocean.in, s4dvar.in</em> and 
        <em class=file>input.nml</em>. Since ROMS has a 
        <em class=program>Bin/subsitute</em> command, it is used to replace temporary
        placeholders with actual values at various parts during the process.</li>
    <li>Advance all the instances of ROMS; each one will produce a restart file and a
        verification observation file.</li>
    <li>Convert all the verification observation files into a single DART observation 
        sequence file with the 
        <a href="../../observations/obs_converters/ROMS/ROMS.html">convert_roms_obs.f90</a>
        program.</li>
    <li>Assimilate. (DART will read and update the ROMS files directly - 
        no conversion is necessary.)</li>
    <li>Update the control files for ROMS in preparation for the next model advance.</li>
</ol>

<P>
The <em class=file>shell_scripts</em> directory has several scripts that are intended to 
provide examples. These scripts <strong>WILL</strong> need to be modified to work on your
system and are heavily internally commented. It will be necessary to read through and 
understand the scripts.  As mentioned before, the ROMS 
<em class=program>Bin/subsitute</em> command is used to replace temporary
placeholders with actual values at various parts during the process.</li>
</P>

<TABLE border=0 cellpadding=5 width=100% summary='script description'>
<THEAD align=left>
<TR><TH> Script </TH>
    <TH> Description </TH> </TR>
</THEAD>
<TBODY valign=top>

<TR><TD><a href="shell_scripts/ensemble.sh">ensemble.sh</a></TD>
    <TD>Was written by Hernan Arango to run an ensemble of ROMS models. 
        It is an appropriate example of what is required from the ROMS perspective.
        It does no data assimilation.
    </TD>
</TR>

<TR><TD><a href="shell_scripts/stage_experiment.csh">stage_experiment.csh</a></TD>
    <TD>prepares a directory for an assimilation experiment.
        The idea is basically that everything you need should be assembled
        by this script and that this should only be run ONCE per experiment.
        After everything is staged in the experiment directory, another script
        can be run to advance the model and perform the assimilation.
        <em class=program>stage_experiment.csh</em> will also modify some of
        the template scripts and copy working versions into the experiment directory.
        This script may be run interactively, i.e. from the UNIX command line.
    </TD>
</TR>

<TR><TD><a href="shell_scripts/submit_multiple_cycles_lsf.csh">submit_multiple_cycles_lsf.csh</a></TD>
    <TD>is an executable script that submits a series of dependent jobs to an LSF
        queuing system. Each job runs <em class=program>cycle.csh</em> in the 
        experiment directory and only runs if the previous dependent job completes 
        successfully. 
    </TD>
</TR>

<TR><TD><a href="shell_scripts/cycle.csh.template">cycle.csh.template</a></TD>
    <TD>is a non-executable template that is modified by 
        <em class=program>stage_experiment.csh</em> and results in an exectuable 
        <em class=program>cycle.csh</em> in the experiment directory. 
        <em class=program>cycle.csh</em> is designed to be run as a batch job 
        and advances the ROMS model states one-by-one for the desired forecast
        length. The assimilation is performed and the control information for
        the next ROMS forecast is updated. Each model execution and 
        <em class=program>filter</em> use the same set of MPI tasks.
    </TD>
</TR>

<TR><TD><a href="shell_scripts/submit_multiple_jobs_slurm.csh">submit_multiple_jobs_slurm.csh</a></TD>
    <TD>is an executable script that submits a series of dependent jobs to an LSF
        queuing system. It is possible to submit <strong>many</strong> jobs the queue,
        but the jobs run one-at-a-time.  Every assimilation cycle is divided into two
        scripts to be able to efficiently set the resources for each phase. 
        <em class=program>advance_ensemble.csh</em> is a job array
        that advances each ROMS instance in separate jobs. When the entire job array
        finishes - and only if they all finish correctly - will the next job start to run.
        <em class=program>run_filter.csh</em> performs the assimilation and prepares
        the experiment directory for another assimilation cycle.
        <em class=program>submit_multiple_jobs_slurm.csh</em> may be run from the command 
        line in the experiment directory. Multiple assimilation cycles can be specified, so
        it is possible to put <strong>many</strong> jobs in the queue.
    </TD>
</TR>

<TR><TD><a href="shell_scripts/advance_ensemble.csh.template">advance_ensemble.csh.template</a></TD>
    <TD>is a non-executable template that is modified by 
        <em class=program>stage_experiment.csh</em> and results in an exectuable 
        <em class=program>advance_ensemble.csh</em> in the experiment directory. 
        <em class=program>advance_ensemble.csh</em> is designed to submit an job array
        to the queueing system (PBS,SLURM, or LSF) to advance the ensemble members in
        separate jobs.
    </TD>
</TR>

<TR><TD><a href="shell_scripts/run_filter.csh.template">run_filter.csh.template</a></TD>
    <TD>is a non-executable template that is modified by 
        <em class=program>stage_experiment.csh</em> and results in an exectuable 
        <em class=program>run_filter.csh</em> in the experiment directory. 
        <em class=program>run_filter.csh</em> is very similar to 
        <em class=program>cycle.csh</em> but does not advance the ROMS model instances.
    </TD>
</TR>

</TABLE>


<P>
The variables from ROMS that are copied into the DART state vector are
controlled by the <em class=file>input.nml</em> <em class=code>model_nml</em> namelist.
See below for the documentation on the &amp;model_nml
entries.  The state vector should include all variables needed to apply the forward 
observation operators as well as the prognostic variables important to restart ROMS.
<br />
<br />
The example <em class=file>input.nml</em> <em class=code>model_nml</em> demonstrates
how to construct the DART state vector. The following table explains in detail each
entry for the <em class=code>variables</em> namelist item:
</P>


<div>
<TABLE border=0 cellpadding=4 width=100% summary='ROMS variables description'>
<THEAD align=left>
<TR><TH>Column 1</TH>
    <TH>Column 2</TH>
    <TH>Column 3</TH>
    <TH>Column 4</TH>
    <TH>Column 5</TH>
</TR>
</THEAD>
<TBODY valign=top>
<TR><TD>Variable name</TD>
    <TD>DART QUANTITY</TD>
    <TD>minimum</TD>
    <TD>maximum</TD>
    <TD>update</TD>
</TR>
</TABLE>
<br />
<br />
<TABLE border=0 cellpadding=4 width=100% summary='ROMS variables description'>
<TBODY valign=top>
<TR><TD>Variable&nbsp;name&nbsp;&nbsp;</TD>
    <TD>This is the ROMS variable name as it appears in the ROMS netCDF file.</TD></TR>
<TR><TD>DART&nbsp;QUANTITY</TD>
    <TD>This is the character string of the corresponding DART QUANTITY. The complete
        list of possible DART QUANTITY values is available in the
        <a href="../../assimilation_code/modules/observations/DEFAULT_obs_kind_mod.html"
        >obs_def_mod</a> that is built by 
        <a href="../../assimilation_code/programs/preprocess/preprocess.html">
        preprocess</a></TD></TR>
<TR><TD>minimum</TD>
    <TD>If the variable is to be updated in the ROMS restart file, this specifies the
        minimum value. If set to 'NA', there is no minimum value.</TD></TR>
<TR><TD>maximum</TD>
    <TD>If the variable is to be updated in the ROMS restart file, this specifies the
        maximum value. If set to 'NA', there is no maximum value.</TD></TR>
<TR><TD>update</TD>
    <TD>The updated variable may or may not be written to the ROMS restart file.<br />
        <em class=code>'UPDATE'</em>&nbsp;&nbsp;means the variable in the restart 
        file is updated. This is case-insensitive.<br />
        <em class=code>'NO_COPY_BACK'</em>&nbsp;&nbsp;(or anything else) means the 
        variable in the restart file remains unchanged.<br /></TD></TR>
</TABLE>
</div>


<!--==================== DESCRIPTION OF A NAMELIST =====================-->

<A NAME="Namelist"></A>
<div class="top">[<a href="#">top</a>]</div><hr />
<H2>NAMELIST</H2>
<P>
This namelist is read from the file <em class=file>input.nml</em>.
Namelists start with an ampersand '&amp;' and terminate with a slash '/'.
Character strings that contain a '/' must be enclosed in quotes to 
prevent them from prematurely terminating the namelist. The default
namelist is presented below, a more realistic namelist is presented at
the end of this section.
</P>

<div class=namelist>
<pre>
&amp;model_nml
   roms_filename               = 'roms_input.nc'
   assimilation_period_days    = 1
   assimilation_period_seconds = 0
   vert_localization_coord     = 3
   debug                       = 0
   variables                   = ''
  /
</pre>
</div>
<div>
<TABLE border=0 cellpadding=10 width=100% summary='namelist description'>
<THEAD align=left>
<TR><TH> Item </TH>
    <TH> Type </TH>
    <TH> Description </TH> </TR>
</THEAD>

<TBODY valign=top>

<TR><TD>roms_filename</TD>
    <TD>character(len=256)</TD>
    <TD>This is the name of the file used to provide information about the ROMS variable
        dimensions, etc.
</TD></TR>

<TR><TD>assimilation_period_days, <br />
        assimilation_period_seconds</TD>
    <TD>integer</TD>
    <TD>Combined, these specify the width of the assimilation window.
        The current model time is used as the center time of the
        assimilation window. All observations in the assimilation window
        are assimilated. BEWARE: if you put observations that occur before
        the beginning of the assimilation_period, DART will error out because
        it cannot move the model 'back in time' to process these observations.</TD>
</TR>

<TR><TD>variables </TD>
    <TD>character(:, 5)</TD>
    <TD>A 2D array of strings, 5 per ROMS variable 
        to be added to the dart state vector.
<OL><LI>ROMS field name - must match netCDF variable name exactly</LI>
    <LI>DART QUANTITY - must match a valid DART QTY_xxx exactly</LI>
    <LI>minimum physical value  - if none, use 'NA'</LI>
    <LI>maximum physical value  - if none, use 'NA'</LI>
    <LI>case-insensitive string describing whether to copy the updated variable
        into the ROMS restart file ('UPDATE') or not (any other value).
        There is generally no point copying diagnostic variables into the restart file.
        Some diagnostic variables may be useful for computing forward operators, however.</LI>
</OL>
</TD></TR>

<TR><TD> vert_localization_coord </TD>
    <TD> integer </TD>
    <TD>Vertical coordinate for vertical localization.
       <UL style="list-style: none;">
           <LI>1 = model level</LI>
           <LI>2 = pressure (in pascals)</LI>
           <LI>3 = height (in meters)</LI>
           <LI>4 = scale height (unitless)</LI>
       </UL>
       Currently, only 3 (height) is supported for ROMS.
</TD></TR>

</TBODY>
</TABLE>
</div>

<br />
<br />

<P>A more realistic ROMS namelist is presented here, along with one of the more
unusual settings that is generally necessary when running ROMS. The
<em class=code>use_precomputed_FOs_these_obs_types</em> variable needs to
list the observation types that are present in the ROMS verification observation
file.   
</P>

<div class=namelist>
<pre>
&amp;model_nml
   roms_filename                = 'roms_input.nc'
   assimilation_period_days     = 1
   assimilation_period_seconds  = 0
   vert_localization_coord      = 3
   debug                        = 1
   variables = 'temp',   'QTY_TEMPERATURE',          'NA', 'NA', 'update',
               'salt',   'QTY_SALINITY',            '0.0', 'NA', 'update',
               'u',      'QTY_U_CURRENT_COMPONENT',  'NA', 'NA', 'update',
               'v',      'QTY_V_CURRENT_COMPONENT',  'NA', 'NA', 'update',
               'zeta',   'QTY_SEA_SURFACE_HEIGHT'    'NA', 'NA', 'update'
  /

&amp;obs_kind_nml
   evaluate_these_obs_types = ''
   assimilate_these_obs_types =          'SATELLITE_SSH',
                                         'SATELLITE_SSS',
                                         'XBT_TEMPERATURE',
                                         'CTD_TEMPERATURE',
                                         'CTD_SALINITY',
                                         'ARGO_TEMPERATURE',
                                         'ARGO_SALINITY',
                                         'GLIDER_TEMPERATURE',
                                         'GLIDER_SALINITY',
                                         'SATELLITE_BLENDED_SST',
                                         'SATELLITE_MICROWAVE_SST',
                                         'SATELLITE_INFRARED_SST'
   use_precomputed_FOs_these_obs_types = 'SATELLITE_SSH',
                                         'SATELLITE_SSS',
                                         'XBT_TEMPERATURE',
                                         'CTD_TEMPERATURE',
                                         'CTD_SALINITY',
                                         'ARGO_TEMPERATURE',
                                         'ARGO_SALINITY',
                                         'GLIDER_TEMPERATURE',
                                         'GLIDER_SALINITY',
                                         'SATELLITE_BLENDED_SST',
                                         'SATELLITE_MICROWAVE_SST',
                                         'SATELLITE_INFRARED_SST'
  /
</pre>
</div>


<!--==================================================================-->
<P><!-- useless paragraph so 'top' aligns correctly --></P>

<A NAME="Interface"></A>
<div class="top">[<a href="#">top</a>]</div><hr />
<H2>OTHER MODULES USED</H2>
<PRE>
types_mod
time_manager_mod
threed_sphere/location_mod
utilities_mod
obs_kind_mod
map_utils
netcdf
typesizes

utilities/default_model_mod.f90
</PRE>

<!--==================================================================-->
<!-- Declare all public entities ...                                  -->
<!-- duplicate public routines template as many times as necessary    -->
<!-- make sure you replace all yyyroutine?? strings                   -->
<!--==================================================================-->
<!--Note to authors. The first row of the table is different.         -->
<!--==================================================================-->

<hr />
<H2>PUBLIC INTERFACES</H2>

<P>
The 18 required public interfaces are standardized for all DART compliant models.
These interfaces allow DART to advance the model, get the model state and metadata 
describing this state, find state variables that are close to a given location,
and do spatial interpolation for a variety of variables required in
observational operators. Some of the interfaces are common to multiple models
and exist in their own modules. Only the interfaces unique to ROMS are described here.
</P>

<TABLE>
<TR><TD><em class=call>use model_mod, only : </em></TD>
                   <TD><A HREF="#get_model_size">get_model_size</A></TD></TR>
<TR><TD>&nbsp;</TD><TD><A HREF="#get_state_meta_data">get_state_meta_data</A></TD></TR>
<TR><TD>&nbsp;</TD><TD><A HREF="#model_interpolate">model_interpolate</A></TD></TR>
<TR><TD>&nbsp;</TD><TD><A HREF="#shortest_time_between_assimilations">shortest_time_between_assimilations</A></TD></TR>
<TR><TD>&nbsp;</TD><TD><A HREF="#static_init_model">static_init_model</A></TD></TR>
<TR><TD>&nbsp;</TD><TD><A HREF="#end_model">end_model</A></TD></TR>
<TR><TD>&nbsp;</TD><TD><A HREF="#nc_write_model_atts">nc_write_model_atts</A></TD></TR>
<TR><TD>&nbsp;</TD><TD><A HREF="#write_model_time">write_model_time</A></TD></TR>
<TR><TD>&nbsp;</TD><TD><A HREF="#read_model_time">read_model_time</A></TD></TR>

<TR><td colspan=2>These routines are not required, but are useful:<td></TR>

<TR><TD>&nbsp;</TD><TD><A HREF="#get_time_information">get_time_information</A></TD></TR>
<TR><TD>&nbsp;</TD><TD><A HREF="#get_location_from_ijk">get_location_from_ijk</A></TD></TR>
</TABLE>

<P>These required interfaces are described in their own module documentation:</P>

<TABLE>
<TR><TD><em class=call>use default_model_mod, only : </em></TD>
                   <TD><A HREF="../utilities/default_model_mod.html#nc_write_model_vars">nc_write_model_vars</A></TD></TR>
<TR><TD>&nbsp;</TD><TD><A HREF="../utilities/default_model_mod.html#pert_model_copies">pert_model_copies</A></TD></TR>
<TR><TD>&nbsp;</TD><TD><A HREF="../utilities/default_model_mod.html#adv_1step">adv_1step</A></TD></TR>
<TR><TD>&nbsp;</TD><TD><A HREF="../utilities/default_model_mod.html#init_time">init_time</A></TD></TR>
<TR><TD>&nbsp;</TD><TD><A HREF="../utilities/default_model_mod.html#init_conditions">init_conditions</A></TD></TR>
<tr><td colspan=2>&nbsp;</td></tr>
<TR><TD><em class=call>use location_model_mod, only : </em></TD>
                   <TD><A HREF="../../assimilation_code/location/threed_sphere/location_mod.html#get_close_obs">get_close_obs</A></TD></TR>
<TR><TD>&nbsp;</TD><TD><A HREF="../../assimilation_code/location/threed_sphere/location_mod.html#get_close_obs">get_close_state</A></TD></TR>
<TR><TD>&nbsp;</TD><TD><A HREF="../../assimilation_code/location/threed_sphere/location_mod.html#convert_vertical_obs">convert_vertical_obs</A></TD></TR>
<TR><TD>&nbsp;</TD><TD><A HREF="../../assimilation_code/location/threed_sphere/location_mod.html#convert_vertical_state">convert_vertical_state</A></TD></TR>
</TABLE>


<P>
The last 4 interfaces are only required for low-order models where advancing
the model can be done by a call to a subroutine. The ROMS model only advances by
executing the program ROMS.exe. Thus the last 4 interfaces only appear as stubs
in the ROMS module.
<br><br>
The interface pert_model_copies is presently not provided for ROMS. The initial
ensemble has to be generated off-line. If coherent structures are not required,
the filter can generate an ensemble with uncorrelated random Gaussian noise of
0.002. This is of course not appropriate for a model like ROMS which has
variables expressed in a wide range of scales. It is thus recommended to
generate the initial ensemble off-line, perhaps with the tools provided in
models/ROMS/PERTURB/3DVAR-COVAR.
</P>

<P>
   A note about documentation style.
   Optional arguments are enclosed in brackets
   <em class=optionalcode>[like this]</em>.
</P>

<!--===================== DESCRIPTION OF A ROUTINE =====================-->

<A NAME="get_model_size"></A>
<br>
<div class=routine>
<em class=call>model_size = get_model_size( )</em>
<pre>
integer :: <em class=code>get_model_size</em>
</pre>
</div>

<div class=indent1>
<!-- Description -->

<P>
Returns the length of the model state vector as an integer.
This includes all nested domains.
</P>

<TABLE width=100% border=0 summary="" cellpadding=3>

<TR><TD valign=top><em class=code>model_size</em></TD>
    <TD>The length of the model state vector.</TD></TR>

</TABLE>

</div>
<br>

<!--===================== DESCRIPTION OF A ROUTINE =====================-->

<A NAME="get_state_meta_data"></A>
<br>
<div class=routine>
<em class=call>call get_state_meta_data (index_in, location
                          <em class=optionalcode>[,&nbsp;var_type]</em>)</em>
<pre>
integer,             intent(in)  :: <em class=code>index_in</em>
type(location_type), intent(out) :: <em class=code>location</em>
integer, optional,   intent(out) :: <em class=optionalcode>var_type</em>
</pre>
</div>

<div class=indent1>
<!-- Description -->

<P>
Returns metadata about a given element in the DART vector, specifically the location and the quantity.
</P>

<TABLE width=100% border=0 summary="" cellpadding=3>

<TR><TD valign=top><em class=code>index_in&nbsp;&nbsp;&nbsp;</em></TD>
    <TD>Index of state vector element about which information is requested.</TD></TR>

<TR><TD valign=top><em class=code>location</em></TD>
    <TD>the location of the indexed state variable.</TD></TR>

<TR><TD valign=top><em class=optionalcode>var_type</em></TD>
    <TD>Returns the DART QUANTITY (QTY_TEMPERATURE, QTY_U_WIND_COMPONENT, etc.)
        of the indexed state variable as an optional argument.</TD></TR>
</TABLE>

</div>
<br>

<!--===================== DESCRIPTION OF A ROUTINE =====================-->

<A NAME="model_interpolate"></A>
<br>
<div class=routine>
<em class=call> call model_interpolate(state_handle, ens_size, location, obs_type, expected_obs, istatus)</em>
<pre>
type(ensemble_type), intent(in)  :: <em class=code>state_handle</em>
integer,             intent(in)  :: <em class=code>ens_size</em>
type(location_type), intent(in)  :: <em class=code>location</em>
integer,             intent(in)  :: <em class=code>obs_type</em>
real(r8),            intent(out) :: <em class=code>expected_obs(:)</em>
integer,             intent(out) :: <em class=code>istatus(:)</em>
</pre>
</div>

<div class=indent1>
<!-- Description -->

<P>
Given a handle to a model state, a physical location, and a desired QUANTITY; 
<em class=code>model_interpolate</em> returns the array of expected observation values
and a status for each. At present, the ROMS 
<em class=code>model_interpolate</em> is under development as the forward operators for
ROMS are precomputed.
</P>

<TABLE width=100% border=0 summary="" cellpadding=3>

<TR><TD valign=top><em class=code>state_handle</em></TD>
    <TD>A model state vector.</TD></TR>

<TR><TD valign=top><em class=code>ens_size</em></TD>
    <TD>The size of the ensemble.</TD></TR>

<TR><TD valign=top><em class=code>location&nbsp;&nbsp;&nbsp;</em></TD>
    <TD>Physical location of interest.</TD></TR>

<TR><TD valign=top><em class=code>obs_type</em></TD>
    <TD>Integer describing the QUANTITY of interest.
    </TD></TR>

<TR><TD valign=top><em class=code>expected_obs</em></TD>
    <TD>The array of interpolated values from the model.
        The length of the array corresponds to the ensemble size.</TD></TR>

<TR><TD valign=top><em class=code>istatus</em></TD>
    <TD>The array of integer flags indicating the status of the interpolation.
        Each ensemble member returns its own status.</TD></TR>

</TABLE>

</div>
<br>

<!--===================== DESCRIPTION OF A ROUTINE =====================-->

<A NAME="shortest_time_between_assimilations"></A>
<br>
<div class=routine>
<em class=call>var = shortest_time_between_assimilations()</em>
<pre>
type(time_type) :: <em class=code>shortest_time_between_assimilations</em>
</pre>
</div>

<div class=indent1>
<!-- Description -->

<P>
Returns the time step (forecast length) of the model;
the smallest increment in time that the model is capable of
advancing the state in a given implementation.
</P>

<TABLE width=100% border=0 summary="" cellpadding=3>

<TR><TD valign=top><em class=code>var&nbsp;&nbsp;&nbsp;</em></TD>
    <TD>Smallest time step of model.</TD></TR>

</TABLE>

</div>
<br>

<!--===================== DESCRIPTION OF A ROUTINE =====================-->

<A NAME="static_init_model"></A>
<br>
<div class=routine>
<em class=call>call static_init_model()</em>
</div>

<div class=indent1>
<!-- Description -->

<P>
Used for runtime initialization of the model. This is the first call
made to the model by any DART compliant assimilation routine. It reads the
model namelist parameters, set the calendar type (the GREGORIAN calendar 
is used with the ROMS model), and determine the dart vector length.
This routine requires that a <em class=file>roms_input.nc</em> is present
in the working directory to retrieve model information 
(grid dimensions and spacing, variable sizes, etc).
</P>

<TABLE width=100% border=0 summary="" cellpadding=3>
</TABLE>

</div>
<br>

<!--===================== DESCRIPTION OF A ROUTINE =====================-->

<A NAME="end_model"></A>
<br>
<div class=routine>
<em class=call>call end_model( )</em>
</div>

<div class=indent1>
<!-- Description -->

<P>
Called when use of a model is completed to clean up storage, etc.
A stub is provided for the ROMS model.
</P>

</div>
<br>

<!--===================== DESCRIPTION OF A ROUTINE =====================-->

<A NAME="nc_write_model_atts"></A>
<br>
<div class=routine>
<em class=call>ierr = nc_write_model_atts(ncid)</em>
<pre>
integer             :: <em class=code>nc_write_model_atts</em>
integer, intent(in) :: <em class=code>ncid</em>
</pre>
</div>

<div class=indent1>
<!-- Description -->

<P>
Function to write model-specific attributes to a netCDF file, usually 
the DART diagnostic files.  This function writes the model metadata 
to a NetCDF file opened to a file identified by ncid.
</P>

<TABLE width=100% border=0 summary="" cellpadding=3>

<TR><TD valign=top><em class=code>ncid&nbsp;&nbsp;&nbsp;</em></TD>
    <TD>Integer file descriptor to previously-opened netCDF file.</TD></TR>

<TR><TD valign=top><em class=code>ierr</em></TD>
    <TD>Returns a 0 for successful completion.</TD></TR>

</TABLE>

</div>
<br>

<!--===================== DESCRIPTION OF A ROUTINE =====================-->

<A NAME="write_model_time"></A>
<br>
<div class=routine>
<em class=call>call write_model_time(ncid, model_time <em class=optionalcode>[, adv_to_time]</em>)</em>
<pre>
integer,         intent(in)           :: <em class=code>ncid</em>
type(time_type), intent(in)           :: <em class=code>model_time</em>
type(time_type), intent(in), optional :: <em class=optionalcode>adv_to_time</em>
</pre>
</div>

<div class=indent1>
<!-- Description -->

<P>
Routine to write the current model time to the requested netCDF file.
This is used for all the DART diagnostic output.
</P>

<TABLE width=100% border=0 summary="" cellpadding=3>

<TR><TD valign=top><em class=code>ncid&nbsp;&nbsp;&nbsp;</em></TD>
    <TD>Integer file descriptor to an open netCDF file.</TD></TR>

<TR><TD valign=top><em class=code>model_time</em></TD>
    <TD>The current time of the model state.</TD></TR>

<TR><TD valign=top><em class=code>adv_to_time</em></TD>
    <TD>The desired time for the next assimilation.</TD></TR>

</TABLE>

</div>
<br>

<!--===================== DESCRIPTION OF A ROUTINE =====================-->

<A NAME="read_model_time"></A>
<br>
<div class=routine>
<em class=call>var = read_model_time(filename)</em>
<pre>
character(len=*),intent(in)  <em class=code>filename</em>
type(time_type)              <em class=code>var</em>
</pre>
</div>

<div class=indent1>
<!-- Description -->

<P>
Routine to read the model time from a netCDF file that has
not been opened. The file is opened, the time is read, 
and the file is closed.
</P>

<TABLE width=100% border=0 summary="" cellpadding=3>
<TR><TD valign=top><em class=code>ncid&nbsp;&nbsp;&nbsp;</em></TD>
    <TD>the name of the netCDF file.</TD></TR>

<TR><TD valign=top><em class=code>var</em></TD>
    <TD>The current time associated with the file.
        Specifically, this is the last 'ocean_time' (normally).</TD></TR>
</TABLE>

</div>


<!--===================== DESCRIPTION OF A ROUTINE =====================-->

<A NAME="get_time_information"></A>
<br>
<div class=routine>
<em class=call>call get_time_information(filename, ncid, var_name, dim_name &amp;
      <em class=optionalcode>[,myvarid, calendar, last_time_index, last_time, 
      origin_time, all_times]</em>)
</em>

<pre>
character(len=*),            intent(in)  :: <em class=code>filename</em>
integer,                     intent(in)  :: <em class=code>ncid</em>
character(len=*),            intent(in)  :: <em class=code>var_name</em>
character(len=*),            intent(in)  :: <em class=code>dim_name</em>
integer,           optional, intent(out) :: <em class=optionalcode>myvarid</em>
character(len=32), optional, intent(out) :: <em class=optionalcode>calendar</em>
integer,           optional, intent(out) :: <em class=optionalcode>last_time_index</em>
type(time_type),   optional, intent(out) :: <em class=optionalcode>last_time</em>
type(time_type),   optional, intent(out) :: <em class=optionalcode>origin_time</em>
type(time_type),   optional, intent(out) :: <em class=optionalcode>all_times(:)</em>
</pre>
</div>


<div class=indent1>
<!-- Description -->

<P>
Routine to determine the variable describing the time in a ROMS netCDF file. 
</P>

<TABLE width=100% border=0 summary="" cellpadding=3>

<TR><TD valign=top><em class=code>filename&nbsp;&nbsp;&nbsp;</em></TD>
   <TD>The name of the netCDF file. This is used for error messages only.</TD></TR>

<TR><TD valign=top><em class=code>ncid</em></TD>
    <TD>The netCDF file ID - the file must be open.</TD></TR>

<TR><TD valign=top><em class=code>var_name</em></TD>
    <TD>The name of the variable containing the current model time.
        This is usually 'ocean_time'.</TD></TR>

<TR><TD valign=top><em class=code>dim_name</em></TD>
    <TD>The dimension specifying the length of the time variable.
        TJH: This should not be needed and should be removed.</TD></TR>

<TR><TD valign=top><em class=optionalcode>myvarid</em></TD>
    <TD>The netCDF variable ID of <em class=code>var_name</em>.</TD></TR>

<TR><TD valign=top><em class=optionalcode>calendar</em></TD>
    <TD>The type of calendar being used by <em class=code>filename</em>.</TD></TR>

<TR><TD valign=top><em class=optionalcode>last_time_index</em></TD>
    <TD>The index of the last time <em class=code>var_name</em>.
        The last time is declared to be the current model time. </TD></TR>

<TR><TD valign=top><em class=optionalcode>last_time</em></TD>
    <TD>The current model time. </TD></TR>

<TR><TD valign=top><em class=optionalcode>origin_time</em></TD>
    <TD>The time defined in the 'units' attributes of <em class=code>var_name</em>.</TD></TR>

<TR><TD valign=top><em class=optionalcode>all_times</em></TD>
    <TD>The entire array of times in <em class=code>var_name</em>.</TD></TR>

</TABLE>

</div>

<!--===================== DESCRIPTION OF A ROUTINE =====================-->

<A NAME="get_location_from_ijk"></A>
<br>
<div class=routine>
<em class=call>var = get_location_from_ijk(filoc, fjloc, fkloc, quantity, location</em>)
</em>

<pre>
real(r8),            intent(in)  :: <em class=code>filoc</em>
real(r8),            intent(in)  :: <em class=code>fjloc</em>
real(r8),            intent(in)  :: <em class=code>fkloc</em>
integer,             intent(in)  :: <em class=code>quantity</em>
type(location_type), intent(out) :: <em class=code>last_time</em>
</pre>
</div>


<div class=indent1>
<!-- Description -->

<P>
Returns the lat,lon,depth given a fractional i,j,k and a specified quantity
as well as a status.
<br />
<br />
Each grid cell is oriented in a counter clockwise direction
for interpolating locations.  First we interpolate in latitude
and longitude, then interpolate in height.  The height/depth of each grid
cell can very on each interpolation, so care is taken when
we interpolate in the horizontal.  Using the 4 different heights
and lat_frac, lon_frac, hgt_frac we can do a simple trilinear
interpolation to find the location given fractional indicies.
<br />
<br />
var = 10 - bad incoming dart_kind<br />
var = 11 - fkloc out of range<br />
var = 12 - filoc or fjloc out of range for u grid<br />
var = 13 - filoc or fjloc out of range for v grid<br />
var = 14 - filoc or fjloc out of range for rho grid<br />
var = 99 - initalized istatus, this should not happen
</P>

<TABLE border=0 summary="" cellpadding=3>

<TR><TD valign=top><em class=code>filoc</em></TD>
    <TD>Fractional x index.<TD></TR>

<TR><TD valign=top><em class=code>fjloc</em></TD>
    <TD>Fractional y index.<TD></TR>

<TR><TD valign=top><em class=code>fkloc</em></TD>
    <TD>Fractional vertical index.<TD></TR>

<TR><TD valign=top><em class=code>quantity</em></TD>
    <TD>The DART quantity of interest.<TD></TR>

<TR><TD valign=top><em class=code>location&nbsp;&nbsp;&nbsp;&nbsp;</em></TD>
    <TD>The latitude, longitude and depth.</TD></TR>

</TABLE>

</div>

<!--==================================================================-->
<!-- Describe the Files Used by this module.                          -->
<!--==================================================================-->
<P><!-- useless paragraph so 'top' aligns correctly --></P>

<A NAME="FilesUsed"></A>
<div class="top">[<a href="#">top</a>]</div><hr />
<H2>FILES</H2>

<P>These are the files used by the DART side of the experiment.
   Additional files necessary to run ROMS are not listed here.
   Depending on the setting in <em class=file>input.nml</em> for
   <em class=code>stages_to_write</em>, there may be more
   DART diagnostic files output. The inflation files are listed here
   because they may be required for an assimilation.
</P>

<UL> <li>input.nml</li>
     <li>varinfo.dat</li>
     <li>ocean.in</li>
     <li>restart_files.txt</li>
     <li>precomputed_files.txt (optionally)</li>
     <li> input_priorinf_[mean,sd].nc (optionally)</li>
     <li>output_priorinf_[mean,sd].nc (optionally)</li>
     <li>  input_postinf_[mean,sd].nc (optionally)</li>
     <li> output_postinf_[mean,sd].nc (optionally)</li>
</UL>

<!--==================================================================-->
<!-- Cite references, if need be.                                     -->
<!--==================================================================-->
<P><!-- useless paragraph so 'top' aligns correctly --></P>

<A NAME="References"></A>
<div class="top">[<a href="#">top</a>]</div><hr />
<H2>REFERENCES</H2>
<a href="https://www.myroms.org">Regional Ocean Modeling System</a>

<!--==================================================================-->
<!-- Describe all the error conditions and codes.                     -->
<!--==================================================================-->
<P><!-- useless paragraph so 'top' aligns correctly --></P>

<A NAME="Errors"></A>
<div class="top">[<a href="#">top</a>]</div><hr />
<H2>ERROR CODES and CONDITIONS</H2>
<div class=errors>
<TABLE border=1 cellspacing=1 cellpadding=10 width=100%>
<TR><TH>Routine</TH><TH>Message</TH><TH>Comment</TH></TR>

<TR><!-- routine --><TD VALIGN=top>static_init_model <BR>
                                   parse_variable_input <BR></TD>
    <!-- message --><TD VALIGN=top>'model_nml:model "variables" not fully specified'</TD>
    <!-- comment --><TD VALIGN=top>There must be 5 items specified for each variable
    intended to be part of the DART vector.  The <em class=code>variables</em> definition 
    in <em class=file>input.nml</em> does not have 5 items per variable.</TD>
</TR>

<TR><!-- routine --><TD VALIGN=top>static_init_model <BR>
                                   parse_variable_input <BR></TD>
    <!-- message --><TD VALIGN=top>'there is no quantity <...> in obs_kind_mod.f90'</TD>
    <!-- comment --><TD VALIGN=top>An unsupported (or misspelled) DART quantity
    is specified in the  <em class=code>variables</em> definition 
    in <em class=file>input.nml</em>.</TD>
</TR>

<TR><!-- routine --><TD VALIGN=top>get_location_from_ijk </TD>
    <!-- message --><TD VALIGN=top>'Routine not finished.'</TD>
    <!-- comment --><TD VALIGN=top>This routine has not been rigorously tested.</TD>
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
Fully support interpolation in addition to relying on the verification observations.
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
