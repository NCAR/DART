<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN"
          "http://www.w3.org/TR/html4/strict.dtd">
<HTML>
<HEAD>
<TITLE>module model_mod for CM1</TITLE>
<link rel="stylesheet" type="text/css" href="../../docs/html/doc.css" />
<link href="../../docs/images/dart.ico" rel="shortcut icon" />
</HEAD>
<BODY>
<A NAME="TOP"></A>

<H1>MODULE <em class=program>model_mod</em> (for CM1)</H1>

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
<A HREF="#PrivateComponents">PRIVATE COMPONENTS</A> /
<A HREF="#Legalese">TERMS OF USE</A>

<H2>Overview</H2>

<P>
Cloud Model 1 (CM1) version 18 (CM1r18) is compatible with the DART.
CM1 is a non-hydrostatic numerical model in Cartesian 3D coordinates 
designed for the study of micro- to mesoscale atmospheric phenomena 
in idealized to semi-idealized simulations.  The CM1 model was 
developed and is maintained by George Bryan at the National Center 
for Atmospheric Research (NCAR) 
Mesoscale and Microscale Meteorology Laboratory (MMM).  
The model code is freely available from the CM1 website: 
<a href=http://www2.mmm.ucar.edu/people/bryan/cm1/>
http://www2.mmm.ucar.edu/people/bryan/cm1</a>
and must be downloaded and compiled outside of DART.
<br />
<br />
This model interface and scripting support were created by Luke Madaus.
<strong>Thanks Luke!</strong>
<br />
<br />
Several modifications to the CM1 namelist <em class=file>namelist.input</em> 
are required to produce model output compatible with DART. The values are described here
and and example is shown below. Using CM1 output files as a prior 
ensemble state in DART requires each ensemble member to produce a 
restart file in netCDF format (<em class=code>&amp;param9:&nbsp;restart_format=2</em>)
only containing output at the analysis time (<em class=code>&amp;param9:&nbsp;restart_filetype=2</em>). 
The only required state variable to be updated is potential temperature (<em class=code>theta</em>) 
and requires the following namelist setting to ensure this is present in the 
CM1 restart files: (<em class=code>&amp;param9:&nbsp;restart_file_theta&nbsp;=&nbsp;.true.,&nbsp;restart_use_theta&nbsp;=&nbsp;.true.</em>).
<br />
<br />
Additional state variables that have been tested within DART include 
<em class=code>ua, va, wa, ppi, u0, v0, u10, v10, t2, th2, tsk, q2, psfc, qv, qc, qr, qi qs, and qg</em>.
At present, observation times are evaluated relative to the date and time specified in section 
<em class=code>&amp;param11</em>.  
Observation locations are specified in meters relative to the domain origin as defined in 
<em class=code>&amp;param2:&nbsp;iorigin</em>.  
</P>

<div class=namelist>
<pre>
&amp;param9 
   restart_format     = 2         <em class=input>restart needs to be netCDF</em>
   restart_filetype   = 2         <em class=input>restart must be the analysis time - ONLY</em>
   restart_file_theta = .true.    <em class=input>make sure <em class=code>theta</em> is in restart file</em>
   restart_use_theta  = .true.
  /
</pre>
</div>

<H2>About testing CM1 and DART</H2>
<P>There are two sets of scripts in the <em class=file>shell_scripts</em> directory.
Luke contributed a set written in python, and the DART team had a set written in csh.
The csh scripts have not been tested in quite some time, so use with the understanding
that they will need work. Those csh scripts and some unfinished python scripts reside
in a <em class=file>shell_scripts/unfinished</em> directory and should be used with the
understanding that they will need work.
</P>
<H3>Strategy and Instructions for using the python scripts.</H3>
<H4>There are some prerequisites:</H4>
<ol><li>CM1 is required to use netCDF restart files.</li>
    <li>A collection of CM1 model states for initial conditions will be available.</li>
    <li>There is a separate observation sequence file for each assimilation time.</li>
    <li>The DART <em class=file>input.nml</em> file has some required values as defined below.</li>
    <li>Each time CM1 is advanced, it will start from the same filename,
    and the restart number in that filename will be 000001 - ALWAYS.
    That filename will be a link to the most current model state.</li>
</ol>

<H4>Testing a cycling experiment.</H4>
<P>The big picture: three scripts (
<em class=program>setup_filter.py</em>, 
<em class=program>run_filter.py</em>, and
<em class=program>advance_ensemble.py</em>
) are alternated to configure an experiment, perform an assimilation on a 
set of restart files, and make the ensemble forecast. 
Time management is controlled through command-line arguments.
<br />
<br />
It is required that you have generated the DART executables before you test.
The term {centraldir} refers to a filesystem and directory that will be used to
run the experiment, the working directory. {centraldir} should have a lot of capacity, 
as ensemble data assimilation will require lots of disk.
The term {dart_dir} will refer to the location of the DART source code.
<br />
<br />
The data referenced in the directories (the initial ensemble, etc.) are provided
as a compressed tar file 
<a href="http://www.image.ucar.edu/pub/DART/CM1/cm1r18_3member_example_data.tar.gz">cm1r18_3member_example_data.tar.gz</a>
You will have to download the tar file, uncompress it, and modify the scripts to 
use these directories instead of the example directories in the scripts.
You will also have to compile your own cm1 executable.
</P>

<ol>
<li>Set some variables in both 
    <em class=program>shell_scripts/setup_filter.py</em> and
    <em class=program>shell_scripts/advance_ensemble.py</em>
    as described below.</li>
    <br />

<li>In the <em class=file>{dart_dir}/models/cm1/shell_scripts</em> directory, run 
    <em>./setup_filter.py -d YYYYmmDDHHMMSS -i</em>  where <em>YYYYmmDDHHMMSS</em>
    is the date and time of the first assimilation cycle
    (the <em>-i</em> option indicates this is the initial setup and 
    extra work will be performed). 
    This will create the working directory {centraldir}, 
    link in required executables, 
    copy in the initial conditions for each member from some predetermined location, 
    copy in the observation sequence file for this assimilation time from some predetermined location,
    modify namelists, and 
    build a queue submission script in the {centraldir}: <em class=program>run_filter.py</em>.</li>
    <br />

<li>Change into {centraldir} and verify the contents of <em class=program>run_filter.py</em>.
    Ensure the assimilation settings in <em class=file>input.nml</em> are correct.
    Once you are satisfied, submit 
    <em class=program>run_filter.py</em> to the queue to perform an assimilation.</li>
    <br />

<li>After the assimilation job completes, check to be sure that the assimilation completed 
    successfully, and the archived files requested in the <em class=program>setup_filter.py</em>
    <em class=code>files_to_archive</em> variable are in
    <em class=file>{centraldir}/archive/</em><em>YYYYmmDDHHMMSS</em></li>
    <br />

<li>Change into <em class=file>{dart_dir}/models/cm1/shell_scripts</em> and advance the 
    ensemble to the next assimilation time by running
    <em>./advance_ensemble.py -d YYYYmmDDHHMMSS -l nnnn</em> where <em>YYYYmmDDHHMMSS</em>
    is the date of the COMPLETED analysis (the start time for the model) and
    <em>nnnn</em> is the length of model integration in seconds (the forecast length).
    (The forecast length option is specified by 'hypen ell' - the lowercase letter L, 
    not the number one.)
    <em class=program>advance_ensemble.py</em> will submit jobs to the queue to advance the
    ensemble.</li>
    <br />

<li>After all ensemble members have successfully completed, run
    <em>./setup_filter.py -d YYYYmmDDHHMMSS</em> where <em>YYYYmmDDHHMMSS</em> is the 
    <strong>new</strong> current analysis time.  Note the -i flag is NOT used here, 
    as we do not need to (should not!) re-initialize the entire directory structure.</li>
    <br />

<li>Change into {centraldir} and submit the <em class=program>run_filter.py</em> script
    to perform the assimilation.</li>
    <br />

<li>Go back to step 4 and repeat steps 4-7 for each assimilation cycle until the end of 
    the experiment.</li>
</ol>

<P>
Within the <em class=program>setup_filter.py</em> and 
<em class=program>advance_ensemble.py</em> scripts, the following
variables need to be set between the "BEGIN USER-DEFINED VARIABLES" and
"END USER-DEFINED VARIABLES" comment blocks:
</P>

<table>

<TABLE border=0 cellpadding=5 cellspacing=3 width=100% summary='python variable setup'>
<THEAD align=left>
<TR><TH>variable </TH>
    <TH>description </TH> </TR>
</THEAD>
<TBODY valign=top>
<TR><TD><em class=code>jobname</em></TD>
    <TD>A name for this experiment, will be included in the working directory path.</TD></TR>

<TR><TD><em class=code>ens_size</em></TD>
    <TD>Number of ensemble members.</TD></TR>

<TR><TD><em class=code>restart_filename</em></TD>
    <TD>The filename for each ensemble member's restart.
        Highly recommended to leave this as 'cm1out_rst_000001.nc'</TD></TR>

<TR><TD><em class=code>window_mins</em></TD>
    <TD>The assimilation window width (in minutes) for each assimilation cycle.</TD></TR>

<TR><TD><em class=code>copy</em></TD>
    <TD>The copy command with desired flags for this system.</TD></TR>

<TR><TD><em class=code>link</em></TD>
    <TD>The link command with desired flags for this system.</TD></TR>

<TR><TD><em class=code>remove</em></TD>
    <TD>The remove command with desired flags for this system.</TD></TR>

<TR><TD><em class=code>files_to_archive</em></TD>
    <TD>A list of DART output files to archive for each assimilation cycle.
        Note that any inflation files generated are automatically carried over.</TD></TR>

<TR><TD><em class=code>centraldir</em></TD>
    <TD>Directory (which will be created if <em class=program>setup_filter.py</em> 
        is run in intialization mode) where the assimilation and model advances 
        will take place.  Should be on a system with enough space to allow for 
        several assimilation cycles of archived output.
</TD></TR>

<TR><TD><em class=code>dart_dir</em></TD>
    <TD>Path to the cm1 subdirectory of DART.</TD></TR>

<TR><TD><em class=code>cm1_dir</em></TD>
    <TD>Path to the cm1 model executable (<em class=program>cm1.exe</em>)</TD></TR>

<TR><TD><em class=code>icdir</em></TD>
    <TD>Path to the ensemble of initial conditions.  It is assumed that
         within this directory, each ensemble member has a subdirectory
         (<em class=file>m1</em>,
         <em class=file>m2</em>, 
         <em class=file>m3</em>, ...) that contains:
        <ul><li>a restart file for cm1 at the desired start time and having the
                filename defined in <em class=code>restart_filename</em> above</li>
            <li>a <em class=file>namelist.input</em> file compatible with the 
                generation of that restart file.</li>
        </ul>
</TD></TR>

<TR><TD><em class=code>obsdir</em></TD>
    <TD>Path to a directory containing observation sequence files to be
        assimilated.  It is assumed that the observation sequence files
        are named following the convention <em class=file>YYYYmmDDHHMMSS_obs_seq.prior</em>,
        where the date of the analysis time whose observations are contained
        in that file is the first part of the file name.
</TD></TR>

<tr><th colspan=2 align=left>
    <em class=program>setup_filter.py</em> and 
    <em class=program>advance_ensemble.py</em> 
    assume that mpi queue submissions are required 
    to run <em class=program>cm1.exe</em> and <em class=program>filter</em>.
    These variables control how that is handled.</td></th>

<TR><TD><em class=code>queue_system</em></TD>
    <TD>The name of the queueing system</TD></TR>

<TR><TD><em class=code>mpi_run_command</em></TD>
    <TD>The command used in a submitted script to execute an mpi task in the queue, 
        including any required flags </TD></TR>

<TR><TD><em class=code>queue_sub_command</em></TD>
    <TD>The command used to submit a script to the queue</TD></TR>

<TR><TD><em class=code>job_sub_info</em></TD>
    <TD>A dictionary of all flags required to execute a job in the
queue, with the key being the flag and the value being the variable.
e.g. {'-P' : 'PROJECT CODE HERE', '-W' : '00:20'}, etc.</TD></TR>

</table>

<H4 class=indent1>
Additionally, <em class=file>{dart_dir}/work/input.nml</em> should be modified
with the desired assimilation settings. Some of the variables listed above will
override the values in <em class=file>{dart_dir}/work/input.nml</em> should be modified
</H4>

<!--=====================================================================-->
<!--===================== DESCRIPTION OF A NAMELIST =====================-->
<!--=====================================================================-->

<P><!-- dumb spacer --></P>
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
&amp;model_nml 
   assimilation_period_days     = 0
   assimilation_period_seconds  = 21600
   model_perturbation_amplitude = 0.2
   cm1_template_file            = 'null'
   calendar                     = 'Gregorian'
   periodic_x                   = .true.
   periodic_y                   = .true.
   periodic_z                   = .false.
   debug                        = 0
   model_variables              = ' '
  /
</pre>
</div>

<div>
<TABLE border=0 cellpadding=3 width=100% summary='namelist description'>
<THEAD align=left>
<TR><TH> Item </TH>
    <TH> Type </TH>
    <TH> Description </TH> </TR>
</THEAD>
<TBODY valign=top>

<TR><TD>assimilation_period_[days,seconds]&nbsp;&nbsp;</TD>
    <TD>integer</TD>
    <TD>This specifies the width of the assimilation window.
        The current model time is used as the center time of the
        assimilation window. All observations in the assimilation window
        are assimilated. BEWARE: if you put observations that occur before
        the beginning of the assimilation_period, DART will error out because
        it cannot move the model 'back in time' to process these observations.
        </TD>
</TR>

<TR><TD>model_perturbation_amplitude</TD>
    <TD>real(r8)</TD>
    <TD>unsupported
        </TD>
</TR>

<TR><TD>cm1_template_file</TD>
    <TD>character(len=256)</TD>
    <TD>filename used to read the variable sizes, location metadata, etc.
        </TD>
</TR>

<TR><TD>calendar</TD>
    <TD>character(len=256)</TD>
    <TD>Character string to specify the calendar in use. Usually 'Gregorian' 
        (since that is what the observations use).
        </TD>
</TR>
<TR><TD>model_variables</TD>
    <TD>character(:,5)</TD>
    <TD>Strings that identify the CM1 variables, their DART
        quantity, the minimum &amp; maximum possible values, and whether
        or not the posterior values should be written to the output file.
        The DART QUANTITY must be one found in the
        <em class=file>DART/obs_kind/obs_kind_mod.f90</em>
        AFTER it gets built by <em class=program>preprocess</em>.
        <TABLE border=0 cellpadding=3 width=100% summary='variable description'>

        <tr><td valign=top><em class=code>model_variables(:,1)</em></td>
            <td valign=top> Specifies the CM1 variable name in the netCDF file.</td></tr>

        <tr><td valign=top><em class=code>model_variables(:,2)</em></td>
            <td valign=top> Specifies the DART quantity for that variable.</td></tr>

        <tr><td valign=top><em class=code>model_variables(:,3)</em></td>
            <td valign=top> Specifies a minimum bound (if any) for that variable.</td></tr>

        <tr><td valign=top><em class=code>model_variables(:,4)</em></td>
            <td valign=top> Specifies a maximum bound (if any) for that variable.</td></tr>

        <tr><td valign=top><em class=code>model_variables(:,5)</em></td>
            <td valign=top> Specifies if the variable should be updated in the restart file.
                 The value may be "UPDATE" or anything else.</td></tr>
        </TABLE></TD>
</TR>

<TR><TD>periodic_x</TD>
    <TD>logical</TD>
    <TD>a value of <em class=code>.true.</em> means the 'X' dimension is periodic.
        </TD>
</TR>

<TR><TD>periodic_y</TD>
    <TD>logical</TD>
    <TD>a value of <em class=code>.true.</em> means the 'Y' dimension is periodic.
        </TD>
</TR>

<TR><TD>periodic_z</TD>
    <TD>logical</TD>
    <TD>unsupported
        </TD>
</TR>

<TR><TD>debug</TD>
    <TD>integer</TD>
    <TD>switch to control the amount of run-time output is produced. Higher values
        produce more output. 0 produces the least.
        </TD>
</TR>

</TABLE>

<br />
<br />

<P>
<strong>Note:</strong> the values above are the default values. A more realistic (useful?)
example is shown below and closely matches the values in the default <em class=file>input.nml</em>. 
</P>

<div class=namelist>
<pre>
&amp;model_nml 
   assimilation_period_days     = 0
   assimilation_period_seconds  = 60
   cm1_template_file            = 'cm1out_rst_000001.nc'
   calendar                     = 'Gregorian'
   periodic_x                   = .true.
   periodic_y                   = .true.
   periodic_z                   = .false.
   debug                        = 0
   model_variables = 'ua'   , 'QTY_U_WIND_COMPONENT'      , 'NULL', 'NULL', 'UPDATE',
                     'va'   , 'QTY_V_WIND_COMPONENT'      , 'NULL', 'NULL', 'UPDATE',
                     'wa'   , 'QTY_VERTICAL_VELOCITY'     , 'NULL', 'NULL', 'UPDATE',
                     'theta', 'QTY_POTENTIAL_TEMPERATURE' , 0.0000, 'NULL', 'UPDATE',
                     'ppi'  , 'QTY_PRESSURE'              , 'NULL', 'NULL', 'UPDATE',
                     'u10'  , 'QTY_10M_U_WIND_COMPONENT'  , 'NULL', 'NULL', 'UPDATE',
                     'v10'  , 'QTY_10M_V_WIND_COMPONENT'  , 'NULL', 'NULL', 'UPDATE',
                     't2'   , 'QTY_2M_TEMPERATURE'        , 0.0000, 'NULL', 'UPDATE',
                     'th2'  , 'QTY_POTENTIAL_TEMPERATURE' , 0.0000, 'NULL', 'UPDATE',
                     'tsk'  , 'QTY_SURFACE_TEMPERATURE'   , 0.0000, 'NULL', 'UPDATE',
                     'q2'   , 'QTY_SPECIFIC_HUMIDITY'     , 0.0000, 'NULL', 'UPDATE',
                     'psfc' , 'QTY_SURFACE_PRESSURE'      , 0.0000, 'NULL', 'UPDATE',
                     'qv'   , 'QTY_VAPOR_MIXING_RATIO'    , 0.0000, 'NULL', 'UPDATE',
                     'qc'   , 'QTY_CLOUD_LIQUID_WATER'    , 0.0000, 'NULL', 'UPDATE',
                     'qr'   , 'QTY_RAINWATER_MIXING_RATIO', 0.0000, 'NULL', 'UPDATE',
                     'qi'   , 'QTY_CLOUD_ICE'             , 0.0000, 'NULL', 'UPDATE',
                     'qs'   , 'QTY_SNOW_MIXING_RATIO'     , 0.0000, 'NULL', 'UPDATE',
                     'qg'   , 'QTY_GRAUPEL_MIXING_RATIO'  , 0.0000, 'NULL', 'UPDATE'

  /
</pre>
</div>

<!--==================================================================-->

<P><!-- dumb spacer --></P>
<A NAME="Interface"></A>
<div class="top">[<a href="#">top</a>]</div><hr />
<H2>OTHER MODULES USED</H2>
<PRE>
types_mod
time_manager_mod
location_mod (multiple choices here)
utilities_mod
POSSIBLY MANY OTHERS DEPENDING ON MODEL DETAILS
</PRE>

<!--==================================================================-->
<!-- Declare all public entities ...                                  -->
<!-- duplicate public routines template as many times as necessary    -->
<!-- make sure you replace all yyyroutine?? strings                   -->
<!--==================================================================-->
<!--Note to authors. The first row of the table is different.         -->
<!--==================================================================-->

<div class="top">[<a href="#">top</a>]</div><hr />
<H2>PUBLIC INTERFACES</H2>

<TABLE>
<TR><TD><em class=call>use model_mod, only : </em></TD>
                   <TD><A HREF="#get_model_size">get_model_size</A></TD></TR>
<TR><TD>&nbsp;</TD><TD><A HREF="#adv_1step">adv_1step</A></TD></TR>
<TR><TD>&nbsp;</TD><TD><A HREF="#get_state_meta_data">get_state_meta_data</A></TD></TR>
<TR><TD>&nbsp;</TD><TD><A HREF="#model_interpolate">model_interpolate</A></TD></TR>
<TR><TD>&nbsp;</TD><TD><A HREF="#get_model_time_step">get_model_time_step</A></TD></TR>
<TR><TD>&nbsp;</TD><TD><A HREF="#static_init_model">static_init_model</A></TD></TR>
<TR><TD>&nbsp;</TD><TD><A HREF="#init_time">init_time</A></TD></TR>
<TR><TD>&nbsp;</TD><TD><A HREF="#init_conditions">init_conditions</A></TD></TR>
<TR><TD>&nbsp;</TD><TD><A HREF="#nc_write_model_atts">nc_write_model_atts</A></TD></TR>
<TR><TD>&nbsp;</TD><TD><A HREF="#nc_write_model_vars">nc_write_model_vars</A></TD></TR>
<TR><TD>&nbsp;</TD><TD><A HREF="#pert_model_copies">pert_model_copies</A></TD></TR>
<TR><TD>&nbsp;</TD><TD><A HREF="#get_close_maxdist_init">get_close_maxdist_init</A></TD></TR>
<TR><TD>&nbsp;</TD><TD><A HREF="#get_close_obs_init">get_close_obs_init</A></TD></TR>
<TR><TD>&nbsp;</TD><TD><A HREF="#get_close_obs">get_close_obs</A></TD></TR>
<TR><TD>&nbsp;</TD><TD><A HREF="#get_close_state_init">get_close_state_init</A></TD></TR>
<TR><TD>&nbsp;</TD><TD><A HREF="#get_close_state">get_close_state</A></TD></TR>
<TR><TD>&nbsp;</TD><TD><A HREF="#query_vert_localization_coord">query_vert_localization_coord</A></TD></TR>
<TR><TD>&nbsp;</TD><TD><A HREF="#vert_convert">vert_convert</A></TD></TR>
<TR><TD>&nbsp;</TD><TD><A HREF="#read_model_time">read_model_time</A></TD></TR>
<TR><TD>&nbsp;</TD><TD><A HREF="#write_model_time">write_model_time</A></TD></TR>
<TR><TD>&nbsp;</TD><TD><A HREF="#end_model">end_model</A></TD></TR>
</TABLE>

<P>
A namelist interface
<a href="#Namelist"><em class=code>&amp;model_nml</em></a>
may be defined by the module, in which case it will be
read from file <em class=file>input.nml</em>.
The details of the namelist are always model-specific 
(there are no generic namelist values).
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
Returns the length of the model state vector.
Required.
</P>

<TABLE width=100% border=0 summary="" cellpadding=3>

<TR><TD valign=top><em class=code>model_size</em></TD>
    <TD>The length of the model state vector.</TD></TR>

</TABLE>

</div>
<br>

<!--===================== DESCRIPTION OF A ROUTINE =====================-->

<A NAME="adv_1step"></A>
<br>
<div class=routine>
<em class=call>call adv_1step(x, time)</em>
<pre>
real(r8), dimension(:), intent(inout) :: <em class=code>x</em>
type(time_type),        intent(in)    :: <em class=code>time</em>
</pre>
</div>

<div class=indent1>
<!-- Description -->

<P>
Does a single timestep advance of the model. The input value of
the vector x is the starting condition and x must be updated to reflect
the changed state after a timestep. The time argument is intent
in and is used for models that need to know the date/time to
compute a timestep, for instance for radiation computations.
This interface is only called if the namelist parameter
async is set to 0 in <em class=program>perfect_model_obs</em> or 
<em class=program>filter</em> or if the program 
<em class=program>integrate_model</em> is to be used to 
advance the model state as a separate executable.
If one of these options is not going to be used 
(the model will <em>only</em> be advanced as a separate 
model-specific executable), this can be a NULL INTERFACE.  
(The subroutine name must still exist, but it can contain no 
code and it will not be called.)
</P>

<TABLE width=100% border=0 summary="" cellpadding=3>

<TR><TD valign=top><em class=code>x</em></TD>
    <TD>State vector of length model_size.</TD></TR>

<TR><TD valign=top><em class=code>time&nbsp;&nbsp;&nbsp;</em></TD>
    <TD>Specifies time of the initial model state.</TD></TR>

</TABLE>

</div>
<br>

<!--===================== DESCRIPTION OF A ROUTINE =====================-->

<A NAME="get_state_meta_data"></A>
<br>
<div class=routine>
<em class=call>call get_state_meta_data (state_handle, index_in, location, 
                          <em class=optionalcode>[,&nbsp;var_type]</em> )</em>
<pre>
type(ensemble_type), intent(in)  :: <em class=code>state_handle</em>
integer,             intent(in)  :: <em class=code>index_in</em>
type(location_type), intent(out) :: <em class=code>location</em>
integer, optional,   intent(out) :: <em class=optionalcode> var_type </em>
</pre>
</div>

<div class=indent1>
<!-- Description -->

<P>
Given a handle to a state structure and an integer index into the state vector, returns the
associated location. A second intent(out) optional argument 
returns the generic quantity of this item, e.g. QTY_TEMPERATURE,
QTY_DENSITY, QTY_SALINITY, QTY_U_WIND_COMPONENT. 
This interface is required to be functional for all applications.
</P>

<TABLE width=100% border=0 summary="" cellpadding=3>

<TR><TD valign=top><em class=code>state_handle&nbsp;&nbsp;&nbsp;</em></TD>
    <TD>The handle to the state structure containing the state vector
        about which information is requested.</TD></TR>

<TR><TD valign=top><em class=code>index_in</em></TD>
    <TD>Index of state vector element about which information is requested.</TD></TR>

<TR><TD valign=top><em class=code>location</em></TD>
    <TD>The location of state variable element.</TD></TR>

<TR><TD valign=top><em class=optionalcode>var_type</em></TD>
    <TD>The generic kind of the state variable element.</TD></TR>

</TABLE>

</div>
<br>

<!--===================== DESCRIPTION OF A ROUTINE =====================-->

<A NAME="model_interpolate"></A>
<br>
<div class=routine>
<em class=call>call model_interpolate(state_handle, x, location, itype, obs_val, istatus)</em>
<pre>
type(ensemble_type),    intent(in)  :: <em class=code>state_handle</em>
real(r8), dimension(:), intent(in)  :: <em class=code>x</em>
type(location_type),    intent(in)  :: <em class=code>location</em>
integer,                intent(in)  :: <em class=code>itype</em>
real(r8),               intent(out) :: <em class=code>obs_val</em>
integer,                intent(out) :: <em class=code>istatus</em>
</pre>
</div>

<div class=indent1>
<!-- Description -->

<P>
Given a handle containing information for a state vector, a state vector, 
a location, and a model state variable kind
interpolates the state variable field to that location and returns 
the value in <em class=code>obs_val</em>. The <em class=code>istatus</em> 
variable should be returned as 0 unless there is some problem in 
computing the interpolation in which case a positive value should be 
returned. The <em class=code>itype</em> variable
is one of the KIND parameters defined in the
<a href="../../assimilation_code/modules/observations/obs_kind_mod.html">obs_kind_mod.f90</a> file
and defines which generic kind of item is being interpolated.
In low-order models that have no notion of kinds of variables this argument may
be ignored. For applications in which only perfect model experiments
with identity observations (i.e. only the value of a particular
state variable is observed), this can be a NULL INTERFACE.
Otherwise it is required (which is the most common case).
</P>

<TABLE width=100% border=0 summary="" cellpadding=3>

<TR><TD valign=top><em class=code>state_handle&nbsp;&nbsp;&nbsp;</em></TD>
    <TD>The handle to the state structure containing information about 
        the state vector about which information is requested.</TD></TR>

<TR><TD valign=top><em class=code>x</em></TD>
    <TD>A model state vector.</TD></TR>

<TR><TD valign=top><em class=code>location&nbsp;&nbsp;&nbsp;</em></TD>
    <TD>Location to which to interpolate.</TD></TR>

<TR><TD valign=top><em class=code>itype</em></TD>
    <TD>Quantity of state field to be interpolated.</TD></TR>

<TR><TD valign=top><em class=code>obs_val</em></TD>
    <TD>The interpolated value from the model.</TD></TR>

<TR><TD valign=top><em class=code>istatus</em></TD>
    <TD>Integer value returning 0 for success.
        Other values can be defined for various failures.</TD></TR>

</TABLE>

</div>
<br>

<!--===================== DESCRIPTION OF A ROUTINE =====================-->

<A NAME="get_model_time_step"></A>
<br>
<div class=routine>
<em class=call>var = get_model_time_step()</em>
<pre>
type(time_type) :: <em class=code>get_model_time_step</em>
</pre>
</div>

<div class=indent1>
<!-- Description -->

<P>
Returns the time step (forecast length) of the model;
the smallest increment in time that the model is capable of 
advancing the state in a given implementation.
The actual value may be set by the model_mod namelist (depends on the model).
This interface is required for all applications. 
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
Called to do one time initialization of the model. As examples,
might define information about the model size or model timestep.
In models that require pre-computed static data, for instance
spherical harmonic weights, these would also be computed here.
Can be a NULL INTERFACE for the simplest models. 
</P>

<TABLE width=100% border=0 summary="" cellpadding=3>
</TABLE>

</div>
<br>

<!--===================== DESCRIPTION OF A ROUTINE =====================-->

<A NAME="init_time"></A>
<br>
<div class=routine>
<em class=call>call init_time(time)</em>
<pre>
type(time_type), intent(out) :: <em class=code>time</em>
</pre>
</div>

<div class=indent1>
<!-- Description -->

<P>
Companion interface to init_conditions. Returns a time that is somehow
appropriate for starting up a long integration of the model.
At present, this is only used if the <em class=program>perfect_model_obs</em> 
namelist parameter <em class=code>read_input_state_from_file&nbsp;=&nbsp;.false.</em>
If this option should not be used in <em class=program>perfect_model_obs</em>, 
calling this routine should issue a fatal error.
</P>

<TABLE width=100% border=0 summary="" cellpadding=3>

<TR><TD valign=top><em class=code>time&nbsp;&nbsp;&nbsp;</em></TD>
    <TD>Initial model time.</TD></TR>

</TABLE>

</div>
<br>

<!--===================== DESCRIPTION OF A ROUTINE =====================-->

<A NAME="init_conditions"></A>
<br>
<div class=routine>
<em class=call>call init_conditions(x)</em>
<pre>
real(r8), dimension(:), intent(out) :: <em class=code>x</em>
</pre>
</div>

<div class=indent1>
<!-- Description -->

<P>
Returns a model state vector, x, that is some sort of appropriate
initial condition for starting up a long integration of the model.
At present, this is only used if the <em class=program>perfect_model_obs</em> 
namelist parameter <em class=code>read_input_state_from_file&nbsp;=&nbsp;.false.</em>
If this option should not be used in <em class=program>perfect_model_obs</em>, 
calling this routine should issue a fatal error.
</P>

<TABLE width=100% border=0 summary="" cellpadding=3>

<TR><TD valign=top><em class=code>x&nbsp;&nbsp;&nbsp;</em></TD>
    <TD>Initial conditions for state vector.</TD></TR>

</TABLE>

</div>
<br>

<!--===================== DESCRIPTION OF A ROUTINE =====================-->

<A NAME="nc_write_model_atts"></A>
<br>
<div class=routine>
<em class=call>ierr = nc_write_model_atts(ncFileID)</em>
<pre>
integer, intent(in)  :: <em class=code>ncFileID</em>
logical, intent(out) :: <em class=code>model_mod_writes_state_variables</em>
integer              :: <em class=code>nc_write_model_atts</em>
</pre>
</div>

<div class=indent1>
<!-- Description -->

<P>
This routine writes the model-specific attributes to netCDF files
that DART creates.  This includes coordinate variables and any 
metadata, but NOT the actual model state vector.
<em class=code>models/template/model_mod.f90</em> contains code that 
can be used for any model as-is.
<br><br>
The typical sequence for adding new dimensions, variables, attributes:  
</P>

<pre>
NF90_OPEN             ! open existing netCDF dataset               
   NF90_redef         ! put into define mode                       
   NF90_def_dim       ! define additional dimensions (if any)     
   NF90_def_var       ! define variables: from name, kind, and dims
   NF90_put_att       ! assign attribute values                    
NF90_ENDDEF           ! end definitions: leave define mode         
   NF90_put_var       ! provide values for variable                
NF90_CLOSE            ! close: save updated netCDF dataset        
</pre>

<TABLE width=100% border=0 summary="" cellpadding=3>

<TR><TD valign=top><em class=code>ncFileID</em></TD>
    <TD>Integer file descriptor to previously-opened netCDF file.</TD></TR>

<TR><TD valign=top><em class=code>model_mod_writes_state_variables&nbsp;&nbsp;&nbsp;</em></TD>
    <TD>logical descriptor to flag whether or not the <em class=code>model_mod</em>
        will write the actual state variables or the DART-intrinsic routines will write them.
    </TD></TR>

<TR><TD valign=top><em class=code>ierr</em></TD>
    <TD>Returns a 0 for successful completion.</TD></TR>

</TABLE>

</div>
<br>

<!--===================== DESCRIPTION OF A ROUTINE =====================-->

<A NAME="nc_write_model_vars"></A>
<br>
<div class=routine>
<em class=call>ierr = nc_write_model_vars(ncFileID, statevec, copyindex, timeindex)</em>
<pre>
integer                            :: <em class=code>nc_write_model_vars</em>
integer,                intent(in) :: <em class=code>ncFileID</em>
real(r8), dimension(:), intent(in) :: <em class=code>statevec</em>
integer,                intent(in) :: <em class=code>copyindex</em>
integer,                intent(in) :: <em class=code>timeindex</em>
</pre>
</div>

<div class=indent1>
<!-- Description -->

<P>
This routine may be used to write the model-specific state vector (data) to a 
netCDF file. Only used if 
<em class=code>model_mod_writes_state_variables&nbsp;=&nbsp;.true.</em>
<br><br>
Typical sequence for adding new dimensions,variables,attributes:   
</P>

<pre>
NF90_OPEN             ! open existing netCDF dataset               
   NF90_redef         ! put into define mode                       
   NF90_def_dim       ! define additional dimensions (if any)      
   NF90_def_var       ! define variables: from name, kind, and dims
   NF90_put_att       ! assign attribute values                    
NF90_ENDDEF           ! end definitions: leave define mode         
   NF90_put_var       ! provide values for variable                
NF90_CLOSE            ! close: save updated netCDF dataset         
</pre>

<TABLE width=100% border=0 summary="" cellpadding=3>

<TR><TD valign=top><em class=code>ncFileID</em></TD>
    <TD>file descriptor to previously-opened netCDF file.</TD></TR>

<TR><TD valign=top><em class=code>statevec</em></TD>
    <TD>A model state vector.</TD></TR>

<TR><TD valign=top><em class=code>copyindex&nbsp;&nbsp;&nbsp;</em></TD>
    <TD> Integer index of copy to be written.</TD></TR>

<TR><TD valign=top><em class=code>timeindex</em></TD>
    <TD>The timestep counter for the given state.</TD></TR>

<TR><TD valign=top><em class=code>ierr</em></TD>
    <TD>Returns 0 for normal completion.</TD></TR>

</TABLE>

</div>
<br>

<!--===================== DESCRIPTION OF A ROUTINE =====================-->

<A NAME="pert_model_copies"></A>
<br>
<div class=routine>
<em class=call>call pert_model_copies(state_ens_handle, ens_size, pert_amp, interf_provided)</em>
<pre>
type(ensemble_type), intent(inout) :: <em class=code>state_ens_handle</em>
integer,             intent(in)    :: <em class=code>ens_size</em>
real(r8),            intent(in)    :: <em class=code>pert_amp</em>
logical,             intent(out)   :: <em class=code>interf_provided</em>
</pre>
</div>

<div class=indent1>
<!-- Description -->

<P>
Given an ensemble handle, the ensemble size, and a perturbation amplitude; 
perturb the ensemble. Used to generate initial conditions for spinning up 
ensembles. If the <em class=code>model_mod</em> does not want to do this, 
instead allowing the default algorithms in <em class=program>filter</em> 
to take effect, <em class=code>interf_provided&nbsp=&nbps;.false.</em> 
and the routine can be trivial.  Otherwise, <em class=code>interf_provided</em> 
must be returned as <em class=code>.true.</em>
</P>

<TABLE width=100% border=0 summary="" cellpadding=3>

<TR><TD valign=top><em class=code>state_ens_handle</em></TD>
    <TD>handle containing and ensemble of state vectors to be perturbed.</TD></TR>

<TR><TD valign=top><em class=code>ens_size</em></TD>
    <TD>the number of ensemble members to perturb.</TD></TR>

<TR><TD valign=top><em class=code>pert_amp</em></TD>
    <TD>the amplitude of the perturbations. The interpretation is based
        on the model-specific implementation.
    </TD></TR>

<TR><TD valign=top><em class=code>interf_provided&nbsp;&nbsp;&nbsp;</em></TD>
    <TD>Returns false if model_mod cannot do this, else true.</TD></TR>

</TABLE>

</div>
<br>

<!--===================== DESCRIPTION OF A ROUTINE =====================-->

<A NAME="get_close_maxdist_init"></A>
<br>
<div class=routine>
<em class=call>call get_close_maxdist_init(gc, maxdist 
   <em class=optionalcode>[,&nbsp;maxdist_list]</em>)</em>
<pre>
type(get_close_type), intent(inout) :: <em class=code>gc</em>
real(r8),             intent(in)    :: <em class=code>maxdist</em>
real(r8), optional,   intent(in)    :: <em class=optionalcode>maxdist_list(:)</em>
</pre>
</div>

<div class=indent1>
<!-- Description -->

<P>
In distance computations any two locations closer than the
given <em class=code>maxdist</em> will be considered close
by the <em class=code>get_close_obs()</em> routine.
In general this is a PASS-THROUGH ROUTINE.
It is listed on the use line for the locations_mod, and in the 
public list for this module, but has no subroutine declaration 
and no other code in this module:
</P>
<pre>
use location_mod, only: get_close_maxdist_init

public :: get_close_maxdist_init
</pre>
<P>
The location module has code which stores maxdist in the gc derived type.
However, if the model needs to alter the value or wants
to supply an alternative implementation it can
intercept the call like so:
</P>
<pre>
use location_mod, only: &amp;
        lm_get_close_maxdist_init =&#62; get_close_maxdist_init

public :: get_close_maxdist_init
</pre>
<P>
In this case a local <em class=code>get_close_maxdist_init()</em> 
routine must be supplied.  To call the original code in the location
module use:
</P>
<pre>
call lm_get_close_maxdist_init(gc, mymaxdist)
</pre>
<P>
This subroutine will be called before <em class=code>get_close_obs_init</em> 
and <em class=code>get_close_obs</em>.
</P>
<P>
In most cases the PASS-THROUGH ROUTINE will be used.
</P>

<TABLE width=100% border=0 summary="" cellpadding=3>

<TR><TD valign=top><em class=code>gc</em></TD>
    <TD>The get_close_type which stores precomputed information
        about the locations to speed up searching</TD></TR>

<TR><TD valign=top><em class=code>maxdist</em></TD>
    <TD>Anything closer than this will be considered close.</TD></TR>

<TR><TD valign=top><em class=code>maxdist_list(:)&nbsp;&nbsp;&nbsp;</em></TD>
    <TD>An array of distances - one for each specific type. This allows
        different types to have different localization behavior.</TD></TR>

</TABLE>

</TABLE>

</div>
<br>

<!--===================== DESCRIPTION OF A ROUTINE =====================-->

<A NAME="get_close_obs_init"></A>
<br>
<div class=routine>
<em class=call>call get_close_obs_init(gc, num, obs)</em>
<pre>
type(get_close_type), intent(inout) :: <em class=code>gc</em>
integer,              intent(in)    :: <em class=code>num</em>
type(location_type),  intent(in)    :: <em class=code>obs(num)</em>
</pre>
</div>

<div class=indent1>
<!-- Description -->

<P>
In general this is a PASS-THROUGH ROUTINE.  The default routine
in the location module precomputes information to accelerate
the distance computations done by <em class=code>get_close_obs()</em>.
Like the other PASS-THROUGH ROUTINES it is listed on
the use line for the locations_mod, and in the public list
for this module, but has no subroutine declaration and
no other code in this module:
</P>
<pre>
use location_mod, only: get_close_obs_init

public :: get_close_obs_init
</pre>  
<P>
The location module code bins the list of locations and
precomputes maximum possible distances between bins.
However, if the model needs to alter the values or wants
to supply an alternative implementation it can 
intercept the call like so:
</P>
<pre>
use location_mod, only: &amp;
        lm_get_close_obs_init =&#62; get_close_obs_init
        
public :: get_close_obs_init
</pre>
<P>
In this case a local <em class=code>get_close_obs_init()</em>
routine must be supplied.  To call the original code in the location
module use:
</P>
<pre>
call lm_get_close_obs_init(gc, num, obs)
</pre>
<P>
This subroutine will be called after <em class=code>get_close_maxdist_init</em> 
and before <em class=code>get_close_obs</em>.
</P>
<P>
In most cases the PASS-THROUGH ROUTINE will be used.
</P>

<TABLE width=100% border=0 summary="" cellpadding=3>

<TR><TD valign=top><em class=code>gc</em></TD>
    <TD>The get_close_type which stores precomputed information
        about the locations to speed up searching</TD></TR>

<TR><TD valign=top><em class=code>num</em></TD>
    <TD>The number of items in the third argument</TD></TR>

<TR><TD valign=top><em class=code>obs&nbsp;&nbsp;&nbsp;</em></TD>
    <TD>A list of locations which will be part
        of the subsequent distance computations</TD></TR>

</TABLE>

</div>
<br>

<!--===================== DESCRIPTION OF A ROUTINE =====================-->

<A NAME="get_close_obs"></A>
<br>
<div class=routine>
<em class=call>call get_close_obs(gc, base_obs_loc, base_obs_kind,
  obs_loc, obs_kind, num_close, close_ind, dist, state_handle</em>) </em>
<pre>
type(get_close_type), intent(in)  :: <em class=code>gc</em>
type(location_type),  intent(in)  :: <em class=code>base_obs_loc</em>
integer,              intent(in)  :: <em class=code>base_obs_kind</em>
type(location_type),  intent(in)  :: <em class=code>obs_loc(:)</em>
integer,              intent(in)  :: <em class=code>obs_kind(:)</em>
integer,              intent(out) :: <em class=code>num_close</em>
integer,              intent(out) :: <em class=code>close_ind(:)</em>
real(r8),             intent(out) :: <em class=code>dist(:)</em>
type(ensemble_type),  intent(in)  :: <em class=code>state_handle</em>
</pre>
</div>

<div class=indent1>
<!-- Description -->

<P>
Given a location and kind, compute the distances to all other locations 
in the <em class=code>obs</em> list.  The return values are the number
of items which are within maxdist of the base, the index numbers in the 
original obs list, and optionally the distances.  The <em class=code>gc</em>
contains precomputed information to speed the computations.
<br><br>
In general this is a PASS-THROUGH ROUTINE.  It is listed on
the use line for the locations_mod, and in the public list
for this module, but has no subroutine declaration and
no other code in this module:
</P>
<pre>
use location_mod, only: get_close_obs

public :: get_close_obs
</pre>  
<P>
However, if the model needs to alter the values or wants
to supply an alternative implementation it can 
intercept the call like so:
</P>
<pre>
use location_mod, only: &amp;
        lm_get_close_obs =&#62; get_close_obs
        
public :: get_close_obs
</pre>
<P>
In this case a local <em class=code>get_close_obs()</em>
routine must be supplied.  To call the original code in the location
module use:
</P>
<pre>
call lm_get_close_obs(gc, base_obs_loc, ...)
</pre>
<P>
This subroutine will be called after <em class=code>get_close_maxdist_init</em> 
and <em class=code>get_close_obs_init</em>.
<br><br>
In most cases the PASS-THROUGH ROUTINE will be used, but some models need
to alter the actual distances depending on the observation or state vector kind,
or based on the observation or state vector location.
It is reasonable in this case to leave 
<em class=code>get_close_maxdist_init()</em>
and <em class=code>get_close_obs_init()</em> as pass-through routines and 
intercept only <em class=code>get_close_obs()</em>.  The local
<em class=code>get_close_obs()</em> can first call the location mod routine
and let it return a list of values, and then inspect the list and alter
or remove any entries as needed.  See the CAM and WRF model_mod files
for examples of this use.
</P>

<TABLE width=100% border=0 summary="" cellpadding=3>

<TR><TD valign=top><em class=code>gc</em></TD>
    <TD>The get_close_type which stores precomputed information
        about the locations to speed up searching</TD></TR>

<TR><TD valign=top><em class=code>base_obs_loc</em></TD>
    <TD>Reference location.  The distances will be computed
        between this location and every other location in the obs list</TD></TR>

<TR><TD valign=top><em class=code>base_obs_kind&nbsp;&nbsp;&nbsp;</em></TD>
    <TD>The kind of base_obs_loc</TD></TR>

<TR><TD valign=top><em class=code>obs_loc(:)</em></TD>
    <TD>Compute the distance between the base_obs_loc and each
        of the locations in this list</TD></TR>

<TR><TD valign=top><em class=code>obs_kind(:)</em></TD>
    <TD>The corresponding kind of each item in the obs list</TD></TR>

<TR><TD valign=top><em class=code>num_close</em></TD>
    <TD>The number of items from the obs list
        which are within maxdist of the base location</TD></TR>

<TR><TD valign=top><em class=code>close_ind(:)</em></TD>
    <TD>The list of index numbers from the obs list 
        which are within maxdist of the base location</TD></TR>

<TR><TD valign=top><em class=code>dist(:)</em></TD>
    <TD>If present, return the distance between each entry
        in the close_ind list and the base location.  If not
        present, all items in the obs list which are closer
        than maxdist will be added to the list but the overhead
        of computing the exact distances will be skipped.</TD></TR>

<TR><TD valign=top><em class=code>state_handle&nbsp;&nbsp;&nbsp;</em></TD>
    <TD>The handle to the state structure containing information about 
        the state vector about which information is requested.</TD></TR>

</TABLE>

</div>
<br>

<!--===================== DESCRIPTION OF A ROUTINE =====================-->

<A NAME="get_close_state_init"></A>
<br>
<div class=routine>
<em class=call>call get_close_state_init(gc, num, obs)</em>
<pre>
type(get_close_type), intent(inout) :: <em class=code>gc</em>
integer,              intent(in)    :: <em class=code>num</em>
type(location_type),  intent(in)    :: <em class=code>obs(num)</em>
</pre>
</div>

<div class=indent1>
<!-- Description -->

<P>
In general this is a PASS-THROUGH ROUTINE.  The default routine
in the location module precomputes information to accelerate
the distance computations done by <em class=code>get_close_state()</em>.
Like the other PASS-THROUGH ROUTINES it is listed on
the use line for the locations_mod, and in the public list
for this module, but has no subroutine declaration and
no other code in this module:
</P>
<pre>
use location_mod, only: get_close_state_init

public :: get_close_state_init
</pre>  
<P>
The location module code bins the list of locations and
precomputes maximum possible distances between bins.
However, if the model needs to alter the values or wants
to supply an alternative implementation it can 
intercept the call like so:
</P>
<pre>
use location_mod, only: &amp;
        loc_get_close_state_init =&#62; get_close_state_init
        
public :: get_close_state_init
</pre>
<P>
In this case a local <em class=code>get_close_state_init()</em>
routine must be supplied.  To call the original code in the location
module use:
</P>
<pre>
call loc_get_close_state_init(gc, num, obs)
</pre>
<P>
This subroutine will be called after <em class=code>get_close_maxdist_init</em> 
and before <em class=code>get_close_state</em>.
</P>
<P>
In most cases the PASS-THROUGH ROUTINE will be used.
</P>

<TABLE width=100% border=0 summary="" cellpadding=3>

<TR><TD valign=top><em class=code>gc</em></TD>
    <TD>The get_close_type which stores precomputed information
        about the locations to speed up searching</TD></TR>

<TR><TD valign=top><em class=code>num</em></TD>
    <TD>The number of items in the third argument</TD></TR>

<TR><TD valign=top><em class=code>obs&nbsp;&nbsp;&nbsp;</em></TD>
    <TD>A list of locations which will be part
        of the subsequent distance computations</TD></TR>

</TABLE>

</div>
<br>

<!--===================== DESCRIPTION OF A ROUTINE =====================-->

<A NAME="get_close_state"></A>
<br>
<div class=routine>
<em class=call>call get_close_state(gc, base_obs_loc, base_obs_kind,
  state_loc, state_kind, num_close, close_ind, dist, state_handle</em>) </em>
<pre>
type(get_close_type), intent(in)  :: <em class=code>gc</em>
type(location_type),  intent(in)  :: <em class=code>base_obs_loc</em>
integer,              intent(in)  :: <em class=code>base_obs_kind</em>
type(location_type),  intent(in)  :: <em class=code>state_loc(:)</em>
integer,              intent(in)  :: <em class=code>state_kind(:)</em>
integer,              intent(out) :: <em class=code>num_close</em>
integer,              intent(out) :: <em class=code>close_ind(:)</em>
real(r8),             intent(out) :: <em class=code>dist(:)</em>
type(ensemble_type),  intent(in)  :: <em class=code>state_handle</em>
</pre>
</div>

<div class=indent1>
<!-- Description -->

<P>
Given a location and kind, compute the distances to all other locations 
in the <em class=code>state_loc</em> list.  The return values are the number
of items which are within maxdist of the base, the index numbers in the 
original state_loc list, and optionally the distances.  The <em class=code>gc</em>
contains precomputed information to speed the computations.
<br><br>
In general this is a PASS-THROUGH ROUTINE.  It is listed on
the use line for the locations_mod, and in the public list
for this module, but has no subroutine declaration and
no other code in this module:
</P>
<pre>
use location_mod, only: get_close_state

public :: get_close_state
</pre>  
<P>
However, if the model needs to alter the values or wants
to supply an alternative implementation it can 
intercept the call like so:
</P>
<pre>
use location_mod, only: &amp;
        lm_get_close_state =&#62; get_close_state
        
public :: get_close_state
</pre>
<P>
In this case a local <em class=code>get_close_state()</em>
routine must be supplied.  To call the original code in the location
module use:
</P>
<pre>
call loc_get_close_state(gc, base_obs_loc, ...)
</pre>
<P>
This subroutine will be called after <em class=code>get_close_maxdist_init</em> 
and <em class=code>get_close_state_init</em>.
<br><br>
In most cases the PASS-THROUGH ROUTINE will be used, but some models need
to alter the actual distances depending on the observation or state vector kind,
or based on the observation or state vector location.
It is reasonable in this case to leave 
<em class=code>get_close_maxdist_init()</em>
and <em class=code>get_close_state_init()</em> as pass-through routines and 
intercept only <em class=code>get_close_state()</em>.  The local
<em class=code>get_close_state()</em> can first call the location mod routine
and let it return a list of values, and then inspect the list and alter
or remove any entries as needed.  See the CAM and WRF model_mod files
for examples of this use.
</P>

<TABLE width=100% border=0 summary="" cellpadding=3>

<TR><TD valign=top><em class=code>gc</em></TD>
    <TD>The get_close_type which stores precomputed information
        about the locations to speed up searching</TD></TR>

<TR><TD valign=top><em class=code>base_obs_loc</em></TD>
    <TD>Reference location.  The distances will be computed
        between this location and every other location in the obs list</TD></TR>

<TR><TD valign=top><em class=code>base_obs_kind&nbsp;&nbsp;&nbsp;</em></TD>
    <TD>The kind of base_obs_loc</TD></TR>

<TR><TD valign=top><em class=code>state_loc(:)</em></TD>
    <TD>Compute the distance between the base_obs_loc and each
        of the locations in this list</TD></TR>

<TR><TD valign=top><em class=code>state_kind(:)</em></TD>
    <TD>The corresponding kind of each item in the state_loc list</TD></TR>

<TR><TD valign=top><em class=code>num_close</em></TD>
    <TD>The number of items from the state_loc list
        which are within maxdist of the base location</TD></TR>

<TR><TD valign=top><em class=code>close_ind(:)</em></TD>
    <TD>The list of index numbers from the state_loc list 
        which are within maxdist of the base location</TD></TR>

<TR><TD valign=top><em class=code>dist(:)</em></TD>
    <TD>If present, return the distance between each entry
        in the close_ind list and the base location.  If not
        present, all items in the state_loc list which are closer
        than maxdist will be added to the list but the overhead
        of computing the exact distances will be skipped.</TD></TR>

<TR><TD valign=top><em class=code>state_handle&nbsp;&nbsp;&nbsp;</em></TD>
    <TD>The handle to the state structure containing information about 
        the state vector about which information is requested.</TD></TR>

</TABLE>

</div>
<br>

<!--===================== DESCRIPTION OF A ROUTINE =====================-->

<A NAME="query_vert_localization_coord"></A>
<br>
<div class=routine>
<em class=call>ivert = query_vert_localization_coord()</em>
</div>

<div class=indent1><!-- Description -->

<P>
Returns the integer code for the vertical coordinate system used for localization. 
</P>

<TABLE width=100% border=0 summary="" cellpadding=3>
<TR><TD valign=top><em class=code>ivert</em></TD>
    <TD>The integer code for the vertical coordinate system. </TD></TR>
</TABLE>
</div>
<br>

<!--===================== DESCRIPTION OF A ROUTINE =====================-->



<A NAME="vert_convert"></A>
<br>
<div class=routine>
<em class=call>call vert_convert(state_handle, location, obs_kind, istatus)</em>
<pre>
type(ensemble_type), intent(in)  :: <em class=code>state_handle</em>
type(location_type), intent(in)  :: <em class=code>location</em>
integer,             intent(in)  :: <em class=code>obs_kind</em>
integer,             intent(out) :: <em class=code>istatus</em>
</pre>
</div>

<div class=indent1><!-- Description -->

<P>
Converts the state to the desired vertical localization coordinate system.
Some models (toy models with no 'real' observations) will not need this. 
Most (real) models have observations in one or more coordinate systems 
(pressure, height) and the model is generally represented in only one 
coordinate system. To be able to interpolate the model state to the 
observation location, or to compute the true distance between the state 
and the observation, it is necessary to convert everything to one 
coodinate system.
</P>

<TABLE width=100% border=0 summary="" cellpadding=3>

<TR><TD valign=top><em class=code>state_handle&nbsp;&nbsp;&nbsp;</em></TD>
    <TD>The handle to that state.</TD></TR>

<TR><TD valign=top><em class=code>location</em></TD>
    <TD>the desired location</TD></TR>

<TR><TD valign=top><em class=code>obs_kind</em></TD>
    <TD>the quantity defining the part of state to convert.</TD></TR>

<TR><TD valign=top><em class=code>istatus</em></TD>
    <TD>Specifies the success or failure of the vertical conversion.
        If <em class=code>istatus&nbsp;=&nbsp;0</em>, the conversion was
        a sucess. Any other value is a failure.
    </TD></TR>

</TABLE>

</div>
<br>

<!--===================== DESCRIPTION OF A ROUTINE =====================-->

<A NAME="read_model_time"></A>
<br>
<div class=routine>
<em class=call>model_time = read_model_time(filename)</em>
<pre>
character(len=*), intent(in) :: <em class=code>filename</em>
type(time_type)              :: <em class=code>model_time</em>
</pre>
</div>

<div class=indent1><!-- Description -->

<P>
Reads the valid time of the model state in a netCDF file. 
There is a default routine in 
<em class=file>assimilation_code/modules/io/dart_time_io_mod.f90</em>
that can be used as a pass-through. That routine will read the <strong>last</strong>
timestep of a 'time' variable - which is the same strategy used for reading
netCDF files that have multiple timesteps in them.
If your model has some other representation of time (i.e. it does not
use a netCDF variable named 'time') - you will have to write this routine.
</P>

<TABLE width=100% border=0 summary="" cellpadding=3>

<TR><TD valign=top><em class=code>ncid</em></TD>
    <TD>handle to an open netCDF file</TD></TR>

<TR><TD valign=top><em class=code>dart_time&nbsp;&nbsp;&nbsp;</em></TD>
    <TD>Specifies the (last) time of the model state.</TD></TR>

</TABLE>

</div>
<br>

<!--===================== DESCRIPTION OF A ROUTINE =====================-->

<A NAME="write_model_time"></A>
<br>
<div class=routine>
<em class=call>call write_model_time(ncid, dart_time)</em>
<pre>
integer,          intent(in) :: <em class=code>ncid</em>
type(time_type),  intent(in) :: <em class=code>dart_time</em>
</pre>
</div>

<div class=indent1><!-- Description -->

<P>
Writes the assimilation time to a netCDF file. 
There is a default routine in 
<em class=file>assimilation_code/modules/io/dart_time_io_mod.f90</em>
that can be used as a pass-through.
If your model has some other representation of time (i.e. it does not
use a netCDF variable named 'time') - you will have to write this routine.
</P>

<TABLE width=100% border=0 summary="" cellpadding=3>

<TR><TD valign=top><em class=code>ncid</em></TD>
    <TD>handle to an open netCDF file</TD></TR>

<TR><TD valign=top><em class=code>dart_time&nbsp;&nbsp;&nbsp;</em></TD>
    <TD>Specifies the time of the assimilation (the current time step).</TD></TR>

</TABLE>

</div>
<br>


<!--===================== DESCRIPTION OF A ROUTINE =====================-->

<A NAME="end_model"></A>
<br>
<div class=routine>
<em class=call>call end_model()</em>
</div>

<div class=indent1>
<!-- Description -->

<P>
Does any shutdown and clean-up needed for model.
Can be a NULL INTERFACE if the model has no need to clean up storage, etc. 
</P>

<TABLE width=100% border=0 summary="" cellpadding=3>
</TABLE>
</div>
<br>

<!--==================================================================-->
<!-- Describe the Files Used by this module.                          -->
<!--==================================================================-->

<P><!-- dumb spacer --></P>
<A NAME="FilesUsed"></A>
<div class="top">[<a href="#">top</a>]</div><hr />
<H2>FILES</H2>
<UL>
   <LI>Models are free to read and write files as they see fit.
</UL>

<!--==================================================================-->
<!-- Cite references, if need be.                                     -->
<!--==================================================================-->

<P><!-- dumb spacer --></P>
<A NAME="References"></A>
<div class="top">[<a href="#">top</a>]</div><hr />
<H2>REFERENCES</H2>
<ol>
<li>Madaus, L. and G. Hakim (2017): "Constraining Ensemble Forecasts of Discrete Convective Initiation with Surface Observations." Mon. Wea. Rev. doi: 10.1175/MWR-D-16-0395.1</li>
</ol>

<!--==================================================================-->
<!-- Describe all the error conditions and codes.                     -->
<!--==================================================================-->

<P><!-- dumb spacer --></P>
<A NAME="Errors"></A>
<div class="top">[<a href="#">top</a>]</div><hr />
<H2>ERROR CODES and CONDITIONS</H2>

<UL><LI>Models are free to issue calls to the error handler as they see fit. 
        No standard error handler calls are mandated.</LI>
</UL>

<H2>KNOWN BUGS</H2>
<P>
none at this time
</P>

<!--==================================================================-->
<!-- Describe Future Plans.                                           -->
<!--==================================================================-->

<P><!-- dumb spacer --></P>
<A NAME="FuturePlans"></A>
<div class="top">[<a href="#">top</a>]</div><hr />
<H2>FUTURE PLANS</H2>
<P>
It is likely that a number of additional optional interfaces will be
added to the model_mod structure. For instance, hints about how to 
divide the state vector into regions for parallel assimilation will
need to be obtained from the model. It is planned that the interf_provided
mechanism used in pert_model_copies will allow those who do not wish
to support enhanced interfaces to add NULL interfaces by simply 
pasting in an interface block.
</P>

<!--==================================================================-->
<!-- PrivateComponents                                                -->
<!--==================================================================-->

<P><!-- dumb spacer --></P>
<A NAME="PrivateComponents"></A>
<div class="top">[<a href="#">top</a>]</div><hr />
<H2>PRIVATE COMPONENTS</H2>
<P>
N/A
</P>

<!--==================================================================-->
<!-- Legalese & Metadata                                              -->
<!--==================================================================-->

<P><!-- dumb spacer --></P>
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
