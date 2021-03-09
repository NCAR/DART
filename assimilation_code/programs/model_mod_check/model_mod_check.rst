<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN"
          "http://www.w3.org/TR/html4/strict.dtd">
<HTML>
<HEAD>
<TITLE>program model_mod_check</TITLE>
<link rel="stylesheet" type="text/css" href="../../../docs/html/doc.css" />
<link href="../../../docs/images/dart.ico" rel="shortcut icon" />
</HEAD>
<BODY>
<A NAME="TOP"></A>

<H1>program <em class=program>model_mod_check</em></H1>

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
<A HREF="#Usage">USAGE </A> / 
<A HREF="#References">REFERENCES</A> /
<A HREF="#Errors">ERRORS</A> /
<A HREF="#FuturePlans">PLANS</A> /
<A HREF="#Legalese">TERMS OF USE</A>

<H2>Overview</H2>

<P>
   <em class="program">model_mod_check</em> tests some of the more
   fundamental routines in any <em class="program">model_mod</em>.
   This is intended to be used when adding a new model to DART - 
   test the pieces as they are written.  As such, this program is
   meant to be hacked up and customized to your own purpose. Right now,
   it reads in model netCDF file(s) - one per domain/nest/whatever -  
   and writes out files, queries the metdata, etc. It also exercises 
   <em class="program">static_init_model()</em>, which is the first routine
   to get right ...
</P>

<!--==================================================================-->
<!--=================== DESCRIPTION OF A NAMELIST  ===================-->
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
&amp;model_mod_check 
   num_ens               = 1
   single_file           = .FALSE.
   input_state_files     = 'null'
   output_state_files    = 'null'
   all_metadata_file     = 'metadata.txt'

   test1thru             = 7
   run_tests             = -1

   x_ind                 = -1
   loc_of_interest       = -1.0, -1.0, -1.0
   quantity_of_interest  = 'NONE'

   interp_test_dlon      = 10.0
   interp_test_dlat      = 10.0
   interp_test_dvert     = 10.0

   interp_test_lonrange  = 0.0, 120.0
   interp_test_latrange  = 0.0, 120.0
   interp_test_vertrange = 0.0, 100.0

   interp_test_dx        = -888888.0
   interp_test_dy        = -888888.0
   interp_test_dz        = -888888.0

   interp_test_xrange    = -888888.0, -888888.0
   interp_test_yrange    = -888888.0, -888888.0
   interp_test_zrange    = -888888.0, -888888.0

   interp_test_vertcoord = 'VERTISHEIGHT'
   verbose               = .FALSE.
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

<TR><TD> num_ens </TD>
    <TD> integer </TD>
    <TD> Provided for future use. Must be 1. Ultimately, The number of 
         ensemble members you would like to read in for testing.
</TD></TR>  

<TR><TD> single_file </TD>
    <TD> logical </TD>
    <TD> If .TRUE. all members are stored in a single restart file.
</TD></TR>  

<TR><TD> input_state_files(:)  </TD>
    <TD> character(len=256) </TD>
    <TD> The name(s) of the NetCDF file(s) containing the model states, one per domain.
         If num_ens &gt; 1 and not single_file, specify a filename for each 
         ensemble member (num_ens). If you have both multiple ensemble members in separate
         files AND multiple domains, specify all the ensemble member filenames for domain 1, 
         then all the ensemble member filenames for domain 2, etc.
</TD></TR>  

<TR><TD> output_state_files(:)  </TD>
    <TD> character(len=256) </TD>
    <TD> The name(s) of the output NetCDF file(s) for testing IO, one per domain.
         If num_ens &gt; 1 and not single_file, specify a filename for each 
         ensemble member (num_ens). If you have both multiple ensemble members in separate
         files AND multiple domains, specify all the ensemble member filenames for domain 1, 
         then all the ensemble member filenames for domain 2, etc.
</TD></TR>  

<TR><TD> all_metadata_file  </TD>
    <TD> character(len=256) </TD>
    <TD> Test 6 produces an exhaustive list of metadata for EVERY element in 
         the DART state vector. The metadata get written to this file name.
</TD></TR>  

<TR><TD> x_ind </TD>
    <TD> integer(i8) </TD>
    <TD> An integer index into the DART state vector.  This will be used to test the
         metadata routines.  Answers questions about location, what variable type is
         stored there, etc.
</TD></TR>  

<TR><TD> loc_of_interest </TD>
    <TD> real(r8), dimension(3) </TD>
    <TD> The lat/lon/level for a <b>particular</b> location.  
         Used in Test 4, the single-point interpolation test.
         Indirectly tests the routine to find the closest gridpoint.
</TD></TR>  

<TR><TD> quantity_of_interest </TD>
    <TD> character(len=32) </TD>
    <TD> Specifies the QUANTITY of the model state to use in Tests 4, 5, and 7.
</TD></TR>

<TR><TD> interp_test_dlon </TD>
    <TD> real(r8) </TD>
    <TD> The distance (measured in degrees) on the longitude interpolation grid.
         Ignored if interpolating with cartesian coordinates. Used in Test 5.
</TD></TR>  

<TR><TD> interp_test_dlat </TD>
    <TD> real(r8) </TD>
    <TD> The distance (measured in degrees) on the latitude interpolation grid.
         Ignored if interpolating with cartesian coordinates. Used in Test 5.
</TD></TR>  

<TR><TD> interp_test_dvert </TD>
    <TD> real(r8) </TD>
    <TD> The distance (measured in interp_vertcoord) on the vertical interpolation grid.
         Ignored if interpolating with cartesian coordinates. Used in Test 5.
</TD></TR>  

<TR><TD> interp_test_lonrange </TD>
    <TD> real(r8) </TD>
    <TD> The range of y to be tested with model_interpolate, with spacing
         <em class=code>interp_test_dlon</em>. Ignored if interpolating with cartesian coordinates. Used in Test 5.
</TD></TR>  

<TR><TD> interp_test_latrange </TD>
    <TD> real(r8) </TD>
    <TD> The range of y to be tested with model_interpolate, with spacing
         <em class=code>interp_test_dlat</em>. Ignored if interpolating with cartesian coordinates. Used in Test 5.
</TD></TR>  

<TR><TD> interp_test_vertrange </TD>
    <TD> real(r8) </TD>
    <TD> The range in the vertical direction to be tested with model_interpolate, with spacing
         <em class=code>interp_test_dvert</em>. Ignored if interpolating with cartesian coordinates. Used in Test 5.
</TD></TR>  

<TR><TD> interp_test_dx </TD>
    <TD> real(r8) </TD>
    <TD> The interval on the x axis of the interpolation grid.
         This is used in Test 5 for models with threed_cartesian coordinates.
</TD></TR>  

<TR><TD> interp_test_dy </TD>
    <TD> real(r8) </TD>
    <TD> The interval on the y axis of the interpolation grid.
         This is used in Test 5 for models with threed_cartesian coordinates.
</TD></TR>  

<TR><TD> interp_test_dz </TD>
    <TD> real(r8) </TD>
    <TD> The interval on the z axis of the interpolation grid.
         This is used in Test 5 for models with threed_cartesian coordinates.
</TD></TR>  

<TR><TD> interp_test_xrange </TD>
    <TD> real(r8) </TD>
    <TD> The range of x to be tested with model_interpolate in Test 5, with spacing <em class=code>interp_test_dx</em>.
</TD></TR>  

<TR><TD> interp_test_yrange </TD>
    <TD> real(r8) </TD>
    <TD> The range of y to be tested with model_interpolate in Test 5, with spacing <em class=code>interp_test_dy</em>.
</TD></TR>  

<TR><TD> interp_test_zrange </TD>
    <TD> real(r8) </TD>
    <TD> The range in the vertical direction to be tested with model_interpolate in Test 5, with spacing <em class=code>interp_test_dz</em>.
</TD></TR>

<TR><TD> interp_test_vertcoord </TD>
    <TD> character(len=32) </TD>
    <TD> Specifies the vertical coordinate system to use during the interpolation tests.
         Valid values are: 'VERTISHEIGHT','VERTISPRESSURE','VERTISLEVEL', and 'VERTISSCALEHEIGHT'.
</TD></TR>  

<TR><TD> test1thru </TD>
    <TD> integer </TD>
    <TD> If <em class=code>test1thru &gt; 0</em>, specifies the last test to be performed.
         All tests get performed sequentially. If <em class=code>test1thru &lt; 0</em>,
         <em class=code>run_tests</em> is used to specify the tests to perform.
        <table><tr><th>test&nbsp;</th><th align=left>summary</th></tr>
        <tr><td>0 </td><td>Mandatory. Tests <em class=program>static_init_model()</em> by calling
                           <em class=program>static_init_assim_model()</em>. Reads
                    <em class=file>input.nml</em> <em class=code>&amp;model_nml</em> </td></tr>
        <tr><td>1 </td><td>Tests <em class=program>get_model_size()</em> and reports on the makeup
                           of the DART state vector.</td></tr>
        <tr><td>2 </td><td>Reads and writes a restart file.</td></tr>
        <tr><td>3 </td><td>Tests <em class=program>get_state_meta_data()</em>
                          for a single index into the DART state.
                           Helps determine if the state vector is constructed correctly.</td></tr>
        <tr><td>4 </td><td>Tests <em class=program>model_interpolate()</em> for a single point.</td></tr>
        <tr><td>5 </td><td>Tests <em class=program>model_interpolate()</em> for a range of interpolation points.</td></tr>
        <tr><td>6 </td><td>Long, expensive test to return the metadata for every element of the state vector.
                            May be useful to decide on known locations for subsequent testing.</td></tr>
        <tr><td>7 </td><td>Find the closest gridpoint to a known location.</td></tr>
         </table>
</TD></TR>

<TR><TD> run_tests(:) </TD>
    <TD> integer </TD>
    <TD> Specifies a list of tests to be performed. Same test numbers as described in test1thru.
         There are some dependencies. Tests 4 and 5 require a valid model state - which is read by Test 2.
         If a required test is not specified, the required test is enabled and run.
         A value of -1 means that <em class=code>test1thru</em> will be used.
</TD></TR>  

<TR><TD>   verbose   </TD>
    <TD>   logical   </TD>
    <TD>Print extra info about the <em class=program>model_mod_check</em> run.
        This is only used for more reporting during Test 5. Be warned - it will generate
        several lines of output for each point in the test!
</TD></TR>


</TBODY> 
</TABLE>
</div>

<P>A more typical namelist for a single ensemble member for a model with an outer grid and a 
   single nested grid is shown below.</P>
<div class=namelist>
<pre>
&amp;model_mod_check_nml
   input_state_files     = 'dart_vector1.nc','dart_vector2.nc'
   output_state_files    = 'check_me1.nc', 'check_me2.nc'
   all_metadata_file     = 'metadata.txt'
   verbose               = .TRUE.
   test1thru             = 5
   run_tests             = -1
   loc_of_interest       = 243.72386169, 52.78578186, 10.0
   x_ind                 = 12666739
   quantity_of_interest  = 'QTY_POTENTIAL_TEMPERATURE'
   interp_test_lonrange  = 144.0, 326.0
   interp_test_dlon      = 1.0
   interp_test_latrange  = -5.0, 80.0
   interp_test_dlat      = 1.0
   interp_test_vertrange = 100.0, 11000.0
   interp_test_dvert     = 200.0
   interp_test_vertcoord = 'VERTISHEIGHT'
  /
</pre>
</div>

<P><!-- simple divider for 'top' --></P>

<!--==================================================================-->

<A NAME="Modules"></A>
<div class="top">[<a href="#">top</a>]</div><hr />
<H2>OTHER MODULES USED</H2>
<PRE>
assimilation_code/location/threed_sphere/location_mod.f90
assimilation_code/location/utilities/default_location_mod.f90
assimilation_code/location/utilities/location_io_mod.f90
assimilation_code/modules/assimilation/adaptive_inflate_mod.f90
assimilation_code/modules/assimilation/assim_model_mod.f90
assimilation_code/modules/assimilation/assim_tools_mod.f90
assimilation_code/modules/assimilation/cov_cutoff_mod.f90
assimilation_code/modules/assimilation/filter_mod.f90
assimilation_code/modules/assimilation/obs_model_mod.f90
assimilation_code/modules/assimilation/quality_control_mod.f90
assimilation_code/modules/assimilation/reg_factor_mod.f90
assimilation_code/modules/assimilation/sampling_error_correction_mod.f90
assimilation_code/modules/assimilation/smoother_mod.f90
assimilation_code/modules/io/dart_time_io_mod.f90
assimilation_code/modules/io/direct_netcdf_mod.f90
assimilation_code/modules/io/io_filenames_mod.f90
assimilation_code/modules/io/state_structure_mod.f90
assimilation_code/modules/io/state_vector_io_mod.f90
assimilation_code/modules/observations/forward_operator_mod.f90
assimilation_code/modules/observations/obs_kind_mod.f90
assimilation_code/modules/observations/obs_sequence_mod.f90
assimilation_code/modules/utilities/distributed_state_mod.f90
assimilation_code/modules/utilities/ensemble_manager_mod.f90
assimilation_code/modules/utilities/netcdf_utilities_mod.f90
assimilation_code/modules/utilities/null_mpi_utilities_mod.f90
assimilation_code/modules/utilities/null_win_mod.f90
assimilation_code/modules/utilities/obs_impact_mod.f90
assimilation_code/modules/utilities/options_mod.f90
assimilation_code/modules/utilities/parse_args_mod.f90
assimilation_code/modules/utilities/random_seq_mod.f90
assimilation_code/modules/utilities/sort_mod.f90
assimilation_code/modules/utilities/time_manager_mod.f90
assimilation_code/modules/utilities/types_mod.f90
assimilation_code/modules/utilities/utilities_mod.f90
assimilation_code/programs/model_mod_check/model_mod_check.f90
models/<em class=input>your_model_here</em>/model_mod.f90
models/model_mod_tools/<em class=input>test_interpolate_threed_sphere.f90</em>
models/model_mod_tools/model_check_utilities_mod.f90
models/utilities/default_model_mod.f90
observations/forward_operators/obs_def_mod.f90
observations/forward_operators/obs_def_utilities_mod.f90
</PRE>
<P>Items highlighted may change based on which model is being tested.</P>

<!--==================================================================-->
<!-- Describe the Files Used by this module.                          -->
<!--==================================================================-->

<A NAME="FilesUsed"></A>
<div class="top">[<a href="#">top</a>]</div><hr />
<H2>FILES</H2>
<UL><LI><em class="file">input.nml</em> is used for 
        <em class="code">model_mod_check_nml</em></LI>

    <LI>The <em class="code">"input_state_files" </em> can either be a
        single file containing multiple restart files, or a single
        NetCDF restart file. One file per domain.</LI>

    <LI>The <em class="code">"output_state_files"</em> is the output netCDF
        files from Test 2. Check the attributes, values, etc.</LI>

    <LI><em class="file">check_me_interptest.nc </em> and
        <em class="file">check_me_interptest.m </em> are the result of
        Test 5.
        </LI>

    <LI><em class=file>"all_metadata_file"</em> is the run-time output of Test 6.</LI>

</UL>

<!--==================================================================-->
<!-- Discuss  typical usage of model_mod_check.                              -->
<!--==================================================================-->

<A NAME="Usage"></A>
<div class="top">[<a href="#">top</a>]</div><hr />
<H2>USAGE</H2>

<P>
Normal circumstances indicate that you are trying to put a new model into
DART, so to be able to build and run <em class="program">model_mod_check</em>,
you will need to create a <em class="file">path_names_model_mod_check</em>
file with the following contents:
</P>
<pre>
assimilation_code/location/threed_sphere/location_mod.f90
assimilation_code/location/utilities/default_location_mod.f90
assimilation_code/location/utilities/location_io_mod.f90
assimilation_code/modules/assimilation/adaptive_inflate_mod.f90
assimilation_code/modules/assimilation/assim_model_mod.f90
assimilation_code/modules/assimilation/assim_tools_mod.f90
assimilation_code/modules/assimilation/cov_cutoff_mod.f90
assimilation_code/modules/assimilation/filter_mod.f90
assimilation_code/modules/assimilation/obs_model_mod.f90
assimilation_code/modules/assimilation/quality_control_mod.f90
assimilation_code/modules/assimilation/reg_factor_mod.f90
assimilation_code/modules/assimilation/sampling_error_correction_mod.f90
assimilation_code/modules/assimilation/smoother_mod.f90
assimilation_code/modules/io/dart_time_io_mod.f90
assimilation_code/modules/io/direct_netcdf_mod.f90
assimilation_code/modules/io/io_filenames_mod.f90
assimilation_code/modules/io/state_structure_mod.f90
assimilation_code/modules/io/state_vector_io_mod.f90
assimilation_code/modules/observations/forward_operator_mod.f90
assimilation_code/modules/observations/obs_kind_mod.f90
assimilation_code/modules/observations/obs_sequence_mod.f90
assimilation_code/modules/utilities/distributed_state_mod.f90
assimilation_code/modules/utilities/ensemble_manager_mod.f90
assimilation_code/modules/utilities/netcdf_utilities_mod.f90
assimilation_code/modules/utilities/null_mpi_utilities_mod.f90
assimilation_code/modules/utilities/null_win_mod.f90
assimilation_code/modules/utilities/obs_impact_mod.f90
assimilation_code/modules/utilities/options_mod.f90
assimilation_code/modules/utilities/parse_args_mod.f90
assimilation_code/modules/utilities/random_seq_mod.f90
assimilation_code/modules/utilities/sort_mod.f90
assimilation_code/modules/utilities/time_manager_mod.f90
assimilation_code/modules/utilities/types_mod.f90
assimilation_code/modules/utilities/utilities_mod.f90
assimilation_code/programs/model_mod_check/model_mod_check.f90
models/<em class=input>your_model_here</em>/model_mod.f90
models/model_mod_tools/<em class=input>test_interpolate_threed_sphere.f90</em>
models/utilities/default_model_mod.f90
observations/forward_operators/obs_def_mod.f90
observations/forward_operators/obs_def_utilities_mod.f90
</pre>
as well as a <em class="file">mkmf_model_mod_check</em> script.
You should be able to look at any other <em class="file">mkmf_xxxx</em> 
script and figure out what to change. Once they exist:
<br />
<br />
<div class="unix">
<pre>
[~/DART/models/yourmodel/work] % <em class="input">csh mkmf_model_mod_check</em>
[~/DART/models/yourmodel/work] % <em class="input">make</em>
[~/DART/models/yourmodel/work] % <em class="input">./model_mod_check</em>
</pre>
</div>

<P>
Unlike other DART components, you are expected
to modify <em class="file">model_mod_check.f90</em> to suit your needs as
you develop your <em class="program">model_mod</em>. The code is roughly 
divided into the following categories:
</P>
<ol><li>Check the geometry information, </li>
    <li>Read/write a restart file, </li>
    <li>Check the construction of the state vector ... i.e. the metadata, </li>
    <li>Interpolate at a single point, </li>
    <li>Interpolate for a range of points.</li>
</ol>

<H3 class=indent1>Test 0. Mandatory.</H3>
<P>
The first test in <em class="program">model_mod_check</em>
reads the namelist and runs
<em class="program">static_init_model</em> - which generally sets the
geometry of the grid, the number of state variables and their shape, etc.
Virtually everything requires knowledge of the grid and state vector,
so this block cannot be skipped.
</P>

<H3 class=indent1>Test 1. Checking the Geometry Information:</H3>
<P>
The first test in <em class="program">model_mod_check</em> exercises
a basic required interface <em class="program">get_model_size()</em>.
This also generates a report on the geometry of the grid, the number of
state variables and their shape, etc. as well as the total number of
elements in the DART state vector.
</P>

<H3 class=indent1>Test 2. Read/writing a restart file:</H3>
<P>
This directly reads and write state variables from the model netCDF file.
This is a nice sanity check to make sure that the DART state vector is
being read in properly.
</P>

<H3 class=indent1>Test 3. Check the construction of the state vector:</H3>
<P>It is critical to return the correct metadata for any given index into
the DART state vector. This code block tests the two most common features of
the metadata. As a bonus, this routine is also quite useful to determine
EXACTLY where to place your first test observation. If you test precisely at
a grid location, you should be able to really get a handle on debugging your
<em class="program">model_interpolate()</em> routine. 
</P>

<H3 class=indent1>Test 4. Test interpolation on a single point.</H3>
<P>
This tests your model's interpolation routine on a single point and
returns the interpolated value.
This requires that Test 2 works - it needs a valid model state with data.
Test 2 is automatically run if this test is selected.
</P>

<H3 class=indent1>Test 5. Test interpolation on a range of values.</H3>
<P>
This tests your model's interpolation routine on a range of values
returns the interpolated grid in <em class=file>check_me_interptest.nc</em>
and <em class=file>check_me_interptest.m</em> which can be read in Matlab
and used to visualize the result.
This requires that Test 2 works - it needs a valid model state with data.
Test 2 is automatically run if this test is selected.
</P>

<H3 class=indent1>Test 6. Exhaustively test the construction of the state vector.</H3>
<P>
This can be a long test, depending on the size of your state vector.
This returns the same data as in Test 3 - but <em class=strong> for every element </em>
in the state vector. The metadata are written to a file specified by <em class=code>all_metadata_file</em>
and <em class=file>check_me_interptest.m</em> which can be read in Matlab
and used to visualize the result.
</P>

<H3 class=indent1>Test 7. Find the closest gridpoint to a test location.</H3>
<P>
This is a good test to verify that <em class=progtram>get_state_meta_data()</em> 
and the grid information are correct. Typically, one would put in a location that 
is actually <strong>on</strong> the grid and see if the correct gridpoint index 
is returned. Repeat the test with slightly different locations until
the next gridpoint is closer. Repeat ...
</P>

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
<P>There are no error conditions to check. This program is intended
to demonstrate simple checks that will allow you to proceed with
improving and testing the <em class="program">model_mod</em>. There
will be plenty of run-time errors, I suggest compiling your code with 
"bounds checking" turned on - at a minimum.
</P>
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
<P>Expanded instructions on how to add a model, and how to methodically
test piece-by-piece.
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
