<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
    "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0" />
<title>module model_mod (MPAS OCN)</title>
<link rel="stylesheet" type="text/css" href=
"../../docs/html/doc.css" />
<link href="../../docs/images/dart.ico" rel="shortcut icon" />
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>MPAS OCN</h1>
<table border="0" summary="" cellpadding="5">
<tr>
<td valign="middle"><img src="../../docs/images/Dartboard7.png"
alt="DART project logo" height="70" /></td>
<td>Jump to <a href="../../docs/index.html">DART Documentation Main
Index</a></td>
</tr>
</table>
<a href="#Namelist">NAMELIST</a> / <a href=
"#Interface">INTERFACES</a> / <a href="#FilesUsed">FILES</a> /
<a href="#References">REFERENCES</a> / <a href="#Errors">ERRORS</a>
/ <a href="#FuturePlans">PLANS</a> / <a href=
"#PrivateComponents">PRIVATE COMPONENTS</a> / <a href=
"#Legalese">TERMS OF USE</a>
<h2>Overview</h2>
<p>The <strong>MPAS OCN</strong> interface for <strong>Data
Assimilation Research Testbed (DART)</strong> is under
development.<br />
<br />
Since MPAS OCN uses netcdf files for their restart mechanism, a
namelist-controlled set of variables is used to build the DART
state vector. Each variable must also correspond to a DART "KIND";
required for the DART interpolate routines. For example:</p>
<pre>&amp;mpas_vars_nml
   mpas_state_variables = 'uReconstructZonal',      'QTY_U_WIND_COMPONENT',
                          'uReconstructMeridional', 'QTY_V_WIND_COMPONENT',
                          'w',                      'QTY_VERTICAL_VELOCITY',
                          'theta',                  'QTY_POTENTIAL_TEMPERATURE',
                          'qv',                     'QTY_VAPOR_MIXING_RATIO',
                          'qc',                     'QTY_CLOUDWATER_MIXING_RATIO',
                          'qr',                     'QTY_RAINWATER_MIXING_RATIO',
                          'qi',                     'QTY_ICE_MIXING_RATIO',
                          'qs',                     'QTY_SNOW_MIXING_RATIO',
                          'qg',                     'QTY_GRAUPEL_MIXING_RATIO',
                          'surface_pressure',       'QTY_SURFACE_PRESSURE'
   /
</pre>
<p>These variables are then adjusted to be consistent with
observations and stuffed back into the same netCDF analysis files.
Since DART is an ensemble algorithm, there are multiple analysis
files for a single analysis time: one for each ensemble member.
Creating the initial ensemble of states is an area of active
research.<br />
<br />
DART reads grid information from the MPAS OCN 'history' file, I
have tried to keep the variable names the same. Internal to the
DART code, the following variables exist:</p>
<table width="100%" cellpadding="3">
<tr>
<td valign="top">integer :: nCells</td>
<td>the number of Cell Centers</td>
</tr>
<tr>
<td valign="top">integer :: nEdges</td>
<td>the number of Cell Edges</td>
</tr>
<tr>
<td valign="top">integer :: nVertices</td>
<td>the number of Cell Vertices</td>
</tr>
<tr>
<td valign="top">integer :: nVertLevels</td>
<td>the number of vertical level midpoints</td>
</tr>
<tr>
<td valign="top">integer :: nVertLevelsP1</td>
<td>the number of vertical level edges</td>
</tr>
<tr>
<td valign="top">integer :: nSoilLevels</td>
<td>the number of soil level ?midpoints?</td>
</tr>
<tr>
<td valign="top">real(r8) :: latCell(:)</td>
<td>the latitudes of the Cell Centers (-90,90)</td>
</tr>
<tr>
<td valign="top">real(r8) :: lonCell(:)</td>
<td>the longitudes of the Cell Centers [0, 360)</td>
</tr>
<tr>
<td valign="top">real(r8) :: zgrid(:,:)</td>
<td>cell center geometric height at cell centers
(ncells,nvert)</td>
</tr>
<tr>
<td valign="top">integer :: CellsOnVertex(:,:)</td>
<td>list of cell centers defining a triangle</td>
</tr>
</table>
<h2>model_mod variable storage</h2>
<p><em class="file">input.nml</em><em class=
"code">&amp;mpas_vars_nml</em> defines the list of MPAS variables
used to build the DART state vector. Combined with an MPAS analysis
file, the information is used to determine the size of the DART
state vector and derive the metadata. To keep track of what
variables are contained in the DART state vector, an array of a
user-defined type called "progvar" is available with the following
components:</p>
<div class="unix">
<pre>
type progvartype
   private
   character(len=NF90_MAX_NAME) :: varname
   character(len=NF90_MAX_NAME) :: long_name
   character(len=NF90_MAX_NAME) :: units
   character(len=NF90_MAX_NAME), dimension(NF90_MAX_VAR_DIMS) :: dimname
   integer, dimension(NF90_MAX_VAR_DIMS) :: dimlens
   integer :: xtype         ! netCDF variable type (NF90_double, etc.) 
   integer :: numdims       ! number of dims - excluding TIME
   integer :: numvertical   ! number of vertical levels in variable
   integer :: numcells      ! number of horizontal locations (typically cell centers)
   logical :: ZonHalf       ! vertical coordinate has dimension nVertLevels
   integer :: varsize       ! prod(dimlens(1:numdims))
   integer :: index1        ! location in dart state vector of first occurrence
   integer :: indexN        ! location in dart state vector of last  occurrence
   integer :: dart_kind
   character(len=paramname_length) :: kind_string
   logical  :: clamping     ! does variable need to be range-restricted before 
   real(r8) :: range(2)     ! being stuffed back into MPAS analysis file.
end type progvartype

type(progvartype), dimension(max_state_variables) :: progvar
</pre></div>
<p>The variables are simply read from the MPAS analysis file and
stored in the DART state vector such that all quantities for one
variable are stored contiguously. Within each variable; they are
stored vertically-contiguous for each horizontal location. From a
storage standpoint, this would be equivalent to a Fortran variable
dimensioned x(nVertical,nHorizontal,nVariables). The
fastest-varying dimension is vertical, then horizontal, then
variable ... naturally, the DART state vector is 1D. Each variable
is also stored this way in the MPAS analysis file.</p>
<div class="indent1">
<h2>The DART interface for MPAS (atm)</h2>
<p>was compiled with the gfortran 4.2.3 compilers and run on a
Mac.<br />
<br />
The DART components were built with the following <em class=
"file">mkmf.template</em> settings:</p>
<pre>
      FC = gfortran
      LD = gfortran
      NETCDF = /Users/thoar/GNU
      INCS = -I${NETCDF}/include
      LIBS = -L${NETCDF}/lib -lnetcdf -lcurl -lhdf5_hl -lhdf5 -lz  -lm
      FFLAGS = -O0 -fbounds-check -frecord-marker=4 -ffpe-trap=invalid $(INCS)
      LDFLAGS = $(FFLAGS) $(LIBS)
   </pre></div>
<a name="conversions" id="conversions"></a>
<div class="indent1">
<h3>Converting between DART files and MPAS analysis files</h3>
<p>is relatively straighforward. Given the namelist mechanism for
determining the state variables and the MPAS history netCDF files
exist, - everything that is needed is readily determined.<br />
<br />
There are two programs - both require the list of MPAS variables to
use in the DART state vector: the <em class=
"code">mpas_vars_nml</em> namelist in the <em class=
"file">input.nml</em> file. The MPAS file name being read and/or
written is - in all instances - specified by the <em class=
"code">model_nml:model_analysis_filename</em> variable in the
<em class="file">input.nml</em> namelist file.</p>
<table width="100%" cellpadding="10">
<tr>
<td valign="top"><a href=
"model_to_dart.html">model_to_dart.f90</a></td>
<td>converts an MPAS analysis file (nominally named <em class=
"file">mpas_analysis.nc</em>) into a DART-compatible file normally
called <em class="file">dart_ics</em> . We usually wind up
linking the actual analysis file to a static name that is used by
DART.</td>
</tr>
<tr>
<td valign="top"><a href=
"dart_to_model.f90">dart_to_model.f90</a></td>
<td>inserts the DART output into an existing MPAS analysis netCDF
file by overwriting the variables in the analysis netCDF file.
There are two different types of DART output files, so there is a
namelist option to specify if the DART file has two time records or
just one (if there are two, the first one is the 'advance_to' time,
followed by the 'valid_time' of the ensuing state). <em class=
"program">dart_to_model</em> updates the MPAS analysis file
specified in <em class="file">input.nml</em><em class=
"code">model_nml:model_analysis_filename</em>. If the DART file
contains an 'advance_to' time, separate control information is
written to an auxiliary file that is used by the <em class=
"file">advance_model.csh</em> script.</td>
</tr>
</table>
</div>
<!-- In namelist.input for the mpas model, we do have start_time, run_duration, and output_interval as below.

  config_start_time = '2010-10-23_03:00:00'  
  config_run_duration = '00_00:15:00'
  config_output_interval = '00_00:15:00'

Here the run_duration and output_interval do not have 'YYYY-MM', but the rest of them are the same as in start_time.
Do you think that is enough for the advance_model script? -->
<div class="indent1">
<p>The header of an MPAS analysis file is presented below - simply
for context. Keep in mind that <strong>many</strong> variables have
been removed for clarity. Also keep in mind that the
multi-dimensional arrays listed below have the dimensions reversed
from the Fortran convention.</p>
</div>
<div class="unix">
<pre>366 mirage2:thoar% ncdump -h mpas_analysis.nc
netcdf mpas_analysis {
dimensions:
        StrLen = 64 ;
        Time = UNLIMITED ; // (1 currently)
        <em>nCells</em> = 10242 ;                                  <em>available in DART</em>
        <em>nEdges</em> = 30720 ;                                  <em>available in DART</em>
        maxEdges = 10 ;
        maxEdges2 = 20 ;
        <em>nVertices</em> = 20480 ;                               <em>available in DART</em>
        TWO = 2 ;
        THREE = 3 ;
        <em>vertexDegree</em> = 3 ;                                <em>available in DART</em>
        FIFTEEN = 15 ;
        TWENTYONE = 21 ;
        R3 = 3 ;
        <em>nVertLevels</em> = 41 ;                                <em>available in DART</em>
        <em>nVertLevelsP1</em> = 42 ;                              <em>available in DART</em>
        nMonths = 12 ;
        nVertLevelsP2 = 43 ;
        <em>nSoilLevels</em> = 4 ;                                 <em>available in DART</em>
variables:
        char <em>xtime(Time, StrLen)</em> ;                        <em>available in DART</em>
        double <em>latCell(nCells)</em> ;                          <em>available in DART</em>
        double <em>lonCell(nCells)</em> ;                          <em>available in DART</em>
        double latEdge(nEdges) ;
        double lonEdge(nEdges) ;
        int indexToEdgeID(nEdges) ;
        double latVertex(nVertices) ;
        double lonVertex(nVertices) ;
        int indexToVertexID(nVertices) ;
        int cellsOnEdge(nEdges, TWO) ;
        int nEdgesOnCell(nCells) ;
        int nEdgesOnEdge(nEdges) ;
        int edgesOnCell(nCells, maxEdges) ;
        int edgesOnEdge(nEdges, maxEdges2) ;
        double weightsOnEdge(nEdges, maxEdges2) ;
        double dvEdge(nEdges) ;
        double dcEdge(nEdges) ;
        double angleEdge(nEdges) ;
        double edgeNormalVectors(nEdges, R3) ;
        double cellTangentPlane(nEdges, TWO, R3) ;
        int cellsOnCell(nCells, maxEdges) ;
        int verticesOnCell(nCells, maxEdges) ;
        int verticesOnEdge(nEdges, TWO) ;
        int edgesOnVertex(nVertices, vertexDegree) ;
        int <em>cellsOnVertex(nVertices, vertexDegree)</em> ;      <em>available in DART</em>
        double kiteAreasOnVertex(nVertices, vertexDegree) ;
        double rainc(Time, nCells) ;
        double cuprec(Time, nCells) ;
        double cutop(Time, nCells) ;
        double cubot(Time, nCells) ;
        double relhum(Time, nCells, nVertLevels) ;
        double qsat(Time, nCells, nVertLevels) ;
        double graupelnc(Time, nCells) ;
        double snownc(Time, nCells) ;
        double rainnc(Time, nCells) ;
        double graupelncv(Time, nCells) ;
        double snowncv(Time, nCells) ;
        double rainncv(Time, nCells) ;
        double sr(Time, nCells) ;
        double surface_temperature(Time, nCells) ;
        double surface_pressure(Time, nCells) ;
        double coeffs_reconstruct(nCells, maxEdges, R3) ;
        double theta_base(Time, nCells, nVertLevels) ;
        double rho_base(Time, nCells, nVertLevels) ;
        double pressure_base(Time, nCells, nVertLevels) ;
        double exner_base(Time, nCells, nVertLevels) ;
        double exner(Time, nCells, nVertLevels) ;
        double h_divergence(Time, nCells, nVertLevels) ;
        double uReconstructMeridional(Time, nCells, nVertLevels) ;
        double uReconstructZonal(Time, nCells, nVertLevels) ;
        double uReconstructZ(Time, nCells, nVertLevels) ;
        double uReconstructY(Time, nCells, nVertLevels) ;
        double uReconstructX(Time, nCells, nVertLevels) ;
        double pv_cell(Time, nCells, nVertLevels) ;
        double pv_vertex(Time, nVertices, nVertLevels) ;
        double ke(Time, nCells, nVertLevels) ;
        double rho_edge(Time, nEdges, nVertLevels) ;
        double pv_edge(Time, nEdges, nVertLevels) ;
        double vorticity(Time, nVertices, nVertLevels) ;
        double divergence(Time, nCells, nVertLevels) ;
        double v(Time, nEdges, nVertLevels) ;
        double rh(Time, nCells, nVertLevels) ;
        double theta(Time, nCells, nVertLevels) ;
        double rho(Time, nCells, nVertLevels) ;
        double qv_init(nVertLevels) ;
        double t_init(nCells, nVertLevels) ;
        double u_init(nVertLevels) ;
        double pressure_p(Time, nCells, nVertLevels) ;
        double tend_theta(Time, nCells, nVertLevels) ;
        double tend_rho(Time, nCells, nVertLevels) ;
        double tend_w(Time, nCells, nVertLevelsP1) ;
        double tend_u(Time, nEdges, nVertLevels) ;
        double qv(Time, nCells, nVertLevels) ;
        double qc(Time, nCells, nVertLevels) ;
        double qr(Time, nCells, nVertLevels) ;
        double qi(Time, nCells, nVertLevels) ;
        double qs(Time, nCells, nVertLevels) ;
        double qg(Time, nCells, nVertLevels) ;
        double tend_qg(Time, nCells, nVertLevels) ;
        double tend_qs(Time, nCells, nVertLevels) ;
        double tend_qi(Time, nCells, nVertLevels) ;
        double tend_qr(Time, nCells, nVertLevels) ;
        double tend_qc(Time, nCells, nVertLevels) ;
        double tend_qv(Time, nCells, nVertLevels) ;
        double qnr(Time, nCells, nVertLevels) ;
        double qni(Time, nCells, nVertLevels) ;
        double tend_qnr(Time, nCells, nVertLevels) ;
        double tend_qni(Time, nCells, nVertLevels) ;

</pre></div>
<!--==================================================================-->
<!--=================== DESCRIPTION OF A NAMELIST  ===================-->
<!--==================================================================-->
<a name="Namelist" id="Namelist"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr />
<h2>NAMELIST</h2>
<p>We adhere to the F90 standard of starting a namelist with an
ampersand '&amp;' and terminating with a slash '/' for all our
namelist input. Consider yourself forewarned that character strings
that contain a '/' must be enclosed in quotes to prevent them from
prematurely terminating the namelist.</p>
<div class="namelist">
<pre>
<em class=
"call">namelist /model_nml/ </em> model_analysis_filename, &amp;
          assimilation_period_days, assimilation_period_seconds, &amp;
          model_perturbation_amplitude, output_state_vector, calendar, debug
</pre></div>
<div class="indent1"><!-- Description -->
<p>This namelist is read in a file called <em class=
"file">input.nml</em>. This namelist provides control over the
assimilation period for the model. All observations within (+/-)
half of the assimilation period are assimilated. The assimilation
period is the minimum amount of time the model can be advanced, and
checks are performed to ensure that the assimilation window is a
multiple of the model dynamical timestep. This also specifies the
MPAS analysis file that will be read and/or written by the
different program units.</p>
<table border="0" cellpadding="3" width="100%">
<tr>
<th align="left">Contents</th>
<th align="left">Type</th>
<th align="left">Description</th>
</tr>
<tr><!--contents-->
<td valign="top">model_analysis_filename</td>
<!--  type  -->
<td valign="top">character(len=256)<br />
<em class="units">[default: 'mpas_analysis.nc']</em></td>
<!--descript-->
<td valign="top">Character string specifying the name of the MPAS
analysis file to be read and/or written by the different program
units.</td>
</tr>
<tr><!--contents-->
<td valign="top">output_state_vector</td>
<!--  type  -->
<td valign="top">logical <em class=
"units">[default: .true.]</em></td>
<!--descript-->
<td valign="top">The switch to determine the form of the state
vector in the output netCDF files. If <em class="code">.true.</em>
the state vector will be output exactly as DART uses it ... one
long array. If <em class="code">.false.</em>, the state vector is
parsed into prognostic variables and output that way -- much easier
to use with 'ncview', for example.</td>
</tr>
<tr><!--contents-->
<td valign="top">assimilation_period_days</td>
<!--  type  -->
<td valign="top">integer <em class=
"units">[default: 1]</em></td>
<!--descript-->
<td valign="top">The number of days to advance the model for each
assimilation.</td>
</tr>
<tr><!--contents-->
<td valign="top">assimilation_period_seconds</td>
<!--  type  -->
<td valign="top">integer <em class=
"units">[default: 0]</em></td>
<!--descript-->
<td valign="top">In addition to <em class=
"code">assimilation_period_days</em>, the number of seconds to
advance the model for each assimilation.</td>
</tr>
<tr><!--contents-->
<td valign="top">model_perturbation_amplitude</td>
<!--  type  -->
<td valign="top">real(r8) <em class=
"units">[default: 0.2]</em></td>
<!--descript-->
<td valign="top">Reserved for future use. 
<!-- The amount of noise to add when trying to perturb a single
                       state vector to create an ensemble. Only used when 
<em class=file>input.nml</em><em class=code>&amp;filter_nml:start_from_restart = .false.</em>
                       See also 
                       <a href="#InitialEnsemble">Generating the initial ensemble</a> 
                       at the start of this document. units: standard deviation 
                       of a gaussian distribution with the mean at the value of 
                       the state vector element. --></td>
</tr>
<tr><!--contents-->
<td valign="top">calendar</td>
<!--  type  -->
<td valign="top">character(len=32)<br />
<em class="units">[default: 'Gregorian']</em></td>
<!--descript-->
<td valign="top">Character string specifying the calendar being
used by MPAS.</td>
</tr>
<tr><!--contents-->
<td valign="top">debug</td>
<!--  type  -->
<td valign="top">integer <em class=
"units">[default: 0]</em></td>
<!--descript-->
<td valign="top">The switch to specify the run-time verbosity.
<em class="code">0</em> is as quiet as it gets. <em class=
"code">&gt; 1</em> provides more run-time messages. <em class=
"code">&gt; 5</em> provides ALL run-time messages.</td>
</tr>
</table>
<h3>Example namelist</h3>
<pre>
&amp;model_nml
   model_analysis_filename      = 'mpas_restart.nc';
   assimilation_period_days     = 0,
   assimilation_period_seconds  = 60,
   model_perturbation_amplitude = 0.2,
   output_state_vector          = .true.,
   calendar                     = 'Gregorian',
   debug                        = 0
   /
</pre></div>
<br />
<!--==================================================================-->
 <a name="mpas_vars_nml" id="mpas_vars_nml"></a><br />
<div class="namelist">
<pre>
<em class="call">namelist /mpas_vars_nml/</em> mpas_state_variables
</pre></div>
<div class="indent1"><!-- Description -->
<p>This namelist is read from <em class="file">input.nml</em> and
contains the list of MPAS variables that make up the DART state
vector.</p>
<table border="0" cellpadding="3" width="100%">
<tr>
<th align="left">Contents</th>
<th align="left">Type</th>
<th align="left">Description</th>
</tr>
<tr><!--contents-->
<td valign="top">mpas_vars_nml</td>
<!--  type  -->
<td valign="top">character(len=NF90_MAX_NAME)::<br />
dimension(160) <em class="units">[default:  see
example]</em></td>
<!--descript-->
<td valign="top">The table that relates the GITM variables to use
to build the DART state vector, and the corresponding DART kinds
for those variables.</td>
</tr>
</table>
<h3 class="indent1">Example</h3>
<p>The following mpas_vars_nml is just for demonstration purposes.
You application will likely involve a different DART state
vector.</p>
<pre>
&amp;mpas_vars_nml
   mpas_state_variables = 'theta',                 'QTY_POTENTIAL_TEMPERATURE',
                          'uReconstructZonal',     'QTY_U_WIND_COMPONENT',
                          'uReconstructMeridional','QTY_V_WIND_COMPONENT',
                          'w',                     'QTY_VERTICAL_VELOCITY',
                          'qv',                    'QTY_VAPOR_MIXING_RATIO',
                          'qc',                    'QTY_CLOUDWATER_MIXING_RATIO',
                          'qr',                    'QTY_RAINWATER_MIXING_RATIO',
                          'qi',                    'QTY_ICE_MIXING_RATIO',
                          'qs',                    'QTY_SNOW_MIXING_RATIO',
                          'qg',                    'QTY_GRAUPEL_MIXING_RATIO'
                          'surface_pressure',      'QTY_SURFACE_PRESSURE'
   /
</pre>
<p>The variables are simply read from the MPAS analysis file and
stored in the DART state vector such that all quantities for one
variable are stored contiguously. Within each variable; they are
stored vertically-contiguous for each horizontal location. From a
storage standpoint, this would be equivalent to a Fortran variable
dimensioned x(nVertical,nHorizontal,nVariables). The
fastest-varying dimension is vertical, then horizontal, then
variable ... naturally, the DART state vector is 1D. Each variable
is also stored this way in the MPAS analysis file.</p>
</div>
<br />
<!--==================================================================-->
 <a name="OtherModulesUsed" id="OtherModulesUsed"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr />
<h2>OTHER MODULES USED</h2>
<pre>
types_mod
time_manager_mod
threed_sphere/location_mod
utilities_mod
obs_kind_mod
mpi_utilities_mod
random_seq_mod
</pre>
<div class="warning">
<div class="title">
<p>Warning</p>
</div>
<p>DAReS staff began creating the MPAS_OCN interface to DART in
preparation for the model's inclusion as the ocean component of the
Community Earth System Model (CESM). The plans for including
MPAS_OCN in CESM were abandoned and the Modular Ocean Model version
6 (MOM6) was included instead. Thus, the documentation on this page
after this point describes an incomplete interface. Please contact
DAReS staff by emailing <a href=
"mailto:dart@ucar.edu">dart@ucar.edu</a> if you want to use DART
with MPAS_OCN.</p>
</div>
<!--==================================================================-->
<!-- Note to authors. The first row of the table is different.        -->
<!--==================================================================-->
<!-- Declare all public entities ...                                  -->
<!-- duplicate public routines template as many times as necessary    -->
<!-- make sure you replace all yyyroutine?? strings                   -->
<!--==================================================================-->
<a name="Interface" id="Interface"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr />
<h2>PUBLIC INTERFACES</h2>
<p>Only a select number of interfaces used are discussed here. Each
module has its own discussion of their routines.</p>
<h3 class="indent1">Required Interface Routines</h3>
<table width="50%">
<tr>
<td><em class="call">use model_mod, only :</em></td>
<td><a href="#get_model_size">get_model_size</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#adv_1step">adv_1step</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_state_meta_data">get_state_meta_data</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#model_interpolate">model_interpolate</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_model_time_step">get_model_time_step</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#static_init_model">static_init_model</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#end_model">end_model</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#init_time">init_time</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#init_conditions">init_conditions</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#nc_write_model_atts">nc_write_model_atts</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#nc_write_model_vars">nc_write_model_vars</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#pert_model_state">pert_model_state</a></td>
</tr>
<tr>
<td> </td>
<td><a href=
"#get_close_maxdist_init">get_close_maxdist_init</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_close_obs_init">get_close_obs_init</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_close_obs">get_close_obs</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#ens_mean_for_model">ens_mean_for_model</a></td>
</tr>
</table>
<h3 class="indent1">Unique Interface Routines</h3>
<table width="50%">
<tr>
<td><em class="call">use model_mod, only :</em></td>
<td><a href="#get_gridsize">get_gridsize</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#restart_file_to_sv">restart_file_to_sv</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#sv_to_restart_file">sv_to_restart_file</a></td>
</tr>
<tr>
<td> </td>
<td><a href=
"#get_gitm_restart_filename">get_gitm_restart_filename</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_base_time">get_base_time</a></td>
</tr>
<tr>
<td> </td>
<td><a href="#get_state_time">get_state_time</a></td>
</tr>
</table>
<table>
<tr>
<td><em class="call">use location_mod, only :</em></td>
<td><a href=
"../../assimilation_code/location/threed_sphere/location_mod.html#get_close_obs">
get_close_obs</a></td>
</tr>
</table>
<p>A note about documentation style. Optional arguments are
enclosed in brackets <em class="optionalcode">[like this]</em>.</p>
<!--==================================================================-->
<h3 class="indent1">Interface routine descriptions</h3>
<!--==================================================================-->
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
<a name="get_model_size" id="get_model_size"></a><br />
<div class="routine"><em class="call">model_size = get_model_size(
)</em>
<pre>
integer :: <em class="code">get_model_size</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Returns the length of the model state vector. Required.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">model_size</em></td>
<td>The length of the model state vector.</td>
</tr>
</table>
</div>
<br />
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="adv_1step" id="adv_1step"></a><br />
<div class="routine"><em class="call">call adv_1step(x, time)</em>
<pre>
real(r8), dimension(:), intent(inout) :: <em class="code">x</em>
type(time_type),        intent(in)    :: <em class="code">time</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p><em class="code">adv_1step</em> is not used for the gitm model.
Advancing the model is done through the <em class=
"program">advance_model</em> script. This is a NULL_INTERFACE,
provided only for compatibility with the DART requirements.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">x</em></td>
<td>State vector of length model_size.</td>
</tr>
<tr>
<td valign="top"><em class="code">time   </em></td>
<td>Specifies time of the initial model state.</td>
</tr>
</table>
</div>
<br />
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_state_meta_data" id="get_state_meta_data"></a><br />
<div class="routine"><em class="call">call get_state_meta_data
(index_in, location, <em class=
"optionalcode">[, var_type]</em> )</em>
<pre>
integer,             intent(in)  :: <em class="code">index_in</em>
type(location_type), intent(out) :: <em class="code">location</em>
integer, optional,   intent(out) :: <em class=
"optionalcode"> var_type </em>
</pre></div>
<div class="indent1"><!-- Description -->
<p><em class="code">get_state_meta_data</em> returns metadata about
a given element of the DART representation of the model state
vector. Since the DART model state vector is a 1D array and the
native model grid is multidimensional, <em class=
"code">get_state_meta_data</em> returns information about the
native model state vector representation. Things like the
<em class="code">location</em>, or the type of the variable (for
instance: temperature, u wind component, ...). The integer values
used to indicate different variable types in <em class=
"code">var_type</em> are themselves defined as public interfaces to
model_mod if required.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class=
"code">index_in   </em></td>
<td>Index of state vector element about which information is
requested.</td>
</tr>
<tr>
<td valign="top"><em class="code">location</em></td>
<td>Returns the 3D location of the indexed state variable. The
<em class="code">location_ type</em> comes from <em class=
"file">DART/assimilation_code/location/threed_sphere/location_mod.f90</em>.
Note that the lat/lon are specified in degrees by the user but are
converted to radians internally.</td>
</tr>
<tr>
<td valign="top"><em class="optionalcode">var_type</em></td>
<td>Returns the type of the indexed state variable as an optional
argument. The type is one of the list of supported observation
types, found in the block of code starting <em class=
"code">! Integer definitions for DART TYPES</em>
in <em class=
"file">DART/assimilation_code/modules/observations/obs_kind_mod.f90</em></td>
</tr>
</table>
<p>The list of supported variables in <em class=
"file">DART/assimilation_code/modules/observations/obs_kind_mod.f90</em>
is created by <em class="program">preprocess</em>.</p>
</div>
<br />
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="model_interpolate" id="model_interpolate"></a><br />
<div class="routine"><em class="call">call model_interpolate(x,
location, itype, obs_val, istatus)</em>
<pre>
real(r8), dimension(:), intent(in)  :: <em class="code">x</em>
type(location_type),    intent(in)  :: <em class=
"code">location</em>
integer,                intent(in)  :: <em class="code">itype</em>
real(r8),               intent(out) :: <em class=
"code">obs_val</em>
integer,                intent(out) :: <em class=
"code">istatus</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Given a model state, <em class="code">model_interpolate</em>
returns the value of the desired observation type (which could be a
state variable) that would be observed at the desired location. The
interpolation method is either completely specified by the model,
or uses some standard 2D or 3D scalar interpolation routines. Put
another way, <em class="code">model_interpolate</em> will apply the
forward operator <strong>H</strong> to the model state to create an
observation at the desired location.<br />
<br />
If the interpolation is valid, <em class="code">istatus = 0</em>.
In the case where the observation operator is not defined at the
given location (e.g. the observation is below the lowest model
level, above the top level, or 'dry'), interp_val is returned as
0.0 and istatus = 1.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">x</em></td>
<td>A model state vector.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">location   </em></td>
<td>Location to which to interpolate.</td>
</tr>
<tr>
<td valign="top"><em class="code">itype</em></td>
<td>Integer indexing which type of observation is desired.</td>
</tr>
<tr>
<td valign="top"><em class="code">obs_val</em></td>
<td>The interpolated value from the model.</td>
</tr>
<tr>
<td valign="top"><em class="code">istatus</em></td>
<td>Integer flag indicating the success of the interpolation.<br />
success == 0, failure == anything else</td>
</tr>
</table>
</div>
<br />
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_model_time_step" id="get_model_time_step"></a><br />
<div class="routine"><em class="call">var =
get_model_time_step()</em>
<pre>
type(time_type) :: <em class="code">get_model_time_step</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p><em class="code">get_model_time_step</em> returns the forecast
length to be used as the "model base time step" in the filter. This
is the minimum amount of time the model can be advanced by
<em class="program">filter</em>. <em class="strong">This is also
the assimilation window</em>. All observations within (+/-) one
half of the forecast length are used for the assimilation. In the
<em class="program">GITM</em> case, this is set from the namelist
values for <em class="file">input.nml</em><em class=
"code">&amp;model_nml:assimilation_period_days,
assimilation_period_seconds</em>.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">var   </em></td>
<td>Smallest time step of model.</td>
</tr>
</table>
</div>
<br />
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="static_init_model" id="static_init_model"></a><br />
<div class="routine"><em class="call">call
static_init_model()</em></div>
<div class="indent1"><!-- Description -->
<p><em class="code">static_init_model</em> is called for runtime
initialization of the model. The namelists are read to determine
runtime configuration of the model, the grid coordinates, etc.
There are no input arguments and no return values. The routine sets
module-local private attributes that can then be queried by the
public interface routines.<br />
<br />
See the GITM documentation for all namelists in <em class=
"file">gitm_in</em> . Be aware that DART reads the GITM
<em class="code">&amp;grid_nml</em> namelist to get the filenames
for the horizontal and vertical grid information as well as the
topography information.<br />
<br />
The namelists (all mandatory) are:<br />
<em class="file">input.nml</em><em class=
"code">&amp;model_mod_nml</em>,<br />
<em class="file">gitm_in</em><em class=
"code">&amp;time_manager_nml</em>,<br />
<em class="file">gitm_in</em><em class=
"code">&amp;io_nml</em>,<br />
<em class="file">gitm_in</em><em class=
"code">&amp;init_ts_nml</em>,<br />
<em class="file">gitm_in</em><em class=
"code">&amp;restart_nml</em>,<br />
<em class="file">gitm_in</em><em class="code">&amp;domain_nml</em>,
and<br />
<em class="file">gitm_in</em><em class=
"code">&amp;grid_nml</em>.</p>
<table width="100%" border="0" summary="" cellpadding="3"></table>
</div>
<br />
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="end_model" id="end_model"></a><br />
<div class="routine"><em class="call">call end_model()</em></div>
<div class="indent1"><!-- Description -->
<p><em class="code">end_model</em> is used to clean up storage for
the model, etc. when the model is no longer needed. There are no
arguments and no return values. The grid variables are
deallocated.</p>
<table width="100%" border="0" summary="" cellpadding="3"></table>
</div>
<br />
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="init_time" id="init_time"></a><br />
<div class="routine"><em class="call">call init_time(time)</em>
<pre>
type(time_type), intent(out) :: <em class="code">time</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p><em class="code">init_time</em> returns the time at which the
model will start if no input initial conditions are to be used.
This is frequently used to spin-up models from rest, but is not
meaningfully supported for the GITM model. The only time this
routine would get called is if the <em class=
"file">input.nml</em><em class=
"code">&amp;perfect_model_obs_nml:start_from_restart</em> is
.false., which is not supported in the GITM model.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">time   </em></td>
<td>the starting time for the model if no initial conditions are to
be supplied. This is hardwired to 0.0</td>
</tr>
</table>
</div>
<br />
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="init_conditions" id="init_conditions"></a><br />
<div class="routine"><em class="call">call init_conditions(x)</em>
<pre>
real(r8), dimension(:), intent(out) :: <em class="code">x</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p><em class="code">init_conditions</em> returns default initial
conditions for model; generally used for spinning up initial model
states. For the GITM model it is just a stub because the initial
state is always provided by the input files.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">x   </em></td>
<td>Initial conditions for state vector. This is hardwired to
0.0</td>
</tr>
</table>
</div>
<br />
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="nc_write_model_atts" id="nc_write_model_atts"></a><br />
<div class="routine"><em class="call">ierr =
nc_write_model_atts(ncFileID)</em>
<pre>
integer             :: <em class="code">nc_write_model_atts</em>
integer, intent(in) :: <em class="code">ncFileID</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p><em class="code">nc_write_model_atts</em> writes model-specific
attributes to an opened netCDF file: In the GITM case, this
includes information like the coordinate variables (the grid
arrays: ULON, ULAT, TLON, TLAT, ZG, ZC, KMT, KMU), information from
some of the namelists, and either the 1D state vector or the
prognostic variables (SALT,TEMP,UVEL,VVEL,PSURF). All the required
information (except for the netCDF file identifier) is obtained
from the scope of the <em class="program">model_mod</em> module.
Both the <em class="file">input.nml</em> and <em class=
"file">gitm_in</em> files are preserved in the netCDF file as
variables <em class="code">inputnml</em> and <em class=
"code">gitm_in</em>, respectively.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class=
"code">ncFileID   </em></td>
<td>Integer file descriptor to previously-opened netCDF file.</td>
</tr>
<tr>
<td valign="top"><em class="code">ierr</em></td>
<td>Returns a 0 for successful completion.</td>
</tr>
</table>
<p><em class="code">nc_write_model_atts</em> is responsible for the
model-specific attributes in the following DART-output netCDF
files: <em class="file">true_state.nc</em>, <em class=
"file">preassim.nc</em>, and <em class="file">analysis.nc</em>.</p>
</div>
<br />
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="nc_write_model_vars" id="nc_write_model_vars"></a><br />
<div class="routine"><em class="call">ierr =
nc_write_model_vars(ncFileID, statevec, copyindex, timeindex)</em>
<pre>
integer,                intent(in) :: <em class=
"code">ncFileID</em>
real(r8), dimension(:), intent(in) :: <em class=
"code">statevec</em>
integer,                intent(in) :: <em class=
"code">copyindex</em>
integer,                intent(in) :: <em class=
"code">timeindex</em>
integer                            :: <em class="code">ierr</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p><em class="code">nc_write_model_vars</em> writes a copy of the
state variables to a NetCDF file. Multiple copies of the state for
a given time are supported, allowing, for instance, a single file
to include multiple ensemble estimates of the state. Whether the
state vector is parsed into prognostic variables (SALT, TEMP, UVEL,
VVEL, PSURF) or simply written as a 1D array is controlled by
<em class="file">input.nml</em><em class=
"code">&amp;model_mod_nml:output_state_vector</em>. If <em class=
"code">output_state_vector = .true.</em> the state vector
is written as a 1D array (the simplest case, but hard to explore
with the diagnostics). If <em class=
"code">output_state_vector = .false.</em> the state
vector is parsed into prognostic variables before being
written.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">ncFileID</em></td>
<td>file descriptor to previously-opened netCDF file.</td>
</tr>
<tr>
<td valign="top"><em class="code">statevec</em></td>
<td>A model state vector.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">copyindex   </em></td>
<td>Integer index of copy to be written.</td>
</tr>
<tr>
<td valign="top"><em class="code">timeindex</em></td>
<td>The timestep counter for the given state.</td>
</tr>
<tr>
<td valign="top"><em class="code">ierr</em></td>
<td>Returns 0 for normal completion.</td>
</tr>
</table>
</div>
<br />
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="pert_model_state" id="pert_model_state"></a><br />
<div class="routine"><em class="call">call pert_model_state(state,
pert_state, interf_provided)</em>
<pre>
real(r8), dimension(:), intent(in)  :: <em class="code">state</em>
real(r8), dimension(:), intent(out) :: <em class=
"code">pert_state</em>
logical,                intent(out) :: <em class=
"code">interf_provided</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Given a model state, <em class="code">pert_model_state</em>
produces a perturbed model state. This is used to generate ensemble
initial conditions perturbed around some control trajectory state
when one is preparing to spin-up ensembles. Since the DART state
vector for the GITM model contains both 'wet' and 'dry' cells, it
is imperative to provide an interface to perturb
<strong>just</strong> the wet cells (<em class=
"code">interf_provided == .true.</em>).<br />
<br />
The magnitude of the perturbation is wholly determined by
<em class="file">input.nml</em><em class=
"code">&amp;model_mod_nml:model_perturbation_amplitude</em> and
<strong>utterly, completely fails</strong>.<br />
<br />
A more robust perturbation mechanism is needed. Until then, avoid
using this routine by using your own ensemble of initial
conditions. This is determined by setting <em class=
"file">input.nml</em><em class=
"code">&amp;filter_nml:start_from_restart = .false.</em></p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">state</em></td>
<td>State vector to be perturbed.</td>
</tr>
<tr>
<td valign="top"><em class="code">pert_state</em></td>
<td>The perturbed state vector.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">interf_provided   </em></td>
<td>Because of the 'wet/dry' issue discussed above, this is always
<em class="code">.true.</em>, indicating a model-specific
perturbation is available.</td>
</tr>
</table>
</div>
<br />
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_close_maxdist_init" id=
"get_close_maxdist_init"></a><br />
<div class="routine"><em class="call">call
get_close_maxdist_init(gc, maxdist)</em>
<pre>
type(get_close_type), intent(inout) :: <em class="code">gc</em>
real(r8),             intent(in)    :: <em class=
"code">maxdist</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Pass-through to the 3-D sphere locations module. See <a href=
"../../assimilation_code/location/threed_sphere/location_mod.html#get_close_maxdist_init">
get_close_maxdist_init()</a> for the documentation of this
subroutine.</p>
</div>
<br />
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_close_obs_init" id="get_close_obs_init"></a><br />
<div class="routine"><em class="call">call get_close_obs_init(gc,
num, obs)</em>
<pre>
type(get_close_type), intent(inout) :: <em class="code">gc</em>
integer,              intent(in)    :: <em class="code">num</em>
type(location_type),  intent(in)    :: <em class=
"code">obs(num)</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Pass-through to the 3-D sphere locations module. See <a href=
"../../assimilation_code/location/threed_sphere/location_mod.html#get_close_obs_init">
get_close_obs_init()</a> for the documentation of this
subroutine.</p>
</div>
<br />
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_close_obs" id="get_close_obs"></a><br />
<div class="routine"><em class="call">call get_close_obs(gc,
base_obs_loc, base_obs_kind, obs, obs_kind, &amp;<br />
          num_close,
close_ind <em class="optionalcode">[, dist]</em>)</em>
<pre>
type(get_close_type),              intent(in ) :: <em class=
"code">gc</em>
type(location_type),               intent(in ) :: <em class=
"code">base_obs_loc</em>
integer,                           intent(in ) :: <em class=
"code">base_obs_kind</em>
type(location_type), dimension(:), intent(in ) :: <em class=
"code">obs</em>
integer,             dimension(:), intent(in ) :: <em class=
"code">obs_kind</em>
integer,                           intent(out) :: <em class=
"code">num_close</em>
integer,             dimension(:), intent(out) :: <em class=
"code">close_ind</em>
real(r8), optional,  dimension(:), intent(out) :: <em class=
"optionalcode">dist</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p>Given a DART location (referred to as "base") and a set of
locations, and a definition of 'close' - return a subset of
locations that are 'close', as well as their distances to the DART
location and their indices. This routine intentionally masks a
routine of the same name in <em class="program">location_mod</em>
because we want to be able to discriminate against selecting 'dry
land' locations.<br />
<br />
Given a single location and a list of other locations, returns the
indices of all the locations close to the single one along with the
number of these and the distances for the close ones. The list of
locations passed in via the <em class="code">obs</em> argument must
be identical to the list of <em class="code">obs</em> passed into
the most recent call to <em class="code">get_close_obs_init()</em>.
If the list of locations of interest changes, <em class=
"code">get_close_obs_destroy()</em> must be called and then the two
initialization routines must be called before using <em class=
"code">get_close_obs()</em> again.<br />
<br />
For vertical distance computations, the general philosophy is to
convert all vertical coordinates to a common coordinate. This
coordinate type is defined in the namelist with the variable
"vert_localization_coord".</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">gc</em></td>
<td>Structure to allow efficient identification of locations
'close' to a given location.</td>
</tr>
<tr>
<td valign="top"><em class="code">base_obs_loc</em></td>
<td>Single given location.</td>
</tr>
<tr>
<td valign="top"><em class="code">base_obs_kind</em></td>
<td>Kind of the single location.</td>
</tr>
<tr>
<td valign="top"><em class="code">obs</em></td>
<td>List of candidate locations.</td>
</tr>
<tr>
<td valign="top"><em class="code">obs_kind</em></td>
<td>Kind associated with candidate locations.</td>
</tr>
<tr>
<td valign="top"><em class="code">num_close</em></td>
<td>Number of locations close to the given location.</td>
</tr>
<tr>
<td valign="top"><em class="code">close_ind</em></td>
<td>Indices of those locations that are close.</td>
</tr>
<tr>
<td valign="top"><em class="optionalcode">dist</em></td>
<td>Distance between given location and the close ones identified
in close_ind.</td>
</tr>
</table>
</div>
<br />
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="ens_mean_for_model" id="ens_mean_for_model"></a><br />
<div class="routine"><em class="call">call
ens_mean_for_model(ens_mean)</em>
<pre>
real(r8), dimension(:), intent(in) :: <em class=
"code">ens_mean</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p><em class="code">ens_mean_for_model</em> normally saves a copy
of the ensemble mean to module-local storage. This is a
NULL_INTERFACE for the GITM model. At present there is no
application which requires module-local storage of the ensemble
mean. No storage is allocated.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">ens_mean</em></td>
<td>State vector containing the ensemble mean.</td>
</tr>
</table>
</div>
<br />
<!--====================================================================-->
<div class="top">[<a href="#">top</a>]</div>
<hr />
<h3 class="indent1">Unique interface routine descriptions</h3>
<!--====================================================================-->
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
<a name="get_gridsize" id="get_gridsize"></a><br />
<div class="routine"><em class="call">call get_gridsize( num_x,
num_y, num_z )</em>
<pre>
integer, intent(out) :: <em class="code">num_x, num_y, num_z</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p><em class="code">get_gridsize</em> returns the dimensions of the
compute domain. The horizontal gridsize is determined from
<em class="file">gitm_restart.nc</em>.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">num_x</em></td>
<td>The number of longitudinal gridpoints.</td>
</tr>
<tr>
<td valign="top"><em class="code">num_y</em></td>
<td>The number of latitudinal gridpoints.</td>
</tr>
<tr>
<td valign="top"><em class="code">num_z</em></td>
<td>The number of vertical gridpoints.</td>
</tr>
</table>
</div>
<br />
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="restart_file_to_sv" id="restart_file_to_sv"></a><br />
<div class="routine"><em class="call">call
restart_file_to_sv(filename, state_vector, model_time)</em>
<pre>
character(len=*),       intent(in)    :: <em class=
"code">filename</em>
real(r8), dimension(:), intent(inout) :: <em class=
"code">state_vector</em>
type(time_type),        intent(out)   :: <em class=
"code">model_time</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p><em class="code">restart_file_to_sv</em> Reads a GITM netCDF
format restart file and packs the desired variables into a DART
state vector. The desired variables are specified in the <em class=
"code">gitm_vars_nml</em> namelist.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">filename</em></td>
<td>The name of the netCDF format GITM restart file.</td>
</tr>
<tr>
<td valign="top"><em class="code">state_vector</em></td>
<td>the 1D array containing the concatenated GITM variables.</td>
</tr>
<tr>
<td valign="top"><em class="code">model_time</em></td>
<td>the time of the model state. The last time in the netCDF
restart file.</td>
</tr>
</table>
</div>
<br />
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="sv_to_restart_file" id="sv_to_restart_file"></a><br />
<div class="routine"><em class="call">call
sv_to_restart_file(state_vector, filename, statedate)</em>
<pre>
real(r8), dimension(:), intent(in) :: <em class=
"code">state_vector</em>
character(len=*),       intent(in) :: <em class=
"code">filename</em>
type(time_type),        intent(in) :: <em class=
"code">statedate</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p><em class="code">sv_to_restart_file</em> updates the variables
in the GITM restart file with values from the DART vector
<em class="code">state_vector</em>. The last time in the file must
match the <em class="code">statedate</em>.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class="code">filename</em></td>
<td>the netCDF-format GITM restart file to be updated.</td>
</tr>
<tr>
<td valign="top"><em class="code">state_vector</em></td>
<td>the 1D array containing the DART state vector.</td>
</tr>
<tr>
<td valign="top"><em class="code">statedate</em></td>
<td>the 'valid_time' of the DART state vector.</td>
</tr>
</table>
</div>
<br />
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_gitm_restart_filename" id=
"get_gitm_restart_filename"></a><br />
<div class="routine"><em class="call">call
get_gitm_restart_filename( filename )</em>
<pre>
character(len=*), intent(out) :: <em class="code">filename</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p><em class="code">get_gitm_restart_filename</em> returns the name
of the gitm restart file - the filename itself is in private module
storage.</p>
<table width="100%" border="0" summary="" cellpadding="3">
<tr>
<td valign="top"><em class=
"code">filename   </em></td>
<td>The name of the GITM restart file.</td>
</tr>
</table>
</div>
<br />
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_base_time" id="get_base_time"></a><br />
<div class="routine"><em class="call">time = get_base_time(
filehandle )</em>
<pre>
integer,          intent(in) :: <em class=
"code">filehandle</em> -OR-
character(len=*), intent(in) :: <em class="code">filehandle</em>
type(time_type),  intent(out) :: <em class="code">time</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p><em class="code">get_base_time</em> extracts the start time of
the experiment as contained in the netCDF restart file. The file
may be specified by either a character string or the integer netCDF
fid.</p>
</div>
<br />
<!--===================== DESCRIPTION OF A ROUTINE =====================-->
 <a name="get_state_time" id="get_state_time"></a><br />
<div class="routine"><em class="call">time = get_state_time(
filehandle )</em>
<pre>
integer,          intent(in) :: <em class=
"code">filehandle</em> -OR-
character(len=*), intent(in) :: <em class="code">filehandle</em>
type(time_type),  intent(out) :: <em class="code">time</em>
</pre></div>
<div class="indent1"><!-- Description -->
<p><em class="code">get_state_time</em> extracts the time of the
model state as contained in the netCDF restart file. In the case of
multiple times in the file, the last time is the time returned. The
file may be specified by either a character string or the integer
netCDF fid.</p>
</div>
<br />
<!--==================================================================-->
<!-- Describe the Files Used by this module.                          -->
<!--==================================================================-->
 <a name="FilesUsed" id="FilesUsed"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr />
<h2>FILES</h2>
<table border="0" width="100%">
<tr>
<th align="left">filename</th>
<th align="left">purpose</th>
</tr>
<tr>
<td>input.nml</td>
<td>to read the model_mod namelist</td>
</tr>
<tr>
<td>gitm_vars.nml</td>
<td>to read the <em class="code">gitm_vars_nml</em> namelist</td>
</tr>
<tr>
<td>gitm_restart.nc</td>
<td>provides grid dimensions, model state, and 'valid_time' of the
model state</td>
</tr>
<tr>
<td>true_state.nc</td>
<td>the time-history of the "true" model state from an OSSE</td>
</tr>
<tr>
<td>preassim.nc</td>
<td>the time-history of the model state before assimilation</td>
</tr>
<tr>
<td>analysis.nc </td>
<td>the time-history of the model state after assimilation</td>
</tr>
<tr>
<td>dart_log.out [default name]</td>
<td>the run-time diagnostic output</td>
</tr>
<tr>
<td>dart_log.nml [default name]</td>
<td>the record of all the namelists actually USED - contains the
default values</td>
</tr>
</table>
<br />
<!--==================================================================-->
<!-- Cite references, if need be.                                     -->
<!--==================================================================-->
 <a name="References" id="References"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr />
<h2>REFERENCES</h2>
<ul>
<li>none</li>
</ul>
<!--==================================================================-->
<!-- Describe all the error conditions and codes.                     -->
<!--==================================================================-->
<a name="Errors" id="Errors"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr />
<h2>ERROR CODES and CONDITIONS</h2>
<div class="errors">
<table border="1" cellspacing="1" cellpadding="10" width="100%">
<tr>
<th>Routine</th>
<th width="50%">Message</th>
<th>Comment</th>
</tr>
<tr><!-- routine -->
<td valign="top">restart_file_to_sv</td>
<!-- message -->
<td valign="top">cannot open file "xxxx" for reading</td>
<!-- comment -->
<td valign="top">The GITM restart file "xxxx" does not exist.</td>
</tr>
<tr><!-- routine -->
<td valign="top">restart_file_to_sv</td>
<!-- message -->
<td valign="top">'WARNING!!! year 0 not supported; setting to year
1</td>
<!-- comment -->
<td valign="top">year 0 ... is not supported in a Gregorian
calendar. Our intent here is to do data assimilation, normally
'real' observations have 'real' dates.</td>
</tr>
<tr><!-- routine -->
<td valign="top">sv_to_restart_file</td>
<!-- message -->
<td valign="top">current time /= model time. FATAL error.</td>
<!-- comment -->
<td valign="top">The DART time does not match the time of the GITM
restart file. This message is preceeded by several lines indicating
the expected times of both DART and GITM.</td>
</tr>
</table>
</div>
<h2>KNOWN BUGS</h2>
<p>none at this time</p>
<!--==================================================================-->
<!-- Describe Future Plans.                                           -->
<!--==================================================================-->
<a name="FuturePlans" id="FuturePlans"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr />
<h2>FUTURE PLANS</h2>
<p>Provide a better mechanism for generating a set of perturbed
initial conditions - <em class="code">pert_model_state()</em></p>
<!--==================================================================-->
<!-- PrivateComponents                                                -->
<!--==================================================================-->
<a name="PrivateComponents" id="PrivateComponents"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr />
<h2>PRIVATE COMPONENTS</h2>
<p>N/A</p>
<!--==================================================================-->
<!-- Legalese & Metadata                                              -->
<!--==================================================================-->
<a name="Legalese" id="Legalese"></a>
<div class="top">[<a href="#">top</a>]</div>
<hr />
<h2>Terms of Use</h2>
<p>DART software - Copyright UCAR. This open source software is
provided by UCAR, "as is", without charge, subject to all terms of
use at <a href=
"http://www.image.ucar.edu/DAReS/DART/DART_download">http://www.image.ucar.edu/DAReS/DART/DART_download</a></p>
<!--==================================================================-->
</body>
</html>
