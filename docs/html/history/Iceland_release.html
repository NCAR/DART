<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
"http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>DART "Iceland release" Documentation</title>
<link rel="stylesheet" type="text/css" href="../../html/doc.css">
<link href="../../images/dart.ico" rel="shortcut icon">
</head>
<body>
<a name="TOP" id="TOP"></a>
<h1>Iceland</h1>
<h2>DART Iceland release documentation</h2>
<div class="attention">
<div class="title">
<p>Attention</p>
</div>
<p>Iceland is a prior release of DART. Its source code is available via the
<a href="https://github.com/NCAR/DART/tree/Iceland">DART repository on Github</a>.
This documentation is preserved merely for reference. See the
<a href="https://dart.ucar.edu/">DART homepage</a> to learn about the latest
release.</p>
</div>
<!--==================================================================-->
<a name="OVERVIEW" id="OVERVIEW"></a>
<h2>Overview of DART</h2>
<p>The Data Assimilation Research Testbed (DART) is designed to
facilitate the combination of assimilation algorithms, models, and
observation sets to allow increased understanding of all three. The
DART programs have been compiled with several Fortran 90 compilers
and run on a linux compute-server and linux clusters. You should
definitely read the <a href="#customizations">Customizations</a>
section.</p>

<p>DART employs a modular programming approach to apply an Ensemble
Kalman Filter which nudges models toward a state that is more
consistent with information from a set of observations. Models may
be swapped in and out, as can different algorithms in the Ensemble
Kalman Filter. The method requires running multiple instances of a
model to generate an ensemble of states. A forward operator
appropriate for the type of observation being used is applied to
each of the states to generate the model's estimate of the
observation. Comparing these estimates and their uncertainty to the
observation and its uncertainty ultimately results in the
adjustments to the model states. Sort of. There's more to it,
described in detail in the tutorial directory of the package.</p>

<p>DART ultimately creates a few netCDF files containing the model
states just before the adjustment <em class=
"file">Prior_Diag.nc</em> and just after the adjustment <em class=
"file">Posterior_Diag.nc</em> as well as a file <em class=
"file">obs_seq.final</em> with the model estimates of the
observations. There is a suite of Matlab® functions that facilitate
exploration of the results.</p>

<p>The Iceland release provides several new models and has a
greatly expanded capability for <b>real</b> observations which
required a fundamentally different implementation of the low-level
routines. It is now required to run a preprocessor on several of
the program units to construct the source code files which will be
compiled by the remaining units. Due to the potentially large
number of observations types possible and for portability reasons,
the preprocessor is actually a F90 program that uses the namelist
mechanism for specifying the observation types to be included. This
also prevents having a gory set of compile flags that is different
for every compiler. One very clever colleague also 'built a better
mousetrap' and figured out how to effectively and robustly read
namelists, detect errors, and generate meaningful error messages.
HURRAY!</p>

<p>The Iceland release has also been tested with more compilers in
an attempt to determine non-portable code elements. It is my
experience that the largest impediment to portable code is the
reliance on the compiler to autopromote <tt>real</tt> variables to
one flavor or another. <a href=
"/DART/PrecisionNotes.html">Different compilers perform this
autopromotion in different ways.</a> Using the F90 "kind" allows
for much more flexible code, in that the use of interface
procedures is possible only when two routines do not have identical
sets of input arguments -- something that happens when the compiler
autopromotes 32bit reals to 64bit reals, for example.</p>
<p>DART programs can require three different types of input. First,
some of the DART programs, those for creating synthetic
observational datasets, require interactive input from the
keyboard. For simple cases, this interactive input can be made
directly from the keyboard. In more complicated cases, a file
containing the appropriate keyboard input can be created and this
file can be directed to the standard input of the DART program.
Second, many DART programs expect one or more input files in DART
specific formats to be available. For instance, <em class=
"program">perfect_model_obs</em>, which creates a synthetic
observation set given a particular model and a description of a
sequence of observations, requires an input file that describes
this observation sequence. At present, the observation files for
DART are in a custom format in either human-readable ascii or more
compact machine-specific binary. Third, many DART modules
(including main programs) make use of the Fortan90 namelist
facility to obtain values of certain parameters at run-time. All
programs look for a namelist input file called <em class=
"file">input.nml</em> in the directory in which the program is
executed. The <em class="file">input.nml</em> file can contain a
sequence of individual Fortran90 namelists which specify values of
particular parameters for modules that compose the executable
program. DART provides a mechanism that automatically generates
namelists with the default values for each program to be run.</p>
<p>DART uses the <a href=
"http://www.unidata.ucar.edu/packages/netcdf/">netCDF</a>
self-describing data format with a particular metadata convention
to describe output that is used to analyze the results of
assimilation experiments. These files have the extension <em class=
"file">.nc</em> and can be read by a number of standard data
analysis tools. A set of <a href=
"http://www.mathworks.com/">Matlab</a> scripts, designed to produce
graphical diagnostics from DART netCDF output files are available.
DART users have also used <a href=
"http://meteora.ucsd.edu/~pierce/ncview_home_page.html">ncview</a>
to create rudimentary graphical displays of output data fields. The
<a href="http://nco.sourceforge.net">NCO</a> tools, produced by
UCAR's Unidata group, are available to do operations like
concatenating, slicing, and dicing of netCDF files.</p>
<!--==================================================================-->
<h2 class="indent1">Document conventions</h2>
<p>Anything underlined is a URL.<br>
<br>
<em class="file">All filenames look like this -- (typewriter font,
green)</em>.<br>
<em class="program">Program names look like this -- (italicized
font, green)</em>.<br>
<em class="input">user input looks like this -- (bold,
magenta)</em>.</p>
<div class="unix">commands to be typed at the command line are
contained in an indented gray box.</div>
<p>And the contents of a file are enclosed in a box with a
border:</p>
<div class="routine">&amp;hypothetical_nml<br>
  obs_seq_in_file_name = "obs_seq.in",<br>
  obs_seq_out_file_name = "obs_seq.out",<br>
  init_time_days = 0,<br>
  init_time_seconds = 0,<br>
  output_interval = 1<br>
&amp;end</div>
<!--==================================================================-->
<a name="Installation" id="Installation"></a>
<hr>
<h1>Installation</h1>
<p>This document outlines the installation of the DART software and
the system requirements. For convenience, some of the original
colloquium exercises are repeated here, mostly just to check the
installation. A few of the <a href=
"/DART/dart_exercise_doc.pdf">exercises from the ASP summer 2003
Colloquium</a> are repeated here, primarily to serve as the
verification of the installation. The entire installation process
is summarized in the following steps:</p>
<ol>
<li><a href="#compilers">Determine which F90 compiler is
available</a>.</li>
<li><a href="#netCDFlib">Determine the location of the <em class=
"code">netCDF</em> library</a>.</li>
<li><a href="#download">Download the DART software bundle and untar
it into the expected source tree</a>.</li>
<li><a href="#customizations">Modify certain DART files to reflect
the available F90 compiler and location of the appropriate
libraries</a>.</li>
<li><a href="#building">Build the executables</a>.</li>
</ol>
<p>We have tried to make the code as portable as possible, but we
do not have access to all compilers on all platforms, so there are
no guarantees. We are interested in your experience building the
system, so please email me (Tim Hoar)
thoar 'at' ucar 'dot' edu (trying to cut down
on the spam).</p>
<p>After the installation, you might want to peruse the
following.</p>
<ul>
<li><a href="#Running">Running the Lorenz_63 Model</a>.</li>
<li><a href="#matlab">Using the Matlab diagnostic scripts</a>.</li>
<li>A short discussion on <a href="#discussion">bias, filter
divergence and covariance inflation.</a></li>
<li>And another one on <a href="#syntheticobservations">synthetic
observations</a>.</li>
</ul>
<!--==================================================================-->
<a name="compilers" id="compilers"></a>
<hr>
<h2>Requirements: an F90 Compiler</h2>
<p>The DART software has been successfully built on several
Linux/x86 platforms with several versions of the <a href=
"http://www.intel.com/software/products/compilers/flin">Intel
Fortran Compiler for Linux</a>, which (at one point) is/was free
for individual scientific use. It has also been built and
successfully run with several versions of each of the following:
<a href="http://www.pgroup.com">Portland Group Fortran
Compiler</a>, <a href="http://www.lahey.com">Lahey Fortran
Compiler</a>, <a href="http://www.pathscale.com">Pathscale Fortran
Compiler</a>, <a href="http://www.absoft.com">Absoft Fortran 90/95
Compiler (Mac OSX)</a>. Since recompiling the code is a necessity
to experiment with different models, there are no binaries to
distribute.</p>
<!--==================================================================-->
<a name="netCDFlib" id="netCDFlib"></a>
<hr>
<h2>Requirements: the <em class="file">netCDF</em> library</h2>
<p>DART uses the <a href=
"http://www.unidata.ucar.edu/packages/netcdf/">netCDF</a>
self-describing data format for the results of assimilation
experiments. These files have the extension <em class=
"file">.nc</em> and can be read by a number of standard data
analysis tools. In particular, DART also makes use of the F90
interface to the library which is available through the <em class=
"file">netcdf.mod</em> and <em class="file">typesizes.mod</em>
modules. <em class="bold">IMPORTANT</em>: different compilers
create these modules with different "case" filenames, and sometimes
they are not <strong>both</strong> installed into the expected
directory. It is required that both modules be present. The normal
place would be in the <tt>netcdf/include</tt> directory, as opposed
to the <tt>netcdf/lib</tt> directory.</p>
<p>If the netCDF library does not exist on your system, you must
build it (as well as the F90 interface modules). The library and
instructions for building the library or installing from an RPM may
be found at the netCDF home page: <a href=
"http://www.unidata.ucar.edu/packages/netcdf/">http://www.unidata.ucar.edu/packages/netcdf/</a>
Pay particular attention to the compiler-specific patches that must
be applied for the Intel Fortran Compiler. (Or the PG compiler, for
that matter.)</p>
<p>The location of the netCDF library, <em class=
"file">libnetcdf.a</em>, and the locations of both <em class=
"file">netcdf.mod</em> and <em class="file">typesizes.mod</em> will
be needed by the makefile template, as described in the <a href=
"#compiling">compiling</a> section.</p>
<!--==================================================================-->
<a name="download" id="download"></a>
<hr>

<h2>Unpacking the distribution.</h2>
<p>This release of the <a href="https://github.com/NCAR/DART/releases/tag/v5.0.0">
DART source code can be downloaded</a> as a compressed zip or tar.gz file.
When extracted, the source tree will begin with a directory
named <em class="file">DART</em> and will be approximately 103.1 Mb.
Compiling the code in this tree (as is usually the case) will
necessitate much more space.</p>
<pre>
<code>
$ gunzip DART-5.0.0.tar.gz
$ tar -xvf DART-5.0.0.tar
</code>
</pre>

<p>The code tree is very "bushy"; there are many directories of
support routines, etc. but only a few directories involved with the
customization and installation of the DART software. If you can
compile and run ONE of the low-order models, you should be able to
compile and run ANY of the low-order models. For this reason, we
can focus on the Lorenz `63 model. Subsequently, the only
directories with files to be modified to check the installation
are:  <em class="file">DART/mkmf</em>,  <em class=
"file">DART/models/lorenz_63/work</em>, and  <em class=
"file">DART/matlab</em> (but only for analysis).</p>
<!--==================================================================-->
<a name="customizations" id="customizations"></a>
<hr>
<h2>Customizing the build scripts -- Overview.</h2>
<p>DART executable programs are constructed using two tools:
<em class="program">make</em> and <em class="program">mkmf</em>.
The <em class="program">make</em> utility is a relatively common
piece of software that requires a user-defined input file that
records dependencies between different source files. <em class=
"program">make</em> then performs a hierarchy of actions when one
or more of the source files is modified. The <em class=
"program">mkmf</em> utility is a custom preprocessor that generates
a <em class="program">make</em> input file (named <em class=
"file">Makefile</em>) and an example namelist <em class=
"file">input.nml.<em class=
"program"><i>program</i>_default</em></em> with the default values.
The <em class="file">Makefile</em> is designed specifically to work
with object-oriented Fortran90 (and other languages) for systems
like DART.</p>
<p><em class="program">mkmf</em> requires two separate input files.
The first is a `template' file which specifies details of the
commands required for a specific Fortran90 compiler and may also
contain pointers to directories containing pre-compiled utilities
required by the DART system. <strong>This template file will need
to be modified to reflect your system</strong>. The second input
file is a `path_names' file which includes a complete list of the
locations (either relative or absolute) of all Fortran90 source
files that are required to produce a particular DART program. Each
'path_names' file must contain a path for exactly one Fortran90
file containing a main program, but may contain any number of
additional paths pointing to files containing Fortran90 modules. An
<em class="program">mkmf</em> command is executed which uses the
'path_names' file and the mkmf template file to produce a
<em class="file">Makefile</em> which is subsequently used by the
standard <em class="program">make</em> utility.</p>
<p>Shell scripts that execute the mkmf command for all standard
DART executables are provided as part of the standard DART
software. For more information on <em class="program">mkmf</em> see
<a href=
"http://www.gfdl.gov/fms/pubrel/j/atm_dycores/doc/dycore_public_manual.html#mkmf">
the FMS mkmf description</a>.</p>

<p>One of the benefits of using <em class="program">mkmf</em> is that
it also creates an example namelist file for each program. The
example namelist is called <em class="file">input.nml.<em class=
"program"><i>program</i>_default</em></em>, so as not to clash with
any exising <em class="file">input.nml</em> that may exist in that
directory.</p>

<a name="template" id="template"></a>
<h3 class="indent1">Building and Customizing the 'mkmf.template'
file</h3>
<p>A series of templates for different compilers/architectures
exists in the <em class="file">DART/mkmf/</em> directory and have
names with extensions that identify either the compiler, the
architecture, or both. This is how you inform the build process of
the specifics of your system. Our intent is that you copy one that
is similar to your system into <em class="file">mkmf.template</em>
and customize it. For the discussion that follows, knowledge of the
contents of one of these templates (i.e. <em class=
"file">mkmf.template.pgf90.ghotiol</em>) is needed: (note that only
the LAST lines are shown here, the head of the file is just a big
comment)</p>

<pre>
<code>
# Makefile template for PGI f90
FC = pgf90
LD = pgf90
CPPFLAGS =
LIST = -Mlist
NETCDF = /contrib/netcdf-3.5.1-cc-c++-pgif90.5.2-4
FFLAGS = -O0 -Ktrap=fp -pc 64 -I$(NETCDF)/include
LIBS = -L$(NETCDF)/lib -lnetcdf
LDFLAGS = $(LIBS)

...
</code>
</pre>

<p>Essentially, each of the lines defines some part of the
resulting <em class="file">Makefile</em>. Since <em class=
"program">make</em> is particularly good at sorting out
dependencies, the order of these lines really doesn't make any
difference. The <em class="code">FC = pgf90</em> line ultimately
defines the Fortran90 compiler to use, etc. The lines which are
most likely to need site-specific changes start with <em class=
"code">FFLAGS</em> and <em class="code">NETCDF</em>, which indicate
where to look for the netCDF F90 modules and the location of the
netCDF library and modules.</p>

<h4 class="indent2">NETCDF</h4>
<p class="indent1">Modifying the <em class="code">NETCDF</em> value
should be relatively straightforward.<br>
Change the string to reflect the location of your netCDF
installation containing <em class="file">netcdf.mod</em> and
<em class="file">typesizes.mod</em>. The value of the <em class=
"code">NETCDF</em> variable will be used by the <em class=
"code">FFLAGS, LIBS,</em> and <em class="code">LDFLAGS</em>
variables.<br></p>

<h4 class="indent2">FFLAGS</h4>
<p class="indent1">Each compiler has different compile flags, so
there is really no way to exhaustively cover this other than to say
the templates as we supply them should work -- depending on the
location of your netCDF. The low-order models can be compiled
without a <em class="code">-r8</em> switch, but the <em class=
"file">bgrid_solo</em> model cannot.</p>
<a name="path_names" id="path_names"></a>
<h3 class="indent1">Customizing the 'path_names_*' file</h3>
<p>Several <em class="file">path_names_*</em> files are provided in
the <em class="file">work</em> directory for each specific model,
in this case: <em class="file">DART/models/lorenz_63/work</em>.</p>
<ol>
<li><em class="file">path_names_preprocess</em></li>
<li><em class="file">path_names_create_obs_sequence</em></li>
<li><em class="file">path_names_create_fixed_network_seq</em></li>
<li><em class="file">path_names_perfect_model_obs</em></li>
<li><em class="file">path_names_filter</em></li>
<li><em class="file">path_names_obs_diag</em></li>
</ol>
<p>Since each model comes with its own set of files, no further
customization is needed.</p>
<!--==================================================================-->
<a name="building" id="building"></a>
<hr>
<h2>Building the Lorenz_63 DART project.</h2>
<p>Currently, DART executables are constructed in a <em class=
"file">work</em> subdirectory under the directory containing code
for the given model. In the top-level DART directory, change to the
L63 work directory and list the contents:</p>
<pre>
<code>
$ cd DART/models/lorenz_63/work
$ ls -1
</code>
</pre>

<p>With the result:</p>
<pre>
<code>
filter_ics
filter_restart
input.nml
mkmf_create_fixed_network_seq
mkmf_create_obs_sequence
mkmf_filter
mkmf_obs_diag
mkmf_perfect_model_obs
mkmf_preprocess
obs_seq.final
obs_seq.in
obs_seq.out
obs_seq.out.average
obs_seq.out.x
obs_seq.out.xy
obs_seq.out.xyz
obs_seq.out.z
path_names_create_fixed_network_seq
path_names_create_obs_sequence
path_names_filter
path_names_obs_diag
path_names_perfect_model_obs
path_names_preprocess
perfect_ics
perfect_restart
Posterior_Diag.nc
Prior_Diag.nc
set_def.out
True_State.nc
workshop_setup.csh
</code>
</pre>

<p>There are six <em class="file">mkmf_</em><em class=
"italic">xxxxxx</em> files for the programs <em class=
"program">preprocess</em>, <em class=
"program">create_obs_sequence</em>, <em class=
"program">create_fixed_network_seq</em>, <em class=
"program">perfect_model_obs</em>, <em class="program">filter</em>,
and <em class="program">obs_diag</em> along with the corresponding
<em class="file">path_names_</em><em class="italic">xxxxxx</em>
files. You can examine the contents of one of the <em class=
"file">path_names_</em><em class="italic">xxxxxx</em> files, for
instance <em class="file">path_names_filter</em>, to see a list of
the relative paths of all files that contain Fortran90 modules
required for the program <em class="program">filter</em> for the
L63 model. All of these paths are relative to your <em class=
"file">DART</em> directory. The first path is the main program
(<em class="file">filter.f90</em>) and is followed by all the
Fortran90 modules used by this program (after preprocessing).</p>
<p>The <em class="program">mkmf_</em><em class="italic">xxxxxx</em>
scripts are cryptic but should not need to be modified -- as long
as you do not restructure the code tree (by moving directories, for
example). The only function of the <em class=
"program">mkmf_</em><em class="italic">xxxxxx</em> script is to
generate a <em class="file">Makefile</em> and an <em class=
"file">input.nml.<em class=
"program"><i>program</i>_default</em></em> file. It is not supposed
to compile anything:</p>
<pre>
<code>
$ csh mkmf_preprocess
$ make
</code>
</pre>
<p>The first command generates an appropriate <em class=
"file">Makefile</em> and the <em class=
"file">input.nml.preprocess_default</em> file. The second command
results in the compilation of a series of Fortran90 modules which
ultimately produces an executable file: <em class=
"program">preprocess</em>. Should you need to make any changes to
the <em class="file">DART/mkmf/mkmf.template</em>, you will need to
regenerate the <em class="file">Makefile</em>.<br>
<br>
The <em class="program">preprocess</em> program actually builds
source code to be used by all the remaining modules. It is
<strong>imperative</strong> to actually <strong>run</strong>
<em class="program">preprocess</em> before building the remaining
executables. This is how the same code can assimilate state vector
'observations' for the Lorenz_63 model and real radar
reflectivities for WRF without needing to specify a set of radar
operators for the Lorenz_63 model!<br>
<br>
<em class="program">preprocess</em> reads the <em class=
"code">&amp;preprocess_nml</em> namelist to determine what
observations and operators to incorporate. For this exercise, we
will use the values in <em class="file">input.nml</em>. <em class=
"program">preprocess</em> is designed to abort if the files it is
supposed to build already exist. For this reason, it is necessary
to remove a couple files (if they exist) before you run the
preprocessor. It is just a good habit to develop.</p>

<pre>
<code>
$ \rm -f ../../../obs_def/obs_def_mod.f90
$ \rm -f ../../../obs_kind/obs_kind_mod.f90
$ ./preprocess
$ ls -l ../../../obs_def/obs_def_mod.f90
$ ls -l ../../../obs_kind/obs_kind_mod.f90
</code>
</pre>

<p>This created <em class=
"file">../../../obs_def/obs_def_mod.f90</em> from <em class=
"file">../../../obs_kind/DEFAULT_obs_kind_mod.F90</em> and several
other modules. <em class=
"file">../../../obs_kind/obs_kind_mod.f90</em> was created
similarly. Now we can build the rest of the project.<br>
<br>
A series of object files for each module compiled will also be left
in the work directory, as some of these are undoubtedly needed by
the build of the other DART components. You can proceed to create
the other five programs needed to work with L63 in DART as
follows:</p>

<pre>
<code>
$ csh mkmf_create_obs_sequence
$ make
$ csh mkmf_create_fixed_network_seq
$ make
$ csh mkmf_perfect_model_obs
$ make
$ csh mkmf_filter
$ make
$ csh mkmf_obs_diag
$ make
</code>
</pre>

<p>The result (hopefully) is that six executables now reside in
your work directory. The most common problem is that the netCDF
libraries and include files (particularly <em class=
"file">typesizes.mod</em>) are not found. Edit the <em class=
"file">DART/mkmf/mkmf.template</em>, recreate the <em class=
"file">Makefile</em>, and try again.</p>
<table border="0" cellpadding="1" width="100%">
<tr>
<th>program</th>
<th>purpose</th>
</tr>
<tr>
<td><em class="program">preprocess</em></td>
<td>creates custom source code for just the observations of
interest</td>
</tr>
<tr>
<td><em class="program">create_obs_sequence</em></td>
<td>specify a (set) of observation characteristics taken by a
particular (set of) instruments</td>
</tr>
<tr>
<td><em class="program">create_fixed_network_seq</em></td>
<td>specify the temporal attributes of the observation sets</td>
</tr>
<tr>
<td><em class="program">perfect_model_obs</em></td>
<td>spinup, generate "true state" for synthetic observation
experiments, ...</td>
</tr>
<tr>
<td><em class="program">filter</em></td>
<td>perform experiments</td>
</tr>
<tr>
<td><em class="program">obs_diag</em></td>
<td>creates observation-space diagnostic files to be explored by
the Matlab® scripts.</td>
</tr>
</table>
<!--==================================================================-->
<!--==================================================================-->
<a name="Running" id="Running"></a>
<hr>
<h2>Running Lorenz_63.</h2>
<p>This initial sequence of exercises includes detailed
instructions on how to work with the DART code and allows
investigation of the basic features of one of the most famous
dynamical systems, the 3-variable Lorenz-63 model. The remarkable
complexity of this simple model will also be used as a case study
to introduce a number of features of a simple ensemble filter data
assimilation system. To perform a synthetic observation
assimilation experiment for the L63 model, the following steps must
be performed (an overview of the process is given first, followed
by detailed procedures for each step):</p>
<h2 class="indent1">Experiment Overview</h2>
<ol>
<li><a href="#integrate">Integrate the L63 model for a long
time</a><br>
starting from arbitrary initial conditions to generate a model
state that lies on the attractor. The ergodic nature of the L63
system means a 'lengthy' integration always converges to some point
on the computer's finite precision representation of the model's
attractor.<br>
<br></li>
<li><a href="#ensemblate">Generate a set of ensemble initial
conditions</a><br>
from which to start an assimilation. Since L63 is ergodic, the
ensemble members can be designed to look like random samples from
the model's 'climatological distribution'. To generate an ensemble
member, very small perturbations can be introduced to the state on
the attractor generated by step 1. This perturbed state can then be
integrated for a very long time until all memory of its initial
condition can be viewed as forgotten. Any number of ensemble
initial conditions can be generated by repeating this
procedure.<br>
<br></li>
<li><a href="#simulate">Simulate a particular observing
system</a><br>
by first creating an 'observation set definition' and then creating
an 'observation sequence'. The 'observation set definition'
describes the instrumental characteristics of the observations and
the 'observation sequence' defines the temporal sequence of the
observations.<br>
<br></li>
<li><a href="#generate">Populate the 'observation sequence' with
'perfect' observations</a><br>
by integrating the model and using the information in the
'observation sequence' file to create simulated observations. This
entails operating on the model state at the time of the observation
with an appropriate forward operator (a function that operates on
the model state vector to produce the expected value of the
particular observation) and then adding a random sample from the
observation error distribution specified in the observation set
definition. At the same time, diagnostic output about the 'true'
state trajectory can be created.<br>
<br></li>
<li><a href="#assimilate">Assimilate the synthetic
observations</a><br>
by running the filter; diagnostic output is generated.</li>
</ol>
<a name="integrate" id="integrate"></a>
<h3 class="indent1">1. Integrate the L63 model for a 'long'
time.</h3>
<em class="program">perfect_model_obs</em> integrates the model for
all the times specified in the 'observation sequence definition'
file. To this end, begin by creating an 'observation sequence
definition' file that spans a long time. Creating an 'observation
sequence definition' file is a two-step procedure involving
<em class="program">create_obs_sequence</em> followed by <em class=
"program">create_fixed_network_seq</em>. After they are both run,
it is necessary to integrate the model with <em class=
"program">perfect_model_obs</em>.

<h4 class="indent1">1.1 Create an observation set definition.</h4>

<p><em class="program">create_obs_sequence</em> creates an
observation set definition, the time-independent part of an
observation sequence. An observation set definition file only
contains the <em class="code">location, type,</em> and <em class=
"code">observational error characteristics</em> (normally just the
diagonal observational error variance) for a related set of
observations. There are no actual observations, nor are there any
times associated with the definition. For spin-up, we are only
interested in integrating the L63 model, not in generating any
particular synthetic observations. Begin by creating a minimal
observation set definition.</p>

<p>In general, for the low-order models, only a single observation set
need be defined. Next, the number of individual scalar observations
(like a single surface pressure observation) in the set is needed.
To spin-up an initial condition for the L63 model, only a single
observation is needed. Next, the error variance for this
observation must be entered. Since we do not need (nor want) this
observation to have any impact on an assimilation (it will only be
used for spinning up the model and the ensemble), enter a very
large value for the error variance. An observation with a very
large error variance has essentially no impact on deterministic
filter assimilations like the default variety implemented in DART.
Finally, the location and type of the observation need to be
defined. For all types of models, the most elementary form of
synthetic observations are called 'identity' observations. These
observations are generated simply by adding a random sample from a
specified observational error distribution directly to the value of
one of the state variables. This defines the observation as being
an identity observation of the first state variable in the L63
model. The program will respond by terminating after generating a
file (generally named <em class="file">set_def.out</em>) that
defines the single identity observation of the first state variable
of the L63 model. The following is a screenshot (much of the
verbose logging has been left off for clarity), the user input
looks <em class="input">like this</em>.</p>

<div class="unix">
<pre>
[unixprompt]$ <em class="input">./create_obs_sequence</em>
 Initializing the utilities module.
 Trying to read from unit           10
 Trying to open file dart_log.out
 
 Registering module :
 $source: /home/dart/CVS.REPOS/DART/utilities/utilities_mod.f90,v $
 $revision: 1.18 $
 $date: 2004/06/29 15:16:40 $
 Registration complete.
 
 &amp;UTILITIES_NML
 TERMLEVEL= 2,LOGFILENAME=dart_log.out                                          
                                                                            
 /
 
 Registering module :
 $source: /home/dart/CVS.REPOS/DART/obs_sequence/create_obs_sequence.f90,v $
 $revision: 1.18 $
 $date: 2004/05/24 15:41:46 $
 Registration complete.

 { ... }

 Input upper bound on number of observations in sequence
<em class="input">10</em>
 
 Input number of copies of data (0 for just a definition)
<em class="input">0</em>

 Input number of quality control values per field (0 or greater)
<em class="input">0</em>

 input a -1 if there are no more obs 
<em class="input">0</em>

</pre>
<pre>
<code>
initialize_module obs_kind_nml values are

-------------- ASSIMILATE_THESE_OBS_TYPES --------------
RAW_STATE_VARIABLE
-------------- EVALUATE_THESE_OBS_TYPES --------------
------------------------------------------------------

    Input -1 * state variable index for identity observations
    OR input the name of the observation kind from table below:
    OR input the integer index, BUT see documentation...
                1 RAW_STATE_VARIABLE
</code>
</pre>
<pre>

<em class="input">-1</em>

 input time in days and seconds
<em class="input">1 0</em>

 Input error variance for this observation definition
<em class="input">1000000</em>

 input a -1 if there are no more obs 
<em class="input">-1</em>

 Input filename for sequence (  set_def.out   usually works well)
<em class="input"> set_def.out </em>
 write_obs_seq  opening formatted file set_def.out
 write_obs_seq  closed file set_def.out
</pre></div>
<h4 class="indent1">1.2 Create an observation sequence
definition.</h4>
<p><em class="program">create_fixed_network_seq</em> creates an
'observation sequence definition' by extending the 'observation set
definition' with the temporal attributes of the observations.<br>
<br>
The first input is the name of the file created in the previous
step, i.e. the name of the observation set definition that you've
just created. It is possible to create sequences in which the
observation sets are observed at regular intervals or irregularly
in time. Here, all we need is a sequence that takes observations
over a long period of time - indicated by entering a 1. Although
the L63 system normally is defined as having a non-dimensional time
step, the DART system arbitrarily defines the model timestep as
being 3600 seconds. If we declare that we have one observation per
day for 1000 days, <!-- 
By declaring we have 1000 observations
taken once per day, --> we create an observation sequence
definition spanning 24000 'model' timesteps; sufficient to spin-up
the model onto the attractor. Finally, enter a name for the
'observation sequence definition' file. Note again: there are no
observation values present in this file. Just an observation type,
location, time and the error characteristics. We are going to
populate the observation sequence with the <em class=
"program">perfect_model_obs</em> program.</p>
<div class="unix">
<pre>
[unixprompt]$ <em class="input">./create_fixed_network_seq</em>

 ...

 Registering module :
 $source: /home/dart/CVS.REPOS/DART/obs_sequence/obs_sequence_mod.f90,v $
 $revision: 1.31 $
 $date: 2004/06/29 15:04:37 $
 Registration complete.
 
 Input filename for network definition sequence (usually  set_def.out  )
<em class="input">set_def.out</em>

 ...

 To input a regularly repeating time sequence enter 1
 To enter an irregular list of times enter 2
<em class="input">1</em>
 Input number of observations in sequence
<em class="input">1000</em>
 Input time of initial ob in sequence in days and seconds
<em class="input">1, 0</em>
 Input period of obs in days and seconds
<em class="input">1, 0</em>
           1
           2
           3
...
         997
         998
         999
        1000
What is output file name for sequence (  obs_seq.in   is recommended )
<em class="input">obs_seq.in</em>
 write_obs_seq  opening formatted file obs_seq.in
 write_obs_seq closed file [blah blah blah]/work/obs_seq.in
</pre></div>
<h4>1.3 Initialize the model onto the attractor.</h4>

<p><em class="program">perfect_model_obs</em> can now advance the
arbitrary initial state for 24,000 timesteps to move it onto the
attractor.</p>

<p><em class="program">perfect_model_obs</em> uses the Fortran90
namelist input mechanism instead of (admittedly gory, but
temporary) interactive input. All of the DART software expects the
namelists to found in a file called <em class=
"file">input.nml</em>. When you built the executable, an example
namelist was created <em class=
"file">input.nml.perfect_model_obs_default</em> that contains all
of the namelist input for the executable. If you followed the
example, each namelist was saved to a unique name. We must now
rename and edit the namelist file for <em class=
"program">perfect_model_obs</em>. Copy <em class=
"file">input.nml.perfect_model_obs_default</em> to <em class=
"file">input.nml</em> and edit it to look like the following:</p>

<pre>
<code>
   &amp;perfect_model_obs_nml
   async = 0,
   adv_ens_command = "./advance_ens.csh",
   obs_seq_in_file_name = "obs_seq.in",
   obs_seq_out_file_name = "obs_seq.out",
   start_from_restart = .false.,
   output_restart = .true.,
   restart_in_file_name = "perfect_ics",
   restart_out_file_name = "perfect_restart",
   init_time_days = 0,
   init_time_seconds = 0,
   output_interval = 1 /

&amp;ensemble_manager_nml
   in_core = .true.,
   single_restart_file_in = .true.,
   single_restart_file_out = .true. /

&amp;assim_tools_nml
   filter_kind = 1,
   cutoff = 0.2,
   sort_obs_inc = .false.,
   cov_inflate = -1.0,
   cov_inflate_sd = 0.05,
   sd_lower_bound = 0.05,
   deterministic_cov_inflate = .true.,
   start_from_assim_restart = .false.,
   assim_restart_in_file_name =
'assim_tools_ics',
   assim_restart_out_file_name =
'assim_tools_restart',
   do_parallel = 0,
   num_domains = 1
   parallel_command = "./assim_filter.csh",
   spread_restoration = .false.,
   cov_inflate_upper_bound = 10000000.0,
   internal_outlier_threshold = -1.0 /

&amp;cov_cutoff_nml
   select_localization = 1 /

&amp;reg_factor_nml
   select_regression = 1,
   input_reg_file = "time_mean_reg"
   save_reg_diagnostics = .false.,
   reg_diagnostics_file = 'reg_diagnostics' /

&amp;obs_sequence_nml
   write_binary_obs_sequence = .false. /

&amp;obs_kind_nml
   assimilate_these_obs_types = 'RAW_STATE_VARIABLE' /

&amp;assim_model_nml
   write_binary_restart_files = .true. /

&amp;model_nml
   sigma = 10.0,
   r = 28.0,
   b = 2.6666666666667,
   deltat = 0.01,
   time_step_days = 0,
   time_step_seconds = 3600 /

&amp;utilities_nml
   TERMLEVEL = 1
   logfilename = 'dart_log.out' /
</code>
</pre>

<p>For the moment, only two namelists warrant explanation. Each
namelists is covered in detail in the html files accompanying the
source code for the module.</p>

<h3 class="indent1">perfect_model_obs_nml</h3>
<table border="1" class="indent1">
<tr>
<th>namelist variable</th>
<th>description</th>
</tr>
<tr>
<td><em class="code">async</em></td>
<td>For the lorenz_63, simply ignore this. Leave it set to '0'</td>
</tr>
<tr>
<td><em class="code">advance_ens_command</em></td>
<td>specifies the shell commands or script to execute when async /=
0</td>
</tr>
<tr>
<td><em class="code">obs_seq_in_file_name</em></td>
<td>specifies the file name that results from running <em class=
"program">create_fixed_network_seq</em>, i.e. the 'observation
sequence definition' file.</td>
</tr>
<tr>
<td><em class="code">obs_seq_out_file_name</em></td>
<td>specifies the output file name containing the 'observation
sequence', finally populated with (perfect?) 'observations'.</td>
</tr>
<tr>
<td><em class="code">start_from_restart</em></td>
<td>When set to 'false', <em class="program">perfect_model_obs</em>
generates an arbitrary initial condition (which cannot be
guaranteed to be on the L63 attractor).</td>
</tr>
<tr>
<td><em class="code">output_restart</em></td>
<td>When set to 'true', <em class="program">perfect_model_obs</em>
will record the model state at the end of this integration in the
file named by <em class="code">restart_out_file_name</em>.</td>
</tr>
<tr>
<td><em class="code">restart_in_file_name</em></td>
<td>is ignored when 'start_from_restart' is 'false'.</td>
</tr>
<tr>
<td><em class="code">restart_out_file_name</em></td>
<td>if <em class="code">output_restart</em> is 'true', this
specifies the name of the file containing the model state at the
end of the integration.</td>
</tr>
<tr>
<td><em class="code">init_time_</em><em class=
"italic">xxxx</em></td>
<td>the start time of the integration.</td>
</tr>
<tr>
<td><em class="code">output_interval</em></td>
<td>interval at which to save the model state.</td>
</tr>
</table>
<h3 class="indent1">utilities_nml</h3>
<table class="indent1" border="1">
<tr>
<th>namelist variable</th>
<th>description</th>
</tr>
<tr>
<td><em class="code">TERMLEVEL</em></td>
<td>When set to '1' the programs terminate when a 'warning' is
generated. When set to '2' the programs terminate only with 'fatal'
errors.</td>
</tr>
<tr>
<td><em class="code">logfilename</em></td>
<td>Run-time diagnostics are saved to this file. This namelist is
used by all programs, so the file is opened in APPEND mode.
Subsequent executions cause this file to grow.</td>
</tr>
</table>
<p>Executing <em class="program">perfect_model_obs</em> will
integrate the model 24,000 steps and output the resulting state in
the file <em class="file">perfect_restart</em>. Interested parties
can check the spinup in the <em class="file">True_State.nc</em>
file.</p>

<pre>
<code>
$ perfect_model_obs
</code>
</pre>

<a name="ensemblate" id="ensemblate"></a>
<h3 class="indent1">2. Generate a set of ensemble initial
conditions.</h3>
<p class="indent1">The set of initial conditions for a 'perfect
model' experiment is created in several steps. 1) Starting from the
spun-up state of the model (available in <em class=
"file">perfect_restart</em>), run <em class=
"program">perfect_model_obs</em> to generate the 'true state' of
the experiment and a corresponding set of observations. 2) Feed the
same initial spun-up state and resulting observations into
<em class="program">filter</em>.<br>
<br>
The first step is achieved by changing a perfect_model_obs namelist
parameter, copying <em class="file">perfect_restart</em> to
<em class="file">perfect_ics</em>, and rerunning <em class=
"program">perfect_model_obs</em>. This execution of <em class=
"program">perfect_model_obs</em> will advance the model state from
the end of the first 24,000 steps to the end of an additional
24,000 steps and place the final state in <em class=
"file">perfect_restart</em>. The rest of the namelists in
<em class="file">input.nml</em> should remain unchanged.</p>

<pre>
<code>
&amp;perfect_model_obs_nml
   async = 0,
   adv_ens_command = "./advance_ens.csh",
   obs_seq_in_file_name = "obs_seq.in",
   obs_seq_out_file_name = "obs_seq.out",
   start_from_restart = .true.,
   output_restart = .true.,
   restart_in_file_name = "perfect_ics",
   restart_out_file_name = "perfect_restart",
   init_time_days = 0,
   init_time_seconds = 0,
   output_interval = 1 /
</code>
</pre>

<pre>
<code>
$ cp perfect_restart perfect_ics
$ perfect_model_obs
</code>
</pre>

<p class="indent1">A <em class="file">True_State.nc</em> file is
also created. It contains the 'true' state of the integration.</p>
<h4 class="indent1">Generating the ensemble</h4>
<p class="indent1">This step (#2 from above) is done with the
program <em class="program">filter</em>, which also uses the
Fortran90 namelist mechanism for input. It is now necessary to copy
the <em class="file">input.nml.filter_default</em> namelist to
<em class="file">input.nml</em> or you may simply insert the
<em class="code">filter_nml</em> namelist block into the existing
<em class="file">input.nml</em>. Having the <em class=
"code">perfect_model_obs</em> namelist in the input.nml does not
hurt anything. In fact, I generally create a single <em class=
"file">input.nml</em> that has all the namelist blocks in it. I
simply copied the filter namelist block from <em class=
"file">input.nml.filter_default</em> and inserted it into our
<em class="file">input.nml</em> for the following example.</p>

<pre>
<code>
&amp;perfect_model_obs_nml
   async = 0,
   adv_ens_command = "./advance_ens.csh",
   obs_seq_in_file_name = "obs_seq.in",
   obs_seq_out_file_name = "obs_seq.out",
   start_from_restart = .true.,
   output_restart = .true.,
   restart_in_file_name = "perfect_ics",
   restart_out_file_name = "perfect_restart",
   init_time_days = 0,
   init_time_seconds = 0,
   output_interval = 1 /

&amp;filter_nml
   async = 0,
   adv_ens_command = "./advance_ens.csh",
   ens_size = 100,
   cov_inflate = 1.0,
   start_from_restart = .false.,
   output_restart = .true.,
   obs_sequence_in_name = "obs_seq.out",
   obs_sequence_out_name = "obs_seq.final",
   restart_in_file_name = "perfect_ics",
   restart_out_file_name = "filter_restart",
   init_time_days = 0,
   init_time_seconds = 0,
   output_state_ens_mean = .true.,
   output_state_ens_spread = .true.,
   output_obs_ens_mean = .true.,
   output_obs_ens_spread = .true.,
   num_output_state_members = 20,
   num_output_obs_members = 20,
   output_interval = 1,
   num_groups = 1,
   outlier_threshold = -1.0 /

&amp;ensemble_manager_nml
   in_core = .true.,
   single_restart_file_in = .true.,
   single_restart_file_out = .true. /

&amp;assim_tools_nml
   filter_kind = 1,
   cutoff = 0.2,
   sort_obs_inc = .false.,
   cov_inflate = -1.0,
   cov_inflate_sd = 0.05,
   sd_lower_bound = 0.05,
   deterministic_cov_inflate = .true.,
   start_from_assim_restart = .false.,
   assim_restart_in_file_name =
'assim_tools_ics',
   assim_restart_out_file_name =
'assim_tools_restart',
   do_parallel = 0,
   num_domains = 1
   parallel_command = "./assim_filter.csh",
   spread_restoration = .false.,
   cov_inflate_upper_bound = 10000000.0,
   internal_outlier_threshold = -1.0 /

&amp;cov_cutoff_nml
   select_localization = 1 /

&amp;reg_factor_nml
   select_regression = 1,
   input_reg_file = "time_mean_reg"
   save_reg_diagnostics = .false.,
   reg_diagnostics_file = 'reg_diagnostics' /

&amp;obs_sequence_nml
   write_binary_obs_sequence = .false. /

&amp;obs_kind_nml
   assimilate_these_obs_types = 'RAW_STATE_VARIABLE'
/

&amp;assim_model_nml
   write_binary_restart_files = .true. /

&amp;model_nml
   sigma = 10.0,
   r = 28.0,
   b = 2.6666666666667,
   deltat = 0.01,
   time_step_days = 0,
   time_step_seconds = 3600 /

&amp;utilities_nml
   TERMLEVEL = 1
   logfilename = 'dart_log.out' /
</code>
</pre>

<p class="indent1">Only the non-obvious(?) entries for <em class=
"code">filter_nml</em> will be discussed.</p>
<table border="1;" class="indent">
<tr>
<th>namelist variable</th>
<th>description</th>
</tr>
<tr>
<td valign="top"><em class="code">ens_size</em></td>
<td>Number of ensemble members. 100 is sufficient for most of the
L63 exercises.</td>
</tr>
<tr>
<td valign="top"><em class="code">cutoff</em></td>
<td>to limit the impact of an observation, set to 0.0 (i.e.
spin-up)</td>
</tr>
<tr>
<td valign="top"><em class="code">cov_inflate</em></td>
<td>A value of 1.0 results in no inflation.(spin-up)</td>
</tr>
<tr>
<td valign="top"><em class="code">start_from_restart</em></td>
<td>when '.false.', <em class="program">filter</em> will generate
its own ensemble of initial conditions. It is important to note
that the filter still makes use of <em class=
"file">perfect_ics</em> by randomly perturbing these state
variables.</td>
</tr>
<tr>
<td valign="top"><em class="code">output_state_ens_mean</em></td>
<td>when '.true.' the mean of all ensemble members is output.</td>
</tr>
<tr>
<td valign="top"><em class="code">output_state_ens_spread</em></td>
<td>when '.true.' the spread of all ensemble members is
output.</td>
</tr>
<tr>
<td valign="top"><em class=
"code">num_output_state_members</em></td>
<td>may be a value from 0 to <em class="code">ens_size</em></td>
</tr>
<tr>
<td valign="top"><em class="code">output_obs_ens_mean</em></td>
<td>when '.true.' Output ensemble mean in observation output
file.</td>
</tr>
<tr>
<td valign="top"><em class="code">output_obs_ens_spread</em></td>
<td>when '.true.' Output ensemble spread in observation output
file.</td>
</tr>
<tr>
<td valign="top"><em class="code">num_output_obs_members</em></td>
<td>may be a value from 0 to <em class="code">ens_size</em></td>
</tr>
<tr>
<td valign="top"><em class="code">output_interval</em></td>
<td>The frequency with which output state diagnostics are written.
Units are in assimilation times. Default value is 1 meaning output
is written at every observation time</td>
</tr>
</table>
<p class="indent1">The filter is told to generate its own ensemble
initial conditions since <em class="code">start_from_restart</em>
is '.false.'. However, it is important to note that the filter
still makes use of <em class="file">perfect_ics</em> which is set
to be the <em class="code">restart_in_file_name</em>. This is the
model state generated from the first 24,000 step model integration
by <em class="program">perfect_model_obs</em>. <em class=
"program">Filter</em> generates its ensemble initial conditions by
randomly perturbing the state variables of this state.</p>
<p class="indent1">The arguments <em class=
"code">output_state_ens_mean</em> and <em class=
"code">output_state_ens_spread</em> are '.true.' so that these
quantities are output at every time for which there are
observations (once a day here) and <em class=
"code">num_output_ens_members</em> means that the same diagnostic
files, <em class="file">Posterior_Diag.nc</em> and <em class=
"file">Prior_Diag.nc</em> also contain values for 20 ensemble
members once a day. Once the namelist is set, execute <em class=
"program">filter</em> to integrate the ensemble forward for 24,000
steps with the final ensemble state written to the <em class=
"file">filter_restart</em>. Copy the <em class=
"program">perfect_model_obs</em> restart file <em class=
"file">perfect_restart</em> (the `true state') to <em class=
"file">perfect_ics</em>, and the <em class="program">filter</em>
restart file <em class="file">filter_restart</em> to <em class=
"file">filter_ics</em> so that future assimilation experiments can
be initialized from these spun-up states.</p>

<pre>
<code>
$ filter
$ cp perfect_restart perfect_ics
$ cp filter_restart filter_ics
</code>
</pre>

<p class="indent1">The spin-up of the ensemble can be viewed by
examining the output in the netCDF files <em class=
"file">True_State.nc</em> generated by <em class=
"program">perfect_model_obs</em> and <em class=
"file">Posterior_Diag.nc</em> and <em class=
"file">Prior_Diag.nc</em> generated by <em class=
"program">filter</em>. To do this, see the detailed discussion of
matlab diagnostics in Appendix I. <a name="simulate" id=
"simulate"></a></p>
<h3 class="indent1">3. Simulate a particular observing system.</h3>
<p class="indent1">Begin by using <em class=
"program">create_obs_sequence</em> to generate an observation set
in which each of the 3 state variables of L63 is observed with an
observational error variance of 1.0 for each observation. To do
this, use the following input sequence (the text including and
after # is a comment and does not need to be entered):</p>
<table class="indent1" bgcolor="#CCCCCC">
<tr>
<td><em class="input">4</em></td>
<td># upper bound on num of observations in sequence</td>
</tr>
<tr>
<td><em class="input">0</em></td>
<td># number of copies of data (0 for just a definition)</td>
</tr>
<tr>
<td><em class="input">0</em></td>
<td># number of quality control values per field (0 or
greater)</td>
</tr>
<tr>
<td><em class="input">0</em></td>
<td># -1 to exit/end observation definitions</td>
</tr>
<!-- observe state variable 1 -->
<tr>
<td><em class="input">-1</em></td>
<td># observe state variable 1</td>
</tr>
<tr>
<td><em class="input">0   0</em></td>
<td># time -- days, seconds</td>
</tr>
<tr>
<td><em class="input">1.0</em></td>
<td># observational variance</td>
</tr>
<tr>
<td><em class="input">0</em></td>
<td># -1 to exit/end observation definitions</td>
</tr>
<!-- observe state variable 2 -->
<tr>
<td><em class="input">-2</em></td>
<td># observe state variable 2</td>
</tr>
<tr>
<td><em class="input">0   0</em></td>
<td># time -- days, seconds</td>
</tr>
<tr>
<td><em class="input">1.0</em></td>
<td># observational variance</td>
</tr>
<tr>
<td><em class="input">0</em></td>
<td># -1 to exit/end observation definitions</td>
</tr>
<!-- observe state variable 3 -->
<tr>
<td><em class="input">-3</em></td>
<td># observe state variable 3</td>
</tr>
<tr>
<td><em class="input">0   0</em></td>
<td># time -- days, seconds</td>
</tr>
<tr>
<td><em class="input">1.0</em></td>
<td># observational variance</td>
</tr>
<tr>
<td><em class="input">-1</em></td>
<td># -1 to exit/end observation definitions</td>
</tr>
<tr>
<td><em class="input">set_def.out</em></td>
<td># Output file name</td>
</tr>
</table>
<p class="indent1">Now, generate an observation sequence definition
by running <em class="program">create_fixed_network_seq</em> with
the following input sequence:</p>
<table class="indent1" bgcolor="#CCCCCC">
<tr>
<td><em class="input">set_def.out</em></td>
<td># Input observation set definition file</td>
</tr>
<tr>
<td><em class="input">1</em></td>
<td># Regular spaced observation interval in time</td>
</tr>
<tr>
<td><em class="input">1000</em></td>
<td># 1000 observation times</td>
</tr>
<tr>
<td><em class="input">0, 43200</em></td>
<td># First observation after 12 hours (0 days, 12 * 3600
seconds)</td>
</tr>
<tr>
<td><em class="input">0, 43200</em></td>
<td># Observations every 12 hours</td>
</tr>
<tr>
<td><em class="input">obs_seq.in</em></td>
<td># Output file for observation sequence definition</td>
</tr>
</table>
<a name="generate" id="generate"></a>
<h3 class="indent1">4. Generate a particular observing system and
true state.</h3>
<p class="indent1">An observation sequence file is now generated by
running <em class="program">perfect_model_obs</em> with the
namelist values (unchanged from step 2):</p>

<pre>
<code>
&amp;perfect_model_obs_nml
   async = 0,
   adv_ens_command = "./advance_ens.csh",
   obs_seq_in_file_name = "obs_seq.in",
   obs_seq_out_file_name = "obs_seq.out",
   start_from_restart = .true.,
   output_restart = .true.,
   restart_in_file_name = "perfect_ics",
   restart_out_file_name = "perfect_restart",
   init_time_days = 0,
   init_time_seconds = 0,
   output_interval = 1 /
</code>
</pre>

<p class="indent1">This integrates the model starting from the
state in <em class="file">perfect_ics</em> for 1000 12-hour
intervals outputting synthetic observations of the three state
variables every 12 hours and producing a netCDF diagnostic file,
<em class="file">True_State.nc</em>.</p>
<a name="assimilate" id="assimilate"></a>
<h3 class="indent1">5. Filtering.</h3>
<p class="indent1">Finally, <em class="program">filter</em> can be
run with its namelist set to:</p>

<pre>
<code>
&amp;filter_nml
   async = 0,
   adv_ens_command = "./advance_ens.csh",
   ens_size = 100,
   cov_inflate = 1.0,
   start_from_restart = .true.,
   output_restart = .true.,
   obs_sequence_in_name = "obs_seq.out",
   obs_sequence_out_name = "obs_seq.final",
   restart_in_file_name = "filter_ics",
   restart_out_file_name = "filter_restart",
   init_time_days = 0,
   init_time_seconds = 0,
   output_state_ens_mean = .true.,
   output_state_ens_spread = .true.,
   output_obs_ens_mean = .true.,
   output_obs_ens_spread = .true.,
   num_output_state_members = 20,
   num_output_obs_members = 20,
   output_interval = 1,
   num_groups = 1,
   outlier_threshold = -1.0 /
</code>
</pre>

<p class="indent1">The large value for the cutoff allows each
observation to impact all other state variables (see Appendix V for
localization). <em class="program">filter</em> produces two output
diagnostic files, <em class="file">Prior_Diag.nc</em> which
contains values of the ensemble mean, ensemble spread, and ensemble
members for 12- hour lead forecasts before assimilation is applied
and <em class="file">Posterior_Diag.nc</em> which contains similar
data for after the assimilation is applied (sometimes referred to
as analysis values).</p>
Now try applying all of the matlab diagnostic functions described
in <a href="#matlab">the Matlab Diagnostics section</a>. 
<!--==================================================================-->
<a name="matlab" id="matlab"></a> 
<!--==================================================================-->
<hr>
<h2>Matlab® Diagnostics</h2>
<p>The output files are netCDF files, and may be examined with many
different software packages. We happen to use Matlab®, and provide
our diagnostic scripts in the hopes that they are useful.</p>
<p>The diagnostic scripts and underlying functions reside in two
places: <em class="file">DART/diagnostics/matlab</em> and
<em class="file">DART/matlab</em>. They are reliant on the
public-domain <a href=
"http://woodshole.er.usgs.gov/staffpages/cdenham/public_html/MexCDF/nc4ml5.html">
netcdf toolbox</a> from <em class=
"file">http://woodshole.er.usgs.gov/staffpages/cdenham/public_html/MexCDF/nc4ml5.html</em>
as well as the public-domain <a href=
"http://www.marine.csiro.au/sw/matlab-netcdf.html">CSIRO
matlab/netCDF interface</a> from <em class=
"file">http://www.marine.csiro.au/sw/matlab-netcdf.html</em>. If
you do not have them installed on your system and want to use
Matlab to peruse netCDF, you must follow their installation
instructions. The 'interested reader' may want to look at the
<em class="code">DART/matlab/startup.m</em> file I use on my
system. If you put it in your <em class="unix">$HOME/matlab</em>
directory, it is invoked every time you start up Matlab.</p>
<p>Once you can access the <em class="program">getnc</em> function
from within Matlab, you can use our diagnostic scripts. It is
necessary to prepend the location of the <em class=
"file">DART/matlab</em> scripts to the <em class=
"code">matlabpath</em>. Keep in mind the location of the netcdf
operators on your system WILL be different from ours ... and that's
OK.<br>
<br></p>
<div class="unix">
<pre>
0[269]0 ghotiol:/&lt;5&gt;models/lorenz_63/work]$ matlab -nojvm

                                             &lt; M A T L A B &gt;
                                 Copyright 1984-2002 The MathWorks, Inc.
                                     Version 6.5.0.180913a Release 13
                                               Jun 18 2002

  Using Toolbox Path Cache.  Type "help toolbox_path_cache" for more info.
 
  To get started, type one of these: helpwin, helpdesk, or demo.
  For product information, visit www.mathworks.com.

&gt;&gt;<em class="input"> which getnc</em>
/contrib/matlab/matlab_netcdf_5_0/getnc.m
&gt;&gt;<em class="input">ls *.nc</em>

ans =

Posterior_Diag.nc  Prior_Diag.nc  True_State.nc


&gt;&gt;<em class="input">path('../../../matlab',path)</em>
&gt;&gt;<em class=
"input">path('../../../diagnostics/matlab',path)</em>
&gt;&gt;<em class="input">which plot_ens_err_spread</em>
../../../matlab/plot_ens_err_spread.m
&gt;&gt;<em class="input">help plot_ens_err_spread</em>

  DART : Plots summary plots of the ensemble error and ensemble spread.
                         Interactively queries for the needed information.
                         Since different models potentially need different 
                         pieces of information ... the model types are 
                         determined and additional user input may be queried.
 
  Ultimately, plot_ens_err_spread will be replaced by a GUI.
  All the heavy lifting is done by PlotEnsErrSpread.
 
  Example 1 (for low-order models)
 
  truth_file = 'True_State.nc';
  diagn_file = 'Prior_Diag.nc';
  plot_ens_err_spread

&gt;&gt;<em class="input">plot_ens_err_spread</em>
</pre></div>
<p>And the matlab graphics window will display the spread of the
ensemble error for each state variable. The scripts are designed to
do the "obvious" thing for the low-order models and will prompt for
additional information if needed. The philosophy of these is that
anything that starts with a lower-case <em class=
"file">plot_<em class="italic">some_specific_task</em></em> is
intended to be user-callable and should handle any of the models.
All the other routines in <em class="file">DART/matlab</em> are
called BY the high-level routines.</p>
<table border="1;" class="indent1">
<tr>
<th>Matlab script</th>
<th>description</th>
</tr>
<tr>
<td><em class="code">plot_bins</em></td>
<td>plots ensemble rank histograms</td>
</tr>
<tr>
<td><em class="code">plot_correl</em></td>
<td>Plots space-time series of correlation between a given variable
at a given time and other variables at all times in a n ensemble
time sequence.</td>
</tr>
<tr>
<td><em class="code">plot_ens_err_spread</em></td>
<td>Plots summary plots of the ensemble error and ensemble spread.
Interactively queries for the needed information. Since different
models potentially need different pieces of information ... the
model types are determined and additional user input may be
queried.</td>
</tr>
<tr>
<td><em class="code">plot_ens_mean_time_series</em></td>
<td>Queries for the state variables to plot.</td>
</tr>
<tr>
<td><em class="code">plot_ens_time_series</em></td>
<td>Queries for the state variables to plot.</td>
</tr>
<tr>
<td><em class="code">plot_phase_space</em></td>
<td>Plots a 3D trajectory of (3 state variables of) a single
ensemble member. Additional trajectories may be superimposed.</td>
</tr>
<tr>
<td><em class="code">plot_total_err</em></td>
<td>Summary plots of global error and spread.</td>
</tr>
<tr>
<td><em class="code">plot_var_var_correl</em></td>
<td>Plots time series of correlation between a given variable at a
given time and another variable at all times in an ensemble time
sequence.</td>
</tr>
</table>
<!--==================================================================-->
<a name="discussion" id="discussion"></a> 
<!--==================================================================-->
<hr>
<h2>Bias, filter divergence and covariance inflation (with the L63
model)</h2>
One of the common problems with ensemble filters is filter
divergence, which can also be an issue with a variety of other
flavors of filters including the classical Kalman filter. In filter
divergence, the prior estimate of the model state becomes too
confident, either by chance or because of errors in the forecast
model, the observational error characteristics, or approximations
in the filter itself. If the filter is inappropriately confident
that its prior estimate is correct, it will then tend to give less
weight to observations than they should be given. The result can be
enhanced overconfidence in the model's state estimate. In severe
cases, this can spiral out of control and the ensemble can wander
entirely away from the truth, confident that it is correct in its
estimate. In less severe cases, the ensemble estimates may not
diverge entirely from the truth but may still be too confident in
their estimate. The result is that the truth ends up being farther
away from the filter estimates than the spread of the filter
ensemble would estimate. This type of behavior is commonly detected
using rank histograms (also known as Talagrand diagrams). You can
see the rank histograms for the L63 initial assimilation by using
the matlab script <em class="program">plot_bins</em>.
<p>A simple, but surprisingly effective way of dealing with filter
divergence is known as covariance inflation. In this method, the
prior ensemble estimate of the state is expanded around its mean by
a constant factor, effectively increasing the prior estimate of
uncertainty while leaving the prior mean estimate unchanged. The
program <em class="program">filter</em> has a namelist parameter
that controls the application of covariance inflation, <em class=
"code">cov_inflate</em>. Up to this point, <em class=
"code">cov_inflate</em> has been set to 1.0 indicating that the
prior ensemble is left unchanged. Increasing <em class=
"code">cov_inflate</em> to values greater than 1.0 inflates the
ensemble before assimilating observations at each time they are
available. Values smaller than 1.0 contract (reduce the spread) of
prior ensembles before assimilating.</p>
<p>You can do this by modifying the value of <em class=
"code">cov_inflate</em> in the namelist, (try 1.05 and 1.10 and
other values at your discretion) and run the filter as above. In
each case, use the diagnostic matlab tools to examine the resulting
changes to the error, the ensemble spread (via rank histogram bins,
too), etc. What kind of relation between spread and error is seen
in this model?</p>
<!--==================================================================-->
<a name="syntheticobservations" id="syntheticobservations"></a> 
<!--==================================================================-->
<hr>
<h2>Synthetic Observations</h2>
<p>Synthetic observations are generated from a `perfect' model
integration, which is often referred to as the `truth' or a `nature
run'. A model is integrated forward from some set of initial
conditions and observations are generated as <em class="equation">y
= H(x) + e</em> where <em class="equation">H</em> is an operator on
the model state vector, <em class="equation">x</em>, that gives the
expected value of a set of observations, <em class=
"equation">y</em>, and <em class="equation">e</em> is a random
variable with a distribution describing the error characteristics
of the observing instrument(s) being simulated. Using synthetic
observations in this way allows students to learn about
assimilation algorithms while being isolated from the additional
(extreme) complexity associated with model error and unknown
observational error characteristics. In other words, for the
real-world assimilation problem, the model has (often substantial)
differences from what happens in the real system and the
observational error distribution may be very complicated and is
certainly not well known. Be careful to keep these issues in mind
while exploring the capabilities of the ensemble filters with
synthetic observations.</p>
<!--==================================================================-->
<!-- Legalese & Metadata                                              -->
<!--==================================================================-->
<a name="Legalese" id="Legalese"></a>
<hr>
<h2>Terms of Use</h2>
<p>DART software - Copyright UCAR. This open source software is
provided by UCAR, "as is", without charge, subject to all terms of
use at <a href=
"http://www.image.ucar.edu/DAReS/DART/DART_download">http://www.image.ucar.edu/DAReS/DART/DART_download</a></p>
<!--==================================================================-->
</body>
</html>
