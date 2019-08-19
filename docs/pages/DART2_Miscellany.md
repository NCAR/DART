---
title: Miscellany
layout: default
---

<span id="useful_software" class="anchor"></span>

### Useful Software

The following free open-source tools have proven to be very useful:

1.  [ncview](http://meteora.ucsd.edu/~pierce/ncview_home_page.html): a
    great visual browser for netCDF files.
2.  [Panoply](http://www.giss.nasa.gov/tools/panoply): another visual
    browser for netCDF, HDF, and GRIB files with many options for map
    projections and data slicing.
3.  [the netCDF Operators (NCO)](http://nco.sourceforge.net/): tools to
    perform operations on netCDF files like concatenating, differencing,
    averaging, etc.
4.  [An MPI
    environment](http://en.wikipedia.org/wiki/Message_Passing_Interface#Overview):
    to run larger jobs in parallel. DART can be used without MPI,
    especially the low order models where the memory use is small. The
    larger models often require MPI so that filter can be run as a
    parallel job, both for speed and memory size reasons. Common options
    are [OpenMPI](http://www.open-mpi.org/) or
    [MPICH](http://www.mpich.org/). See the DART MPI introduction in
    [mpi_intro.html](dart_mpi.html).
5.  [Observation Processing And Wind Synthesis
    (OPAWS)](http://code.google.com/p/opaws/): OPAWS can process NCAR
    Dorade (sweep) and NCAR EOL Foray (netcdf) radar data. It analyzes
    (grids) data in either two-dimensions (on the conical surface of
    each sweep) or three-dimensions (Cartesian). Analyses are output in
    netcdf, Vis5d, and/or DART (Data Assimilation Research Testbed)
    formats.

The following licensed (commercial) tools have proven to be very useful:

1.  [MATLAB®](http://www.mathworks.com/products/matlab/): An interactive
    and programming language for computation and visualization. We
    supply our diagnostic and plotting routines as MATLAB® scripts.  
    Free alternatives to MATLAB® (for which we unfortunately do not have
    the resources to support, but would happily accept user
    contributions) include:
    1.  [Octave](http://www.gnu.org/software/octave)
    2.  [SciPy](http://www.scipy.org/) plus
    3.  [matplotlib](http://matplotlib.org/)
    4.  The [R](http://www.r-project.org/) programming language has
        similiar functionality but a different enough syntax that the
        diagnostic and plotting routines we supply which work with
        MATLAB® are unlikely to be easy to port.
    5.  The [NCAR Command Language
        (NCL)](http://www.ncl.ucar.edu/Document/Manuals/Getting_Started/introduction.shtml)
        package also does computation and plotting, but again with a
        different enough syntax that the routines supplied with DART
        would need major rewriting to work.

<span id="RegressionTestResults" class="anchor"></span>

\[[top](#)\]

-----

## DART platforms/compilers/batch systems

We work to keep the DART code highly portable. We avoid
compiler-specific constructs, require no system-specific functions, and
try as much as possible to be easy to build on new platforms.

DART has been compiled and run on Apple laptops and workstations, Linux
clusters small and large, SGI Altix systems, IBM Power systems, IBM
Intel systems, Cray systems.

DART has been compiled with compilers from Intel, PGI, Cray, GNU, IBM,
Pathscale.

MPI versions of DART have run under batch systems including LSF, PBS,
Moab/Torque, and Sun Grid Engine.

<span id="platform_notes" class="anchor"></span> 

\[[top](#)\]

-----

## Platform-specific notes.

Most of the platform-specific notes are in the appropriate
*mkmf.template.xxxx.yyyy* file. There are very few situations that
require making additional changes.

### gfortran

For some reason, the *gfortran* compiler does not require an *interface*
to the *system()* routine while all the other compilers we have tested
**do** need the interface. This makes it impossible to have a module
that is compiler-independent. The interface is needed in the
*null_mpi_utilities_mod.f90*, and/or *mpi_utilities_mod.f90*. The
problem surfaces at link time :

~~~
null_mpi_utilities_mod.o(.text+0x160): In function `__mpi_utilities_mod__shell_execute':
: undefined reference to `system_'
null_mpi_utilities_mod.o(.text+0x7c8): In function `__mpi_utilities_mod__destroy_pipe':
: undefined reference to `system_'
null_mpi_utilities_mod.o(.text+0xbb9): In function `__mpi_utilities_mod__make_pipe':
: undefined reference to `system_'
collect2: ld returned 1 exit status
make: *** [preprocess] Error 1
~~~

There is a script to facilitate making the appropriate change to
*null_mpi_utilities_mod.f90* and *mpi_utilities_mod.f90*. Run the
shell script *DART/mpi_utilities/fixsystem* with no arguments to simply
'flip' the state of these files (i.e. if the system block is defined, it
will undefine the block by commenting it out; if the block is commented
out, it will define it by uncommenting the block). If you want to
hand-edit *null_mpi_utilities_mod.f90* and *mpi_utilities_mod.f90*
- look for the comment block that starts *\! BUILD TIP* and follow the
directions in the comment block.  

### module mismatch errors

Compilers create modules in their own particular manner ... a module
built by one compiler may not (will usually not) be useable by another
compiler. Sometimes it happens that the Fortran90 modules for the netCDF
interface compiled by compiler *A* is trying to be used by compiler *B*.
This generally results in an error message like:

~~~ 
Fatal Error: File 'netcdf.mod' opened at (1) is not a <pick_your_compiler> module file
make: *** [utilities_mod.o] Error 1
~~~

The only solution here is to make sure the *mkmf.template* file is
referencing the appropriate netCDF installation.

### *endian*-ness errors

The *endian*-ness of the binary files is specific to the chipset, not
the compiler or the code (normally). There are some models that require
a specific *endian* binary file. Most compilers have some sort of
ability to read and/or write binary files of a specific (or non-native)
endianness by throwing some compile flags. It is generally an
'all-or-nothing' approach in that trying to micromanage which files are
opened with native endianness and which files are openened with the
non-native endianness is generally too time-consuming and fraught with
error to be of much use. If the compile flags exist and are known to us,
we try to include them in the comment section of the individual
*mkmf.template.xxxx.yyyy* file.  
  
These errors most often manifest themselves as 'time' errors in the DART
execution. The restart/initial conditions files have the valid time of
the ensuing model state as the first bit of information in the header,
and if these files are 'wrong'-endian, the encoded times are
nonsensical.

### MPI

If you want to use MPI and are interested in testing something simple
before total immersion: try running the MPI test routines in the
*DART/doc/mpi* directory. This directory contains some small test
programs which use both MPI and the netCDF libraries. It may be simpler
to debug any build problems here, and if you need to submit a problem
report to your system admin people these single executables are much
simpler than the entire DART build tree.

<span id="FAQ" class="anchor"></span> 

\[[top](#)\]

-----

<span id="FAQ"></span>

## Frequently Asked Questions for DART

  - [General Info about DART](#General)
  - [Installation Questions](#Install)
  - [Questions about Using DART](#Using)
  - [The WRF Weather Model and DART](#WRFDART)
  - [The CESM Climate Model and DART](#CESMDART)

### General Info about DART
  
<span id="General"></span>

#### What kind of data assimilation does DART do?

There are two main techniques for doing data assimilation: variational
and ensemble methods. DART uses a variety of ensemble Kalman filter
techniques.

#### What parts of the DART source should I feel free to alter?

We distribute the full source code for the system so you're free to edit
anything you please. However, the system was designed so that you should
be able to add code in a few specific places to add a new model, work
with new observation types, or change the assimilation algorithm.  
  
To add a new model you should be able to add a new
`DART/models/XXX/model_mod.f90` file to interface between your model and
DART. We expect that you should not have to alter any code in your model
to make it work with DART.  
  
To add new observation types you should be able to add a new
`DART/obs_def/obs_def_XXX_mod.f90` file. If there is not already a
converter for this observation type you can add a converter in
`DART/observations/XXX`.  
  
If you are doing data assimilation algorithm research you may be
altering some of the core DART routines in the
`DART/assim_tools/assim_tools_mod.f90` or `DART/filter/filter.f90`
files. Please feel free to email DART support (dart at ucar.edu) for
help with how to do these modifications so they work with the parallel
version of DART correctly.  
  
If you add support for a new observation type, a new model, or filter
kind, we'd love for you to send a copy of it back to us for inclusion in
the DART
distribution.

#### What systems and compilers do you support? What other tools do I need?

We run on almost any Linux-based system including laptops, clusters, and
supercomputers. This includes IBMs, Crays, SGIs, Macs. We discourage
trying to use Windows but it has been done using the CygWin package.  
  
We require a Fortran 90 compiler. Common ones in use are from GNU
(gfortran), Intel, PGI, PathScale, IBM, and g95.  
  
We need a compatible NetCDF library, which means compiled with the same
compiler you build DART with, and built with the Fortran interfaces.  
  
You can run DART as a single program without any additional software. To
run in parallel on a cluster or other multicore platform you will need a
working MPI library and runtime system. If one doesn't come with your
system already OpenMPI is a good open-source option.  
  
Our diagnostic routines are MATLAB® scripts, which is a commercial
math/visualization package. Some users use IDL, NCL, or R but they have
to adapt our scripts themselves.

### Installation Questions

<span id="Install"></span>

#### How do I get started?

Go to the extensive [DART web pages](../index.html) where there are
detailed instructions on checking the source out of our subversion
server, compiling, running the tutorials, and examples of other users'
applications of DART.  
  
If you really hate reading instructions you can try looking at the
README in the top level directory. But if you run into problems please
read the [full setup instructions](DART2_Starting.html#installing)
before contacting us for help. We will start out suggesting you read
those web pages first anyway.

#### I'm trying to build with MPI and getting errors.

The MPI compiler commands are usually scripts or programs which add
additional arguments and then call the standard Fortran compiler. If
there is more than one type of compiler on a system you must find the
version of MPI which was compiled to wrap around the compiler you are
using.  
  
In the `DART/developer_tests/mpi_utilities/tests` directory are some
small programs which can be used to test compiling and running with
MPI.  
  
If you are using version 1.10.0 of OpenMPI and getting compiler errors
about being unable to find a matching routine for calls to MPI_Get()
and/or MPI_Reduce(), please update to version 1.10.1 or later. There
were missing interfaces in the 1.10.0 release which are fixed in the
1.10.1
release.

#### I'm getting errors from subversion when I'm trying to check out a copy of DART.

If you are behind some kind of firewall they may not allow the ports
needed to talk to the subversion server. If you can, try from a machine
outside the firewall, or talk to your system support people about how to
access a subversion server. Sometimes there are machines which allow
subversion access which share filesystems with machines that are not
allowed. It's always better if you can go back later and update your
copy of DART with subversion to keep it in sync with the server; making
a tar of the checked out source and moving it to another machine won't
let you do
this.

#### I'm getting errors related to NetCDF when I try to build the executables.

Any application that uses the NetCDF data libraries must be compiled
with exactly the same compiler as the libraries were built with. On
systems which have either multiple compilers, or multiple versions of
the same compiler, there is the possibilty that the libraries don't
match the compiler you're using to compile DART. Options here are:

  - If there are multiple versions of the NetCDF libraries installed,
    find a method to select the right version; e.g. specify the exact
    path to the include files and libraries in your mkmf.template file,
    or load the right module if your system uses the 'module' command to
    select software options.
  - Change the version of the compiler you are using to build DART to
    match the one used to build NetCDF.
  - Build your own version of the NetCDF libraries with the compiler you
    prefer to use. See [this web
    page](http://www.unidata.ucar.edu/software/netcdf/docs/building_netcdf_fortran.html)
    for help in building the libraries. DART requires only the basic
    library with the NetCDF 3 interfaces, but will work with NetCDF 4
    versions. Building NetCDF 4 does require additional libraries such
    as HDF, libz, etc.

If you believe you are using the right version of the compiler, then
check to see if the Fortran interfaces have been compiled into a single
library with the C code, or if there are two libraries, libnetcdf.a and
libnetcdff.a (note the 2 f's in the second library). The library lines
in your mkmf.template must reference either one or both libraries,
depending on what exists. This is a choice that is made by the person
who built the NetCDF libraries and cannot be predicted
beforehand.

#### I'm getting errors about undefined symbol "_system_" when I try to compile.

If you're running the Lanai release or code from the trunk later than
2013, the DART Makefiles should automatically call a script in the
`DART/mpi_utilities` directory named `fixsystem`. This script tries to
alter the MPI source code in that directory to work with your compiler.
If you still get a compiler error look at this script and see if you
have to add a case for the name of your compiler.  
  
If you're running the Kodiak release or earlier, you have to run
`fixsystem` yourself before compiling. We distributed the code so it
would work without change for the `gfortran` compiler, but all other
compilers require that you run `fixsystem` before trying to
compile.  

#### I have a netCDF library but I'm getting errors about unrecognized module format.

The netCDF libraries need to be built by the same version of the same
compiler as you are building DART. If your system has more than one
compiler on it (e.g. intel ifort and gfortran) or multiple versions of
the same compiler (e.g. gfortran 4.1 and 4.5) you must have a version of
the netCDF libraries which was built with the same version of the same
compiler as you're using to build
DART.

#### I'm getting errors about `-lnetcdff` not found, or I'm getting errors about undefined symbols for things in the netCDF libraries.

There are several important options when the netCDF libraries are built
that change what libraries you get and whether you have what you need.
The problems we run into most frequently are:

  - If the netCDF installation includes only the C library routines and
    not the Fortran interfaces. We require both.
  - The C routines are always in `-lnetcdf`, but the Fortran interfaces
    can either be included in that single library or placed in a
    separate `-lnetcdff` library (note 2 f's).
  - If HDF support is included, additional libraries are required to
    link an executable. Most of our `mkmf` template files have comments
    about the usual list of required libraries that you need to include.

Bottom line: What you need to set for the library list in your
`DART/mkmf/mkmf.template` file depends on how your netCDF was
built.
  
#### My model runs in single precision and I want to compile DART the same way.

We recommend that you run an assimilation with Fortran double precision
reals (e.g. all real values are real\*8 or 64 bits). However if your
model is compiled in single precision (real\*4 or 32 bits) then there is
an option to build DART the same way. Edit `DART/common/types_mod.f90`
and change the definition of `R8` to equal `R4` (comment out the
existing line and comment in the following line). Rebuild all DART
executables and it will run with single precision reals. We declare
every real variable inside DART with an explicit size, so we do not
recommend using compiler flags to try to change the default real
variable precision because it will not affect the DART
code.

### Questions about Using DART

<span id="Using"></span>

#### I'm trying to run an MPI filter and I'm getting N copies of every message.

Look in the log or in the standard output for the message:
'initialize_mpi_utilities: Running with N MPI processes.' Instead of
this message, if you see: `initialize_mpi_utilities: Running single
process` then you have NOT successfully compiled with MPI; you are
running N duplicate copies of a single-task program. Rerun the
quickbuild.csh script with the -mpi flag to force it to build filter
with mpif90 or whatever the mpi compiler wrapper is called on your
system.

#### How does DART interact with running my model?

If you are running one of the "low-order" models (e.g. one of the Lorenz
models, the null model, the pe2lyr model, etc), the easiest way to run
is to let DART control advancing the model when necessary. You run the
"filter" executable and it runs both the assimilation and model advances
until all observations in the input observation sequence file have been
assimilated. See the "async" setting in the [filter namelist
documentation](https://ncar.github.io/DART/api/v2.1.10/program/filter.html)
for more information.  
  
If you are running a large model with a complicated configuration and/or
run script, you will probably want to run the assimilation separately
from the model advances. To do this, you will need to script the
execution, and break up the observations into single timestep chunks per
file. The scripting will need to create filter input files from the
model files, link the current observation file to the input filename in
the namelist, copy or rename any inflation files from the previous
assimilation step, run filter, convert the filter output to model input
files, and then run the model. There are example scripts which do this
in the WRF shell_scripts directory, also the MPAS shell_scripts
directory. These scripts are both highly model-dependent as well as
computing system dependent.  
  
If you are running any of the CESM models (e.g. CAM, POP, CLM) then the
scripts to set up a CESM case with assimilation are provided in the DART
distribution. At run time, the run script provided by CESM is used.
After the model advance a DART script is called to do the assimilation.
The "multi-instance" capability of CESM is used to manage the multiple
copies of the components which are needed for assimilation, and to run
them all as part of a single job.

#### After assimilating, my model variables are out of range.

One of the assumptions of the Kalman filter is that the model states and
the observation values have gaussian distributions. The assimilation can
work successfully even if this is not actually true but there are
certain cases where this leads to problems.  
  
If any of the model state values must remain bounded, for example values
which must remain positive, or must remain between 0 and 1, you may have
to add some additional code to ensure the posterior values obey these
constraints. It is not an indication of an error if after the
assimilation some values are outside the required range.  
  
Most users deal with this, successfully, by letting the assimilation
update the values as it will, and then during the step where the model
data is converted from DART format to the model native format, any
out-of-range values are changed at that time. For example, the WRF model
has a namelist item in the &model_nml namelist which can be set at
run-time to list which variables have minimum and/or maximum values and
the conversion code will enforce the given limits.  
  
Generally this works successfully, but if the observations or the model
are biased and the assimilation is continuously trying to move model
state out of range, the distribution can become seriously unbalanced. In
this case another solution, which requires more coding, might be to
convert the values to a log scale on import to DART, do the assimilation
with the log of the observation values, and then convert back to the
original scale at export time. This ensures the values stay positive,
which is common requirement for legal
values.

#### After assimilating, my job finished successfully but no values changed.

This is a common problem, especially when adding a new observation type
or trying to assimilate with a new model. But it can happen at any time
and can be confusing about why nothing is changing. See [this web
page](DART2_Diagnostics.html#DidItWork) for a list of common causes of
the assimilation output state being the same as the input state, and how
to determine which one is responsible.

#### You have lots of namelists. How can I tell what to set?

Each module in DART has an html web page which describes the namelists
in detail. Start with `DART/index.html` and follow the links to all the
other modules and namelists in the system. If you want help with setting
up an experiment the `DART/filter/filter.html` page has some
introductory advice for some of the more important namelist settings.

#### I'm not getting an error but I am getting MPI timeouts

If your job is getting killed for no discernable reason but is usually
during computing prior or posterior forward operators, or during writing
the diagnostics file, the problem may be caused by the MPI timeout
limit. This usually happens only when the number of MPI tasks is much
larger than the number of ensemble members, and there are very slow
forward operator computations or very large states to write into the
diagnostics files. In the standard DART distribution only the first N
tasks (where N is the number of ensemble members) are doing work during
the forward operators, or only 1 task for writing diagnostic files. All
the other tasks will be waiting at an MPI barrier. If they wait there
long enough they reach the timeout threshold which assumes that at least
one or more other tasks have failed and so they exit.  
  
The solutions are either to set an environment variable that lengthens
the timeout threshold, run with fewer MPI tasks, or ask the DART team to
be a Beta user of a newer version of DART which does not have such large
time differentials between different MPI tasks.

#### filter is finishing but my job is hanging at exit

If filter finishes running, including the final timestamp message to the
log file, but then the MPI job does not exit (the next line in the job
script is not reached), and you have set the MPI timeout to be large to
avoid the job being killed by MPI timeouts, then you have run into a bug
we also have seen. We believe this to be an MPI library bug which only
happens under a specific set of circumstances. We can reproduce it but
cannot find a solution. The apparent bug happens more frequently with
larger processor counts (usually larger than about 4000 MPI tasks), so
if you run into this situation try running with a smaller MPI task count
if possible, and not setting the MPI debug flags. We have seen this
happen on the NCAR supercomputer Yellowstone with both the MPICH2 and
PEMPI MPI libraries.

### The WRF Weather Model and DART

<span id="WRFDART"></span>

#### How does DART interact with running WRF?

Most users with large WRF domains run a single cycle of filter to do
assimilation, and then advance each ensemble member of WRF from a
script, possibly submitting them in a batch to the job queues.  
  
For smaller WRF runs, if WRF can be compiled without MPI (the 'serial'
configuration) then filter can cycle inside the same program, advancing
multiple ensemble members in parallel. See the WRF documentation pages
for more
details.

#### I have completed running *filter* and I have the *filter_restart.\#\#\#\#* files. Can you refer me to the utility to convert them back to a set of *wrfinput_d01* files?

If you are using the *advance_model.csh* script that is distributed
with DART, it will take care of converting the filter output files back
to the WRF input files for the next model advance.  
  
If you are setting up a free run or doing something different than what
the basic script supports, read on to see what must be done.  
  
When you finish running DART it will have created a set of
*sssss.\#\#\#\#* restart files, where the *sssss* part of the filename
comes from the setting of *&filter_nml :: restart_out_file_name*
(and is frequently *filter_restart*). The *.\#\#\#\#* is a 4 digit
number appended by filter based on the ensemble number. These files
contain the WRF state vector data that was used in the assimilation,
which is usually a subset of all the fields in a *wrfinput_d01* file.  
  
*dart_to_wrf* is the standard utility to insert the DART state
information into a WRF input file, e.g. *wrfinput_d01*. For multiple
WRF domains, a single run of the converter program will update the
*_d02*, *_d03*, ..., files at the same time as the *_d01* file.  
  
In the *input.nml* file, set the following:

~~~
&dart_to_wrf_nml
   model_advance_file = .false.
   dart_restart_name  = 'filter_restart.####',
/
~~~

where '\#\#\#\#' is the ensemble member number. There is no option to
alter the input/output WRF filename. Run *dart_to_wrf*. Remember to
preserve each *wrfinput_d01* file or you will simply keep overwriting
the information in the same output file. Repeat for each ensemble member
and you will be ready to run WRF to make ensemble forecasts.  
  
If filter is advancing the WRF model, and you want to spawn forecasts
from intermediate assimilation steps:  
Use the *assim_model_state_ic.\#\#\#\#* files instead of the
*filter_restart.\#\#\#\#* files, and set the *model_advance_file*
namelist item to be *.true.* .  

### The CESM Climate Model and DART

<span id="CESMDART"></span>

#### How does DART interact with CESM?

The CESM climate model comes with its own configuration, build, run, and
archive scripts already. The DART distribution supplies a 'setup' script
that calls the CESM scripts to build a new case, and then add a section
to the CESM run script so the DART code will be run after each CESM
model advance. The DART setup scripts are needed only when building a
new case. At run-time the CESM run scripts are used to start the job.
The CESM "multi-instance" capability is used to run multiple ensemble
members as part of a single job.

#### I want to assimilate with one of the CESM models. Where do I start?

We use the CESM framework to execute the CESM model components, and then
call the DART assimilation via an addition to the standard CESM run
script. We provide a set of setup scripts in our
`DART/models/XXX/shell_scripts` directories, where `XXX` is currently
one of: `cam`, `POP`, `clm`, or `CESM`. Start with the shell script, set
the options you want there, and then run the script. It calls the
standard CESM 'build_case' scripts, and stages the files that will be
needed for assimilation. See comments in the appropriate setup script
for more details of how to
proceed.

#### I'm getting a mysterious run-time error from CESM about box rearranging.

Certain versions of CESM (including CESM1_5_alpha02d) won't run with 3
instances (ensemble members). We are unsure what other instance sizes
fail. The error message is about box rearranging from
box_rearrange.F90. This is a problem in CESM and should be reported via
their Bugzilla process. 4 instances works
fine.

#### I'm getting 'update_reg_list' errors trying to assimilate with POP.

If you are trying to assimilate with POP and you get this error:

~~~
ERROR FROM:  
 routine: update_reg_list  
 message: max_reg_list_num ( 80) is too small ... increase  
~~~

The most likely cause is that the POP-DART model interface code is
trying to read the POP grid information and the default file is in the
wrong kind of binary for this system (big-endian and not little-endian).
At this point the easiest solution is to rebuild the DART executables
with a flag to swap the bytes as it is reading binary files. For the
Intel compiler, see the comments at the top of the mkmf file about
adding '-convert bigendian' to the FFLAGS line.

#### I'm getting asked to confirm removing files when CESM is built.

If you have the `rm` (remove) command aliased to require you to confirm
removing files, the CESM build process will stop and wait for you to
confirm removing the files. You should reply yes when prompted.  
  
If you have questions about the DART setup scripts and how they interact
with CESM it is a good idea to set up a standalone CESM case without any
DART scripts or commands to be sure you have a good CESM environment
before trying to add DART. The DART setup script uses CESM scripts and
commands and cannot change how those scripts behave in your environment.

#### I'm getting module errors when CESM is built.

DART only uses the plain NetCDF libraries for I/O. CESM can be
configured to use several versions of NetCDF including PIO, parallel
netCDF, and plain netCDF. Be sure you have the correct modules loaded
before you build CESM. If there are questions, try setting up a CESM
case without DART and resolve any build errors or warnings there before
using the DART scripts.

\[[top](#)\]

-----

<span id="widgets"></span>

### DART tools/preprocessors/widgets - under construction
