---
title: Getting Started
layout: default
---

## Getting Started with DART

The entire process is summarized in the following steps:

1.  [Determine which F90 compiler is available](#fortran90).
2.  [Determine the location of (or build) the *netCDF* library](#netCDFlib).
3.  [Download the DART software into the expected source tree](#installing).
4.  [Modify certain DART files to reflect the available F90 compiler and location of the appropriate libraries](#customizations).
5.  [Build the executables](#building).
6.  [Verifying the result](#matlab).

If you can compile and run ONE of the low-order models, you should be
able to compile and run ANY of the low-order models. For this reason, we
can focus on the Lorenz_63 model. Consequently, the only directories
with files to be modified to check the installation are usually:
`DART/build_templates` and `DART/models/lorenz_63/work`.  

We have tried to make the code as portable as possible, but we do not
have access to all compilers on all platforms, so there are no
guarantees. We are interested in your experience building the system, so
please send us a note at dart @ ucar .edu

<span id="requirements" class="anchor"></span> [](#requirements)   

----

## System Requirements

DART is intended to be highly portable among Unix/Linux operating
systems. At this point we have no plans to port DART to Windows
machines.  

Minimally, you will need:

1.  [a Fortran90 compiler](#fortran90),
2.  the [netCDF libraries](http://www.unidata.ucar.edu/software/netcdf/)
    built with the F90 interface,
3.  *perl* (just about any version),
4.  an environment that understands *csh* or *tcsh*, *sh*, *ksh*, and
5.  the old unix standby ... *make*
6.  and more than 1Gb of disk space for the DART distribution.

History has shown that it is a very good idea to make sure your run-time
environment has the following:  

~~~
limit stacksize unlimited  
limit datasize unlimited
~~~

Additionally, what has proven to be *nice* (but not required) is:

1.  [ncview](http://meteora.ucsd.edu/~pierce/ncview_home_page.html): a
    great visual browser for netCDF files.
2.  [the netCDF Operators (NCO)](http://nco.sourceforge.net/): tools to
    perform operations on netCDF files like concatenating, slicing, and
    dicing
3.  Some sort of MPI environment. Put another way, DART does not come
    with MPICH, LAM-MPI, or OpenMPI; but we use them all the time. You
    should become familiar with [the DART MPI introduction](dart_mpi.md).
4.  If you want to use the DART diagnostic scripts, you will need a
    basic MATLAB®. No additional toolboxes are required. No third-party
    toolboxes are required.

<span id="fortran90" class="anchor"></span> [](#fortran90)  

### Requirements: a Fortran90 compiler

The DART software is written in standard Fortran 90, with no compiler-specific
extensions. It has been compiled with and run with several versions of each of
the following:
[GNU Fortran Compiler ("gfortran")](http://gcc.gnu.org/fortran) (free),
[Intel Fortran Compiler for Linux and OSX](http://software.intel.com/en-us/intel-composer-xe),
[IBM XL Fortran Compiler](http://www-01.ibm.com/software/awdtools/fortran/),
[Portland Group Fortran Compiler](http://www.pgroup.com/), and
[Lahey Fortran Compiler](http://www.lahey.com/).
Since recompiling the code is a necessity to experiment with different models,
there are no binaries to distribute.

<span id="netCDFlib" class="anchor"></span> [](#netCDFlib)  

### Requirements: the netCDF library

DART uses the [netCDF](https://www.unidata.ucar.edu/software/netcdf/)
self-describing data format for the results of assimilation experiments.
These files have the extension *.nc* and can be read by a number of
standard data analysis tools. In particular, DART also makes use of the
F90 interface to the library which is available through the `netcdf.mod`
and `typesizes.mod` modules. *IMPORTANT*: different compilers create
these modules with different "case" filenames, and sometimes they are
not **both** installed into the expected directory. It is required that
both modules be present. The normal place would be in the
`netcdf/include` directory, as opposed to the `netcdf/lib` directory.  

If the netCDF library does not exist on your system, you must build it
(as well as the F90 interface modules). The library and instructions for
building the library or installing from an RPM may be found at the
netCDF home page: <https://www.unidata.ucar.edu/software/netcdf/>  

NOTE: The location of the netCDF library, `libnetcdf.a`, and the
locations of both `netcdf.mod` and `typesizes.mod` will be needed later.
Depending on the version of netCDF and the build options selected, the
fortran interface routines may be in a separate library named
`libnetcdff.a` (note the 2 "f"'s). In this case both libraries are
required to build executables.

<span id="download" class="anchor"></span> [](#download)  

\[[top](#)\]

-----

## Download DART

The DART source code is distributed on the GitHub repository
[NCAR/DART](https://github.com/NCAR/DART)
with the documentation served through GitHub Pages at
[https://ncar.github.io/DART](https://ncar.github.io/DART).

If you ever intend to contribute your work back to DART, we ask that
you fork the repository to facilitate issuing pull requests. *:thumbsup:*!

If you follow the instructions on the download site, you should wind up
with a directory named **my_path_to/DART**, which we call
_**DARTHOME**_. Compiling the code in this tree (as is usually the case)
will necessitate much more space.  

<span id="installing" class="anchor"></span> [](#installing)  

\[[top](#)\]

-----

## Testing DART

### document conventions

> commands to be typed at the command line appear as a blockquote

And the contents of a file are enclosed in a box:

~~~
&hypothetical_nml  
  obs_seq_in_file_name = "obs_seq.in",  
  obs_seq_out_file_name = "obs_seq.out",  
  init_time_days = 0,  
  init_time_seconds = 0,  
  output_interval = 1  
&end
~~~

<span id="customizations" class="anchor"></span> [](#customizations)  

### Customizing the build scripts -- Overview.

DART executable programs are constructed using two tools: *mkmf*, and
*make*. The *make* utility is a very common piece of software that
requires a user-defined input file that records dependencies between
different source files. *make* then performs a hierarchy of actions when
one or more of the source files is modified. *mkmf* is a *perl* script
that generates a *make* input file (named *Makefile*) and an example
namelist `input.nml._**program**_default` with the default values. The
*Makefile* is designed specifically to work with object-oriented
Fortran90 (and other languages) for systems like DART.  

*mkmf* (think *"make makefile"*) requires two separate input files. The
first is a 'template' file which specifies details of the commands
required for a specific Fortran90 compiler and may also contain pointers
to directories containing pre-compiled utilities required by the DART
system. **This template file will need to be modified to reflect your system**.

The second input file is a 'path_names' file which are
supplied by DART and can be used without modification. An *mkmf* command
is executed which uses the 'path_names' file and the mkmf template file
to produce a `Makefile` which is subsequently used by the standard
*make* utility.  

Shell scripts that execute the mkmf command for all standard DART
executables are provided as part of the standard DART software. For more
information on the FMS [mkmf](https://github.com/NOAA-GFDL/mkmf) please
see the [mkmf documentation](https://extranet.gfdl.noaa.gov/~vb/mkmf.html).  
Be aware that we have slightly modified *mkmf* such that it also creates
an example namelist file for each program. The example namelist is
called *input.nml._program_default*, so as not to clash with any
existing `input.nml` that may exist in that directory.

<span id="template" class="anchor"></span> [](#template)  

#### Building and Customizing the 'mkmf.template' file

A series of templates for different compilers/architectures exists in
the *DART/build_templates* directory and have names with extensions
that identify the compiler, the architecture, or both. This is how you
inform the build process of the specifics of your system. **Our intent
is that you copy one that is similar to your system into
`DART/build_templates/mkmf.template` and customize it.** For the
discussion that follows, knowledge of the contents of one of these
templates (i.e. `DART/build_templates/mkmf.template.intel.linux`) is
needed. Note that only the LAST lines are shown here, the head of the
file is just a big comment (worth reading).

~~~
...  
MPIFC = mpif90  
MPILD = mpif90  
FC = *ifort*  
LD = *ifort*  
NETCDF = */usr/local*  
INCS = -I${NETCDF}/include  
LIBS = -L${NETCDF}/lib -lnetcdf  
FFLAGS = -O2 $(INCS)  
LDFLAGS = $(FFLAGS) $(LIBS)  
~~~

| variable | value |
| :------- | :---- |
| FC       | the Fortran compiler |
| LD       | the name of the loader; typically, the same as the Fortran compiler |
| NETCDF   | the location of your netCDF installation containing `netcdf.mod` and `typesizes.mod`. Note that the value of the *NETCDF* variable will be used by the *FFLAGS, LIBS,* and *LDFLAGS* variables. |

<span id="path_names" class="anchor"></span> [](#path_names)  

#### Customizing the 'path_names_*' file

Several *path_names_&ast;*  files are provided in the *work* directory for
each specific model, in this case: `DART/models/lorenz_63/work`. Since
each model comes with its own set of files, the *path_names_&ast;* files
need no customization.

<span id="building" class="anchor"></span> [](#building)  

### Building the Lorenz_63 DART project.

All DART programs are compiled the same way. Each model directory has a
directory called `work` that has the components to build the
executables. This is an example of how to build two programs for the
lorenz_63 model: *preprocess* and *obs_diag*. *preprocess* needs to be
built and run to create the source codes that support observations -
which are (potentially) different for every model. Once that has been
done, any other DART program may be built in the same was as *obs_diag*
in this example.

> cd DART/models/lorenz_63/work  
> ./mkmf_preprocess  
> make  
> ./preprocess  
> ./mkmf_obs_diag  
> make  

Currently, DART executables are built in a `work` subdirectory under the
directory containing code for the given model. The Lorenz_63 model has
seven *mkmf_xxxxxx* files (some models have many more) for the
following programs:  

| Program | Purpose |
| :------ | :------ |
| [preprocess](https://ncar.github.io/DART/api/v2.1.10/program/preprocess.html) | creates custom source code for just the observations of interest |
| [create_obs_sequence](https://ncar.github.io/DART/api/v2.1.10/program/create_obs_sequence.html) | specify a (set) of observation characteristics taken by a particular (set of) instruments |
| [create_fixed_network_seq](https://ncar.github.io/DART/api/v2.1.10/program/create_fixed_network_seq.html) | specify the temporal attributes of the observation sets |
| [perfect_model_obs](https://ncar.github.io/DART/api/v2.1.10/program/perfect_model_obs.html) | spinup, generate "true state" for synthetic observation experiments, ... |
| [filter](https://ncar.github.io/DART/api/v2.1.10/program/filter.html) | perform experiments |
| *obs_diag* | creates observation-space diagnostic files to be explored by the MATLAB® scripts. |
| [obs_sequence_tool](https://ncar.github.io/DART/api/v2.1.10/program/obs_sequence_tool.html) | manipulates observation sequence files. It is not generally needed (particularly for low-order models) but can be used to combine observation sequences or convert from ASCII to binary or vice-versa. Since this is a specialty routine - we will not cover its use in this document. |

*quickbuild.csh* is a script that will build every executable in the
directory. There is an optional argument that will additionally build
the mpi-enabled versions - which is not the intent of this set of
instructions. Running *quickbuild.csh* will compile all the executables.

> cd DART/models/lorenz_63/work  
> ./quickbuild.csh -nompi

The result (hopefully) is that seven executables now reside in your work directory.
*The most common problem* is that the netCDF libraries and include files
(particularly `typesizes.mod`) are not found.
Find them, edit the `DART/build_templates/mkmf.template` to point to their
location, recreate the `Makefile`, and try again.

### Checking the build -- running something.

This section is not intended to provide any details of why we are doing
what we are doing - this is sort of a 'black-box' test.
The `DART/models/lorenz_63/work` directory is distributed with input files
ready to run a simple experiment: use 20 ensemble members to assimilate
observations 'every 6 hours' for 50 days.
Simply run the programs *perfect_model_obs* and *filter* to generate
the results to compare against known results.  

The initial conditions files and observations sequences are in ASCII, so
there is no portability issue, but there may be some roundoff error in
the conversion from ASCII to machine binary. With such a highly
nonlinear model, small differences in the initial conditions will result
in a different model trajectory. Your results should start out looking
VERY SIMILAR and may diverge with time.  

The Manhattan release uses netCDF files for the input file format.
Creating the netCDF files from their ASCII representation is a trivial
operation - simply running a command that comes with any netCDF
installation: *ncgen*. This is done automatically by the *lorenz_63 quickbuild.csh*,
but is repeated here for clarity. Once the netCDF input files
are created, simply running *perfect_model_obs* and *filter* is
easy:

> ncgen -o perfect_input.nc perfect_input.cdl  
> ncgen -o filter_input.nc filter_input.cdl  
> ./perfect_model_obs  
> ./filter

There should now be the following output files:

|                     |                   |
| :------             | :------           |
| **from executable "perfect_model_obs"** |       |
| `perfect_output.nc` | a netCDF file containing the model trajectory ... the **truth** |
| `obs_seq.out`       | The observations (harvested as the true model was advanced) that were assimilated. |
| **from executable "filter"** |      |
| `preassim.nc`       | A netCDF file of the ensemble model states just before assimilation. This is the **prior**. |
| `filter_output.nc`  | A netCDF file of the ensemble model states after assimilation. |
| `obs_seq.final`     | The observations and ensemble estimates of the 'observations'. |
| **from both**       |      |
| `dart_log.out`      | The run-time log of the experiment.  This grows with each execution and may safely be deleted at any time. |
| `dart_log.nml`      | A record of the input settings of the experiment.  This file may safely be deleted at any time. |

Note that if you change the `input.nml` namelist values controlling
inflation and file output, several (perhaps many) more files are created.  

The [DART/docs/tutorial](Tutorial.md)
documents are an excellent way to kick the tires on DART and learn about
ensemble data assimilation. If you've been able to build the Lorenz 63
model, you have correctly configured your `mkmf.template` and you can
run anything in the tutorial.

<span id="matlab" class="anchor"></span> [](#matlab)

\[[top](#)\]

-----

## Configuring MATLAB®

The Manhattan release of DART uses native MATLAB netCDF support and no
longer requires any third-party toolboxes. Furthermore, no additional
MATLAB toolboxes are required. To allow your environment to seamlessly
use the DART MATLAB functions, your MATLABPATH must be set such that you
have access to a couple of DART directories. Do something like the
following at the MATLAB prompt, using the real path to your DART
installation:

> \>\> addpath('path_to_dart/diagnostics/matlab','-BEGIN')  
> \>\> addpath('path_to_dart/docs/DART_LAB/matlab','-BEGIN')

It's very convenient to put these it in your `startup.m` so
they get run every time MATLAB starts up. You will have to copy
`startup.m` to whatever directory you have specified as your Matlab `userpath`.

If you don't want the DART environment every time you use Matlab,
DART provides an example `diagnostics/matlab/startup.m` that checks
to see if you are working in a DART repository and will modify your
matlabpath accordingly. It is internally documented. Again, you will have to copy
`diagnostics/matlab/startup.m` to whatever directory you have specified
as your Matlab `userpath`.

<span id="verify" class="anchor"></span> [](#verify)  

\[[top](#)\]

-----

## Verify: Are the results correct? (requires MATLAB®)

The initial conditions files and observations sequences are in ASCII, so there
is no portability issue, but there may be some roundoff error in the conversion
from ASCII to machine binary. With such a highly nonlinear model, small
differences in the initial conditions will result in a different model
trajectory. Even compiler options will cause tiny differences that ultimately
result in noticeably different trajectories.
Your results should start out looking VERY SIMILAR and may diverge with time.  

The simplest way to determine if the installation is successful is to run some
of the functions we have available in `DART/diagnostics/matlab/`. Usually, we
launch MATLAB from the `DART/models/lorenz_63/work` directory and use the
MATLAB *addpath* command to make the `DART/matlab/` functions available. In
this case, we know the true state of the model that is consistent with the
observations. The following MATLAB scripts compare the ensemble members with
the truth and can calculate an error.  
  
<table>
<tr>
<td width="50%">
<pre>
<code>
[unix prompt] cd DART/models/lorenz_63/work
[unix prompt] matlab -nodesktop
(lots of startup messages I'm skipping)

    >> addpath ../../../diagnostics/matlab
    >> plot_total_err
    Input name of true model trajectory file;
    <cr> for perfect_output.nc
    perfect_output.nc
    Input name of ensemble trajectory file;
    <cr> for preassim.nc
    preassim.nc
    Comparing true_state.nc and
              preassim.nc
    >> plot_ens_time_series
    Input name of ensemble trajectory file;
    <cr> for preassim.nc

    Comparing true_state.nc and  
              preassim.nc  
    Using Variable state IDs 1  2  3  

    pinfo =  

      struct with fields:  

                     model: 'Lorenz_63'  
                   def_var: 'state'  
            num_state_vars: 1  
                num_copies: 20  
           num_ens_members: 20  
          ensemble_indices: [1 2 3 ... 18 19 20]  
             min_state_var: 1  
             max_state_var: 3  
            def_state_vars: [1 2 3]  
                     fname: 'preassim.nc'  
                truth_file: 'true_state.nc'  
                diagn_file: 'preassim.nc'  
                truth_time: [1 200]  
                diagn_time: [1 200]   
                      vars: {'state'}   
                      time: [200x1 double]  
        time_series_length: 200   
                       var: 'state'  
                  var_inds: [1 2 3]
</code>
</pre>
</td>
<td width="50%">
<img src="../images/lorenz_63_total_err.png"       width="500" alt="xxxx" /><br />  
<img src="../images/lorenz_63_ens_time_series.png" width="500" alt="xxxx" />
</td>
</tr>
</table>  
   
From the *plot_ens_time_series* graphic, you can see the individual
green ensemble members getting more constrained as time evolves. If your
figures look similar to these, that's pretty much what you're looking
for and you should feel pretty confident that everything is
working. Don't miss the opportunity to rotate the 'butterfly' plot.  
  
\[[top](#)\]
  
-----
   
### What to do if things **do not** look correct
  
FIXME ... tell them what to do here

