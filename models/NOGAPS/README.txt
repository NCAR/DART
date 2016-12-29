# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$


TJH: 30 July 2013
	The DART interfaces for NOGAPS have been approved for distribution.
	The NOGAPS code itself is NOT available for public download.
	As such, we have tried to keep the interfaces current as far as
	DART requirements have evolved, but we have been unable to test
	them in a real environment. These routines are from 2010 or 2011.


Here is my take on the experiment setup for a perfect model assimilation 
experiment with DART and NOGAPS. I'm not sure of the version of NOGAPS 
(or I'd record that here). NOGAPS requires the LAPACK library. 
The system used to derive the NOGAPS/DART code used the Intel 10.1 
compiler suite with the Intel Math Kernel Libraries on an 
Intel x86-64 chipset (little-endian). The following compiler 
arugments were required in the mkmf.template: 

INCS = -I$(NETCDF)/include
LIBS = -L$(NETCDF)/lib -lnetcdf -lmkl -lmkl_lapack -lguide -lpthread
FFLAGS  = -g -O2 -fpe0 -convert big_endian -assume byterecl -132 $(INCS)
LDFLAGS = $(FFLAGS) $(LIBS)


NOGAPS, filter, wakeup_filter, dart_to_nogaps, and nogaps_to_dart are 
necessarily compiled with MPI support. perfect_model_obs and trans_time 
are single-threaded.

There must be a 'src' directory in the DART/models/NOGAPS directory
with a libnogaps.a, imp.h, param.h, and the  following modules (not all
of which may be needed):
basediab.mod
cubic.mod
cupsum.mod
diagnos.mod
dyncons.mod
fields.mod
gcvars.mod
mpinog.mod
navdas.mod
symem.mod
times.mod

These are expected to be in the DART/models/NOGAPS/src
directory by the mkmf_ scripts.


This is the BIG PICTURE:

0) [DART/models/]NOGAPS/shell_script/config.csh must be edited to reflect the 
   locations of the executables, etc.

1) A set of NOGAPS files are converted to DART initial conditions 
   files by the script NOGAPS/shell_scripts/run_nogapsIC_to_dart.csh.
   run_nogapsIC_to_dart.csh is a batch script that uses Job Array
   syntax to launch N "identical" jobs - one for each ensemble member
   to convert to a DART initial conditions file.  This script also
   creates the "experiment directory", the run-time filesystem for all
   future assimilation runs of this experiment. 

2) NOGAPS/shell_scripts/run_perfect.csh   will select ONE of these
   initial conditions files 
   (input.nml:&perfect_model_obs_nml:restart_in_file_name)
   as though it is the true model trajectory - and harvest 
   synthetic observations from this trajectory at the times/locations 
   specified in obs_seq.in The synthetic observations will be in obs_seq.out

3) NOGAPS/shell_scripts/run_filter.csh   will launch multiple MPI-aware
   programs that successively use the entire processor set. At this time,
   it is not possible to have some of the executables simultaneously use 
   subsets of the processor pool.

Now for the detail:

------------------------------------------------------------
Section 0: editing config.csh
------------------------------------------------------------

	This is pretty self-explanatory (I hope). 

The idea is that the NOGAPS/shell_scripts/config.csh will be copied
to the experiment directory to be used by the subsequent processes.
Fundamentally, the combination of $scratch_dir and $experiment_name
define the variable $experiment_dir (aka CENTRALDIR). This is 
where all the run-time action takes place.

------------------------------------------------------------
Section 1: run_nogapsIC_to_dart.csh 
Generating the DART initial conditions
------------------------------------------------------------

&nogaps_to_dart_nml
   nogaps_restart_filename    = 'specfiles/shist000000',
   nogaps_to_dart_output_file = 'dart_new_vector',
   nogaps_data_time_filename  = 'CRDATE.dat'
   /

	The run_nogapsIC_to_dart.csh script is written to exploit
the Job Array syntax of LSF and PBS and can easily be modified to
incorporate others. The idea is simple. Multiple copies of the job
are spawned when the job is submitted ONCE. Each copy of the job
has a unique Array ID or Task Identifier or ... I translate all the
queueing-system specific variables to generic ones and use the generic
ones throughout the rest of the script. This one script will work on 
multiple platforms.

	The number of copies of jobs spawned is controlled through
the job array syntax. We use one copy for each desired ensemble member.
This must be hand-entered as part of the LSF/PBS directives. Each has
their own separate syntax.

	Each ensemble member gets a unique directory and the files and
executables are copied to these directories. The weak spot here is that
ALL of the ensemble members are trying to link and populate the SAME
climo directory. Since this is readonly anyway, I'd prefer to just
point everything to the souce of the link and be done with it. 
(i.e. in config.csh just 'set climo = "${climo_dir}/climo${resolution}"')
Dan created tarfiles with the specfiles for multiple ensemble members. 
The appropriate tarfile (which remains in-situ and) is unpacked. 
The contents are then moved to the expected locations. i.e. those 
specified in input.nml:&nogaps_to_dart_nml:nogaps_restart_filename
('specfiles/shist000000').

	The namelists for NOGAPS are created. All directories use relative
pathnames (i.e. '.') for the shortest possible names. There were commments
in the original scripts that the length of the strings was a concern.
Since the namelists used to use 'temp' - which was a relative link - or
'$HOME' which was the same place as 'temp' and '.', there is no reason
to make things more complicated than necessary.

	FIXME: nogaps_to_dart() uses the following time processing logic.
The NOGAPS pieces need to know the 'current' time of the NOGAPS state - 
this is contained in CRDATE.dat. The DART code needs to know the time
of the NOGAPS state - currently from a different file ('nogaps_data.time')
containing exactly the same information. This is crazy..
this should be the same as CRDATE.dat  ...

	All the bits and pieces necessary to run nogaps_to_dart are
assembled in the unique run directory and nogaps_to_dart is run.
The output file is namelist specified as 'dart_new_vector' and is 
renamed to be consistent with what 'filter' will expect from the
input.nml:$filter_nml:restart_in_file_name (usually "filter_ic").

------------------------------------------------------------
Section 2: (optional) run_perfect.csh 
Running a perfect model experiment.
------------------------------------------------------------

	I am going to assume that the target observation sequence
file is already created somewhere and is called 'obs_seq.in'.
Furthermore, 'perfect_model_obs' is a single-threaded application
while the model is MPI-aware. This means that only one MPI-aware
application is running at one time - a pretty simple scenario.

	There are several comment blocks for PBS or LSF directives
that make it possible to use the same script for both batch queueing
systems. The first executable portion of the script simply translates
the queueing-system-specific variables to generic names that can be
used throughout the remainder of the script. The 'experiment_dir'
(i.e. CENTRALDIR) is known from the original config.csh, so the first
that happens is to 'cd' to CENTRALDIR.

	All the executables and input control files are copied
to CENTRALDIR. The ensemble member to be used as THE TRUTH is defined 
by input.nml:&perfect_model_obs_nml:restart_in_file_name
(right now it is set to = "filter_ic.0001") which must be a pre-existing
file in CENTRALDIR (created by run_nogapsIC_to_dart.csh).

	Really all that is left is to set the value of the MPI 
command that is needed by the model executable. If you are using 
a queueing system, the MPI command is already known 
(from config.csh, actually); if not, there is some work to be done.
The block with the comment 
"# WARNING: This block is untested ..." 
is, well, untested and unlikely to work without modification.

	This is a great way to test changes to the advance_model.csh
script. The same advance_model.csh script can be (should be?) used
by both perfect_model_obs and by filter.

------------------------------------------------------------
Section 3: run_filter.csh 
Running an assimilation experiment.
------------------------------------------------------------

	run_filter.csh has the same strategy as run_perfect.csh
as pertains to the submission directives and variable-name translation.
All the input files/executables are copied to CENTRALDIR.
There is some shell trickery to extract bits from the input.nml - namely;
the ensemble size, the filter 'async' variable, and the string
containing the model advance command.  All of these have bearing
on the logic of the script.

	Essentially, if the model is an MPI-aware program and filter
is an MPI-aware program ... getting the O/S to run both of these at
the same time has been tricky. Filter runs in the background and
quite literally goes to sleep while the model executes. When the model
advance is complete, 'wakeup_filter' is executed to wake filter and
continue. The communication for this is through named pipes - which 
are like files but don't have the delay of the filesystem. The one
problem with this is that sometimes filter fails and doesn't exit
cleanly - causing the job to hang. system-dependent, but we're working
on a more reliable mechanism.

	Since NOGAPS and the convert routines ARE MPI-aware,
the $parallel_model variable must be TRUE. The logic for the
parallel_model = .false. is included for completeness only.

	In the current implementation, the 'filter_control00000' file
created by filter contains all the information to advance all of
the ensemble members - one after another - sequentially.

------------------------------------------------------------
Explanation of advance_model.csh - how one member advances
------------------------------------------------------------

	advance_model.csh gets spawned by either perfect_model_obs
or filter. The input.nml:&filter_nml:adv_ens_command
(or &perfect_model_obs_nml:adv_ens_command) specify which
shell script gets invoked. For the scenario when all of the MPI
tasks available to the job are used to advance a single ensemble member,
'advance_model.csh' is the right choice. 'advance_model_batch.csh' is
under development for the scenario when the model advance is a 
self-contained job for the queueing system.  The remainder of the
discussion is for 'advance_model.csh' and is intended to be the
'cleanest' (i.e. simplest?) example upon which to build future scripts.  

	advance_model.csh is called with precisely three arguments.
The process number of the caller, the number of ensemble members that
must be advanced by the process, and the name of the control file for
the process. These are used to control the iterations (one for each
ensemble member) inside advance_model.csh.

	A fundamental tenet is that all the files needed by
advance_model.csh are available in the run-time directory of 'filter'
(i.e. CENTRALDIR).  'filter' will fork the advance_model.csh script - 
which will cd to a local directory in which to advance the model. 
All the files are moved/copied into this local directory, the work 
is done, and the output file get moved/copied back to CENTRALDIR. 

1) Some error-checking is performed to ensure the directories
   required by NOGAPS exist. I have no idea what is supposed to
   be in those directories, so ...

   The (private, local) sub-directory is created and populated with
   generic bits from CENTRALDIR. The DART state vector file is
   queried (by trans_time) to extract the current/valid date of 
   the state vector, the target or "advance_to" date, and the 
   forecast length - in hours. trans_time is based on a little DART
   utility and is customized for NOGAPS - so it is the assembla
   repository as opposed to the general DART repository. trans_time
   expects the input filename to be 'temp_ic' and the output
   file containing the time information to be 'time_info'. These
   are hardwired. 'time_info' has three things - one per line:
   <dtg>
   <dtgnext>
   <endtau - the forecast length>

   The NOGAPS namelists are created. There had been a circular
   dependence on some environment variables to specify both
   absolute and relative pathname information to the same directory;
   all pathnames are now relative to the current private run-time
   directory by using the "./" methodology.

2) The dart state is converted to a nogaps specfile. nogaps_to_dart
   is MPI-aware and requires some of the nogaps code - which requires 
   a file named 'CRDATE.dat' containing a valid "$dtg". nogaps_to_dart 
   also creates a file that is not needed - the default name is
   'dart_data.time' - This is the same format as from the 'trans_time'
   routine - containing the same data - so it's redundant.

3) The model is advanced using 'CRDATE.dat' (i.e. $dtg) for a forecast
   length of $endtau (via namelist) to the time specified by $dtgnext 
   FIXME: CRDATE.dat is manually updated with the new time/date. 
   If possible, this is the number one thing I would fix. If the 
   model advance fails but does not cause an error exit, the whole 
   machine can march blindly on ...
   
4) Convert the nogaps state to a dart state - which requires
   some of the nogaps code - which requires a CRDATE.dat file.
   The updated DART state is moved back to CENTRALDIR with the
   required filename.

5) The indices to extract information from the 'filter_control00000' file
   are updated and the current working directory is moved back to 
   CENTRALDIR. 

	After all ensemble members have been advanced, the 
'filter_control00000' file is removed. IMPORTANT - if this file still 
exists when the advance_model.csh script has finished - IT IS AN ERROR - 
and filter will die a very theatrical death.

------------------------------------------------------------
------------------------------------------------------------

	All of this is predicated on the ability to assimilate 
as many cycles as you want in a single job - which is unrealistic
and not very smart. When DART finishes, it writes out restart files
that can be used as input for subsequent assimilations. Coming
up with a naming scheme to archive these files is left as an
exercise ... as are the scripts that manipulate the observation 
sequences and/or times that appear in the input.nml:&filter_nml namelist. 


