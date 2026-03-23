Aether Cube Sphere 
==================

This document describes the DART interface to the Aether ionosphere-thermosphere model in its cube
sphere implementation.

In addition to the standard DART programs associated with a given model, this interface creates
two additional programs for the Aether cube sphere: ``aether_to_dart``, and ``dart_to_aether``.
It also has a program in the ``developer_tests directory``, ``test_aether_grid`` 
that tests the geometry related aspects of the grid,
verifies model interpolation assumptions, and compares to a template filter_input_file. 

aether_to_dart
--------------

This program reads Aether restart files and combines them into a single filter input file for
DART. This program reads entries from two namelists in ``input.nml``:

In ``&aether_to_dart_nml``:

- ``aether_file_directory`` Path to the Aether restart files
- ``dart_file_directory``  Path to the where ``filter_input`` files will be created

In ``&transform_state_nml``:

- ``np`` Number of grid points across each cube sphere face (no halos)
- ``nblocks`` Total number of Aether grid files 
- ``nhalos`` Number of Aether halo rows
- ``scalar_f10_7`` True if F10.7 is a scalar in the DART state vector, false means it is has a value at each column.

Using aether_to_dart
~~~~~~~~~~~~~~~~~~~~

When executing ``aether_to_dart`` a range of ensemble members must be specified as 
command line arguments:

.. code-block::

    aether_to_dart 1 10

This example would translate ensemble members 1 to 10. There are three types of Aether netcdf
restart files required in the Aether input directory. The first are grid files that include the 
metadata that defines the location of Aether grid points and halos for each block. These file
are named ``grid_g####.nc`` where #### is the number of the block, starting from 0000. The second
are neutrals files that contain the data for neutral quantities. There is one file for each
block for each ensemble member. These are named ``neutrals_m&&&&_g####`` where &&&& is the ensemble
member starting from 0000 and #### is the block. The third are ions files named 
``ions_m&&&&_g####`` with the same numbering as the neutrals files. A single filter input file is 
created for each ensemble with the name ``filter_input_&&&&`` where &&&& is the ensemble number
starting at 0001 (note the offset in numbering from the Aether files). 

dart_to_aether
--------------

This program reads filter output files and inserts the data into Aether block restart files.
This program reads entries from two namelists in ``input.nml``"

In ``&dart_to_aether_nml``:

- ``dart_file_directory``  Path to the ``filter_output`` files that are read
- ``aether_file_directory`` Path to the Aether restart files that will be modified

In ``&transform_state_nml``:

- ``np`` Number of grid points across each cube sphere face (no halos)
- ``nblocks`` Total number of Aether grid files 
- ``nhalos`` Number of Aether halo rows
- ``scalar_f10_7`` True if F10.7 is a scalar in the DART state vector, false means it is has a

Using dart_to_aether
~~~~~~~~~~~~~~~~~~~~

When executing ``dart_to_aether`` a range of ensemble members must be specified as 
command line arguments:

    ``./dart_to_aether 1 10``

This example would translate ensemble members 1 to 10. A single filter input file is 
read for each ensemble with the name ``filter_output_&&&&`` where &&&& is the ensemble number
starting at 0001. The results are written into previously created aether neutrals and ions files
that are named as described above for ``aether_to_dart``. The Aether directory must also contain
grid files that contain the metadata for the Aether restarts. It is possible to overwrite the
Aether files in the same directory that was used for input to ``aether_to_dart``, or to copy the files
from that directory to a new directory to be updated by ``dart_to_aether``.

Using perfect_model_obs or filter with the aether model_mod
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Simple example workflow
~~~~~~~~~~~~~~~~~~~~~~~
This section describes a simple workflow that demonstrates and tests the capabilities of
the DART/Aether system.  To do this simple example, you will need to contact 
Jeff Anderson, jla@ucar.edu, or dart@ucar.edu to request a copy of the file ``TEST_INPUT.tar.gz``. 
Place this file in the ``models/aether_cube_sphere/`` directory and execute the command
``tar -xzvf TEST_INPUT.tar.gz`` to create the ``TEST_INPUT`` directory.
that contains a small set of Aether grid input files. 
There are 6 blocks, one covering each face of the cubed sphere. 
The neutrals files contain two variables, O2 and Temperature. The ions files contain two variables, 
O2+ and N2+. The ``TEST_INPUT`` directory only contains a single ensemble member for the ions and 
neutrals files. 

Demonstration steps:

1. In the models/aether_cube_sphere directory, execute the matlab script 
``perturb_aether_ensemble.m``. This generates 10 ensemble members in the ``TEST_INPUT``
directory. All variables are perturbed
so that the prior correlations between any two variables are one. 

2. In the ``aether_cube_sphere`` direcctory, copy the ``TEST_INPUT`` directory to the
``TEST_OUTPUT`` directory, ``cp -rf TEST_INPUT TEST_OUTPUT``.

3. Build all of the DART programs by executing ``quickbuild.sh nompi`` in the 
``aether_cube_sphere/work`` directory. 

4. Run ``aether_to_dart 1 10`` in the ``aether_cube_sphere/work`` directory.

5. Run ``perfect_model_obs`` in the ``aether_cube_sphere/work`` directory. This creates 
synthetic observations using the file ``obs_seq.in`` for metadata and creating the file
``obs_seq.out``.

6. Run ``filter`` in the ``work`` directory to do a single step ensemble assimilation.

7. Run ``dart_to_aether 1 10`` in the ``work`` directory to create updated aether restart
files in the ``aether_cube_sphere/TEST_OUTPUT`` directory.

8. Use matlab script ``plot_filter_lat_lon.m`` in the ``aether_cube_sphere`` directory to
interactively view the increments for different variables for the DART
``filter_input_&&&&.nc`` and ``filter_output_&&&&.nc`` files

9. Use matlab script ``plot_aether_lat_lon.m`` in the ``aether_cube_sphere`` directory to
view increments between the input aether restart files in the ``TEST_INPUT`` directory
and the updated aether restart filtes in the ``TEST_OUTPUT`` directory.

The ``obs_seq.in`` file defines the observations that are created by ``perfect_model_obs`` and
then assimilated by ``filter``. In this case, there are 6 observations, one each of 
temperature, density of O2, density of O2P, density of N2P, a ground station GPS
vertical total electron content, and a slant GPS total electron content. Each is at
a different horizontal location and the first 4 are at different vertical locations.

The file ``create_obs_seq.input`` in the ``aether_cube_sphere`` directory contains input that
can be read by the program ``create_obs_sequence`` to create the default ``obs_seq.in`` file

Work in Progress
~~~~~~~~~~~~~~~~

**Time:**
The method by which model time is read into DART has not been finalized at this time. All tests
to date use time that is manually inserted into the ``perfect_model_obs`` and ``filter`` namelist entries
``init_time_days`` and ``init_time_seconds``. The specifics of the how time is included in Aether input 
files needs to be clarified so that the model_mod can read this directly from the filter restart
files. Aether is not currently using time that is consistent with any calendar supported by DART,
so this may require code in ``aether_to_dart.f90`` that translates the aether time to a time that 
DART understands.

**F10.7:**
Aether restart netcdf files do not currently include parameter values like F10.7. For now, 
the ``aether_to_dart`` and ``dart_to_aether`` programs do not do not do input/output with Aether.
Default (meaningless) F10.7 values are generated by transform_state_mod.f90. Once Aether
modelers decide where the F10.7 value will come from, code will need to be added to read
the F10.7 value in ``model_to_dart`` and write it in ``dart_to_model``,
but obvious hooks are available in ``transform_state_mod.f90``. This module implements the
basics of two ways to do F10.7 estimation. The first is to have a single scalar value of 
F10.7 in the DART state. Subroutine ``get_state_meta_data`` provides some initial suggestions for
the location associated with a scalar F10.7 that are taken from Alexey Morozov's work in 
GITM. Because this requires the time, which is not yet available from Aether, this requires
additional implementation. Aether scientists also need to confirm that the subsolar point
is the right choice for a location. Alexey also implemented a different localization 
algorithm for F10.7 in GITM. Aether scientists should work with DART experts to determine
if and how this would be implemented in Aether. Under namelist control, ``aether_to_dart``
can also treat F10.7 as a horizontally distributed variable, basically copying the same value
of F10.7 to each horizontal column. The value at each column is updated and ``dart_to_aether``
currently just averages the posterior values. Other choices for weighted averages are
scientifically interesting and could be explored by aether/DART collaborations.


**VTEC:**
The established forward operator for vertically integrated electron content in DART is found in 
the ``observations/forward_operators/obs_def_upper_atm_mod.f90``. It assumes that the DART state 
includes a 3D field with quantity ``QTY_DENSITY_ION_E`` and that the state also includes the 
geometric height of each grid point in ``QTY_GEOMETRIC_HEIGHT``. The subroutine 
``get_expected_gnd_gps_vtec`` integrates the density in a column. This subroutine was originally 
developed for GITM and then extended for TIEGCM. Unlike GITM, Aether does not include the
ION density in its restart netcdf files. The ``aether_to_dart.f90`` sums up the density of all 
variables in the ions files that have units of /m3 and puts this into the filter_input file that
is created. Aether model experts should verify both the creation of the density field and the
way that a vertical integral is computed to confirm that these are consistent with the model
and the available observations. Note that there are other electron content forward operators
that may also need to be evaluated by model experts before use.

**Slant VTEC:**
There is a subroutine called ``get_expected_slant_gps_vtec`` in
``/observations/forward_operators/obs_def_upper_atm_mod.f90``. It does exaclty the same thing
as the vtec described above. However, it includes extended metadata in the obs_seq files. 
These are two locations descriptions, one for the satellite postion (lon, lat, height), 
and one for a ground point (lon, lat, height). One way to implement a slant vtec forward
operator would be to trace a ray between the satellite and the ground and get the density
at each level along the ray. Other ways of describing the geometry of the ray may be more
appropriate. Aether developers and observation experts should be able to use the example
code to easily implement the forward operator once the exact method for tracing the ray
from the satellite is implemented. 

Testing the grid computations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The program ``test_aether_grid`` in ``developer_tests/aether_grid`` can be run 
with the namelist setting used for a ``filter`` run to 
verify the geometry in the ``model_mod`` and to confirm consistency with the aether template file
selected by the ``template_file`` entry in the ``model_nml`` namelist. Note that an aether template
filter file must have been created in the ``aether_cube_sphere/work`` directory before this test
is run.

