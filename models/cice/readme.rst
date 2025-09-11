.. index:: sea ice

CICE
====

Overview
--------

The `Community Ice CodE (CICE) <https://github.com/CICE-Consortium/CICE>`_ is a
sea ice model that was first developed by Elizabeth Hunke as the Los Alamos Sea
Ice Model. Its code base and capabilities have grown as a result of continued
development by the broader geosciences community, an effort organized by the
CICE Consortium.

Dr. Cecilia Bitz implemented support for the CICE (v5) model (as part of CESM)
in DART; this was later updated for compatibility with CICE v6. The DART model 
interface was developed to work with CICE's dynamical core on an Arakawa B-grid.
[1]_ When CICE is coupled to POP in CESM, the ocean and sea ice
grids are identical.

According to the CICE manual:

  The spatial discretization is specialized for a generalized orthogonal B-grid
  as in Murray (1996) [2]_ or Smith et al. (1995). [3]_ The ice and snow area,
  volume and energy are given at the center of the cell, velocity is defined at
  the corners, and the internal ice stress tensor takes four different values
  within a grid cell; bilinear approximations are used for the stress tensor
  and the ice velocity across the cell, as described in Hunke and Dukowicz
  (2002). [4]_ This tends to avoid the grid decoupling problems associated with
  the B-grid.

Hence, in the DART interface:

- ``U``, ``V`` are at grid cell corners
- ``T``, ``h``, ``hs``, and the various scalar quantities are at grid cell centers

CICE is under development to work with other grids, such as the
unstructured grid in MPAS and the C-grid in MOM.

Namelist
--------

.. code-block:: fortran

  &model_nml
     assimilation_period_days     = 1
     assimilation_period_seconds  = 0
     model_perturbation_amplitude = 0.00002
     update_dry_cell_walls        = .false.
     binary_grid_file_format      = 'big_endian'
     debug                        = 1
     model_state_variables        = 'aicen', 'QTY_SEAICE_CONCENTR',   'UPDATE',
                                    'vicen', 'QTY_SEAICE_VOLUME',     'UPDATE',
                                    ...
                                    'vsnon', 'QTY_SEAICE_SNOWVOLUME', 'UPDATE',
  /

Description of each namelist entry
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

+------------------------------+---------------+---------------------------------+
| Item                         | Type          | Description                     |
+==============================+===============+=================================+
| time_step_days               | integer       | Number of days for dimensional  |
|                              |               | timestep, mapped to deltat.     |
+------------------------------+---------------+---------------------------------+
| time_step_seconds            | integer       | Number of seconds for           |
|                              |               | dimensional timestep, mapped to |
|                              |               | deltat.                         |
+------------------------------+---------------+---------------------------------+
| model_perturbation_amplitude | real(r8)      | Perturbation amplitude          |
+------------------------------+---------------+---------------------------------+
| update_dry_cell_walls        | logical       | Currently does nothing.         |
|                              |               | Additional code is needed to    |
|                              |               | detect the cells which are      |
|                              |               | wet but within 1 cell of        |
|                              |               | the bottom/sides/etc.           |
+------------------------------+---------------+---------------------------------+
| binary_grid_file_format      | character(64) | Byte sequence for the binary    |
|                              |               | grid. Valid values are native,  |
|                              |               | big_endian & little_endian.     |
+------------------------------+---------------+---------------------------------+
| debug                        | integer       | When set to 0, debug statements |
|                              |               | are not printed. Higher numbers |
|                              |               | mean more debug reporting.      |
+------------------------------+---------------+---------------------------------+
| model_state_variables        | character(*)  | List of model state variables   |
+------------------------------+---------------+---------------------------------+


Postprocessing
---------------

The CICE model interface includes a postprocessing program ``dart_to_cice`` that 
modifies the CICE state to be consistent with the state output from DART assimilation.

This postprocessing is required when the assimilation configuration cannot appropriately 
account for the bounds on sea ice variables in the model. Prior to the introduction of 
the QCEF framework, the use of Gaussian distributions in the EAKF, in particular, often 
results in updates to sea ice concentration (SIC) that exceeded 0 or 1; it was also 
possible to produce adjustments to sea ice volume or thickness that were negative. [6]_ 
While these errors do not hamper the DA filtering, they result in various CICE and 
Icepack model crashes when the adjusted restart files are used to initialize subsequent 
forecasts. 

The utilization of bounded distributions in the QCEF framework reduces erroneous updates 
in the observation space, removing a source of the errors postprocessing seeks to address. 
[7]_ [8]_ However, the need for postprocessing remains, since the model’s state variables are 
represented using an ice thickness distribution, [5]_ which divides sea ice area and volume 
in each grid cell into discrete ice thickness categories. Observable sea ice variables like 
SIC and SIT are aggregated from the individual category fractions of sea ice area and volume.
During the filtering process, observation increments are regressed onto categorized state 
variables, and the final updated aggregate variable is recalculated as the sum (or weighted
sum) of the categories. In some cases, it is still not possible to appropriately account for
bounds on each category of a state variable and the bounds on their sum. [8]_ This is the case
for categorized sea ice area and total SIC, since SIC itself cannot be less than 0 or greater
than 1, but each area category can only berestricted to the interval [0, 1] in the DA process.
As such, it is still currently possible to assimilate observations that result in SIC greater
than 1 and postprocessing steps are still required for CICE-DART. 

There are currently 3 CICE/Icepack postprocessing options available in DART. The first two
(``aicen_simple_squeeze`` or ``vice_simple_squeeze``) are mass-aware rescaling approaches 
that conserve sea ice area or volume post-DA, respectively. The third and recommended default 
option is ``cice_rebalancing``, which adapts a series of functions originally developed to 
internally “rebalance” the CICE ice thickness distribution under conditions of unusually 
rapid or unexpected ice change. These routines are housed in the ``ice_postprocessing`` 
module, but called and deployed in ``dart_to_cice``. 

The postprocessing options for ``dart_to_cice`` are controlled by the namelist
``dart_to_cice_nml``. 


.. code-block:: fortran

  &dart_to_cice_nml
     dart_to_cice_input_file    = 'dart_restart.nc'
     original_cice_restart_file = 'cice_restart.nc'
     postprocessed_output_file  = 'postprocessed_restart.nc'
     postprocess                = 'cice'
  /


+------------------------------+---------------+---------------------------------+
| Item                         | Type          | Description                     |
+==============================+===============+=================================+
| dart_to_cice_input_file      | character(*)  | Output from filter              |
+------------------------------+---------------+---------------------------------+
| original_cice_restart_file   | character(*)  | Original CICE restart file that |
|                              |               | was input to filter             |
+------------------------------+---------------+---------------------------------+
| postprocessed_output_file    | character(*)  | Postprocessed CICE restart file |
+------------------------------+---------------+---------------------------------+
| postprocess                  | character(*)  | Postprocessing method. Must be  |
|                              |               | one of the following:           |
|                              |               | 'cice' (default),               |
|                              |               | 'aice',                         |
|                              |               | 'vice'                          |
+------------------------------+---------------+---------------------------------+

References
~~~~~~~~~~

.. [1] Arakawa, Akio and Vivian R. Lamb, 1977: Computational Design of the
       Basic Dynamical Processes of the UCLA General Circulation Model.
       *Methods in Computational Physics: Advances in Research and
       Applications*, **17**, 173–265, `doi:10.1016/B978-0-12-460817-7.50009-4
       <https://doi.org/10.1016/B978-0-12-460817-7.50009-4>`__

.. [2] Murray, Ross J., 1996: Explicit Generation of Orthogonal Grids for Ocean
       Models. *Journal of Computational Physics*, **126**, 251–273, 
       `doi:10.1006/jcph.1996.0136 <https://doi.org/10.1006/jcph.1996.0136>`__

.. [3] Smith, Richard D., Samuel Kortas and Bertrand Meltz, 1995: Curvilinear
       Coordinates for Global Ocean Models. Technical Report LA-UR95-1146, Los
       Alamos National Laboratory.

.. [4] Hunke, Elizabeth C., and John K. Dukowicz, 2002: The
       Elastic–Viscous–Plastic Sea Ice Dynamics Model in General Orthogonal
       Curvilinear Coordinates on a Sphere—Incorporation of Metric Terms.
       *Monthly Weather Review*, **130**, 1848–1865, 
       `doi:10.1175/1520-0493(2002)130%3C1848:TEVPSI%3E2.0.CO;2
       <https://doi.org/10.1175/1520-0493(2002)130%3C1848:TEVPSI%3E2.0.CO;2>`__

.. [5] Thorndike, A. S., D. A. Rothrock, G. A. Maykut, and R. Colony, 1975: 
       The thickness distribution of sea ice. *Journal of Geophysical Research*, 
       **80(33)**, 4501–4513, `doi:10.1029/JC080i033p04501 
       <https://doi.org/10.1029/JC080i033p04501>`__

.. [6] Zhang, Y., C. M. Bitz, J. L. Anderson, N. Collins, J. Hendricks, T. Hoar, 
       K. Raeder, and F. Massonnet, 2018: Insights on Sea Ice Data Assimilation 
       from Perfect Model Observing System Simulation Experiments. 
       *Journal of Climate*, **31**, 5911–5926, `doi:10.1175/JCLI-D-17-0904.1 
       <https://doi.org/10.1175/JCLI-D-17-0904.1>`__

.. [7] Riedel, C. P., M. M. Wieringa, and J. L. Anderson, 2025: Exploring Bounded 
       Nonparametric Ensemble Filter Impacts on Sea Ice Data Assimilation. 
       *Monthly Weather Review*, **153**, 637–654, `doi:10.1175/MWR-D-24-0096.1
       <https://doi.org/10.1175/MWR-D-24-0096.1>`__

.. [8] Wieringa, M. M., Riedel, C., Anderson, J. L., and Bitz, C. M., 2024: Bounded 
       and categorized: targeting data assimilation for sea ice fractional coverage 
       and nonnegative quantities in a single-column multi-category sea ice model.
       *The Cryosphere*, **18**, 5365–5382, `doi:10.5194/tc-18-5365-2024
       <https://doi.org/10.5194/tc-18-5365-2024>`__
