CICE
====

Overview
--------

The `Community Ice CodE (CICE) <https://github.com/CICE-Consortium/CICE>`_ is a
sea ice model that was first developed by Elizabeth Hunke as the Los Alamos Sea
Ice Model. Its code base and capabilities have grown as a result of continued
development by the broader geosciences community, an effort organized by the
CICE Consortium.

Dr. Cecilia Bitz implemented support for the CICE model (as part of CESM) in DART.
The DART model interface was developed to work with CICE's dynamical core on an
Arakawa B-grid. [1]_ When CICE is coupled to POP in CESM, the ocean and sea ice
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
