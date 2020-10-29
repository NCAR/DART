###########
CICE README
###########

Contents
========

#. `Overview`_
#. `Development Notes`_
#. `Namelist`_
#. `Terms of Use`_
#. `References`_

Overview
========

The `Community Ice CodE (CICE) <https://github.com/CICE-Consortium/CICE>`_ is a
sea ice model that was first developed by Elizabeth Hunke as the Los Alamos Sea
Ice Model. Its code base and capabilities have grown as a result of continued
development by the broader geosciences community, an effort organized by the
CICE Consortium.

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

- U, V are at grid cell corners
- T, h, hs, and the various scalar quantities are at grid cell centers

Development Notes
=================

CICE is under active development to work with other grids, such as the
unstructured grid in MPAS and the C-grid in MOM. Due to this activity, this
README contains development notes chronicling the development of the model
interface. We ask that developers continue to document the development both in
this README and with descriptive comments in the source code.

- Possible bug found in ``model_mod.f90`` for POP where set_date is sent sec
  this day. The routine wants sec this min.

Notes from Cecilia M. Bitz on 14 May 2016
-----------------------------------------

- Created ``../../obs_def/obs_def_cice_mod.f90`` to make new obs_kinds used in
  ``model_mod.f90`` and ``input.nml``.
- Not sure about ``QTY_TRACERARRAY_CATS``
- Model mod assumes that the grid is identical to POP
- Leaving this part but it may be unneeded in CICE ``INTERFACE
  vector_to_prog_var MODULE PROCEDURE vector_to_2d_prog_var
  MODULE PROCEDURE vector_to_3d_prog_var END INTERFACE``
- Not used in pop so of course not used now in CICE either, why? ``subroutine
  vector_to_3d_prog_var(x, varindex,
  data_3d_array) subroutine get_gridsize(num_x, num_y, num_z)``
- Come back here, some changes made below but need to look line-by-line still
  ``subroutine get_state_meta_data(state_handle, index_in, location,
  var_type)``

Fortran Files
~~~~~~~~~~~~~

- ``dart_to_cice.f90`` Think it is done
- ``cice_to_dart.f90`` Is trivial so hope it's done too
- ``model_mod_check.f90``
- ``dart_cice_mod.f90`` Should it have a get_cat_dim?
- ``model_mod.f90`` I do not understand this part, but appears in clm too:
  ``INTERFACE vector_to_prog_var MODULE PROCEDURE vector_to_1d_prog_var ! this
  is in clm MODULE PROCEDURE vector_to_2d_prog_var ! this is in pop MODULE
  PROCEDURE vector_to_3d_prog_var ! this is in pop END INTERFACE``
- ``test_dipole_interp.f90`` Also trivial, nothing to change?
- ``rma-cice/location/`` Has a bunch of subdirs each with a location_mod.f90

.. code-block:: bash

  -rw-r--r--  1 bitz  staff  142664 May 26 17:12 model_mod.f90
  -rw-r--r--  1 bitz  staff    4439 May 21 07:55 cice_to_dart.f90
  -rw-r--r--  1 bitz  staff    5676 May 21 07:49 dart_to_cice.f90
  -rw-r--r--  1 bitz  staff   24008 May 18 21:55 dart_cice_mod.f90
  -rw-r--r--  1 bitz  staff   24294 May 14 16:30 model_mod_check.f90
  -rw-r--r--  1 bitz  staff    2270 May 14 16:30 test_dipole_interp.f90

Code Snippet from model_mod.f90
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: fortran

  ! these routines must be public and you cannot change
  ! the arguments - they will be called *from* the DART code.
  public :: get_model_size,                &
            adv_1step,                     &
            get_state_meta_data,           &
            model_interpolate,             &
            get_model_time_step,           &
            static_init_model,             &
            end_model,                     &
            init_time,                     &
            init_conditions,               &
            nc_write_model_atts,           &
            nc_write_model_vars,           &
            pert_model_copies,             &
            get_close_maxdist_init,        &
            get_close_obs_init,            &
            get_close_obs,                 &
            query_vert_localization_coord, &
            vert_convert,                  &
            construct_file_name_in,        &
            read_model_time,               &
            write_model_time

Namelist
========

.. code-block:: fortran

  &model_nml
     assimilation_period_days     = 1
     assimilation_period_seconds  = 0
     model_perturbation_amplitude = 0.00002
     binary_grid_file_format      = 'big_endian'
     debug                        = 1
     model_state_variables        = 'aicen', 'QTY_SEAICE_CONCENTR',   'UPDATE',
                                    'vicen', 'QTY_SEAICE_VOLUME',     'UPDATE',
                                    ...
                                    'vsnon', 'QTY_SEAICE_SNOWVOLUME', 'UPDATE',
  /

Description of each namelist entry
----------------------------------

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

Terms of Use
============

|Copyright| University Corporation for Atmospheric Research

Licensed under the `Apache License, Version 2.0
<http://www.apache.org/licenses/LICENSE-2.0>`__. Unless required by applicable
law or agreed to in writing, software distributed under this license is
distributed on an "as is" basis, without warranties or conditions of any kind,
either express or implied.

.. |Copyright| unicode:: 0xA9 .. copyright sign

References
==========

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
