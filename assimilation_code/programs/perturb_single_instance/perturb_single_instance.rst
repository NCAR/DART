PROGRAM ``perturb_single_instance``
===================================

Overview
--------

Utility program to generate an ensemble of perturbed ensemble member restart files. This program can be run in parallel
and used as a stand alone program.

Namelist
--------

This namelist is read from the file ``input.nml``. Namelists start with an ampersand '&' and terminate with a slash '/'.
Character strings that contain a '/' must be enclosed in quotes to prevent them from prematurely terminating the
namelist.

::

   &perturb_single_instance
      ens_size               = ''
      input_files            = ''      
      output_files           = ''
      output_file_list       = ''
      perturbation_amplitude = 0.0     
      single_restart_file_in = .false.      
     /

.. container::

   +------------------------+-------------------------------------------+---------------------------------------------+
   | Item                   | Type                                      | Description                                 |
   +========================+===========================================+=============================================+
   | ens_size               | integer                                   | Total number of ensemble members.           |
   +------------------------+-------------------------------------------+---------------------------------------------+
   | input_files            | character(len=256),dimension(num_domains) | The restart file you would like to perturb  |
   |                        |                                           | from.                                       |
   +------------------------+-------------------------------------------+---------------------------------------------+
   | output_file_list       | character(len=256)                        | A file containing a list of the desired     |
   |                        |                                           | output names.                               |
   +------------------------+-------------------------------------------+---------------------------------------------+
   | output_files           | character(len=256)                        | An array of filenames                       |
   +------------------------+-------------------------------------------+---------------------------------------------+
   | perturbation_amplitude | real(r8)                                  | The desired perturbation amplitude. If the  |
   |                        |                                           | model provides an interface then it will    |
   |                        |                                           | use that subroutine, otherwise it will      |
   |                        |                                           | simply add gaussian noise to the entire     |
   |                        |                                           | state, and this is the standard deviation.  |
   +------------------------+-------------------------------------------+---------------------------------------------+
   | single_restart_file_in | logical                                   | A boolean, specifying if you have a single  |
   |                        |                                           | file restart, such as the case for lower    |
   |                        |                                           | order models.                               |
   +------------------------+-------------------------------------------+---------------------------------------------+

Below is an example of a typical namelist for the perturb_single_instance.

::

   &perturb_single_instance_nml
      ens_size         = 3
      input_files      = 'caminput.nc'
      output_files     = 'cam_pert1.nc','cam_pert2.nc','cam_pert3.nc'
   /

| 

Files
-----

-  inputfile.nc (description file that will be perturbed)
-  output_file_list.txt (a file containing a list of restart files) and,
-  perturb_single_instance.nml

References
----------

-  none
