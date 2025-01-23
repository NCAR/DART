Pangu-Weather
==============


Overview
--------

Pangu-Weather model is an AI model for global weather forecasting. The Pangu-Weather model is described and available at `Pangu-Weather <https://github.com/198808xc/Pangu-Weather?tab=readme-ov-file>`_.

Dr. Nuo Chen implemented the support for the Pangu-Weather model in DART based on the CAM-FV DART interface.

Pangu-Weather model was trained with 0.25 degree ERA5 reanalysis data, making the horizontal dimensions fixed at (721, 1440) and the a total 13 vertical levels fixed at (1000hPa, 925hPa, 850hPa, 700hPa, 600hPa, 500hPa, 400hPa, 300hPa, 250hPa, 200hPa, 150hPa, 100hPa and 50hPa **in the exact order**). All state variables are assumed to be on the mass point.

The Pangu DART Interface
=========================

How to perform the assimilation?

Files required:

* Initial ensemble anlysis or reanalysis files in your data storage directory
* Conversion from initial ensemble files to npy files
* DART obs_seq format observation files
* ``landmask.npy`` quater degree landmask created by WPS geogrid.exe 
* ``terrain.npy`` quater degree landmask created by WPS geogrid.exe 
* ``pangu_weather_6.onnx`` Pangu-Weather model
* In your working directory:

  * ``inference_cpu.py`` or ``inference_gpu.py`` depending on whether you have access to GPU to run Pangu-Weather model, modify ``model_path=$dir_to_pangu_weather_model/pangu_weather_6.onnx``
  * ``input.nml``
  * ``sampling_error_correction_table.nc``
  * ``convert_pgout_to_nc.py``; modify ``landmask = $dir_to_landmask/landmask.npy``, ``terrainmask = $dir_to_terrain/terrain.npy``, and ``work_dir = $dir_to_working_directory/``
  * ``convert_dartout_to_npy.py``
  * ``assimilate.sh``
  * ``run_filter.csh``

Steps: 

#. Download ``pangu_weather_6.onnx`` from `Pangu-Weather <https://github.com/198808xc/Pangu-Weather?tab=readme-ov-file>`_.
#. Create an virtual environment with ```pangu-dart-cpu.yml``` or ``pangu-dart-gpu.yml``
#. Prepare the initial ensemble anlysis or reanalysis files and convert them using ``convert_initial_conditions.py``
#. In ``$DART/build_template/`` choose a ``mkmf.template.xx`` that suits your system, then ``cp mkmf.template.xx mkmf.template``
#. In ``$DART/models/pangu/work/``, run ``./quichbuild.sh``. It compiles the necessary DART executables based on the ``mkmf.template`` that you set. Among all the executables generated, ``./filter`` performs the actual assimilation. After the compilation, move `./filter` to the working directory.
#. Modify file locations in ``inference_cpu.py`` (or ``inference_gpu.py``), ``convert_pgout_to_nc.py`` as suggested above
#. Modify the user defined section in ``assimilate.sh`` 
#. Modify the hpc settings in ``run_filter.csh`` and ``assimilate.sh`` 
#. Modify ``input.nml``. See `Namelists`_ and for more details.
#. run ``qsub ./assimilate.sh`` or ``./assimilate.sh`` to perform the assimilation cycle



.. list-table:: User Defined variables
   :widths: 20 10 50
   :header-rows: 1

   * - Item
     - Type 
     - Description     
   * - num_instances
     - integer
     - ensemble size
   * - old_date
     - string
     - cycle start date "yyyy-mm-dd-HH" 
   * - output_dir
     - string
     - directory to store the output files (e.g. output_mean, output_sd, obs_seq.final and ensemble output files )
   * - obs_dir
     - string
     - where the observations are located. Observation should be in DART obs_seq format, see :ref:`observations`.


Namelists 
---------

DART assembles the namelists for all the relevant modules into a single namelist file ``$DART/models/pangu/input.nml``.
Namelists star with an ampersand ``&`` and terminate with a slash ``/``. 
Character strings that contain a ``/`` must be enclosed in quotes to prevent them from interfering with the namelist structure.
Text outside of the ``&`` and ``/`` pair is ignored.

Here is a list of the model_nml variales and default values.

.. code-block:: fortran

    &model_nml
        cam_template_filename               = 'pginput_0001.nc'
        vertical_localization_coord         = 'PRESSURE'
        use_log_vertical_scale              = .false.
        state_variables  =
            'T',     'QTY_TEMPERATURE',         'NA', 'NA', 'UPDATE'
            'U',     'QTY_U_WIND_COMPONENT',    'NA', 'NA', 'UPDATE'
            'V',     'QTY_V_WIND_COMPONENT',    'NA', 'NA', 'UPDATE'
            'Q',     'QTY_SPECIFIC_HUMIDITY',   'NA', 'NA', 'UPDATE'
        assimilation_period_days            = 0
        assimilation_period_seconds         = 21600
        debug_level                         = 0
   /  

utilities like `pert_copies`, `fields_to_perturb`, `perturbation_amplitude`
and options like `no_obs_assim_above_level`, `model_damping_ends_at_level`, `no_normalization_of_scale_heights`, `use_log_vertical_scale` are not currently supported.

.. list-table:: model_nml in input.nml
   :widths: 20 10 50
   :header-rows: 1

   * - Item
     - Type 
     - Description     
   * - cam_template_filename
     - character(len=128)
     - Pangu input template file used to provide configuration information, such as the longitude, latitude, land mask, etc.
   * - vertical_localization_coord
     - character(len=128)
     - The vertical coordinate to which all vertical locations are converted in model_mod. Valid options is "PRESSURE".
   * - use_log_vertical_scale
     - logical
     - Use the log of the vertical distances when interpolating. This is only used for locations having which_vert = VERTISPRESSURE. It should be .true. when vertical_localization_coord = "scaleheight" or "height".  
   * - state_variables 
     - character (len=64) dimension(100)
     - Character string table that includes: column 1. Pangu variable names to be read into the state vector; column 2, the corresponding DART QTY (quantity); cloumn 3 and 4, if a bounded quantity, the minimum and maximum valid values, Column 5. the string 'UPDATE' indicates that the updated values should be written back to the output file. 'NOUPDATE' will skip writing this field at the end of the assimilation.
   * - assimilation_period_days
     - integer 
     - With assimilation_period_seconds, sets the assimilation cycle length. They should match the model forecast step. The common global assimilation window is 0 days, 21600 seconds (6 hours). They also set the assimilation window width.
   * - assimilation_period_seconds
     - integer 
     - See assimilation_period_days   
   * - debug_level
     - integer 
     - Set this to increasingly larger values to print out more debugging information. Note that this can be very verbose. Use with care.



Features not implemented and future development plan
-----------------------------------------------------

* Allow ensemble generation from single intial condition files.
* Implement the ability to discard observations at too high or low levels, including damping options. 
* Build the support for vertical localization in HEIGHT, SCALEHEIGHT, and LEVEL is the PRESSURE coordinate.
* Assimilation of the surface variables (MSLP, U10, V10, T2M)
* Ability to specify the model pressure level in the namelist or read from the input file
