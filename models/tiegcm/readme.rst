.. _tiegcm:

TIEGCM
======


Overview
--------

The Thermosphere Ionosphere Electrodynamic General Circulation Model 
(`TIEGCM <http://www.hao.ucar.edu/modeling/tgcm/tie.php>`__) is developed by the NSF NCAR
High Altitude Observatory (`HAO <https://www2.hao.ucar.edu/>`__).


DART-TIEGCM has been used to assimilate neutral mass density
retrieved from satellite-borne accelerometers and electron density obtained from ground-based and space-based GNSS
signals. Unlike other ionospheric data assimilation applications, this approach allows simultaneous assimilation of
thermospheric and ionospheric parameters by taking advantage of the coupling of plasma and neutral constituents
described in TIEGCM. DART/TIEGCM's demonstrated capability to infer under-observed thermospheric parameters from
abundant electron density observations has important implications for the future of upper atmosphere research.

DART is designed so that the TIEGCM source code can be used with no modifications.  TIEGCM and DART run as separate
executables.
The TIEGCM 2.0 source code and User's Guide is available from HAO:

- `TIEGCM 2.0 source code <http://www.hao.ucar.edu/modeling/tgcm/download.php>`__

- `TIEGCM 2.0 User's Guide <https://www.hao.ucar.edu/modeling/tgcm/tiegcm2.0/userguide/html/>`__


DART-TIEGCM namelist options
----------------------------

The ``model_nml`` namelist contains the TIEGCM specific options for DART.
``model_nml`` is read from the file ``input.nml``.
Namelists start with an ampersand '&' and terminate with a slash '/'.
Character strings that contain a '/' must be enclosed in quotes to prevent them from prematurely terminating the
namelist.

.. code-block:: text

   &model_nml 
      tiegcm_restart_file_name    = 'tiegcm_restart_p.nc'
      tiegcm_secondary_file_name  = 'tiegcm_s.nc'
      model_res                   = 5.0
      assimilation_period_seconds = 3600
      estimate_f10_7              = .false.
      f10_7_file_name             = 'f10_7.nc'
      debug                       = 0
      variables = 'NE',    'QTY_ELECTRON_DENSITY',          '1000.0',  'NA',      'restart',    'UPDATE'
                  'OP',    'QTY_DENSITY_ION_OP',            'NA',      'NA',      'restart',    'UPDATE',
                  'TI',    'QTY_TEMPERATURE_ION',           'NA',      'NA',      'restart',    'UPDATE',
                  'TE',    'QTY_TEMPERATURE_ELECTRON',      'NA',      'NA',      'restart',    'UPDATE',
                  'OP_NM', 'QTY_DENSITY_ION_OP',            'NA',      'NA',      'restart',    'UPDATE',
                  'O1',    'QTY_ATOMIC_OXYGEN_MIXING_RATIO','0.00001', '0.99999', 'secondary',  'NO_COPY_BACK',
                  'O2',    'QTY_MOLEC_OXYGEN_MIXING_RATIO', '0.00001', '0.99999', 'secondary',  'NO_COPY_BACK',
                  'TN',    'QTY_TEMPERATURE',               '0.0',     '6000.0',  'secondary',  'NO_COPY_BACK',
                  'ZG',    'QTY_GEOMETRIC_HEIGHT',          'NA',      'NA',      'secondary',  'NO_COPY_BACK',
      /



+-----------------------------+----------------------+---------------------------------------+
| Namelist entry              | Type                 | Description                           |
+=============================+======================+=======================================+
| tiegcm_restart_file_name    | character(len=256)   | The TIEGCM restart template           |
+-----------------------------+----------------------+---------------------------------------+
| tiegcm_secondary_file_name  | character(len=256)   | The TIEGCM secondary template         |
+-----------------------------+----------------------+---------------------------------------+
| model_res                   | real(r8)             | TIEGCM model resolution 5.0 or 2.5    |
|                             |                      | degrees                               |
+-----------------------------+----------------------+---------------------------------------+
| assimilation_period_seconds | integer              | This specifies the width of the       |
|                             |                      | assimilation window. The current      |
|                             |                      | model time is used as the center time |
|                             |                      | of the assimilation window. All       |
|                             |                      | observations in the assimilation      |
|                             |                      | window are assimilated. BEWARE: if    |
|                             |                      | you put observations that occur       |
|                             |                      | before the beginning of the           |
|                             |                      | assimilation_period, DART will error  |
|                             |                      | out because it cannot move the model  |
|                             |                      | 'back in time' to process these       |
|                             |                      | observations.                         |
|                             |                      | ``assimilation_period_seconds`` must  |
|                             |                      | be an integer number of TIEGCM        |
|                             |                      | dynamical timesteps (as specified by  |
|                             |                      | tiegcm.nml:STEP) AND be able to be    |
|                             |                      | expressed by tiegcm.nml:STOP. Since   |
|                             |                      | STOP has three components:            |
|                             |                      | day-of-year, hour, and minute, the    |
|                             |                      | ``assimilation_period_seconds`` must  |
|                             |                      | be an integer number of minutes.      |
+-----------------------------+----------------------+---------------------------------------+
| estimate_f10_7              | logical              | Switch to specify that the f10.7      |
|                             |                      | index should be estimated by          |
|                             |                      | augmenting the DART state vector with |
|                             |                      | a scalar. The location of the f10.7   |
|                             |                      | index is taken to be longitude of     |
|                             |                      | local noon and latitude zero.         |
+-----------------------------+----------------------+---------------------------------------+
| f10_7_file_name             | character(len=256)   | If ``estimate_f107=.true.``           |
|                             |                      | f10.7 will be part of the dart state. |
|                             |                      | The variable f10_7 is read from       |
|                             |                      | ``f10_7_file_name``                   |
|                             |                      | An example f10_7.cdl file is          |
|                             |                      | given in the work directory           |
+-----------------------------+----------------------+---------------------------------------+
| debug                       | integer              | Set to 0 (zero) for minimal output.   |
|                             |                      | Successively larger values generate   |
|                             |                      | successively more output.             |
+-----------------------------+----------------------+---------------------------------------+
| variables                   | character            | Six strings to describe the TIEGCM    |
|                             | (MAX_NUM_VARIABLES * | variables to be used in DART.         |
|                             | 6)                   | A description of the six strings is   |
|                             |                      | given below.                          |
+-----------------------------+----------------------+---------------------------------------+


::

      variables = 'NAME', 'QTY', 'MIN', 'MAX', 'FILE', 'UPDATE'


``NAME`` The variable name in the TIEGCM netCDF file. 

``QTY`` The DART quantity for the variable.

``MIN`` The minimum bound (if any) for the variable. Enter 'NA' for no minimum.

``MAX`` The a maximum bound (if any) for the variable.  

``FILE`` The tiegcm netcdf file containing the variable. 'restart' or 'secondary'

``UPDATE`` filter will update the variable in the TIEGCM netcdf file. Use ``NO_COPY_BACK`` to prevent
filter from updating the variable.

Below is an example showing the namelist options necessary to add f10.7 to the DART state

:: 

      &filter_nml
         input_state_file_list        = 'restart_p_files.txt', 'secondary_files.txt', 'f10.7.txt'
         output_state_file_list       = 'out_restart_p_files.txt', 'out_secondary_files.txt', 'out_f10.7.txt' 
      
      &model_nml
        estimate_f10_7 = .true.
        f10_7_file_name = 'f10_7.nc'
        variables =  'NE',    'QTY_ELECTRON_DENSITY',          '1000.0',  'NA',      'restart',    'UPDATE'
                     ...
                     'ZG',    'QTY_GEOMETRIC_HEIGHT',          'NA',      'NA',      'secondary',  'NO_COPY_BACK',
                     'f10_7'  'QTY_1D_PARAMETER'               'NA',      'NA',      'calculate', 'UPDATE'


References
----------

-  Matsuo, T., and E. A. Araujo-Pradere (2011),
   Role of thermosphere-ionosphere coupling in a global ionosphere specification,
   *Radio Science*, **46**, RS0D23, `doi:10.1029/2010RS004576 <http://dx.doi.org/doi:10.1029/2010RS004576>`__
  
-  Lee, I. T., T, Matsuo, A. D. Richmond, J. Y. Liu, W. Wang, C. H. Lin, J. L. Anderson, and M. Q. Chen (2012),
   Assimilation of FORMOSAT-3/COSMIC electron density profiles into thermosphere/Ionosphere coupling model by using
   ensemble Kalman filter,
   *Journal of Geophysical Research*, **117**, A10318,
   `doi:10.1029/2012JA017700 <http://dx.doi.org/doi:10.1029/2012JA017700>`__
  
-  Matsuo, T., I. T. Lee, and J. L. Anderson (2013),
   Thermospheric mass density specification using an ensemble Kalman filter,
   *Journal of Geophysical Research*, **118**, 1339-1350,
   `doi:10.1002/jgra.50162 <http://dx.doi.org/doi:10.1002/jgra.50162>`__
  
-  Lee, I. T., H. F. Tsai, J. Y. Liu, Matsuo, T., and L. C. Chang (2013),
   Modeling impact of FORMOSAT-7/COSMIC-2 mission on ionospheric space weather monitoring,
   *Journal of Geophysical Research*, **118**, 6518-6523,
   `doi:10.1002/jgra.50538 <http://dx.doi.org/doi:10.1002/jgra.50538>`__
  
-  Matsuo, T. (2014),
   Upper atmosphere data assimilation with an ensemble Kalman filter, in Modeling the Ionosphere-Thermosphere System,
   *Geophys. Monogr. Ser.*, vol. 201, edited by J. Huba, R. Schunk, and G. Khazanov, pp. 273-282, John Wiley & Sons,
   Ltd, Chichester, UK, `doi:10.1002/9781118704417 <http://dx.doi.org/doi:10.1002/9781118704417>`__
  
-  Hsu, C.-H., T. Matsuo, W. Wang, and J. Y. Liu (2014),
   Effects of inferring unobserved thermospheric and ionospheric state variables by using an ensemble Kalman filter on
   global ionospheric specification and forecasting,
   *Journal of Geophysical Research*, **119**, 9256-9267,
   `doi:10.1002/2014JA020390 <http://dx.doi.org/doi:10.1002/2014JA020390>`__
  
-  Chartier, A., T. Matsuo, J. L. Anderson, G. Lu, T. Hoar, N. Collins, A. Coster, C. Mitchell, L. Paxton, G. Bust
   (2015),
   Ionospheric Data Assimilation and Forecasting During Storms,
   *Journal of Geophysical Research*, `doi:10.1002/2014JA020799 <https://doi.org/10.1002/2014JA020799>`__
  
