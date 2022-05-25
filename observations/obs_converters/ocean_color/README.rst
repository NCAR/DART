VIIRS/AQUA Satellite Ocean Color
================================

Overview
--------

The `OCEAN COLOR <https://oceandata.sci.gsfc.nasa.gov/>`__ (NASA's ocean color web) supports the collection,
processing, calibration, archive and distribution of ocean-related products from a number of missions
which are supported within the framework and facilities of the NASA Ocean Data Processing System. This
converter support ocean color data from two different missions:

1. `AQUA/MODIS <https://oceancolor.gsfc.nasa.gov/data/aqua/>`__ (or Moderate Resolution Imaging
Spectroradiometer) is a key instrument aboard the Terra (EOS AM) and Aqua (EOS PM).
Terra MODIS and Aqua MODIS view the entire Earth's surface every 2 days.

2. `VIIRS <https://oceancolor.gsfc.nasa.gov/data/viirs-snpp/>`__ (or Visible and Infrared Imager/Radiometer
Suite) is a multi-disciplinary instrument that is being flown on the Joint Polar Satellite System. VIIRS
is the successor to MODIS for Earth Science data product generation.

Data source
-----------

In order to have access to the data, you will need a log in on: `https://earthdata.nasa.gov/ <https://earthdata.nasa.gov/>`__

Once a username and a password have been created, the ocean color data can be downloaded from:
`data-access-webpage <https://oceandata.sci.gsfc.nasa.gov/api/file_search>`__ where a manual file
search can be performed or using the provided script: ``shell_scripts/get_ocdata.sh``

Programs and Scripts
--------------------

The programs and scripts in the ``observations/obs_converters/ocean_color`` directory download the data
and create DART observation sequence (obs_seq) files. All programs are built in the ``ocean_color/work``
directory and can be done by running ``./quickbuild.sh`` In addition to the main converter (i.e.,
``convert_sat_chl``), other programs such ``advance_time`` and ``preprocess`` utilities will be built.

**Converter namelist** ``convert_sat_chl_nml``:
This namelist is added to the rest of DART program namelists in file ``input.nml``. Namelists start
with an ampersand '&' and terminate with a slash '/'.

::

   &convert_sat_chl_nml
      file_in         = '../data/V2020350.L3m_DAY_SNPP_CHL_chlor_a_4km.nc'
      file_out        = '../data/obs_seq_chl'
      chl_thresh      = 0.03
      subsample_intv  = 1
      special_mask    = .true.
      debug           = .false.
      /

.. container::

  +-----------------+-----------+-----------------------------------------------------------------------------------------+
  | Contents        | Type      | Description                                                                             |
  +=================+===========+=========================================================================================+
  | file_in         | character | Name of the input netcdf data file. Example of absolute path:                           |
  |                 |           | ``$DART/observations/obs_converters/ocean_color/data/V2020336.L3m_DAY_SNPP_CHL_4km.nc`` |
  +-----------------+-----------+-----------------------------------------------------------------------------------------+
  | file_out        | character | Partial filename for the output file.  The date and time are appended to ``file_out``   |
  |                 |           | to construct a unique filename reflecting the time of the observations in the file.     |
  +-----------------+-----------+-----------------------------------------------------------------------------------------+
  | chl_thresh      | real(r8)  | When the observed chlorophyll values are small, a threshold value is used for the obs.  |
  |                 |           | Example: ``chl_thresh = 0.03``                                                          |
  +-----------------+-----------+-----------------------------------------------------------------------------------------+
  | subsample_intv  | integer   | It is possible to 'thin' the observations. ``subsample_intv``                           |
  |                 |           | allows one to take every nth observation.                                               |
  +-----------------+-----------+-----------------------------------------------------------------------------------------+
  | special_mask    | logical   | A simple procedure to ignore data in certain areas.                                     |
  |                 |           | Users can the code in convert_sat_chl.f90 to change the mask location                   |
  +-----------------+-----------+-----------------------------------------------------------------------------------------+
  | debug           | logical   | Print extra information during the ``convert_sat_chl`` execution.                       |
  +-----------------+-----------+-----------------------------------------------------------------------------------------+

It's widely known that chlorophyll follows a log-normal distribution. After reading the observed chlorophyll data
from the input netcdf file, a **LOG10** transformation is applied such that the resulting distribution is normal.
The idea behind this step is to prepare the data for the EnKF which assumes both the state and observations to
be Gaussian.

The ``get_ocdata.sh`` script is placed inside the ``shell_scripts`` directory. Technically, this is the only the script that the
user needs to run. Prior to running the script, one should edit some parameters such as: date range, the resolution
of the data, the instrument, username and password to access the data site, frequency, domain area coordinates, etc.
These parameters are hard coded in the script. The tasks that this script do are:

#. It Downloads ocean color data for the requested period and frequency one file at a time.
#. It retains the chlorophyll variable for the requested domain and gets rid of unnecessary information in the netcdf data files.
#. It modifies the netcdf files by adding a time dimension with the necessary attributes.
#. It runs the obs converter ``convert_sat_chl`` for each data file.

So, each data file will have a corresponding obs_seq file. The observation sequence files are named such that
the observation time is part of the file name. For example: ``obs_seq_chl.2020-12-15-00000``

This script is *NOT* intended to be turnkey script; it will certainly needs to be customized for your use. This script
will make use of some DART utilities and namelists such as ``advance_time`` and ``input.nml``. If an existing observation
sequence file of the same output name is found when the converter is run again,
it will replace that file.
