SSUSI F16 EDR-DSK format to observation sequence converters
===========================================================

Overview
--------

The Special Sensor Ultraviolet Spectrographic Imager `SSUSI <http://http://ssusi.jhuapl.edu/>`__ is designed to remotely
sense the ionosphere and thermosphere. The following is repeated from the SSUSI home page:

   *Overview
   Beginning in 2003, the Defense Meteorological Satellite Program (DMSP) satellites began carrying the SSUSI instrument
   - a combination of spectrographic imaging and photometric systems designed to remotely sense the ionosphere and
   thermosphere.
   The long term focus of the SSUSI program is to provide data concerning the upper atmospheric response to the sun over
   the changing conditions of the solar cycle. Data collected by SSUSI instrument can help identify structure in the
   equatorial and polar regions.
   Mission
   SSUSI was designed for the DMSP Block 5D-3 satellites. These satellites are placed into nearly polar, sun-synchronous
   orbits at an altitude of about 850 km. SSUSI is a remote-sensing instrument which measures ultraviolet (UV) emissions
   in five different wavelength bands from the Earth's upper atmosphere. SSUSI is mounted on a nadir-looking panel of
   the satellite. The multicolor images from SSUSI cover the visible Earth disk from horizon to horizon and the
   anti-sunward limb up to an altitude of approximately 520 km.
   The UV images and the derived environmental data provide the Air Force Weather Agency (Offutt Air Force Base,
   Bellevue, NE) with near real-time information that can be utilized in a number of applications, such as maintenance
   of high frequency (HF) communication links and related systems and assessment of the environmental hazard to
   astronauts on the Space Station.*

| ``convert_f16_edr_dsk.f90`` will extract the ON2 observations from the F16 "edr-dsk" format files and create DART
  observation sequence files. There is one additional preprocessing step before the edr-dsk files may be converted.
| The ON2_UNCERTAINTY variable in the netcdf files have IEEE NaN values, but none of the required metadata to interpret
  them correctly. These 2 lines will add the required attributes so that NaNs are replaced with a fill value that can be
  queried and manipulated. Since the ON2_UNCERTAINTY is a standard deviation, it is sufficient to make the fill value
  negative. See the section on Known Bugs

.. container:: unix

   ::

      ncatted -a _FillValue,ON2_UNCERTAINTY,o,f,NaN        input_file.nc
      ncatted -a _FillValue,ON2_UNCERTAINTY,m,f,-1.0       input_file.nc

Data sources
------------

http://ssusi.jhuapl.edu/data_products

Please read their `data usage <http://ssusi.jhuapl.edu/home_data_usage>`__ policy.

Programs
--------

``DART/observations/SSUSI/convert_f16_edr_dsk.f90`` will extract ON2 data from the distribution files and create DART
observation sequence (obs_seq) files. Build it in the ``SSUSI/work`` directory by running the ``./quickbuild.sh``
script located there. In addition to the converters, the ``advance_time`` and ``obs_sequence_tool`` utilities will be
built.

An example data file is in the ``data`` directory. An example scripts for adding the required metadata to the
ON2_UNCERTAINTY variable in the ``shell_scripts`` directory. These are *NOT* intended to be turnkey scripts; they will
certainly need to be customized for your use. There are comments at the top of the scripts saying what options they
include, and should be commented enough to indicate where changes will be likely to need to be made.

Errors
------

The code for setting observation error variances is using fixed values, and we are not certain if they are correct.
Incoming QC values larger than 0 are suspect, but it is not clear if they really signal unusable values or whether there
are some codes we should accept.

Known Bugs
----------

The netCDF files - as distributed - have NaN values to indicate "MISSING".
This makes it exceptionally hard to read or work with, as almost everything
will core dump when trying to perform any math with NaNs.
``convert_f16_edr_dsk.f90`` tries to count how many values are missing. If the
NaN has not been replaced with a numerically valid MISSING value, the following
FATAL ERROR is generated (by the Intel compiler, with debug and traceback enabled):


.. container:: unix 

  ::

     set_nml_output Echo NML values to log file only
     Trying to open namelist log dart_log.nml
     forrtl: error (65): floating invalid
     Image              PC                Routine            Line        Source             
     convert_f16_edr_d  000000000051717D  MAIN__                    143  convert_f16_edr_dsk.f90
     convert_f16_edr_d  0000000000409B3C  Unknown               Unknown  Unknown
     libc.so.6          0000003101E1ED5D  Unknown               Unknown  Unknown
     convert_f16_edr_d  0000000000409A39  Unknown               Unknown  Unknown
     Abort (core dumped)


The solution is to replace the NaN values with a viable MISSING value using
the ``shell_scripts/netcdf_manip.csh`` script.
It relies on the netCDF Operators, freely available 
http://nco.sourceforge.net

