NOAA GOES-R Series Advanced Baseline Imager (ABI) Level 1b Radiances
====================================================================

The data are available from
`NOAA-NCEI <https://data.nodc.noaa.gov/cgi-bin/iso?id=gov.noaa.ncdc:C01501>`__

The **convert_goes_ABI_L1b** program converts ABI Level 1b Radiances in
netCDF format to a DART observation sequence file with
``GOES_16_ABI_RADIANCE`` observations (there is a namelist option to
select other GOES satellites, which will have the appropriate
observation type).

   The Advanced Baseline Imager (ABI) instrument samples the radiance of
   the Earth in sixteen spectral bands using several arrays of detectors
   in the instrument’s focal plane. Single reflective band ABI Level 1b
   Radiance Products (channels 1 - 6 with approximate center wavelengths
   0.47, 0.64, 0.865, 1.378, 1.61, 2.25 microns, respectively) are
   digital maps of outgoing radiance values at the top of the atmosphere
   for visible and near-infrared (IR) bands. Single emissive band ABI
   L1b Radiance Products (channels 7 - 16 with approximate center
   wavelengths 3.9, 6.185, 6.95, 7.34, 8.5, 9.61, 10.35, 11.2, 12.3,
   13.3 microns, respectively) are digital maps of outgoing radiance
   values at the top of the atmosphere for IR bands. Detector samples
   are compressed, packetized and down-linked to the ground station as
   Level 0 data for conversion to calibrated, geo-located pixels (Level
   1b Radiance data). The detector samples are decompressed,
   radiometrically corrected, navigated and resampled onto an invariant
   output grid, referred to as the ABI fixed grid.

   Cite as: GOES-R Calibration Working Group and GOES-R Series Program,
   (2017): NOAA GOES-R Series Advanced Baseline Imager (ABI) Level 1b
   Radiances. [indicate subset used]. NOAA National Centers for
   Environmental Information. doi:10.7289/V5BV7DSR. [access date].

Specifying a vertical location
------------------------------

Top of the atmosphere radiance observations are sensitive to the
atmospheric constituents (e.g. water vapor) that reside in the vertical
profile of the atmosphere. Given this is an integrated quantity and does
not depend on a single vertical location, it may be appropriate to leave
the vertical location undefined (i.e. VERTISUNDEF) within the ``obs_seq.out``
file, by setting the ``vloc_pres_hPa = -1`` (See namelist options below). This approach, 
however, limits the application of vertical localization during the assimilation step.

Alternatively, for some applications it may be appropriate to assign 
a vertical location to the radiance observation, by setting the ``vloc_pres_hPa``
to a vertical pressure level (hPa). This is an ongoing area
of observation-space localization research, and is the standard
workaround pioneered by Lili Lei and Jeff Whittaker.

Radiance versus Brightness Temperature
--------------------------------------

This converter assigns the observation type as ``GOES[16-19]_ABI_RADIANCE``.
The default setup in DART is that this radiance observation type is assigned
the quantity ``QTY_RADIANCE``.  Radiance observations are commonly expressed 
in spectrally resolved units (mW/cm/m^2/sr). 

Alternatively, radiances can also be expressed as brightness temperatures
(units: Kelvin) and the DART code also supports observation quantities of 
``QTY_BRIGHTNESS_TEMPERATURE``. Both the spectral and temperature units
quantify the same physical properties of the atmosphere 
(emitted and reflected radiation), however, in some applications it may
be advantageous to use brightness temperatures given their Gaussian 
distribution.  This is an ongoing area of research.


An overview of the namelist options
-----------------------------------

A description of the most important namelist options. Note that supplying
an observation error value is mandatory. This is not the full list of namelist
variables.

+-------------------------+------------+-----------------------------+
| Variable                | Default    |      Description            |
+=========================+============+=============================+
| x_thin                  | 0          | Skip this many observations |
|                         |            | per X scan                  |
+-------------------------+------------+-----------------------------+
| y_thin                  | 0          | Skip this many observations |
|                         |            | per Y scan                  |
+-------------------------+------------+-----------------------------+
| goes_num                | 16         | GOES Satellite number       |
+-------------------------+------------+-----------------------------+
| reject_dqf_1            | .true.     | Bad scan rejection critera. |
|                         |            | If .true. and DQF /= 0, the |
|                         |            | scan is rejected. If        |
|                         |            | .false. any DQF > 1         |
|                         |            | rejects the scan.           |
+-------------------------+------------+-----------------------------+
| verbose                 | .false.    | Run-time output verbosity   |
+-------------------------+------------+-----------------------------+
| obs_err                 | MISSING_R8 | The observation error       |
|                         |            | standard deviation in units |
|                         |            | of radiance. IMPORTANT:     |
|                         |            | the user must supply a      |
|                         |            | value other than MISSING_R8.|
|                         |            | Be aware that the           |
|                         |            | observation sequence files  |
|                         |            | convert this to a variance. |
+-------------------------+------------+-----------------------------+
| vloc_pres_hPa           | -1.0       | If a positive value, the    |
|                         |            | vertical location of the    |
|                         |            | observation (hPa) is        |
|                         |            | assigned with a vertical    |
|                         |            | coordinate of               |
|                         |            | VERTISPRESSURE. If negative |
|                         |            | value there is no vertical  |
|                         |            | location and the coordinate |
|                         |            | is VERTISUNDEFINED.         |
+-------------------------+------------+-----------------------------+

Radiance metadata supplied to obs_seq.out file
----------------------------------------------

This converter is designed to supply metadata to the radiative transfer
model (RTTOV) :doc:`../../../observations/forward_operators/obs_def_rttov_mod`
that supports the calculation of the expected radiance 
observation during the assimilation step.  Below is a description
of this metadata information.

The integer values that describe the platform/satellite/sensor/channel
combination for the GOES satellite are required by the RTTOV radiance
model to assign the appropriate coefficent file during the radiance
calculation. For more details refer to the 
`RTTOV user guide. <https://www.nwpsaf.eu/site/software/rttov/documentation/>`__


+-------------------------+------------+-----------------------------+
| Variable                | Value      |      Description            |
+=========================+============+=============================+
| sat az/el               | Supplied   | GOES satellite azimuth and  |
|                         | by data    | elevation angle             |
|                         | file       |                             |
+-------------------------+------------+-----------------------------+
| sun az/el               | -888888    | Sun azimuth and elevation   |
|                         |            | angle. The default is for   |
|                         |            | missing values. For IR      |
|                         |            | channels reflected solar    |
|                         |            | has no impact. For NIR/     |
|                         |            | VIS/UV providing sun az/el  |
|                         |            | may be desirable but not    |
|                         |            | required by RTTOV.          |
+-------------------------+------------+-----------------------------+
| platform                | 4          | Platform number supplied    |
|                         |            | by converter code.          |
+-------------------------+------------+-----------------------------+
| sat_id                  | 16         | Satellite ID number         |
|                         |            | supplied by converter code. |
+-------------------------+------------+-----------------------------+
| instrument              | 44         | Sensor number supplied by   |
|                         |            | converter code              |
+-------------------------+------------+-----------------------------+
| channel                 | Supplied   | The wavelength channel      |
|                         | by data    | (1-16) assigned from band_id|
|                         | file       | variable in data file.      |
+-------------------------+------------+-----------------------------+





