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

Jeff Steward added (PR 48) the capability to specify a vertical location
if desired. This allows for localization in the vertical. 

  It’s sometimes helpful, even though definitely wrong from a theoretical
  standpoint, to give a vertical location to satellite observations
  (which are integrated quantities). This has been an issue with
  observation-space localization for some time, and this is the standard
  workaround pioneered by Lili Lei and Jeff Whittaker.

A short description of the namelist options
-------------------------------------------

This table is meant to familiarize you with some of the options
available. Until we fully implement automatic documentation generation,
you would be well advised to familiarize yourself with the code. This is
not the full list of namelist variables …

+-------------------------+---------------------+---------------------+
| variable                | default             | meaning             |
+=========================+=====================+=====================+
| x_thin                  | 0                   | Skip this many per  |
|                         |                     | X scan.             |
+-------------------------+---------------------+---------------------+
| y_thin                  | 0                   | Skip this many per  |
|                         |                     | Y scan.             |
+-------------------------+---------------------+---------------------+
| goes_num                | 16                  | GOES Satellite      |
|                         |                     | number.             |
+-------------------------+---------------------+---------------------+
| reject_dqf_1            | .true.              | Bad scan rejection  |
|                         |                     | criteria. If .true. |
|                         |                     | and DQF /= 0, the   |
|                         |                     | scan is rejected.   |
|                         |                     | If .false. any DQF  |
|                         |                     | > 1 rejects the     |
|                         |                     | scan.               |
+-------------------------+---------------------+---------------------+
| verbose                 | .false.             | Run-time output     |
|                         |                     | verbosity           |
+-------------------------+---------------------+---------------------+
| obs_err                 | MISSING_R8          | The observation     |
|                         |                     | error standard      |
|                         |                     | deviation (std dev, |
|                         |                     | in radiance units)  |
|                         |                     | TODO: make this     |
|                         |                     | more sophisticated. |
|                         |                     | You must supply a   |
|                         |                     | value other than    |
|                         |                     | MISSING_R8. Be      |
|                         |                     | aware that the      |
|                         |                     | observation         |
|                         |                     | sequence files      |
|                         |                     | convert this to a   |
|                         |                     | variance.           |
+-------------------------+---------------------+---------------------+
| vloc_pres_hPa           | -1.0                | The vertical        |
|                         |                     | location of this    |
|                         |                     | observation (hPa).  |
|                         |                     | A negative means    |
|                         |                     | there is no         |
|                         |                     | vertical location   |
|                         |                     | (which is typical   |
|                         |                     | for a               |
|                         |                     | ve                  |
|                         |                     | rtically-integrated |
|                         |                     | quantity).          |
+-------------------------+---------------------+---------------------+
