MPD
===

The Micro Pulse Differential Absorption Lidar
(MPD) data were collected during field campaigns and testing periods
by the Earth Observing Laboratory (EOL).

The differential absorption lidar (DIAL) technique uses two separate
laser wavelengths: an absorbing wavelength (online) and a non-absorbing
wavelength (offline). The ratio of the range-resolved backscattered
signals between the online and offline wavelengths is proportional to
the amount of water vapor in the atmosphere, which allows the retrieval
of absolute humidity profiles above the lidar site.

This observation converter takes absolute humidity (g/m3) profiles
retrieved from the MPD data and converts them to the format used by
DART. The ``obs_converter/MPD/work/convert_to_text.py`` script reads the
netCDF files from each MPD site and combines them into text files, one
for each date and time. The ``obs_converter/MPD/work/MPD_text_to_obs``
program translates the text files to the DART ``obs_seq.out`` format.

Test data for a single site and an example output can be downloaded from
https://www.image.ucar.edu/pub/DART/MPD/MPD.tar.gz

For more details of the retrieval and quality control process, and
inquire about data availability for your research project, please
contact Tammy Weckwerth at EOL, NSF NCAR.
