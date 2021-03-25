GMI Brightness Temperatures
============================

This directory contains the code to convert the GMI Brightness
Temperatures in HDF5 format to the DART observation sequence file
format.

The dataset of interest is: “GPM GMI Common Calibrated Brightness
Temperatures Collocated L1C 1.5 hours 13 km V05 (GPM_1CGPMGMI) at GES
DISC” **not** the _R set! The *short name* for this dataset is
‘GPM_1CGPMGMI’.

The introductory paragraph for the dataset is:

   Version 5 is the current version of the data set. Version 4 is no
   longer available and has been superseded by Version 5. All 1C
   products have a common L1C data structure, simple and generic. Each
   L1C swath includes scan time, latitude and longitude, scan status,
   quality, incidence angle, Sun glint angle, and the intercalibrated
   brightness temperature (Tc). One or more swaths are included in a
   product. The radiometer data are recalibrated to a common basis so
   that precipitation products derived from them are consistent. 1CGMI
   contains common calibrated brightness temperatures from the GMI
   passive microwave instrument flown on the GPM satellite. 1C-R GMI is
   a remapped version of 1CGMI which is explained at the end of this
   section. Swath S1 has 9 channels which are similar to TRMM TMI (10V
   10H 19V 19H 23V 37V 37H 89V 89H). Swath S2 has 4 channels similar to
   AMSU-B (166V 166H 183+/-3V 183+/-8V). Data for both swaths is
   observed in the same revolution of the instrument.

The citation information for this dataset is:
   Title: GPM GMI Common
   Calibrated Brightness Temperatures Collocated L1C 1.5 hours 13 km V05 
   Version: 05 
   Creator: Wesley Berg 
   Publisher: Goddard Earth Sciences
   Data and Information Services Center (GES DISC) > Release Date:
   2016-03-03T00:00:00.000Z 

Linkage:
https://disc.gsfc.nasa.gov/datacollection/GPM_1CGPMGMI_05.html

Instructions to download the GPM_1CGPMGMI dataset for the GMI converter
------------------------------------------------------------------------

1. Go to https://earthdata.nasa.gov
2. Log in (or create an account if necessary)
3. Search for GMI L1C (the “c” here is for cross-calibrated with other
   satellites)
4. Scroll down past datasets to “Matching results.”

-  Follow the link to the GMI common calibrated data set: “GPM GMI
   Common Calibrated Brightness Temperatures Collocated L1C 1.5 hours 13
   km V05 (GPM_1CGPMGMI) at GES DISC” dataset (**NOT** the _R set)

5. You should now be at the
   https://cmr.earthdata.nasa.gov/search/concepts/C1383813813-GES_DISC.html
   page.

-  Select the ‘Download data’ tab
-  Select ‘Earthdata search’
-  Select the GPM link under ‘Matching datasets’

6. You can now select ‘Granule filters’ to choose your start and end
   dates.
7. Select the granules you want, then click ‘download all’ and ‘download
   data’
8. Click download access script
9. Follow the instructions on that page to download the data.

| Each granule is about 28M and has names like:
| ``1C.GPM.GMI.XCAL2016-C.20160621-S001235-E014508.013137.V05A.HDF5``

Guidelines for converting the observations, thinning, superobbing, etc.
are forthcoming. For more background on assimilating radiances in DART,
please read https://dart.ucar.edu/pages/Radiance_support.html

When running the DART converter, two swaths (S1, S2) are converted to
observations. S1 and S2 have different channels and different
“postings,” meaning actual observation locations. They are more or less
right next to each other …

https://disc.gsfc.nasa.gov/datasets/GPM_1CGPMGMI_05/summary

| Swath S1 has 9 channels which are similar to TRMM TMI (10V 10H 19V 19H
  23V 37V 37H 89V 89H).
| Swath S2 has 4 channels similar to AMSU-B (166V 166H 183+/-3V
  183+/-8V).
| Data for both swaths is observed in the same revolution of the
  instrument.

Partial run-time output for one file (no thinning, whole globe,
i.e. about 8 million observations):

::

     ...
      Data Metadata: observation
        QC Metadata: GMI QC
    First timestamp: day=151747, sec=6309
      calendar Date: 2016 Jun 21 01:45:09
     Last timestamp: day=151747, sec=11863
      calendar Date: 2016 Jun 21 03:17:43
      Number of obs processed  :               5734296
      ---------------------------------------------------------
                         GPM_1_GMI_TB 5734296 obs
    
     add_swath_observations:  Converted      5734296  obs for swath /S1; total GMI obs =      5734296
    
      Data Metadata: observation
        QC Metadata: GMI QC
    First timestamp: day=151747, sec=6309
      calendar Date: 2016 Jun 21 01:45:09
     Last timestamp: day=151747, sec=11863
      calendar Date: 2016 Jun 21 03:17:43
      Number of obs processed  :               8279480
      ---------------------------------------------------------
                         GPM_1_GMI_TB 8279480 obs
    
     add_swath_observations:  Converted      2545184  obs for swath /S2; total GMI obs =      8279480

     write_obs_seq  opening unformatted observation sequence file "obs_seq.gmi"
     write_obs_seq  closed observation sequence file "obs_seq.gmi"
     convert_gmi_L1.f90 Finished successfully.
     ...
