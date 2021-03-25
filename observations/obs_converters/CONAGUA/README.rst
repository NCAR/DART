CONAGUA
=======

The streamflow observations from CONAGUA are naturally in a Microsoft 
database format. Mirce converts these one-at-a-time to a csv format.
The filenames have a gage identifier in them, there is another file
that has the lat/lon of the gage.

.. code-block:: text

   /glade/scratch/mirce/LaSierra/Observations/

The existing DART csv readers are:

.. code-block:: text

   vi -R Ameriflux/level4_to_obs.f90 \
   CHAMP/CHAMP_density_text_to_obs.f90 \
   CNOFS/CNOFS_text_to_obs.f90 \
   COSMOS/COSMOS_development.f90 \
   COSMOS/COSMOS_to_obs.f90 \
   MODIS/MOD15A2_to_obs.f90 \
   ROMS/convert_roms_obs.f90 \
   gnd_gps_vtec/gnd_gps_vtec_text_to_obs.f90 \
   gps/convert_cosmic_gps_cdf.f90 \
   gps/convert_cosmic_ionosphere.f90 \
   quikscat/quikscat_JPL_mod.f90 \
   snow/snow_to_obs.f90 \
   text/text_to_obs.f90 \
   text_GITM/text_to_obs.f90

One of these should be close enough. Some are more sophisticated in that
they try to determine which column contains the string that identifies the year, mondy, day, etc. - 
as opposed to hardcoding the knowledge about which column is which.

These are the meanings for each of the column headers in the daily observation files:
pk_anio = Year
pk_mes = Month
ngasto_d01, d02 ...and so on up to d31 = Streamflow in day 01, day 02 ...day 31
The streamflow is in cms

