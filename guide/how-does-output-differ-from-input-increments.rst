.. role:: bolditalic
  :class: bolditalic



Computing filter increments using a complex model
===================================================

The innovations to the model state are easy to derive. Use the `NCO
Operator <http://nco.sourceforge.net/>`__ *ncdiff* to difference the two DART
diagnostic netCDF files to create the innovations. Be sure to check the
*CopyMetaData* variable to figure out what *copy* is of interest. Then, use
*ncview* to explore the innovations or the inflation values or …

If the assimilation used state-space inflation, the inflation fields will be
added as additional ‘copies’. A sure sign of trouble is if the inflation fields
grow without bound. As the observation network changes, expect the inflation
values to change.

The only other thing I look for in state-space is that the increments are
‘reasonable’. As the assimilation ‘burns in’, the increments are generally
larger than increments from an assimilation that has been cycling for a long
time. If the increments keep getting bigger, the ensemble is continually
drifting away from the observation. Not good. In *ncview*, it is useful to
navigate to the copy/level of interest and re-range the data to values
appropriate to the current data and then hit the ‘>>’ button to animate the
image. It should be possible to get a sense of the magnitude of the innovations
as a function of time.

Example from a model of intermediate complexity: the bgrid model
----------------------------------------------------------------

I ran a perfect model experiment with the bgrid model in the DART-default
configuration and turned on some adaptive inflation for this example. To fully
demonstrate the adaptive inflation, it is useful to have an observation
network that changes through time. I created two observation sequence files:
one that had a single ‘RADIOSONDE_TEMPERATURE’ observation at the surface with
an observation error variance of 1.5 degrees Kelvin - repeated every 6 hours
for 6 days (24 timesteps); and one that had 9 observations locations clustered
in about the same location that repeated every 6 hours for 1.5 days (6
timesteps). I merged the two observation sequences into one using
``obs_sequence_tool`` and ran them through ``perfect_model_obs`` to derive the
observation values and create an ``obs_seq.out`` file to run through
``filter``.

.. note::

   Other models may have their ensemble means and spreads and inflation values
   in separate files. `See the table of possible filenames. <#FilenameTable>`__

.. code-block:: bash

   $ cd ${DARTROOT}/models/bgrid_solo/work
   $ ncdiff analysis.nc preassim.nc Innov.nc
   $ ncview preassim.nc &
   $ ncview Innov.nc &
   $ ncdump -v MemberMetadata preassim.nc
   netcdf preassim {
   dimensions:
           metadatalength = 64 ;
           member = 20 ;
           time = UNLIMITED ; // (24 currently)
           NMLlinelen = 129 ;
           NMLnlines = 303 ;
           StateVariable = 28200 ;
           TmpI = 60 ;
           TmpJ = 30 ;
           lev = 5 ;
           VelI = 60 ;
           VelJ = 29 ;
   variables:
           char MemberMetadata(member, metadatalength) ;
                   MemberMetadata:long_name = "Metadata for each copy/member" ;
           ...
           double ps(time, member, TmpJ, TmpI) ;
                   ps:long_name = "surface pressure" ;
                   ps:units = "Pa" ;
                   ps:units_long_name = "pascals" ;
           double t(time, member, lev, TmpJ, TmpI) ;
                   t:long_name = "temperature" ;
                   t:units = "degrees Kelvin" ;
           double u(time, member, lev, VelJ, VelI) ;
                   u:long_name = "zonal wind component" ;
                   u:units = "m/s" ;
           double v(time, member, lev, VelJ, VelI) ;
                   v:long_name = "meridional wind component" ;
                   v:units = "m/s" ;
           double ps_mean(time, TmpJ, TmpI) ;        The ensemble mean   is now a separate variable.
           double t_mean(time, lev, TmpJ, TmpI) ;    The ensemble spread is now a separate variable.
           double u_mean(time, lev, VelJ, VelI) ;    If I was using inflation, they would also be separate variables.
           double v_mean(time, lev, VelJ, VelI) ;
           double ps_sd(time, TmpJ, TmpI) ;
           double t_sd(time, lev, TmpJ, TmpI) ;
           double u_sd(time, lev, VelJ, VelI) ;
           double v_sd(time, lev, VelJ, VelI) ;

   data:
     MemberMetadata =
     "ensemble member      1 ",
     "ensemble member      2 ",
     "ensemble member      3 ",
     "ensemble member      4 ",
     "ensemble member      5 ",
     "ensemble member      6 ",
     "ensemble member      7 ",
     "ensemble member      8 ",
     "ensemble member      9 ",
     "ensemble member     10 ",
     "ensemble member     11 ",
     "ensemble member     12 ",
     "ensemble member     13 ",
     "ensemble member     14 ",
     "ensemble member     15 ",
     "ensemble member     16 ",
     "ensemble member     17 ",
     "ensemble member     18 ",
     "ensemble member     19 ",
     "ensemble member     20 " ;
   }

This is an exploration of the ``preassim.nc`` file. Note that I selected the
‘**t**’ field, turned the coastlines ‘off’ under the ‘Opts’ button, used the
‘Repl’ instead of ‘Bi-lin’ (to more faithfully represent the model resolution),
*navigated to copy 23 of 24 (in this case, the* :bolditalic:`inflation mean` *)* select
the **inflation mean variable of your choice** and advanced to the last
timestep. The image plot is pretty boring, but does indicate that the inflation
values are restricted to where I put the observations. Right-clicking on the
‘Range’ button automatically re-ranges the colorbar to the min/max of the
current data. Clicking on any location generates a time series figure.

This is an exploration of the ``Innov.nc`` file as created by *ncdiff*. Note
that the titles are somewhat misleading because they reflect information from
the first file given to *ncdiff*. This time I left the rendering as ‘Bi-lin’
(which obfuscates the model resolution), *navigated to copy 1 of 24 (in this
case, the* :bolditalic:`ensemble mean` *)* selected the **t_mean** variable and advanced
to the 6th timestep. Right-click on the ‘Range’ button to reset the colorbar.
The image plot confirms that the innovations are restricted to a local region.
Clicking on any location generates a time series.

This is fundamentally the same as the previous panel except that I have now
selected the ‘**u**’ **u_mean** variable. Despite the fact the observations were
only of ‘**t**’, the assimilation has generated (rightly so) increments to the
‘**u**’ state variable.
