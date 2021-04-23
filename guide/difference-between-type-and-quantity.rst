The difference between observation TYPE and QUANTITY
====================================================

Since DART is designed to assimilate data from any data source into any model,
the assimilation algorithms need a way to define how observational data sources
relate to model state variables.

DART does this by defining a single generic observation QUANTITY, such as zonal
wind, and mapping many specific observation TYPEs, corresponding to source 
observations, to the single QUANTITY.

For example, QuikSCAT and radiosondes are both capable of measuring zonal wind.
DART defines two observation TYPEs:

- ``QKSWND_U_WIND_COMPONENT`` for the QuikSCAT observations of zonal wind
- ``RADIOSONDE_U_WIND_COMPONENT`` for the radiosonde observations of zonal wind

and relates both of these TYPES to a single QUANTITY: QTY_U_WIND_COMPONENT.

Thus TYPE and QUANTITY have a many-to-one relationship. This distinction
enables you to assimilate or evaluate observation platforms independently 
of one another with a single observation sequence file; 
reducing the possibility of error.

The forward observation operators are implemented based on observation
QUANTITY. When requested, the model generates a QTY_U_WIND_COMPONENT, it
doesn't need to know that it will be compared to an observation from QuikSCAT
or one from a radiosonde.

.. tip::

   It is usually scientifically very interesting to be able to compare the
   assimilations one TYPE of observation verus another. An observation
   sequence file can have many types of observations. DART has the capability
   to assimilate (or evaluate) any combination of observation types without
   getting bogged down in dataset management. The same observation sequence can
   be used for experiments that include or exclude certain observation types.
   This procedure can ensure that you are actually performing the experiment
   that you think you are performing.
