Advice for models with multiple vertical coordinate options
===========================================================

DART vertical types for the 3D sphere locations type (threed_sphere)
--------------------------------------------------------------------

A location when using the
:doc:`../assimilation_code/location/threed_sphere/location_mod` 
location module consists of a Latitude (-90 to 90), a Longitude (0 to 360), and
a vertical value and type. The value is a real number. Possible types are:

- Height (in meters)
- Pressure (in Pascals)
- Model Level (index number)
- Scale Height (unitless)
- Surface (if value used, elevation of surface in meters)
- Undefined (entire vertical column)

If the model grid locations, all observation locations, and the choice of
localization coordinate are all using the identical vertical type then no
vertical conversion routines are needed. However, this is seldom the case.

Multiple vertical coordinate types
----------------------------------

Most Earth System models and observations use latitude and longitude for
horizontal coordinates or can generate them if needed (e.g. spectral models can
transform their state into Lat/Lon coords).  But often vertical coordinates
pose additional complications.

Some models use terrain-following vertical coordinates, or a mix of pressure
and terrain coordinates. Observation vertical locations are frequently reported 
in height or in pressure.

Additionally, if vertical localization is to be done in a different coordinate
than the model or observations (e.g. scale height), then conversion routines
are needed.

Vertical conversion routines typically take a DART ``location_type`` derived
type and a desired output vertical coordinate type as inputs, and either update
the location derived type or return a separate location type with the value
converted to the requested type.

The conversion code may require additional auxilliary arrays from the model in
order to convert the vertical coordinates accurately.

Varying vertical levels
-----------------------

If the computation of the vertical location depends on any of the fields in the
state (e.g. pressure), then different ensemble members may compute different
vertical locations.

Forward operators
~~~~~~~~~~~~~~~~~

During computation of expected values (Forward Operators), each ensemble member
should compute the most accurate value regardless of whether the location in
the model grid is consistent from member to member.

Localization
~~~~~~~~~~~~

During assimilation the distance between model state values and the observation
must be computed and only a single value can be returned, not an ensemble of
distances. If part of the state is needed to compute the vertical location the
ensemble mean is available to compute a single value which is representative of
the entire state.

Choice of when conversion is done 
---------------------------------

During assimilation there is a choice for when vertical
conversion is done: all at the start, before the sequential assimilation loop;
or on demand, as each observation is processed. 

The default behavior is to convert all observation locations before the sequential
assimilation loop, and convert all model state locations on demand. 
The default setting for vertical conversion in assim_tools_nml:

.. code-block:: text

    &assim_tools_nml
       convert_all_obs_verticals_first   = .true.,
       convert_all_state_verticals_first = .false.,
    /

When deciding to when to convert verticals, consider the following:

- How many state elements are updated, e.g. if you have a large :ref:`cutoff <localization>`
  and most of the state is being updated  during assimilation, then converting all state verticals first 
  may be more efficient.
- If you are :ref:`distributing the mean <data-distribution>` across MPI processes this will increase the cost of
  each conversion. You may not want to convert all state verticals first if the mean is distributed and 
  few state elements are updated.
- How many observations are being assimilated, e.g. if you have a large number of observations, then converting
  all observation verticals first may be more efficient as the conversion is done in parallel.

Before doing large scale experiments, we recommend profiling the assimilation to determine the best setting 
for your model and assimilation setup.