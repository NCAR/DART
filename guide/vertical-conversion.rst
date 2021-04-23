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

During the assimilation phase of filter there are two options for when vertical
conversion is done: all at the start, or on demand.  If the observations to be
assimilated are expected to impact all or almost all of the state, doing all
vertical conversion at the start is more efficient. If the observations are
expected to impact only a small percentage of the state variables then doing it
on demand is more efficient.

The options here are namelist selectable at runtime and the impact on total
runtime can be easily measured and compared.
