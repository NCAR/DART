CHAMP
=====

This is a modification of a standard ``text_to_obs`` converter that comes that
comes with DART.

This observation converter reads CHAMP and GRACE density files, as described in
Sutton (2011) [1]_ and outputs ``obs_seq`` files that can be assimilated using
DART.

.. warning::

   If an ``obs_seq.out`` file already exists, this converter  automatically adds
   new observations to that file without deleting it. This is done to allow the
   wrapper script (work/convert.sh) to process sequentially numbered
   Density_*.ascii files, as documented in the comments in convert.sh. If you
   don't want this behavior, comment out lines 129-132 in text_to_obs.f90 and
   rebuild.

Namelist
--------

Please inspect the ``text_to_obs_nml`` namelist in ``work/input.nml`` to ensure
the input and output filenames are specified properly.

.. note::

   The `work/Density_3deg_02_335.ascii` file is truncated to 2 datapoints
   merely to demonstrate the format. It isn't meant to be used for real
   experiments.

Author
------

Thank you to Alexey Morozov for contributing this observation converter.

References
----------

.. [1]  Sutton, Erik K., 2011: Accelerometer-Derived Atmospheric Density
        from the CHAMP and GRACE Satellites.

