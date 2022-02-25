Lorenz 84
=========

Overview
--------

This model was described in Lorenz (1984). [1]_ In Lorenz 84, DART advances the
model, gets the model state and metadata describing this state, find states
variables that are close to a given location, and does spatial interpolation
for model state variables. The distinctive part of the model interfaces is the
`namelist`_.

The system of equations is:

.. math::

   \frac{dx}{dt} = -y^2-z^2-ax+aF
   \frac{dy}{dt} = xy-bxz-y+G
   \frac{dz}{dt} = bxy+xz-z

and, within DART, the model parameters have default values of:

.. math::

   a=\frac{1}{4}, b=4, F=8, G=\frac{5}{4}

that can be altered by editing the ``&model_nml`` `namelist`_ in the
``input.nml`` file.

The Lorenz 84 model has a ``work/workshop_setup.csh`` script that compiles and runs 
an example.  This example is referenced specifically in Section 7 of the 
:doc:`DART tutorial <../../theory/readme>`
and is intended to provide insight into model/assimilation behavior.
The example **may or may not** result in good (*or even decent!*) results!

The Lorenz 84 model may be used instead of the Lorenz 63 model in many sections
of the Tutorial. It has a more complex attractor, is not as periodic as Lorenz 63
and may be more challenging for certain filter variants.

Namelist
--------

The ``&model_nml`` namelist is read from the ``input.nml`` file. Namelists
start with an ampersand ``&`` and terminate with a slash ``/``. Character
strings that contain a ``/`` must be enclosed in quotes to prevent them from
prematurely terminating the namelist.

.. code-block:: fortran

  &model_nml
    a      = 0.25,
    b      = 4.00,
    f      = 8.00,
    g      = 1.25,
    deltat = 0.01,
    time_step_days    = 0,
    time_step_seconds = 3600
  /

Description of each namelist entry
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

+-------------------+----------+-------------------------------------+
| Item              | Type     | Description                         |
+===================+==========+=====================================+
| a                 | real(r8) | Model parameter.                    |
+-------------------+----------+-------------------------------------+
| b                 | real(r8) | Model parameter.                    |
+-------------------+----------+-------------------------------------+
| f                 | real(r8) | Model parameter.                    |
+-------------------+----------+-------------------------------------+
| g                 | real(r8) | Model parameter.                    |
+-------------------+----------+-------------------------------------+
| deltat            | real(r8) | Non-dimensional timestep. This is   |
|                   |          | mapped to the dimensional timestep  |
|                   |          | specified by time_step_days and     |
|                   |          | time_step_seconds.                  |
+-------------------+----------+-------------------------------------+
| time_step_days    | integer  | Number of days for dimensional      |
|                   |          | timestep, mapped to deltat.         |
+-------------------+----------+-------------------------------------+
| time_step_seconds | integer  | Number of seconds for dimensional   |
|                   |          | timestep, mapped to deltat.         |
+-------------------+----------+-------------------------------------+

References
~~~~~~~~~~

.. [1] Lorenz, Edward N., 1984: Irregularity: A Fundamental Property of the
       Atmosphere. *Tellus*, **36A**, 98-110, 
       `doi:10.1111/j.1600-0870.1984.tb00230.x
       <https://doi.org/10.1111/j.1600-0870.1984.tb00230.x>`__
