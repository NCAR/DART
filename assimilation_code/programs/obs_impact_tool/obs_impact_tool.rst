PROGRAM ``obs_impact_tool``
===========================

Overview
--------

The standard DART algorithms compute increments for an observation and then compute corresponding increments for each
model state variable due to that observation. To do this, DART computes a sample regression coefficient using the prior
ensemble distributions of a state variable and the observation. The increments for each member of the observation are
multiplied by this regression coefficient and then added to the corresponding prior ensemble member for the state
variable. However, in many cases, it is appropriate to reduce the impact of an observation on a state variable; this is
called localization. The standard DART algorithms allow users to specify a localization that is a function of the
horizontal (and optionally vertical) distance between the observation and the state variable. The localization is a
value between 0 and 1 and multiplies the regression coefficient when updating state ensemble members.

Sometimes, it may be desirable to do an additional localization that is a function of the 
type of observation and the
state vector quantity. This program allows users to construct a table that is read by 
filter at run-time to localize the
impact of sets of observation types on sets of state vector quantities. Users can create 
named sets of observation types
and sets of state vector quantities and specify a localization for the impact of the 
specified observation types on the state vector quantities.

An example would be to create a subset of observations of tracer concentration for a variety of tracers, and a subset of
dynamic state variable quantities like temperatures and wind components. It has been common to set this localization
value to 0 so that tracer observations have no impact on dynamic state quantities, however, the tool allows values
between 0 and 1 to be specified.

This tool allows related collections of observation types and state vector quantities to be named and then express the
relationship of the named groups to each other in a concise way. It can also define relationships by exceptions.

All the listed observation types and state vector quantities must be known by the system.
If they are not, look at the
&preprocess_nml :: input_items namelist which specifies which *obs_def_xxx_mod.f90* files 
are included, which is where observation types are defined.
Quantities for different regimes (atmosphere, ocean, land, etc.) are defined in
``assimilation_code/modules/observations/xxx_quantities_mod.f90`` and explained in
:doc:`../../modules/observations/obs_kind_mod`

Format of the input file can be any combination of these types of sections:

.. container::

   ::



      # hash mark starts a comment.

      # the GROUP keyword starts a group and must be followed
      # by a name.  All types or quantities listed before the END
      # line becomes members of this group.

      # GROUPs cannot contain nested groups.

      GROUP groupname1
       QTY_xxx  QTY_xxx  QTY_xxx
       QTY_xxx                          # comments can be here
      END GROUP

      GROUP groupname2
       QTY_xxx  
       QTY_xxx  
       QTY_xxx
       QTY_xxx
      END GROUP

      # GROUPs can also be defined by specifying ALL, ALLQTYS,
      # or ALLTYPES and then EXCEPT and listing the types or
      # quantities which should be removed from this group.
      # ALL EXCEPT must be the first line in a group, and all
      # subsequent items are removed from the list.
      # The items listed after EXCEPT can include the names
      # of other groups.

      GROUP groupnameM
      ALL EXCEPT QTY_xxx QTY_xxx
      QTY_xxx
      END GROUP

      GROUP groupnameN
      ALL EXCEPT groupnameY
      END GROUP


      # once any groups have been defined, a single instance
      # of the IMPACT table is specified by listing a TYPE,
      # QTY, or group in column 1, then a QTY or GROUP
      # in column 2 (the second name cannot be a specific type).
      # column 3 must be 0.0 or 1.0.  subsequent entries
      # that overlap previous entries have precedence
      # (last entry wins).

      IMPACT
       QTY_xxx     QTY_xxx      0.0
       QTY_xxx     groupname1   0.0
       groupname1  QTY_xxx      0.0
       groupname1  groupname1   0.0
      END IMPACT

Namelist interface ``&obs_impact_tool_nml`` must be read from file ``input.nml``.

Namelist
--------

This namelist is read from the file ``input.nml``. Namelists start with an ampersand '&' and terminate with a slash '/'.
Character strings that contain a '/' must be enclosed in quotes to prevent them from prematurely terminating the
namelist.

::

   &obs_impact_tool_nml
     input_filename          = 'cross_correlations.txt'
     output_filename         = 'control_impact_runtime.txt'
     debug                   = .false.
     /

| 

.. container::

   +-----------------+--------------------+-----------------------------------------------------------------------------+
   | Item            | Type               | Description                                                                 |
   +=================+====================+=============================================================================+
   | input_filename  | character(len=512) | Name of an ascii text file which describes how the interaction of           |
   |                 |                    | observations to state vector values and observations to other observations  |
   |                 |                    | should be controlled. See the Overview section for details about the format |
   |                 |                    | of the input file entries.                                                  |
   +-----------------+--------------------+-----------------------------------------------------------------------------+
   | output_filename | character(len=512) | Name of an ascii text file which created by this tool. It can be read at    |
   |                 |                    | filter run time to control the impact of observations on state vector items |
   |                 |                    | and other observation values. The format of this file is set by this tool   |
   |                 |                    | and should not be modified by hand. Rerun this tool to recreate the file.   |
   +-----------------+--------------------+-----------------------------------------------------------------------------+
   | debug           | logical            | If true print out debugging info.                                           |
   +-----------------+--------------------+-----------------------------------------------------------------------------+

| 

Examples
--------

To prevent chemistry species from impacting the meterological variables in the model state, and vice versa:

.. container::

   ::

      GROUP chem
       QTY_CO QTY_NO QTY_C2H4
      END GROUP

      GROUP met
       ALLQTYS EXCEPT chem
      END GROUP

      IMPACT
       chem   met    0.0
       met    chem   0.0
      END IMPACT

Modules used
------------

::

   types_mod
   utilities_mod
   parse_args_mod

Files
-----

-  two text files, one input and one output.
-  obs_impact_tool.nml

References
----------

-  none
