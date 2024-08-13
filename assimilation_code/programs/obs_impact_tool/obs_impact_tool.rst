PROGRAM ``obs_impact_tool``
===========================

Overview
--------

The standard DART algorithms work by calculating increments for an observation and then determining corresponding
increments for each variable in the state due to that observation. This is done by computing a sample regression 
coefficient using the prior ensemble distributions of a state variable and the observation. The increments for each member 
of the ensemble are multiplied by this coefficient and then added to the corresponding prior ensemble member for the variable.

However, in many cases it is necessary to limit the influence of an observation on a variable; this is known as localization.
DART provides a way to specify a localization, known as cutoff, based on the horizontal and vertical distance between the observation 
and the state variable.

In some situations, you may want additional localization based on the type of observation and the state quanity. 
``obs_impact_tool`` allows you to create a table that filter reads during runtime to localize the impact of certain types of 
observations on specific state vector quantities. You can define sets of observation types and state vector quantities, and 
specify localization for the impact of those observation types on the state vector quantities.

For example, you can create a subset of observations related to tracer concentration for various tracers, and a subset of 
dynamic state variables like temperatures and wind components. Typically, it is common practice to set this localization value 
to 0 to prevent tracer observations from affecting dynamic state quantities. However, ``obs_impact_tool`` allows you to specify values
between 0 and 1.


#. Build ``obs_sequence_tool`` by adding ``obs_impact_tool`` to the list of serial_programs in the quickbuild.sh script for the model you are using.
   Run ./quickbuild.sh to build all the DART programs.
#. Create an input file for ``obs_sequence_tool`` to define the impacts of observations. In the examples on this page, the input file
   is called `cross_correlations.txt`.  
   The format of the input file can be any combination of the following types of sections:

   .. code:: bash

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

   The following is an example of an input file to prevent chemistry species from impacting the meterological variables in the model state, and vice versa:

   .. code:: bash

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


#. Run ``obs_impact_tool`` using your `cross_correlations.txt` as input. ``obs_impact_tool`` will create an output file,
   named `control_impact_runtime.txt` in this example.

   .. code:: text
   
      &obs_impact_tool_nml
        input_filename          = 'cross_correlations.txt'
        output_filename         = 'control_impact_runtime.txt'
        /
   

#. Set the following namelist options in :ref:`&assim_tools_nml<assim_tools>` to use `control_impact_runtime.txt` in filter. 
   Filter will apply your selected observation impacts during assimilation.

   .. code:: text
   
      &assim_tools_nml
        adjust_obs_impact               = .true.
        obs_impact_filename             = 'control_impact_runtime.txt'
        /


obs_impact_tool Namelist
------------------------

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
