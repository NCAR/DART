module tracer_manager_mod
=========================

Overview
--------

Code to manage the simple addition of tracers to the FMS code. This code keeps track of the numbers and names of tracers
included in a tracer table.

.. container::

   This code is a grouping of calls which will allow the simple introduction of tracers into the FMS framework. It is
   designed to allow users of a variety of component models interact easily with the dynamical core of the model.
   In calling the tracer manager routines the user must provide a parameter identifying the model that the user is
   working with. This parameter is defined within field_manager as MODEL_X where X is one of [ATMOS, OCEAN, LAND, ICE].
   In many of these calls the argument list includes model and tracer_index. These are the parameter corresponding to
   the component model and the tracer_index N is the Nth tracer within the component model. Therefore a call with
   MODEL_ATMOS and 5 is different from a call with MODEL_OCEAN and 5.

| 

Other modules used
------------------

.. container::

   ::

      mpp_mod
             mpp_io_mod
                fms_mod
      field_manager_mod

Public interface
----------------

.. container::

   ::

      use tracer_manager_mod [, only:  tracer_manager_init,
                                       register_tracers,
                                       get_number_tracers,
                                       get_tracer_indices,
                                       get_tracer_index,
                                       assign_tracer_field,
                                       tracer_manager_end,
                                       get_tracer_field,
                                       get_tracer_tlevels,
                                       get_tracer_tendency,
                                       get_tracer_names,
                                       get_family_name,
                                       check_if_prognostic,
                                       find_family_members,
                                       add_members_to_family,
                                       split_family_into_members,
                                       set_tracer_profile,
                                       query_method,
                                       query_combined,
                                       set_tracer_atts ]

   tracer_manager_init:
      Routine to initialize the tracer manager
   register_tracers:
      A routine to register the tracers included in a component model.
   get_number_tracers:
      A routine to return the number of tracers included in a component model.
   get_tracer_indices:
      Routine to return the component model tracer indices as defined within the tracer manager.
   get_tracer_index:
      Function which returns the number assigned to the tracer name.
   assign_tracer_field:
      Routine to point the appropriate field within the tracer_type to the appropriate field within the component model.
   tracer_manager_end:
      Routine to write to the log file that the tracer manager is ending.
   get_tracer_field:
      A function to retrieve the present timestep data.
   get_tracer_tlevels:
      A function to retrieve the three time levels data.
   get_tracer_tendency:
      A function to retrieve the tendency data.
   get_tracer_names:
      Routine to find the names associated with a tracer number.
   get_family_name:
      Routine to return the family name for tracer n.
   check_if_prognostic:
      Function to see if a tracer is prognostic or diagnostic.
   find_family_members:
      Subroutine to find which tracers are members of family family_name.
   add_members_to_family:
      Routine to sum up the members of a family of tracers so that they may be advected and diffused as one tracer.
   split_family_into_members:
      Subroutine that sets the present value of the member of a tracer family according to the fraction of the family
      that it was in the previous step.
   set_tracer_profile:
      Subroutine to set the tracer field to the wanted profile.
   query_method:
      A function to query the "methods" associated with each tracer.
   query_combined:
      A function to query whether families of tracers have been combined already.
   set_tracer_atts:
      A subroutine to allow the user set the tracer longname and units from the tracer initialization routine.

| 

Public data
-----------

.. container::

   None.

Public routines
---------------

a. .. rubric:: Tracer_manager_init
      :name: tracer_manager_init

   ::

      call tracer_manager_init 

   **DESCRIPTION**
      This routine writes the version and tagname to the logfile and sets the module initialization flag.

b. .. rubric:: Register_tracers
      :name: register_tracers

   ::

      call register_tracers (model, num_tracers,num_prog,num_diag,num_family)

   **DESCRIPTION**
      This routine returns the total number of valid tracers, the number of prognostic and diagnostic tracers and the
      number of families of tracers.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``model``                                                 | A parameter to identify which model is being used.        |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``num_tracers``                                           | The total number of valid tracers within the component    |
      |                                                           | model.                                                    |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``num_prog``                                              | The number of prognostic tracers within the component     |
      |                                                           | model.                                                    |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``num_diag``                                              | The number of diagnostic tracers within the component     |
      |                                                           | model.                                                    |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``num_family``                                            | The number of family tracers within the component model.  |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

c. .. rubric:: Get_number_tracers
      :name: get_number_tracers

   ::

      call get_number_tracers (model, num_tracers,num_prog,num_diag,num_family)

   **DESCRIPTION**
      This routine returns the total number of valid tracers, the number of prognostic and diagnostic tracers and the
      number of families of tracers.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``model``                                                 | A parameter to identify which model is being used.        |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``num_tracers``                                           | The total number of valid tracers within the component    |
      |                                                           | model.                                                    |
      |                                                           | [integer, optional]                                       |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``num_prog``                                              | The number of prognostic tracers within the component     |
      |                                                           | model.                                                    |
      |                                                           | [integer, optional]                                       |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``num_diag``                                              | The number of diagnostic tracers within the component     |
      |                                                           | model.                                                    |
      |                                                           | [integer, optional]                                       |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``num_family``                                            | The number of family tracers within the component model.  |
      |                                                           | [integer, optional]                                       |
      +-----------------------------------------------------------+-----------------------------------------------------------+

d. .. rubric:: Get_tracer_indices
      :name: get_tracer_indices

   ::

      call get_tracer_indices (model, ind, prog_ind, diag_ind, fam_ind)

   **DESCRIPTION**
      If several models are being used or redundant tracers have been written to the tracer_table, then the indices in
      the component model and the tracer manager may not have a one to one correspondence. Therefore the component model
      needs to know what index to pass to calls to tracer_manager routines in order that the correct tracer information
      be accessed.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``model``                                                 | A parameter to identify which model is being used.        |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``ind``                                                   | An array containing the tracer manager defined indices    |
      |                                                           | for all the tracers within the component model.           |
      |                                                           | [integer, optional, dimension(:)]                         |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``prog_ind``                                              | An array containing the tracer manager defined indices    |
      |                                                           | for the prognostic tracers within the component model.    |
      |                                                           | [integer, optional, dimension(:)]                         |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``diag_ind``                                              | An array containing the tracer manager defined indices    |
      |                                                           | for the diagnostic tracers within the component model.    |
      |                                                           | [integer, optional, dimension(:)]                         |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``fam_ind``                                               | An array containing the tracer manager defined indices    |
      |                                                           | for the family tracers within the component model.        |
      |                                                           | [integer, optional, dimension(:)]                         |
      +-----------------------------------------------------------+-----------------------------------------------------------+

e. .. rubric:: Get_tracer_index
      :name: get_tracer_index

   ::

      value= get_tracer_index (model, name, indices, verbose)

   **DESCRIPTION**
      This is a function which returns the index, as implied within the component model.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``model``                                                 | A parameter to identify which model is being used.        |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``name``                                                  | The name of the tracer (as assigned in the field table).  |
      |                                                           | [character]                                               |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``indices``                                               | An array of the component model indices. This array can   |
      |                                                           | be found by calling get_tracer_indices.                   |
      |                                                           | [integer, optional, dimension(:)]                         |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``verbose``                                               | A flag to allow the message saying that a tracer with     |
      |                                                           | this name has not been found. This should only be used    |
      |                                                           | for debugging purposes.                                   |
      |                                                           | [logical, optional]                                       |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``get_tracer_index``                                      | The index of the tracer named "name". If indices is       |
      |                                                           | passed then the result is the array index which           |
      |                                                           | corresponds to tracer named "name".                       |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

f. .. rubric:: Assign_tracer_field
      :name: assign_tracer_field

   ::

      call assign_tracer_field (model,index, data, data_tlevels, tendency)

   **DESCRIPTION**
      The generality provided here is that one can point the three dimensional tracer field at either a two time level
      scheme [data and tendency] or a three time level scheme [data_tlevels]. The tracer manager points the appropriate
      tracer_type field at the data supplied from the component model.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``model``                                                 | A parameter representing the component model in use.      |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``index``                                                 | The tracer number that you wish to assign a tracer field  |
      |                                                           | for.                                                      |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``data``                                                  | The 3D field that is associated with the present time     |
      |                                                           | step in the component model.                              |
      |                                                           | [real, target, optional, dimension(:,:,:)]                |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``tendency``                                              | The 3D field that is associated with the tendency time    |
      |                                                           | step in the component model.                              |
      |                                                           | [real, target, optional, dimension(:,:,:)]                |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``data_tlevels``                                          | The 4D field that is associated with the tracer field in  |
      |                                                           | the component model.                                      |
      |                                                           | [real, target, optional, dimension(:,:,:,:)]              |
      +-----------------------------------------------------------+-----------------------------------------------------------+

g. .. rubric:: Tracer_manager_end
      :name: tracer_manager_end

   ::

      call tracer_manager_end 

   **DESCRIPTION**
      Routine to write to the log file that the tracer manager is ending.

h. .. rubric:: Get_tracer_field
      :name: get_tracer_field

   ::

      array= get_tracer_field (model, tracer_index)

   **DESCRIPTION**
      Function to point to the 3D field associated with a tracer.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``model``                                                 | A parameter representing the component model in use.      |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``tracer_index``                                          | The tracer number within the component model.             |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``data``                                                  | The tracer field is returned in this array.               |
      |                                                           | [real, pointer, dimension(:,:,:)]                         |
      +-----------------------------------------------------------+-----------------------------------------------------------+

i. .. rubric:: Get_tracer_tlevels
      :name: get_tracer_tlevels

   ::

      array= get_tracer_tlevels (model, tracer_index)

   **DESCRIPTION**
      Function to point to the 4D field associated with a tracer.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``model``                                                 | A parameter representing the component model in use.      |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``tracer_index``                                          | The tracer number within the component model.             |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``data``                                                  | The tracer field is returned in this array.               |
      |                                                           | [real, pointer, dimension(:,:,:,:)]                       |
      +-----------------------------------------------------------+-----------------------------------------------------------+

j. .. rubric:: Get_tracer_tendency
      :name: get_tracer_tendency

   ::

      array= get_tracer_tendency (model, tracer_index)

   **DESCRIPTION**
      Function to point to the 3D field associated with a tracer.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``model``                                                 | A parameter representing the component model in use.      |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``tracer_index``                                          | The tracer number within the component model.             |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``data``                                                  | The tracer tendency field is returned in this array.      |
      |                                                           | [real, pointer, dimension(:,:,:)]                         |
      +-----------------------------------------------------------+-----------------------------------------------------------+

k. .. rubric:: Get_tracer_names
      :name: get_tracer_names

   ::

      call get_tracer_names (model,n,name,longname, units)

   **DESCRIPTION**
      This routine can return the name, long name and units associated with a tracer.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``model``                                                 | A parameter representing the component model in use.      |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``n``                                                     | Tracer number.                                            |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``name``                                                  | Field name associated with tracer number.                 |
      |                                                           | [character]                                               |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``longname``                                              | The long name associated with tracer number.              |
      |                                                           | [character, optional]                                     |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``units``                                                 | The units associated with tracer number.                  |
      |                                                           | [character, optional]                                     |
      +-----------------------------------------------------------+-----------------------------------------------------------+

l. .. rubric:: Get_family_name
      :name: get_family_name

   ::

      call get_family_name (model,n,name)

   **DESCRIPTION**
      You may wish to use this routine to retrieve the name of the family that a tracer belongs to.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``model``                                                 | A parameter representing the component model in use.      |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``n``                                                     | Tracer number that you want the family name for.          |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``name``                                                  | The family name.                                          |
      |                                                           | [character]                                               |
      +-----------------------------------------------------------+-----------------------------------------------------------+

m. .. rubric:: Check_if_prognostic
      :name: check_if_prognostic

   ::

      logical = check_if_prognostic (model, n)

   **DESCRIPTION**
      All tracers are assumed to be prognostic when read in from the field_table However a tracer can be changed to a
      diagnostic tracer by adding the line "tracer_type","diagnostic" to the tracer description in field_table.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``model``                                                 | A parameter representing the component model in use.      |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``n``                                                     | Tracer number that you want the family name for.          |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``check_if_prognostic``                                   | A logical flag set TRUE if the tracer is prognostic.      |
      |                                                           | [logical]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

n. .. rubric:: Find_family_members
      :name: find_family_members

   ::

      call find_family_members (model, family_name,is_family_member)

   **DESCRIPTION**
      Subroutine to find which tracers are members of family family_name. This will return a logical array where the
      array positions corresponding to the tracer numbers for family members are set .TRUE.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``model``                                                 | A parameter representing the component model in use.      |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``family_name``                                           | The family name of the members one is seeking.            |
      |                                                           | [character]                                               |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``is_family_member``                                      | A logical array where the tracer number is used as the    |
      |                                                           | index to signify which tracer is part of the family. i.e. |
      |                                                           | If tracers 1, 3, and 7 are part of the same family then   |
      |                                                           | is_family_member(1), is_family_member(3), and             |
      |                                                           | is_family_member(7) are set TRUE.                         |
      |                                                           | [logical, dimension(:)]                                   |
      +-----------------------------------------------------------+-----------------------------------------------------------+

o. .. rubric:: Add_members_to_family
      :name: add_members_to_family

   ::

      call add_members_to_family (model,family_name, cur, prev, next)

   **DESCRIPTION**
      Routine to sum up the members of a family of tracers so that they may be advected and diffused as one tracer. This
      should only be used in conjunction with split_family_into_members and should be placed before the advection scheme
      is called.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``model``                                                 | A parameter representing the component model in use.      |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``n``                                                     | Tracer number.                                            |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``cur``                                                   | Array index for the current time step. This is only of    |
      |                                                           | use with a three timestep model.                          |
      |                                                           | [integer, optional]                                       |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``prev``                                                  | Array index for the previous time step. This is only of   |
      |                                                           | use with a three timestep model.                          |
      |                                                           | [integer, optional]                                       |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``next``                                                  | Array index for the next time step. This is only of use   |
      |                                                           | with a three timestep model.                              |
      |                                                           | [integer, optional]                                       |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **NOTE**
      This should be used with extreme caution. Unless the family member distributions are similar to each other
      spatially, advection as one tracer and subsequent splitting will result in a different result to advecting each
      tracer separately. The user should understand the possible repercussions of this before using it.

p. .. rubric:: Split_family_into_members
      :name: split_family_into_members

   ::

      call split_family_into_members (model,family_name,cur,prev,next)

   **DESCRIPTION**
      Subroutine that sets the present value of the member of a tracer family according to the fraction of the family
      that it was in the previous step.
      This splits the transported family into the constituent members. This should only be used in conjunction with
      *add_members_to_family* and should be placed after the advection scheme is called.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``model``                                                 | A parameter representing the component model in use.      |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``family_name``                                           | The name of the family of tracers that you would like to  |
      |                                                           | split up.                                                 |
      |                                                           | [character]                                               |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``cur``                                                   | Array index for the current time step. This is only of    |
      |                                                           | use with a three timestep model.                          |
      |                                                           | [integer, optional]                                       |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``prev``                                                  | Array index for the previous time step. This is only of   |
      |                                                           | use with a three timestep model.                          |
      |                                                           | [integer, optional]                                       |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``next``                                                  | Array index for the next time step. This is only of use   |
      |                                                           | with a three timestep model.                              |
      |                                                           | [integer, optional]                                       |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **NOTE**
      This should be used with extreme caution. Unless the family member distributions are similar to each other
      spatially, advection as one tracer and subsequent splitting will result in a different result to advecting each
      tracer separately. The user should understand the possible repercussions of this before using it.

q. .. rubric:: Set_tracer_profile
      :name: set_tracer_profile

   ::

      call set_tracer_profile (model, n, surf_value, multiplier)

   **DESCRIPTION**
      If the profile type is 'fixed' then the tracer field values are set equal to the surface value. If the profile
      type is 'profile' then the top/bottom of model and surface values are read and an exponential profile is
      calculated, with the profile being dependent on the number of levels in the component model. This should be called
      from the part of the dynamical core where tracer restarts are called in the event that a tracer restart file does
      not exist.
      This can be activated by adding a method to the field_table e.g. "profile_type","fixed","surface_value = 1e-12"
      would return values of surf_value = 1e-12 and a multiplier of 1.0 One can use these to initialize the entire field
      with a value of 1e-12.
      "profile_type","profile","surface_value = 1e-12, top_value = 1e-15" In a 15 layer model this would return values
      of surf_value = 1e-12 and multiplier = 0.6309573 i.e 1e-15 = 1e-12*(0.6309573^15) In this case the model should be
      MODEL_ATMOS as you have a "top" value.
      If you wish to initialize the ocean model, one can use bottom_value instead of top_value.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``model``                                                 | A parameter representing the component model in use.      |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``n``                                                     | Tracer number.                                            |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``surf_value``                                            | The surface value that will be initialized for the tracer |
      |                                                           | [real]                                                    |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``multiplier``                                            | The vertical multiplier for the tracer Level(k-1) =       |
      |                                                           | multiplier \* Level(k)                                    |
      |                                                           | [real]                                                    |
      +-----------------------------------------------------------+-----------------------------------------------------------+

r. .. rubric:: Query_method
      :name: query_method

   ::

      logical = query_method (method_type, model, n, name, control)

   **DESCRIPTION**
      A function to query the "methods" associated with each tracer. The "methods" are the parameters of the component
      model that can be adjusted by user by placing formatted strings, associated with a particular tracer, within the
      field table. These methods can control the advection, wet deposition, dry deposition or initial profile of the
      tracer in question. Any parametrization can use this function as long as a routine for parsing the name and
      control strings are provided by that routine.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``method_type``                                           | The method that is being requested.                       |
      |                                                           | [character]                                               |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``model``                                                 | A parameter representing the component model in use.      |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``n``                                                     | Tracer number that you want the family name for.          |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``name``                                                  | A string containing the modified name to be used with     |
      |                                                           | method_type. i.e. "2nd_order" might be the default for    |
      |                                                           | advection. One could use "4th_order" here to modify that  |
      |                                                           | behaviour.                                                |
      |                                                           | [character]                                               |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``control``                                               | A string containing the modified parameters that are      |
      |                                                           | associated with the method_type and name.                 |
      |                                                           | [character, optional]                                     |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``query_method``                                          | A flag to show whether method_type exists with regard to  |
      |                                                           | tracer n. If method_type is not present then one must     |
      |                                                           | have default values.                                      |
      |                                                           | [logical]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **NOTE**
      | At present the tracer manager module allows the initialization of a tracer profile if a restart does not exist
        for that tracer. Options for this routine are as follows
      | Tracer profile setup

      ::


         ==================================================================
         |method_type |method_name |method_control |
         ==================================================================
         |profile_type |fixed |surface_value = X | |profile_type |profile
         |surface_value = X, top_value = Y |(atmosphere) |profile_type
         |profile |surface_value = X, bottom_value = Y |(ocean)
         ==================================================================

      | 

s. .. rubric:: Query_combined
      :name: query_combined

   ::

      logical = query_combined (model, index)

   **DESCRIPTION**
      A function to query whether families of tracers have been combined already. This function should only be used in
      conjunction with add_members_to_family and split_family_into_members.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``model``                                                 | A parameter representing the component model in use.      |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``index``                                                 | Tracer number.                                            |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``query_combined``                                        | A flag to show whether the tracer family has been         |
      |                                                           | combined.                                                 |
      |                                                           | [logical]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

t. .. rubric:: Set_tracer_atts
      :name: set_tracer_atts

   ::

      call set_tracer_atts (model, name, longname, units)

   **DESCRIPTION**
      A function to allow the user set the tracer longname and units from the tracer initialization routine. It seems
      sensible that the user who is coding the tracer code will know what units they are working in and it is probably
      safer to set the value in the tracer code rather than in the field table.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``model``                                                 | A parameter representing the component model in use.      |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``name``                                                  | Tracer name.                                              |
      |                                                           | [character]                                               |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``longname``                                              | A string describing the longname of the tracer for output |
      |                                                           | to NetCDF files                                           |
      |                                                           | [character, optional]                                     |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``units``                                                 | A string describing the units of the tracer for output to |
      |                                                           | NetCDF files                                              |
      |                                                           | [character, optional]                                     |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``set_tracer_atts``                                       | A flag to show that                                       |
      |                                                           | [character, optional]                                     |
      +-----------------------------------------------------------+-----------------------------------------------------------+

Data sets
---------

.. container::

   None.

Error messages
--------------

.. container::

   **FATAL in register_tracers**
      invalid model type
      The index for the model type is invalid.
   **NOTE in register_tracers**
      No tracers are available to be registered.
      No tracers are available to be registered. This means that the field table does not exist or is empty.
   **FATAL in register_tracers**
      MAX_TRACER_FIELDS exceeded
      The maximum number of tracer fields has been exceeded.
   **NOTE in register_tracers**
      There is only 1 tracer for tracer family X. Making an orphan.
      A tracer has been given a family name but that family has only this member. Therefore it should be an orphan.
   **FATL in register_tracers**
      MAX_TRACER_FIELDS needs to be increased
      The number of tracer fields has exceeded the maximum allowed. The parameter MAX_TRACER_FIELDS needs to be
      increased.
   **FATAL in get_number_tracers**
      Model number is invalid.
      The index of the component model is invalid.
   **Fatal in get_tracer_indices**
      index array size too small in get_tracer_indices
      The global index array is too small and cannot contain all the tracer numbers.
   **FATAL in get_tracer_indices**
      family array size too small in get_tracer_indices
      The family index array is too small and cannot contain all the tracer numbers.
   **FATAL in get_tracer_indices**
      prognostic array size too small in get_tracer_indices
      The prognostic index array is too small and cannot contain all the tracer numbers.
   **FATAL in get_tracer_indices**
      diagnostic array size too small in get_tracer_indices
      The diagnostic index array is too small and cannot contain all the tracer numbers.
   **NOTE in get_tracer_index**
      tracer with this name not found: X
   **FATAL in assign_tracer_field**
      invalid index
      The index that has been passed to this routine is invalid.
   **FATAL in assign_tracer_field**
      At least one of data, data_tlevels or tendency must be passed in here.
      At least one of data, data_tlevels or tendency must be passed to assign_tracer_field Otherwise there is not much
      point in calling this routine.
   **FATAL in get_tracer_field**
      invalid index
      The index that has been passed to this routine is invalid. Check the index that is being passed corresponds to a
      valid tracer name.
   **FATAL in get_tracer_field**
      invalid index
      The index that has been passed to this routine is invalid. Check the index that is being passed corresponds to a
      valid tracer name.
   **FATAL in get_tracer_field**
      tracer field array not allocated
      The tracer array has not been allocated. This means that a call to assign_tracer_field is absent in the code.
   **FATAL in get_tracer_tlevels**
      invalid index
      The index that has been passed to this routine is invalid. Check the index that is being passed corresponds to a
      valid tracer name.
   **FATAL in get_tracer_tlevels**
      invalid index
      The index that has been passed to this routine is invalid. Check the index that is being passed corresponds to a
      valid tracer name.
   **FATAL in get_tracer_tlevels**
      tracer field array not allocated
      The tracer array has not been allocated. This means that a call to assign_tracer_field is absent in the code.
   **FATAL in get_tracer_tendency**
      invalid index
      The index that has been passed to this routine is invalid. Check the index that is being passed corresponds to a
      valid tracer name.
   **FATAL in get_tracer_tendency**
      invalid index
      The index that has been passed to this routine is invalid. Check the index that is being passed corresponds to a
      valid tracer name.
   **FATAL in get_tracer_tendency**
      tracer tendency field array not allocated
      The tracer array has not been allocated. This means that a call to assign_tracer_field is absent in the code.
