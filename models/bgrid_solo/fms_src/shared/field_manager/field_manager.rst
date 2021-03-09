module field_manager_mod
========================

Overview
--------

The field manager reads entries from a field table and stores this information along with the type of field it belongs
to. This allows the component models to query the field manager to see if non-default methods of operation are desired.
In essence the field table is a powerful type of namelist. Default values can be provided for all the fields through a
namelist, individual fields can be modified through the field table however.

.. container::

   The field table consists of entries in the following format.
   The first line of an entry should consist of three quoted strings. The first quoted string will tell the field
   manager what type of field it is. At present the supported types of fields are "tracer" for tracers, "xland_mix" for
   cross-land mixing, and, "checkerboard" for checkerboard null mode.
   The second quoted string will tell the field manager which model the field is being applied to. The supported types
   at present are "atmos_mod" for the atmosphere model, "ocean_mod" for the ocean model, "land_mod" for the land model,
   and, "ice_mod" for the ice model.
   The third quoted string should be a unique name that can be used as a query.
   The second and following lines of each entry are called methods in this context. Methods can be developed within any
   module and these modules can query the field manager to find any methods that are supplied in the field table.
   These lines can consist of two or three quoted strings. The first string will be an identifier that the querying
   module will ask for. The second string will be a name that the querying module can use to set up values for the
   module. The third string, if present, can supply parameters to the calling module that can be parsed and used to
   further modify values.
   An entry is ended with a backslash (/) as the final character in a row.
   Comments can be inserted in the field table by having a # as the first character in the line.
   An example of a field table entry could be
   ::

      "tracer","atmos_mod","sphum"/
      "tracer","atmos_mod","sf6"
      "longname","sulf_hex"
      "advection_scheme_horiz","2nd_order"
      "Profile_type","Fixed","surface_value = 0.0E+00"/

   In this example we have two field entries.
   The first is a simple declaration of a tracer called "sphum".
   The second is for a tracer called "sf6". Methods that are being applied to this tracer include initiating the long
   name of the tracer to be "sulf_hex", changing the horizontal advection scheme to be second order, and initiating the
   tracer with a profile with fixed values, in this example all zero.

| 

Other modules used
------------------

.. container::

   ::

         mpp_mod
      mpp_io_mod
         fms_mod

Public interface
----------------

.. container::

   ::

      use field_manager_mod [, only:  field_manager_init,
                                      field_manager_end,
                                      find_field_index,
                                      get_field_info,
                                      get_field_method,
                                      get_field_methods,
                                      parse ]

   field_manager_init:
      Routine to initialize the field manager.
   field_manager_end:
      Destructor for field manager.
   find_field_index:
      Function to return the index of the field.
   get_field_info:
      This routine allows access to field information given an index.
   get_field_method:
      A routine to get a specified method.
   get_field_methods:
      A routine to obtain all the methods associated with a field.
   parse:
      A function to parse an integer or an array of integers, a real or an array of reals, a string or an array of
      strings.

| 

Public data
-----------

.. container::

   ===================== ================== ======= ===== =============================
   Name                  Type               Value   Units Description
   ===================== ================== ======= ===== =============================
   NUM_MODELS            integer, parameter 5       ---   Number of models.
   module_is_initialized logical            .false. ---   Field manager is initialized.
   MODEL_ATMOS           integer, parameter 1       ---   Atmospheric model.
   MODEL_OCEAN           integer, parameter 2       ---   Ocean model.
   MODEL_LAND            integer, parameter 3       ---   Land model.
   MODEL_ICE             integer, parameter 4       ---   Ice model.
   ===================== ================== ======= ===== =============================

Public routines
---------------

a. .. rubric:: Field_manager_init
      :name: field_manager_init

   ::

      call field_manager_init (nfields, table_name)

   **DESCRIPTION**
      This routine reads from a file containing formatted strings. These formatted strings contain information on which
      schemes are needed within various modules. The field manager does not initialize any of those schemes however. It
      simply holds the information and is queried by the appropriate module.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``table_name``                                            | The name of the field table. The default name is          |
      |                                                           | field_table.                                              |
      |                                                           | [character, optional, dimension(len=128)]                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``nfields``                                               | The number of fields.                                     |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

b. .. rubric:: Field_manager_end
      :name: field_manager_end

   ::

      call field_manager_end 

   **DESCRIPTION**
      This subroutine writes to the logfile that the user is exiting field_manager and changes the initialized flag to
      false.

c. .. rubric:: Find_field_index
      :name: find_field_index

   ::

      value= find_field_index ( model, field_name )

   **DESCRIPTION**
      This function when passed a model number and a field name will return the index of the field within the field
      manager. This index can be used to access other information from the field manager.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``model``                                                 | The number indicating which model is used.                |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

d. .. rubric:: Get_field_info
      :name: get_field_info

   ::

      call get_field_info ( n,fld_type,fld_name,model,num_methods )

   **DESCRIPTION**
      When passed an index, this routine will return the type of field, the name of the field, the model which the field
      is associated and the number of methods associated with the field.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``n``                                                     | The field index.                                          |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``fld_type``                                              | The field type.                                           |
      |                                                           | [character, dimension(*)]                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``fld_name``                                              | The name of the field.                                    |
      |                                                           | [character, dimension(*)]                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``model``                                                 | The number indicating which model is used.                |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``num_methods``                                           | The number of methods.                                    |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

e. .. rubric:: Get_field_method
      :name: get_field_method

   ::

      call get_field_method ( n,m,method )

   **DESCRIPTION**
      This routine, when passed a field index and a method index will return the method text associated with the
      field(n) method(m).
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``n``                                                     | The field index.                                          |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``m``                                                     | The method index.                                         |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

f. .. rubric:: Get_field_methods
      :name: get_field_methods

   ::

      call get_field_methods ( n,methods )

   **DESCRIPTION**
      When passed a field index, this routine will return the text associated with all the methods attached to the
      field.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``n``                                                     | The field index.                                          |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

g. .. rubric:: Parse
      :name: parse

   ::

      number = parse (text, label, value)

   **DESCRIPTION**
      Parse is an integer function that decodes values from a text string. The text string has the form: "label=list"
      where "label" is an arbitrary user defined label describing the values being decoded, and "list" is a list of one
      or more values separated by commas. The values may be integer, real, or character. Parse returns the number of
      values decoded.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``text``                                                  | The text string from which the values will be parsed.     |
      |                                                           | [character(len=*)]                                        |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``label``                                                 | A label which describes the values being decoded.         |
      |                                                           | [character(len=*)]                                        |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``value``                                                 | The value or values that have been decoded.               |
      |                                                           | [integer, real, character(len=*)]                         |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``parse``                                                 | The number of values that have been decoded. This allows  |
      |                                                           | a user to define a large array and fill it partially with |
      |                                                           | values from a list. This should be the size of the value  |
      |                                                           | array.                                                    |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

Data sets
---------

.. container::

   None.

Error messages
--------------

.. container::

   **NOTE in field_manager_init**
      No field table available, so no fields are being registered.
      The field table does not exist.
   **FATAL in field_manager_init**
      max fields exceeded
      Maximum number of fields for this module has been exceeded.
   **FATAL in field_manager_init**
      Too many fields in tracer entry.
      There are more that 3 fields in the tracer entry. This is probably due to separating the parameters entry into
      multiple strings. The entry should look like
      "Type","Name","Control1=XXX,Control2=YYY"
      and not like
      "Type","Name","Control1=XXX","Control2=YYY"
   **FATAL in field_manager_init**
      Maximum number of methods for field exceeded
      Maximum number of methods allowed for entries in the field table has been exceeded.
   **NOTE in field_manager_init**
      field with identical name and model name duplicate found, skipping
      The name of the field and the model name are identical. Skipping that field.
   **FATAL in field_manager_init**
      error reading field table
      There is an error in reading the field table.
   **FATAL in get_field_info**
      invalid field index
      The field index is invalid because it is less than 1 or greater than the number of fields.
   **FATAL in get_field_method**
      invalid field index
      The field index is invalid because it is less than 1 or greater than the number of fields.
   **FATAL in get_field_method**
      invalid method index
      The method index is invalid because it is less than 1 or greater than the number of methods.
   **FATAL in get_field_methods**
      invalid field index
      The field index is invalid because it is less than 1 or greater than the number of fields.
   **FATAL in get_field_methods**
      method array too small
      The method array is smaller than the number of methods.

.. container::

   top
