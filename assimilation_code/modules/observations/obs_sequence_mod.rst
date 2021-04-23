MODULE obs_sequence_mod
=======================

Overview
--------

Provides interfaces to the observation type and observation sequence type. An observation contains everything there is
to know about an observation including all metadata contained in the observation definition and any number of copies of
data associated with the observation (for instance an actual observation, an ensemble of first guess values, etc). An
observation sequence is a time-ordered set of observations that is defined by a linked list so that observations can be
easily added or deleted. A number of commands to extract observations depending on the times at which they were taken
are provided. For now, the observations are only ordered by time, but the ability to add extra sort keys could be added.

These routines are commonly used in conversion programs which read observation data from various formats and create a
DART observation sequence in memory, and then write it out to a file. See the
`observations <../../../observations/obs_converters/README.md>`__ directory for examples of programs which create and
manipulate observations using this routines.

| 

Other modules used
------------------

::

   types_mod
   location_mod (depends on model_choice)
   obs_def_mod
   time_manager_mod
   utilities_mod
   obs_kind_mod

Public interfaces
-----------------

============================== ===================================
*use obs_sequence_mod, only :* obs_sequence_type
\                              init_obs_sequence
\                              interactive_obs_sequence
\                              get_num_copies
\                              get_num_qc
\                              get_num_obs
\                              get_max_num_obs
\                              get_copy_meta_data
\                              get_qc_meta_data
\                              get_next_obs
\                              get_prev_obs
\                              get_next_obs_from_key
\                              get_prev_obs_from_key
\                              insert_obs_in_seq
\                              delete_obs_from_seq
\                              set_copy_meta_data
\                              set_qc_meta_data
\                              get_first_obs
\                              get_last_obs
\                              add_copies
\                              add_qc
\                              write_obs_seq
\                              read_obs_seq
\                              append_obs_to_seq
\                              get_obs_from_key
\                              get_obs_time_range
\                              set_obs
\                              get_time_range_keys
\                              get_num_times
\                              static_init_obs_sequence
\                              destroy_obs_sequence
\                              read_obs_seq_header
\                              get_expected_obs
\                              delete_seq_head
\                              delete_seq_tail
\                              
\                              LINKS BELOW FOR OBS_TYPE INTERFACES
\                              
\                              obs_type
\                              init_obs
\                              destroy_obs
\                              get_obs_def
\                              set_obs_def
\                              get_obs_values
\                              set_obs_values
\                              replace_obs_values
\                              get_qc
\                              set_qc
\                              replace_qc
\                              write_obs
\                              read_obs
\                              interactive_obs
\                              copy_obs
\                              assignment(=)
============================== ===================================

| 

.. container:: type

   ::

      type obs_sequence_type
         private
         integer                       :: num_copies
         integer                       :: num_qc
         integer                       :: num_obs
         integer                       :: max_num_obs
         character(len=64), pointer    :: copy_meta_data(:)
         character(len=64), pointer    :: qc_meta_data(:)
         integer                       :: first_time
         integer                       :: last_time
         type(obs_type), pointer       :: obs(:)
      end type obs_sequence_type

.. container:: indent1

   The obs_sequence type represents a series of observations including multiple copies of data and quality control
   fields and complete metadata about the observations. The sequence is organized as an integer pointer linked list
   using a fixed array of storage for obs (type obs_type). Each observation points to the previous and next observation
   in time order (additional sort keys could be added if needed) and has a unique integer key (see obs_type below). The
   maximum number of observations in the sequence is represented in the type as max_num_obs, the current number of
   observations is in num_obs. The number of quality control (qc) fields per observation is num_qc and the number of
   data values associated with each observation is num_copies. Metadata for each copy of the data is in copy_meta_data
   and metadata for the qc fields is in qc_meta_data. The first and last pointers into the time linked list are in
   first_time and last_time. A capability to write and read an obs_sequence structure to disk is available. At present,
   the entire observation sequence is read in to core memory. An on-disk implementation may be necessary for very large
   observational datasets.

   ============== ===============================================================
   Component      Description
   ============== ===============================================================
   num_copies     Number of data values associated with each observation.
   num_qc         Number of qc fields associated with each observation.
   num_obs        Number of observations currently in sequence.
   max_num_obs    Upper bounds on number of observations in sequence.
   copy_meta_data Text describing each copy of data associated with observations.
   qc_meta_data   Text describing each quality control field.
   first_time     Location of first observation in sequence.
   last_time      Location of last observation in sequence.
   obs            Storage for all of the observations in the sequence.
   ============== ===============================================================

| 

.. container:: type

   ::

      type obs_type
         private
         integer            :: key
         type(obs_def_type) :: def
         real(r8), pointer  :: values(:)
         real(r8), pointer  :: qc(:)
         integer            :: prev_time
         integer            :: next_time
         integer            :: cov_group
      end type obs_type

.. container:: indent1

   Structure to represent everything known about a given observation and to help with storing the observation in the
   observation sequence structure (see above). The prev_time and next_time are integer pointers that allow a linked list
   sorted on time to be constructed. If needed, other sort keys could be introduced (for instance by time available?).
   Each observation in a sequence has a unique key and each observation has an obs_def_type that contains all the
   definition and metadata for the observation. A set of values is associated with the observation along with a set of
   qc fields. The cov_group is not yet implemented but will allow non-diagonal observation error covariances in a future
   release.

   ========= ====================================================================
   Component Description
   ========= ====================================================================
   key       Unique integer key when in an obs_sequence.
   def       The definition of the observation (see obs_def_mod).
   values    Values associated with the observation.
   qc        Quality control fields associated with the observation.
   prev_time When in an obs_sequence, points to previous time sorted observation.
   next_time When in an obs_sequence, points to next time sorted observation.
   cov_group Not currently implemented.
   ========= ====================================================================

| 

.. container:: routine

   *call init_obs_sequence(seq, num_copies, num_qc, expected_max_num_obs)*
   ::

      type(obs_sequence_type), intent(out) :: seq
      integer,                 intent(in)  :: num_copies
      integer,                 intent(in)  :: num_qc
      integer,                 intent(in)  :: expected_max_num_obs

.. container:: indent1

   Constructor to create a variable of obs_sequence_type. This routine must be called before using an obs_sequence_type.
   The number of copies of the data to be associated with each observation (for instance the observation from an
   instrument, an ensemble of prior guesses, etc.) and the number of quality control fields associated with each
   observation must be specified. Also, an estimated upper bound on the number of observations to be stored in the
   sequence is helpful in making creation of the sequence efficient.

   ======================== ============================================================================
   ``seq``                  The observation sequence being constructed
   ``num_copies``           Number of copies of data to be associated with each observation
   ``num_qc``               Number of quality control fields associated with each observation
   ``expected_max_num_obs`` An estimate of the largest number of observations the sequence might contain
   ======================== ============================================================================

| 

.. container:: routine

   *var = interactive_obs_sequence()*
   ::

      type(obs_sequence_type) :: interactive_obs_sequence

.. container:: indent1

   Uses input from standard in to create an observation sequence. Initialization of the sequence is handled by the
   function.

   ======= ===================================================
   ``var`` An observation sequence created from standard input
   ======= ===================================================

| 

.. container:: routine

   *var = get_num_copies(seq)*
   ::

      integer                             :: get_num_copies
      type(obs_sequence_type), intent(in) :: seq

.. container:: indent1

   Returns number of copies of data associated with each observation in an observation sequence.

   ======= =============================================================================
   ``var`` Returns number of copies of data associated with each observation in sequence
   ``seq`` An observation sequence
   ======= =============================================================================

| 

.. container:: routine

   *var = get_num_qc(seq)*
   ::

      integer                             :: get_num_qc
      type(obs_sequence_type), intent(in) :: seq

.. container:: indent1

   Returns number of quality control fields associated with each observation in an observation sequence.

   ======= =====================================================================================
   ``var`` Returns number of quality control fields associated with each observation in sequence
   ``seq`` An observation sequence
   ======= =====================================================================================

| 

.. container:: routine

   *var = get_num_obs(seq)*
   ::

      integer                             :: get_num_obs
      type(obs_sequence_type), intent(in) :: seq

.. container:: indent1

   Returns number of observations currently in an observation sequence.

   ======= ===================================================================
   ``var`` Returns number of observations currently in an observation sequence
   ``seq`` An observation sequence
   ======= ===================================================================

| 

.. container:: routine

   *var = get_max_num_obs(seq)*
   ::

      integer                             :: get_max_num_obs
      type(obs_sequence_type), intent(in) :: seq

.. container:: indent1

   Returns maximum number of observations an observation sequence can hold.

   ======= =======================================================================
   ``var`` Returns maximum number of observations an observation sequence can hold
   ``seq`` An observation sequence
   ======= =======================================================================

| 

.. container:: routine

   *var = get_copy_meta_data(seq, copy_num)*
   ::

      character(len=64)                   :: get_copy_meta_data
      type(obs_sequence_type), intent(in) :: seq
      integer,                 intent(in) :: copy_num

.. container:: indent1

   Returns metadata associated with a given copy of data in an observation sequence.

   ============ =======================================================================
   ``var``      Returns metadata associated with a copy of data in observation sequence
   ``seq``      An observation sequence
   ``copy_num`` Return metadata for this copy
   ============ =======================================================================

| 

.. container:: routine

   *var = get_qc_meta_data(seq,qc_num)*
   ::

      character(len=64)                   :: get_qc_meta_data
      type(obs_sequence_type), intent(in) :: seq
      integer,                 intent(in) :: qc_num

.. container:: indent1

   Returns metadata associated with a given copy of quality control fields associated with observations in an
   observation sequence.

   ========== ================================================
   ``var``    Returns metadata associated with a given qc copy
   ``seq``    An observation sequence
   ``qc_num`` Return metadata for this copy
   ========== ================================================

| 

.. container:: routine

   *call get_next_obs(seq, obs, next_obs, is_this_last)*
   ::

      type(obs_sequence_type), intent(in)  :: seq
      type(obs_type),          intent(in)  :: obs
      type(obs_type),          intent(out) :: next_obs
      logical,                 intent(out) :: is_this_last

.. container:: indent1

   Given an observation in a sequence, returns the next observation in the sequence. If there is no next observation,
   is_this_last is set to true.

   ================ ========================================
   ``seq``          An observation sequence
   ``obs``          Find the next observation after this one
   ``next_obs``     Return the next observation here
   ``is_this_last`` True if obs is the last obs in sequence
   ================ ========================================

| 

.. container:: routine

   *call get_prev_obs(seq, obs, prev_obs, is_this_first)*
   ::

      type(obs_sequence_type), intent(in)  :: seq
      type(obs_type),          intent(in)  :: obs
      type(obs_type),          intent(out) :: prev_obs
      logical,                 intent(out) :: is_this_first

.. container:: indent1

   Given an observation in a sequence, returns the previous observation in the sequence. If there is no previous
   observation, is_this_first is set to true.

   ================= =============================================
   ``seq``           An observation sequence
   ``obs``           Find the previous observation before this one
   ``prev_obs``      Return the previous observation here
   ``is_this_first`` True if obs is the first obs in sequence
   ================= =============================================

| 

.. container:: routine

   *call get_next_obs_from_key(seq, last_key_used, next_obs, is_this_last)*
   ::

      type(obs_sequence_type), intent(in)  :: seq
      integer,                 intent(in)  :: last_key_used
      type(obs_type),          intent(out) :: next_obs
      logical,                 intent(out) :: is_this_last

.. container:: indent1

   Given the last key used in a sequence, returns the next observation in the sequence. If there is no next observation,
   is_this_last is set to true.

   ================= ========================================
   ``seq``           An observation sequence
   ``last_key_used`` Find the next observation after this key
   ``next_obs``      Return the next observation here
   ``is_this_last``  True if obs is the last obs in sequence
   ================= ========================================

| 

.. container:: routine

   *call get_prev_obs_from_key(seq, last_key_used, prev_obs, is_this_first)*
   ::

      type(obs_sequence_type), intent(in)  :: seq
      integer,                 intent(in)  :: last_key_used
      type(obs_type),          intent(out) :: prev_obs
      logical,                 intent(out) :: is_this_first

.. container:: indent1

   Given the last key used in a sequence, returns the previous observation in the sequence. If there is no previous
   observation, is_this_first is set to true.

   ================= =============================================
   ``seq``           An observation sequence
   ``last_key_used`` Find the previous observation before this key
   ``prev_obs``      Return the previous observation here
   ``is_this_first`` True if obs is the first obs in sequence
   ================= =============================================

| 

.. container:: routine

   *call get_obs_from_key(seq, key, obs)*
   ::

      type(obs_sequence_type), intent(in)  :: seq
      integer,                 intent(in)  :: key
      type(obs_type),          intent(out) :: obs

.. container:: indent1

   Each entry in an observation sequence has a unique integer key. This subroutine returns the observation given an
   integer key.

   ======= ====================================
   ``seq`` An observation sequence
   ``key`` Return the observation with this key
   ``obs`` The returned observation
   ======= ====================================

| 

.. container:: routine

   *call insert_obs_in_seq(seq, obs [, prev_obs])*
   ::

      type(obs_sequence_type),  intent(inout) :: seq
      type(obs_type),           intent(inout) :: obs
      type(obs_type), optional, intent(in)    :: prev_obs

.. container:: indent1

   Inserts an observation in a sequence in appropriate time order. If the optional argument prev_obs is present, the new
   observation is inserted directly after the prev_obs. If an incorrect prev_obs is provided so that the sequence is no
   longer time ordered, bad things will happen.

   ========== =======================================================================
   ``seq``    An observation sequence
   ``obs``    An observation to be inserted in the sequence
   *prev_obs* If present, says the new observation belongs immediately after this one
   ========== =======================================================================

| 

.. container:: routine

   *call delete_obs_from_seq(seq, obs)*
   ::

      type(obs_sequence_type), intent(inout) :: seq
      type(obs_type),          intent(inout) :: obs

.. container:: indent1

   Given an observation and a sequence, removes the observation with the same key from the observation sequence.

   ======= ===============================================
   ``seq`` An observation sequence
   ``obs`` The observation to be deleted from the sequence
   ======= ===============================================

| 

.. container:: routine

   *call set_copy_meta_data(seq, copy_num, meta_data)*
   ::

      type(obs_sequence_type), intent(inout) :: seq
      integer,                 intent(in)    :: copy_num
      character(len=64),       intent(in)    :: meta_data

.. container:: indent1

   Sets the copy metadata for this copy of the observations in an observation sequence.

   ============= ==================================
   ``seq``       An observation sequence
   ``copy_num``  Set metadata for this copy of data
   ``meta_data`` The metadata
   ============= ==================================

| 

.. container:: routine

   *call set_qc_meta_data(seq, qc_num, meta_data)*
   ::

      type(obs_sequence_type), intent(inout) :: seq
      integer,                 intent(in)    :: qc_num
      character(len=64),       intent(in)    :: meta_data

.. container:: indent1

   Sets the quality control metadata for this copy of the qc in an observation sequence.

   ============= ===========================================
   ``seq``       An observation sequence
   ``qc_num``    Set metadata for this quality control field
   ``meta_data`` The metadata
   ============= ===========================================

| 

.. container:: routine

   *var = get_first_obs(seq, obs)*
   ::

      logical                              :: get_first_obs
      type(obs_sequence_type), intent(in)  :: seq
      type(obs_type),          intent(out) :: obs

.. container:: indent1

   Returns the first observation in a sequence. If there are no observations in the sequence, the function returns
   false, else true.

   ======= =============================================
   ``var`` Returns false if there are no obs in sequence
   ``seq`` An observation sequence
   ``obs`` The first observation in the sequence
   ======= =============================================

| 

.. container:: routine

   *var = get_last_obs(seq, obs)*
   ::

      logical                              :: get_last_obs
      type(obs_sequence_type), intent(in)  :: seq
      type(obs_type),          intent(out) :: obs

.. container:: indent1

   Returns the last observation in a sequence. If there are no observations in the sequence, the function returns false,
   else true.

   ======= =============================================
   ``var`` Returns false if there are no obs in sequence
   ``seq`` An observation sequence
   ``obs`` The last observation in the sequence
   ======= =============================================

| 

.. container:: routine

   *call add_copies(seq, num_to_add)*
   ::

      type(obs_sequence_type), intent(inout) :: seq
      integer,                 intent(in)    :: num_to_add

.. container:: indent1

   Increases the number of copies of data associated with each observation by num_to_add. The current implementation
   re-creates the entire observation sequence by deallocating and reallocating each entry with a larger size.

   ============== ===============================
   ``seq``        An observation sequence
   ``num_to_add`` Number of copies of data to add
   ============== ===============================

| 

.. container:: routine

   *call add_qc(seq, num_to_add)*
   ::

      type(obs_sequence_type), intent(inout) :: seq
      integer,                 intent(in)    :: num_to_add

.. container:: indent1

   Increases the number of quality control fields associated with each observation by num_to_add. The current
   implementation re-creates the entire observation sequence by deallocating and reallocating each entry with a larger
   size.

   ============== =======================================
   ``seq``        An observation sequence
   ``num_to_add`` Number of quality control fields to add
   ============== =======================================

| 

.. container:: routine

   *call read_obs_seq(file_name, add_copies, add_qc, add_obs, seq)*
   ::

      character(len=*),        intent(in)  :: file_name
      integer,                 intent(in)  :: add_copies
      integer,                 intent(in)  :: add_qc
      integer,                 intent(in)  :: add_obs
      type(obs_sequence_type), intent(out) :: seq

.. container:: indent1

   Read an observation sequence from ``file_name``. The sequence will have enough space for the number of observations
   in the file plus any additional space requested by the "add_xx" args. It is more efficient to allocate the additional
   space at create time rather than try to add it in later. The arguments can specify that the caller wants to add
   additional data copies associated with each observation, or to add additional quality control fields, or to add space
   for additional observations. The format of the file (``formatted`` vs. ``unformatted``) has been automatically
   detected since the I release. The obs_sequence file format with I and later releases has a header that associates
   observation type strings with an integer which was not present in previous versions. I format files are no longer
   supported.

   ============== ================================================================================
   ``file_name``  Read from this file
   ``add_copies`` Add this number of copies of data to the obs_sequence on file
   ``add_qc``     Add this number of qc fields to the obs_sequence on file
   ``add_obs``    Add space for this number of additional observations to the obs_sequence on file
   ``seq``        The observation sequence read in with any additional space
   ============== ================================================================================

| 

.. container:: routine

   *call write_obs_seq(seq, file_name)*
   ::

      type(obs_sequence_type), intent(in) :: seq
      character(len=*),        intent(in) :: file_name

.. container:: indent1

   Write the observation sequence to file file_name. The format is controlled by the namelist parameter
   write_binary_obs_sequence.

   ============= ===============================
   ``seq``       An observation sequence
   ``file_name`` Write the sequence to this file
   ============= ===============================

| 

.. container:: routine

   *call set_obs(seq,obs [, key_in])*
   ::

      type(obs_sequence_type), intent(inout) :: seq
      type(obs_type),          intent(in)    :: obs
      integer, optional,       intent(in)    :: key_in

.. container:: indent1

   Given an observation, copies this observation into the observation sequence using the key specified in the
   observation. If the optional key_in argument is present, the observation is instead copied into this element of the
   observation sequence (and the key is changed to be key_in).

   ======== ===========================================================
   ``seq``  An observation sequence
   ``obs``  Observation to be put in sequence
   *key_in* If present, the obs is copied into this key of the sequence
   ======== ===========================================================

| 

.. container:: routine

   *call append_obs_to_seq(seq, obs)*
   ::

      type(obs_sequence_type), intent(inout) :: seq
      type(obs_type),          intent(inout) :: obs

.. container:: indent1

   Append an observation to an observation sequence. An error results if the time of the observation is not equal to or
   later than the time of the last observation currently in the sequence.

   ======= =======================================
   ``seq`` An observation sequence
   ``obs`` Append this observation to the sequence
   ======= =======================================

| 

.. container:: routine

   *call get_obs_time_range(seq, time1, time2, key_bounds, num_keys, out_of_range [, obs])*
   ::

      type(obs_sequence_type),  intent(in)  :: seq
      type(time_type),          intent(in)  :: time1
      type(time_type),          intent(in)  :: time2
      integer, dimension(2),    intent(out) :: key_bounds
      integer,                  intent(out) :: num_keys
      logical,                  intent(out) :: out_of_range
      type(obs_type), optional, intent(in)  :: obs

.. container:: indent1

   Given a time range specified by a beginning and ending time, find the keys that bound all observations in this time
   range and the number of observations in the time range. The routine get_time_range_keys can then be used to get a
   list of all the keys in the range if desired. The logical out_of_range is returned as true if the beginning time of
   the time range is after the time of the latest observation in the sequence. The optional argument obs can increase
   the efficiency of the search through the sequence by indicating that all observations before obs are definitely at
   times before the start of the time range.

   ================ ====================================================================================
   ``seq``          An observation sequence
   ``time1``        Lower time bound
   ``time2``        Upper time bound
   ``key_bounds``   Lower and upper bounds on keys that are in the time range
   ``num_keys``     Number of keys in the time range
   ``out_of_range`` Returns true if the time range is entirely past the time of the last obs in sequence
   *obs*            If present, can start search for time range from this observation
   ================ ====================================================================================

| 

.. container:: routine

   *call get_time_range_keys(seq, key_bounds, num_keys, keys)*
   ::

      type(obs_sequence_type),      intent(in)  :: seq
      integer, dimension(2),        intent(in)  :: key_bounds
      integer,                      intent(in)  :: num_keys
      integer, dimension(num_keys), intent(out) :: keys

.. container:: indent1

   Given the keys of the observations at the start and end of a time range and the number of observations in the time
   range (these are returned by ``get_obs_time_range()``), return a list of the keys of all observations in the time
   range. Combining the two routines allows one to get a list of all observations in any time range by key. The ``keys``
   array must be at least ``num_keys`` long to hold the return values.

   ============== ==================================================
   ``seq``        An observation sequence
   ``key_bounds`` Keys of first and last observation in a time range
   ``num_keys``   Number of obs in the time range
   ``keys``       Output list of keys of all obs in the time range
   ============== ==================================================

| 

.. container:: routine

   *var = get_num_times(seq)*
   ::

      integer                             :: get_num_times
      type(obs_sequence_type), intent(in) :: seq

.. container:: indent1

   Returns the number of unique times associated with observations in an observation sequence.

   ======= =====================================================
   ``var`` Number of unique times for observations in a sequence
   ``seq`` An observation sequence
   ======= =====================================================

| 

.. container:: routine

   *var = get_num_key_range(seq, key1, key2)*
   ::

      integer                             :: get_num_key_range
      type(obs_sequence_type), intent(in) :: seq
      integer, optional,       intent(in) :: key1, key2

.. container:: indent1

   Returns the number of observations between the two given keys. The default key numbers are the first and last in the
   sequence file. This routine can be used to count the actual number of observations in a sequence and will be accurate
   even if the sequence has been trimmed with delete_seq_head() or delete_seq_tail().

   ======== ===========================================================================
   ``var``  Number of unique times for observations in a sequence
   ``seq``  An observation sequence
   ``key1`` The starting key number. Defaults to the first observation in the sequence.
   ``key2`` The ending key number. Defaults to the last observation in the sequence.
   ======== ===========================================================================

| 

.. container:: routine

   *call static_init_obs_sequence()*

.. container:: indent1

   Initializes the obs_sequence module and reads namelists. This MUST BE CALLED BEFORE USING ANY OTHER INTERFACES.

| 

.. container:: routine

   *call destroy_obs_sequence(seq)*
   ::

      type(obs_sequence_type), intent(inout) :: seq

.. container:: indent1

   Releases all allocated storage associated with an observation sequence.

   ======= =======================
   ``seq`` An observation sequence
   ======= =======================

| 

.. container:: routine

   *call read_obs_seq_header(file_name, num_copies, num_qc, num_obs, max_num_obs, file_id, read_format, pre_I_format [,
   close_the_file])*
   ::

      character(len=*),   intent(in)  :: file_name
      integer,            intent(out) :: num_copies
      integer,            intent(out) :: num_qc
      integer,            intent(out) :: num_obs
      integer,            intent(out) :: max_num_obs
      integer,            intent(out) :: file_id
      character(len=*),   intent(out) :: read_format
      logical,            intent(out) :: pre_I_format
      logical, optional,  intent(in)  :: close_the_file

.. container:: indent1

   Allows one to see the global metadata associated with an observation sequence that has been written to a file without
   reading the whole file.

   +------------------+--------------------------------------------------------------------------------------------------+
   | ``file_name``    | File contatining an obs_sequence                                                                 |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``num_copies``   | Number of copies of data associated with each observation                                        |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``num_qc``       | Number of quality control fields associated with each observation                                |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``num_obs``      | Number of observations in sequence                                                               |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``max_num_obs``  | Maximum number of observations sequence could hold                                               |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``file_id``      | File channel/descriptor returned from opening the file                                           |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``read_format``  | Either the string ``'unformatted'`` or ``'formatted'``                                           |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``pre_I_format`` | Returns .true. if the file was written before the observation type string/index number table was |
   |                  | added to the standard header starting with the I release.                                        |
   +------------------+--------------------------------------------------------------------------------------------------+
   | *close_the_file* | If specified and .TRUE. close the file after the header has been read. The default is to leave   |
   |                  | the file open.                                                                                   |
   +------------------+--------------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *call init_obs(obs, num_copies, num_qc)*
   ::

      type(obs_type), intent(out) :: obs
      integer,        intent(in)  :: num_copies
      integer,        intent(in)  :: num_qc

.. container:: indent1

   Initializes an obs_type variable. This allocates storage for the observation type and creates the appropriate
   obs_def_type and related structures. IT IS ESSENTIAL THAT OBS_TYPE VARIABLES BE INITIALIZED BEFORE USE.

   ============== ====================================================
   ``obs``        An obs_type data structure to be initialized
   ``num_copies`` Number of copies of data associated with observation
   ``num_qc``     Number of qc fields associated with observation
   ============== ====================================================

| 

.. container:: routine

   *call destroy_obs(obs)*
   ::

      type(obs_type), intent(inout) :: obs

.. container:: indent1

   Destroys an observation variable by releasing all associated storage.

   ======= =======================================
   ``obs`` An observation variable to be destroyed
   ======= =======================================

| 

.. container:: routine

   *call get_obs_def(obs, obs_def)*
   ::

      type(obs_type),     intent(in)  :: obs
      type(obs_def_type), intent(out) :: obs_def

.. container:: indent1

   Extracts the definition portion of an observation.

   =========== =========================================
   ``obs``     An observation
   ``obs_def`` The definition portion of the observation
   =========== =========================================

| 

.. container:: routine

   *call set_obs_def(obs, obs_def)*
   ::

      type(obs_type),     intent(out) :: obs
      type(obs_def_type), intent(in)  :: obs_def

.. container:: indent1

   Given an observation and an observation definition, insert the definition in the observation structure.

   =========== =======================================================
   ``obs``     An observation whose definition portion will be updated
   ``obs_def`` The observation definition that will be inserted in obs
   =========== =======================================================

| 

.. container:: routine

   *call get_obs_values(obs, values [, copy_indx])*
   ::

      type(obs_type),         intent(in)  :: obs
      real(r8), dimension(:), intent(out) :: values
      integer, optional,      intent(in)  :: copy_indx

.. container:: indent1

   Extract copies of the data from an observation. If *copy_indx* is present extract a single value indexed by
   *copy_indx* into ``values(1)``. *copy_indx* must be between 1 and ``num_copies``, inclusive. If *copy_indx* is not
   present extract all copies of data into the ``values`` array which must be ``num_copies`` long (See
   ``get_num_copies``.)

   =========== ===============================================================
   ``obs``     Observation from which to extract values
   ``values``  The values extracted
   *copy_indx* If present extract only this copy, otherwise extract all copies
   =========== ===============================================================

| 

.. container:: routine

   *call get_qc(obs, qc [, qc_indx])*
   ::

      type(obs_type),         intent(in)  :: obs
      real(r8), dimension(:), intent(out) :: qc
      integer, optional,      intent(in)  :: qc_indx

.. container:: indent1

   Extract quality control fields from an observation. If *qc_indx* is present extract a single field indexed by
   *qc_indx* into ``qc(1)``. *qc_indx* must be between 1 and ``num_qc``, inclusive. If *qc_indx* is not present extract
   all quality control fields into the ``qc`` array which must be ``num_qc`` long (See ``get_num_qc``.)

   ========= ===================================================================
   ``obs``   Observation from which to extract qc field(s)
   ``qc``    Extracted qc fields
   *qc_indx* If present extract only this field, otherwise extract all qc fields
   ========= ===================================================================

| 

.. container:: routine

   *call set_obs_values(obs, values [, copy_indx])*
   ::

      type(obs_type),         intent(out) :: obs
      real(r8), dimension(:), intent(in)  :: values
      integer, optional,      intent(in)  :: copy_indx

.. container:: indent1

   Set value(s) of data in this observation. If *copy_indx* is present set the single value indexed by *copy_indx* to
   ``values(1)``. *copy_indx* must be between 1 and ``num_copies``, inclusive. If *copy_indx* is not present set all
   copies of data from the ``values`` array which must be ``num_copies`` long (See ``get_num_copies``.)

   =========== ===============================================================
   ``obs``     Observation whose values are being set
   ``values``  Array of value(s) to be set
   *copy_indx* If present set only this copy of data, otherwise set all copies
   =========== ===============================================================

| 

.. container:: routine

   *call replace_obs_values(seq, key, values [, copy_indx])*
   ::

      type(obs_sequence_type), intent(inout) :: seq
      integer,                 intent(in)    :: key
      real(r8), dimension(:),  intent(in)    :: values
      integer, optional,       intent(in)    :: copy_indx

.. container:: indent1

   Set value(s) of data in the observation from a sequence with the given ``key``. If *copy_indx* is present set the
   single value indexed by *copy_indx* to ``values(1)``. *copy_indx* must be between 1 and ``num_copies``, inclusive. If
   *copy_indx* is not present set all copies of data from the ``values`` array which must be ``num_copies`` long (See
   ``get_num_copies``.)

   =========== ===============================================================
   ``seq``     Sequence which contains observation to update
   ``key``     Key to select which observation
   ``values``  Array of value(s) to be set
   *copy_indx* If present set only this copy of data, otherwise set all copies
   =========== ===============================================================

| 

.. container:: routine

   *call set_qc(obs, qc [, qc_indx])*
   ::

      type(obs_type),         intent(out) :: obs
      real(r8), dimension(:), intent(in)  :: qc
      integer, optional,      intent(in)  :: qc_indx

.. container:: indent1

   Sets the quality control fields in an observation. If *qc_indx* is present set a single field indexed by *qc_indx* to
   ``qc(1)``. *qc_indx* must be between 1 and ``num_qc``, inclusive. If *qc_indx* is not present set all quality control
   fields from the ``qc`` array which must be ``num_qc`` long (See ``get_num_qc``.)

   ========= =================================================================
   ``obs``   Observation having its qc fields set
   ``qc``    Input values of qc fields
   *qc_indx* If present update only this field, otherwise update all qc fields
   ========= =================================================================

| 

.. container:: routine

   *call replace_qc(seq, key, qc [, qc_indx])*
   ::

      type(obs_sequence_type), intent(inout) :: seq
      integer,                 intent(in)    :: key
      real(r8), dimension(:),  intent(in)    :: qc
      integer, optional,       intent(in)    :: qc_indx

.. container:: indent1

   Set value(s) of the quality control fields in the observation from a sequence with the given ``key``. If *qc_indx* is
   present set the single value indexed by *qc_indx* to ``qc(1)``. *qc_indx* must be between 1 and ``num_qc``,
   inclusive. If *qc_indx* is not present set all quality control fields from the ``qc`` array which must be ``num_qc``
   long (See ``get_num_qc``.)

   ========= ==================================================================
   ``seq``   Observation sequence containing observation to update
   ``key``   Key to select which observation
   ``qc``    Input values of qc fields
   *qc_indx* If present, only update single qc field, else update all qc fields
   ========= ==================================================================

| 

.. container:: routine

   *call write_obs(obs, file_id, num_copies, num_qc)*
   ::

      type(obs_type), intent(in) :: obs
      integer,        intent(in) :: file_id
      integer,        intent(in) :: num_copies
      integer,        intent(in) :: num_qc

.. container:: indent1

   Writes an observation and all its associated metadata to a disk file that has been opened with a format consistent
   with the namelist parameter ``write_binary_obs_sequence``.

   ============== =========================================================================
   ``obs``        Observation to be written to file
   ``file_id``    Channel open to file for writing
   ``num_copies`` The number of copies of data associated with the observation to be output
   ``num_qc``     The number of qc fields associated with the observation to be output
   ============== =========================================================================

| 

.. container:: routine

   *call read_obs(file_id, num_copies, add_copies, num_qc, add_qc, key, obs, read_format [, max_obs])*
   ::

      integer,            intent(in)    :: file_id
      integer,            intent(in)    :: num_copies
      integer,            intent(in)    :: add_copies
      integer,            intent(in)    :: num_qc
      integer,            intent(in)    :: add_qc
      integer,            intent(in)    :: key
      type(obs_type),     intent(inout) :: obs
      character(len=*),   intent(in)    :: read_format
      integer, optional,  intent(in)    :: max_obs

.. container:: indent1

   Reads an observation from an obs_sequence file. The number of copies of data and the number of qc values associated
   with each observation must be provided. If additional copies of data or additional qc fields are needed, arguments
   allow them to be added. WARNING: The key argument is no longer used and should be removed.

   +-----------------+---------------------------------------------------------------------------------------------------+
   | ``file_id``     | Channel open to file from which to read                                                           |
   +-----------------+---------------------------------------------------------------------------------------------------+
   | ``num_copies``  | Number of copies of data associated with observation in file                                      |
   +-----------------+---------------------------------------------------------------------------------------------------+
   | ``add_copies``  | Number of additional copies of observation to be added                                            |
   +-----------------+---------------------------------------------------------------------------------------------------+
   | ``num_qc``      | Number of qc fields associated with observation in file                                           |
   +-----------------+---------------------------------------------------------------------------------------------------+
   | ``add_qc``      | Number of additional qc fields to be added                                                        |
   +-----------------+---------------------------------------------------------------------------------------------------+
   | ``key``         | No longer used, should be deleted                                                                 |
   +-----------------+---------------------------------------------------------------------------------------------------+
   | ``obs``         | The observation being read in                                                                     |
   +-----------------+---------------------------------------------------------------------------------------------------+
   | ``read_format`` | Either the string ``'formatted'`` or ``'unformatted'``                                            |
   +-----------------+---------------------------------------------------------------------------------------------------+
   | *max_obs*       | If present, specifies the largest observation key number in the sequence. This is used only for   |
   |                 | additional error checks on the next and previous obs linked list values.                          |
   +-----------------+---------------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *call interactive_obs(num_copies, num_qc, obs, key)*
   ::

      integer,        intent(in)    :: num_copies
      integer,        intent(in)    :: num_qc
      type(obs_type), intent(inout) :: obs
      integer,        intent(in)    :: key

.. container:: indent1

   Use standard input to create an observation. The number of values, number of qc fields, and an observation
   type-specific key associated with the observation are input. (Note that the key here is not the same as the key in an
   observation sequence.)

   ============== =====================================================================================================
   ``num_copies`` Number of copies of data to be associated with observation
   ``num_qc``     Number of qc fields to be associated with observation
   ``obs``        Observation created via standard input
   ``key``        An observation type-specific key can be associated with each observation for use by the obs_def code.
   ============== =====================================================================================================

| 

.. container:: routine

   *call copy_obs(obs1, obs2)*
   ::

      type(obs_type), intent(out) :: obs1
      type(obs_type), intent(in)  :: obs2

.. container:: indent1

   Copies the observation type obs2 to obs1. If the sizes of obs fields are not compatible, the space in obs1 is
   deallocated and reallocated with the appropriate size. This is overloaded to assignment(=).

   ======== ===============================
   ``obs1`` Copy obs2 to here (destination)
   ``obs2`` Copy into obs1 (source)
   ======== ===============================

| 

.. container:: routine

   *call get_expected_obs_from_def_distrib_state(state_handle, ens_size, copy_indices, key, & obs_def, obs_kind_ind,
   state_time, isprior, assimilate_this_ob, evaluate_this_ob, expected_obs, & istatus)*
   ::

      type(ensemble_type), intent(in)  :: state_handle
      integer,             intent(in)  :: ens_size
      integer,             intent(in)  :: copy_indices(ens_size)
      integer,             intent(in)  :: key
      type(obs_def_type),  intent(in)  :: obs_def
      integer,             intent(in)  :: obs_kind_ind
      type(time_type),     intent(in)  :: state_time
      logical,             intent(in)  :: isprior
      integer,             intent(out) :: istatus(ens_size)
      logical,             intent(out) :: assimilate_this_ob, evaluate_this_ob
      real(r8),            intent(out) :: expected_obs(ens_size)

.. container:: indent1

   Used to compute the expected value of a set of observations in an observation sequence given a model state vector.
   Also returns a status variable that reports on problems taking forward operators. This version returns forward
   operator values for the entire ensemble in a single call.

   ====================== ============================================================================
   ``state_handle``       An observation sequence
   ``keys``               List of integer keys that specify observations in seq
   ``ens_index``          The ensemble number for this state vector
   ``state``              Model state vector
   ``state_time``         The time of the state data
   ``obs_vals``           Returned expected values of the observations
   ``istatus``            Integer error code for use in quality control (0 means no error)
   ``assimilate_this_ob`` Returns true if this observation type is being assimilated
   ``evaluate_this_ob``   Returns true if this observation type is being evaluated but not assimilated
   ====================== ============================================================================

| 

.. container:: routine

   *call delete_seq_head(first_time, seq, all_gone)*
   ::

      type(time_type),         intent(in)    :: first_time
      type(obs_sequence_type), intent(inout) :: seq
      logical,                 intent(out)   :: all_gone

.. container:: indent1

   Deletes all observations in the sequence with times before first_time. If no observations remain, return all_gone as
   .true. If no observations fall into the time window (e.g. all before first_time or empty sequence to begin with), no
   deletions are done and all_gone is simply returned as .true.

   ============== ==========================================================================================
   ``first_time`` Delete all observations with times before this
   ``seq``        An observation sequence
   ``all_gone``   Returns true if there are no valid observations remaining in the sequence after first_time
   ============== ==========================================================================================

| 

.. container:: routine

   *call delete_seq_tail(last_time, seq, all_gone)*
   ::

      type(time_type),         intent(in)    :: last_time
      type(obs_sequence_type), intent(inout) :: seq
      logical,                 intent(out)   :: all_gone

.. container:: indent1

   Deletes all observations in the sequence with times after last_time. If no observations remain, return all_gone as
   .true. If no observations fall into the time window (e.g. all after last_time or empty sequence to begin with), no
   deletions are done and all_gone is simply returned as .true.

   ============= ==========================================================================================
   ``last_time`` Delete all observations with times after this
   ``seq``       An observation sequence
   ``all_gone``  Returns true if there are no valid observations remaining in the sequence before last_time
   ============= ==========================================================================================

| 

Namelist
--------

This namelist is read from the file ``input.nml``. Namelists start with an ampersand '&' and terminate with a slash '/'.
Character strings that contain a '/' must be enclosed in quotes to prevent them from prematurely terminating the
namelist.

::

   &obs_sequence_nml
      write_binary_obs_sequence = .false.
      read_binary_file_format   = 'native'
     /

| 

.. container::

   +---------------------------+-------------------+--------------------------------------------------------------------+
   | Item                      | Type              | Description                                                        |
   +===========================+===================+====================================================================+
   | write_binary_obs_sequence | logical           | If true, write binary obs_sequence files. If false, write ascii    |
   |                           |                   | obs_sequence files.                                                |
   +---------------------------+-------------------+--------------------------------------------------------------------+
   | read_binary_file_format   | character(len=32) | The 'endian'ness of binary obs_sequence files. May be 'native'     |
   |                           |                   | (endianness matches hardware default), 'big-endian',               |
   |                           |                   | 'little-endian', and possibly 'cray'. Ignored if observation       |
   |                           |                   | sequence files are ASCII.                                          |
   +---------------------------+-------------------+--------------------------------------------------------------------+

| 

Files
-----

-  obs_sequence_mod.nml in input.nml
-  Files for reading and writing obs_sequences and obs specified in filter_nml.

References
----------

-  none

Private components
------------------

N/A
