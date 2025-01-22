MODULE ensemble_manager_mod
===========================

Overview
--------

Manages storage and a number of operations for multiple copies of a vector. The most obvious use is to manage ensembles
of model state vectors. In this case, the number of copies stored for each state vector element is the ensemble size
plus one or more additional copies like the mean, variance, associated inflation values, etc. The ensemble_manager
provides routines to compute the mean and variance of a subset of the copies, to track the time associated with the
copies, and to write and read restart files. Most importantly, it provides a capability to do transposes between two
storage representations of an ensemble. In one representation, each process stores all copies of a subset of the state
variables while in the other, each process stores all of the state variables for a subset of copies. The ensemble
manager is also used to manage ensembles of observation priors and quality control and ensembles of forward observation
operator error status.

The ensemble manager interacts strongly with the multiple process capability of the Message Passing Interface (MPI)
libraries. It is used to partition the data so each MPI process stores only a subset of the copies and variables,
dividing the data as evenly as possible across the processes. At no time during the execution does any one process have
to store the entire dataset for all ensemble members (unless running in serial mode without MPI, or if running with 1
MPI task).

The ensemble manager is set of general purpose data management routines. For run-time efficiency, the derived type
information is not marked private which means other modules can directly manipulate the data arrays. However it means
much care must be taken to access the most recently updated representation of the data, either the copies or variables
arrays.

A set of sanity check routines have been added to track the last modified version of the data: the copies array or the
vars array. Before directly reading or writing these arrays call one of the 'prepare' routines to indicate what kind of
data access you are about to make. If the most recently updated data is not as expected an error message will occur.
After the direct access if the following operations detect that the data they are operating on is not the most recently
updated they will print an error message. Routines inside the ensemble manager that alter the copies or vars will set
the state automatically so these routines are only necessary to call if you are directly accessing the copies or vars
arrays from outside the ensemble manager.

Namelist
--------

This namelist is read from the file ``input.nml``. Namelists start with an ampersand '&' and terminate with a slash '/'.
Character strings that contain a '/' must be enclosed in quotes to prevent them from prematurely terminating the
namelist.

::

   &ensemble_manager_nml
      layout                      = 1
      tasks_per_node              = 1
      communication_configuration = 1
      debug                       = .false.
     /

| 

.. container::

   +-----------------------------+---------+----------------------------------------------------------------------------+
   | Item                        | Type    | Description                                                                |
   +=============================+=========+============================================================================+
   | layout                      | integer | Determines the logical process (PE) layout across MPI tasks. 1 is PE = MPI |
   |                             |         | task. 2 is a round-robin layout around the nodes. Layout 2 results in a    |
   |                             |         | more even usage of memory across nodes. This may allow you to run with a   |
   |                             |         | larger state vector without hitting the memory limit of the node. It may   |
   |                             |         | give a slight (5%) increase in performance, but this is machine dependent. |
   |                             |         | It has no effect on serial runs.                                           |
   +-----------------------------+---------+----------------------------------------------------------------------------+
   | tasks_per_node              | integer | The number of MPI tasks per hardware node is generally fixed when a batch  |
   |                             |         | job is submitted. This namelist item tells the ensemble manager what the   |
   |                             |         | user selected at that time. Once a program is running the code has no      |
   |                             |         | control to change how MPI tasks are assigned to physical CPUs. This number |
   |                             |         | is used only if layout = 2, and it allows the code spread high-memory-use  |
   |                             |         | PEs to different hardware nodes by assigning them in a round-robin order.  |
   |                             |         | The job will still run if this number does not match the real              |
   |                             |         | "tasks_per_node" at the hardware level, but it may run out of memory if    |
   |                             |         | the mismatch causes multiple high-memory-use tasks to be run on the same   |
   |                             |         | node.                                                                      |
   +-----------------------------+---------+----------------------------------------------------------------------------+
   | communication_configuration | integer | For most users, the default value of 1 is the best choice. However there   |
   |                             |         | are multiple strategies for the internal MPI communication patterns (see   |
   |                             |         | \*Note below). Values from 1 to 4 select different options; try the        |
   |                             |         | various options to see if one might be faster than the others.             |
   +-----------------------------+---------+----------------------------------------------------------------------------+
   | debug                       | logical | If true print debugging information.                                       |
   +-----------------------------+---------+----------------------------------------------------------------------------+

| 

| *\*Note about MPI communication flags:*
| The communication_configuration flags select various combinations of the internal settings for use_copy2var_send_loop
  and use_var2copy_rec_loop. These flags change the order of the MPI send and MPI receives in the the routines
  all_copies_to_all_vars and all_vars_to_all_copies. The figures below show the data transferred between tasks for an 80
  member ensemble. The left figure is using 96 tasks, the right figure is using 512 tasks. As the number of tasks
  increases, the 'all to all' data transfer becomes a 'some to all, all to some' transfer and the order of MPI send and
  MPI receives becomes increasingly important. The default values give a performance advantage as the number of tasks
  becomes much greater than the the ensemble size. However, for small numbers of tasks, i.e. less than the ensemble
  size, changing the default values may improve performance.

.. container::

   ======================= =========================
   |communication pattern| |communication pattern 2|
   ======================= =========================

| 

Other modules used
------------------

::

   types_mod
   utilities_mod
   assim_model_mod
   time_manager_mod
   random_seq_mod
   mpi_utilities_mod
   sort_mod

Public interfaces
-----------------

================================== ===========================
*use ensemble_manager_mod, only :* init_ensemble_manager
\                                  read_ensemble_restart
\                                  write_ensemble_restart
\                                  get_copy
\                                  put_copy
\                                  broadcast_copy
\                                  set_ensemble_time
\                                  get_ensemble_time
\                                  end_ensemble_manager
\                                  duplicate_ens
\                                  get_my_num_copies
\                                  get_my_copies
\                                  get_my_num_vars
\                                  get_my_vars
\                                  get_copy_owner_index
\                                  get_var_owner_index
\                                  all_vars_to_all_copies
\                                  all_copies_to_all_vars
\                                  compute_copy_mean
\                                  compute_copy_mean_sd
\                                  compute_copy_mean_var
\                                  print_ens_handle
\                                  map_pe_to_task
\                                  map_task_to_pe
================================== ===========================

A note about documentation style. Optional arguments are enclosed in brackets *[like this]*.

| 

.. container:: type

   ::

      type ensemble_type
         !DIRECT ACCESS INTO STORAGE IS ALLOWED; BE CAREFUL
         integer :: num_copies
         integer :: num_vars
         integer :: my_num_copies
         integer :: my_num_vars
         integer, pointer :: my_copies(:)
         integer, pointer :: my_vars(:)
         ! Storage in next line is to be used when each PE has all copies of subset of vars
         real(r8), pointer :: copies(:, :)  ! Dimensioned (num_copies, my_num_vars)
         ! Storage on next line is used when each PE has subset of copies of all vars
         real(r8), pointer :: vars(:, :)    ! Dimensioned (num_vars, my_num_copies)
         ! Time is only related to var complete
         type(time_type), pointer :: time(:)
         integer :: distribution_type
         integer :: id_num
         integer, allocatable :: task_to_pe_list(:) ! List of tasks
         integer, allocatable :: pe_to_task_list(:) ! List of tasks
         ! Flexible my_pe, layout_type which allows different task layouts for different ensemble handles
         integer :: my_pe
         integer :: layout_type
      end type ensemble_type

.. container:: indent1

   Provides a handle for an ensemble that manages copies of a vector. For efficiency, the type internals are not private
   and direct access to the storage arrays is used throughout DART.

   +-------------------+-------------------------------------------------------------------------------------------------+
   | Component         | Description                                                                                     |
   +===================+=================================================================================================+
   | num_copies        | Global number of copies of the vector.                                                          |
   +-------------------+-------------------------------------------------------------------------------------------------+
   | num_vars          | Global number of elements (variables) in the vector.                                            |
   +-------------------+-------------------------------------------------------------------------------------------------+
   | my_num_copies     | Number of copies stored by this process.                                                        |
   +-------------------+-------------------------------------------------------------------------------------------------+
   | my_num_vars       | Number of variables stored by this process.                                                     |
   +-------------------+-------------------------------------------------------------------------------------------------+
   | my_copies         | Dimensioned to size my_num_copies. Contains a list of the global indices of copies stored by    |
   |                   | this process.                                                                                   |
   +-------------------+-------------------------------------------------------------------------------------------------+
   | my_vars           | Dimensioned to size my_num_vars. Contains a list of the global indices of variables stored by   |
   |                   | this process.                                                                                   |
   +-------------------+-------------------------------------------------------------------------------------------------+
   | copies            | Dimensioned (num_copies, my_num_vars). Storage for all copies of variables stored by this       |
   |                   | process.                                                                                        |
   +-------------------+-------------------------------------------------------------------------------------------------+
   | vars              | Dimensioned (num_vars, my_num_copies). Storage for all variables of copies stored by this       |
   |                   | process.                                                                                        |
   +-------------------+-------------------------------------------------------------------------------------------------+
   | time              | Dimensioned my_num_copies. A time_type that stores time associated with a given copy of the     |
   |                   | vector.                                                                                         |
   +-------------------+-------------------------------------------------------------------------------------------------+
   | distribution_type | Does nothing at present. Can be used for future releases to control the layout of different     |
   |                   | copies and variables in storage.                                                                |
   +-------------------+-------------------------------------------------------------------------------------------------+
   | valid             | Flag to track whether the copies array has the most recently updated data, the vars array is    |
   |                   | most recently modified, or if both the arrays have identical data, like after a transpose.      |
   +-------------------+-------------------------------------------------------------------------------------------------+
   | id_num            | Internal number unique to each ensemble handle, used for debugging purposes.                    |
   +-------------------+-------------------------------------------------------------------------------------------------+
   | task_to_pe_list   | Mapping from MPI task number to logical Processing Element (PE) number. Enables different       |
   |                   | assignment of MPI tasks to PEs. If the number of MPI tasks is larger than the number of copies  |
   |                   | of the vector, when the ensemble is var complete then the first N MPI tasks have allocated      |
   |                   | 'vars' arrays and the remaining ones do not. Assigning the MPI tasks round-robin to             |
   |                   | multi-processor nodes can make the memory usage more uniform across nodes, which may allow more |
   |                   | MPI tasks per node than the standard layout.                                                    |
   +-------------------+-------------------------------------------------------------------------------------------------+
   | pe_to_task_list   | Logical PE to MPI task mapping. See above for more description.                                 |
   +-------------------+-------------------------------------------------------------------------------------------------+
   | my_pe             | The logical PE number for the MPI task.                                                         |
   +-------------------+-------------------------------------------------------------------------------------------------+
   | layout_type       | Controls the mapping type between MPI tasks and PEs. Currently type 1 is the standard layout    |
   |                   | (one-to-one mapping) and type 2 is a round-robin mapping where each node gets a task in turn    |
   |                   | before assigning a second task to each node, until all tasks are assigned.                      |
   +-------------------+-------------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *call init_ensemble_manager(ens_handle, num_copies, num_vars [, distribution_type_in] [, layout_type])*
   ::

      type(ensemble_type), intent(out) :: ens_handle
      integer,             intent(in)  :: num_copies
      integer,             intent(in)  :: num_vars
      integer, optional,   intent(in)  :: distribution_type_in
      integer, optional,   intent(in)  :: layout_type

.. container:: indent1

   Initializes an instance of an ensemble. Storage is allocated and the size descriptions in the ensemble_type are
   initialized.

   +------------------------+--------------------------------------------------------------------------------------------+
   | ``ens_handle``         | Handle for the ensemble being initialized                                                  |
   +------------------------+--------------------------------------------------------------------------------------------+
   | ``num_copies``         | Number of copies of vector.                                                                |
   +------------------------+--------------------------------------------------------------------------------------------+
   | ``num_vars``           | Number of variables in the vector.                                                         |
   +------------------------+--------------------------------------------------------------------------------------------+
   | *distribution_type_in* | Controls layout of storage on PEs. Currently only option 1 is supported.                   |
   +------------------------+--------------------------------------------------------------------------------------------+
   | *layout_type*          | Controls layout of MPI tasks on PEs. Type 1 is the default, where MPI tasks are assigned   |
   |                        | to PEs on a one-to-one basis. Type 2 is a round-robin assignment where each node gets one  |
   |                        | task before the nodes are assigned a second task. If running with more MPI tasks than      |
   |                        | ``num_copies``, this can result in a more uniform usage of memory across the nodes.        |
   +------------------------+--------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *call read_ensemble_restart(ens_handle, start_copy, end_copy, start_from_restart, file_name [, init_time] [,
   force_single_file])*
   ::

      type(ensemble_type),       intent(inout) :: ens_handle
      integer,                   intent(in)    :: start_copy
      integer,                   intent(in)    :: end_copy
      logical,                   intent(in)    :: start_from_restart
      character(len=*),          intent(in)    :: file_name
      type(time_type), optional, intent(in)    :: init_time
      logical, optional,         intent(in)    :: force_single_file

.. container:: indent1

   Read in a set of copies of a vector from file ``file_name``. The copies read are place into global copies
   start_copy:end_copy in the ens_handle. If start_from_restart is false, then only a single copy of the vector is read
   from the file and then it is perturbed using routines in assim_model_mod to generate the required number of copies.
   The read can be from a single file that contains all needed copies or from a different file for each copy. This
   choice is controlled by the namelist entry single_restart_file_in. However, the optional argument force_single_file
   forces the read to be from a single file if it is present and true. This is used for ensembles that contain the
   inflation values for state space inflation. If multiple files are to be read, the file names are generated by
   appending integers to the input file_name. If the input is a single file all reads are done sequentially by process 0
   and then shipped to the PE that stores that copy. If the input is multiple files each MPI task reads the copies it
   stores directly and independently.

   ====================== ===============================================================================================
   ``ens_handle``         Handle of ensemble.
   ``start_copy``         Global index of first of continguous set of copies to be read.
   ``end_copy``           Global index of last of contiguous set of copies to be read, copies(start_copy:end_copy).
   ``start_from_restart`` If true, read all copies from file. If false, read one copy and perturb to get required number.
   ``file_name``          Name of file from which to read.
   *init_time*            If present, set time of all copies read to this value.
   *force_single_file*    If present and true, force the read to be from a single file which contains all copies.
   ====================== ===============================================================================================

| 

.. container:: routine

   *call write_ensemble_restart(ens_handle, file_name, start_copy, end_copy [, force_single_file])*
   ::

      type(ensemble_type), intent(inout) :: ens_handle
      character(len=*),    intent(in)    :: file_name
      integer,             intent(in)    :: start_copy
      integer,             intent(in)    :: end_copy
      logical, optional,   intent(in)    :: force_single_file

.. container:: indent1

   Writes a set of copies of a vector to file file_name. The copies written are from global copies start_copy:end_copy
   in the ens_handle. The write can be to a single file or to a different file for each copy. This choice is controlled
   by the namelist entry single_restart_file_out. However, the optional argument force_single_file forces the write to
   be to a single file if it is present and true. This is used for ensembles that contain the inflation values for state
   space inflation. If multiple files are to be written, the file names are generated by appending integers to the input
   file_name. If the output is a single file all copies are shipped from the PE that stores that copy to process 0, and
   then written out sequentially. If the output is to multiple files each MPI task writes the copies it stores directly
   and independently.

   =================== ============================================================================================
   ``ens_handle``      Handle of ensemble.
   ``file_name``       Name of file from which to read.
   ``start_copy``      Global index of first of continguous set of copies to be written.
   ``end_copy``        Global index of last of contiguous set of copies to be written, copies(start_copy:end_copy).
   *force_single_file* If present and true, force the write to be to a single file which contains all copies.
   =================== ============================================================================================

| 

.. container:: routine

   *call get_copy(receiving_pe, ens_handle, copy, vars [, mtime])*
   ::

      integer,                   intent(in)  :: receiving_pe
      type(ensemble_type),       intent(in)  :: ens_handle
      integer,                   intent(in)  :: copy
      real(r8), dimension(:),    intent(out) :: vars
      type(time_type), optional, intent(out) :: mtime

.. container:: indent1

   Retrieves a copy of the state vector, indexed by the global index copy. The process that is to receive the copy is
   receiving_pe and the copy is returned in the one dimensional array vars. The time of the copy is also returned if
   mtime is present. This is generally used for operations, like IO, that require a single processor to do things with
   the entire state vector. Data is only returned in vars on the receiving PE; vars on all other PEs is unset.

   +------------------+--------------------------------------------------------------------------------------------------+
   | ``receiving_pe`` | This process ends up with the requested copy of the state vector.                                |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``ens_handle``   | Handle for ensemble.                                                                             |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``copy``         | The global index of the copy of the state vector that is to be retrieved.                        |
   +------------------+--------------------------------------------------------------------------------------------------+
   | ``vars``         | One dimensional array in which the requested copy of the state vector is returned. Data is only  |
   |                  | returned in vars on the receiving PE; vars on all other PEs is unset.                            |
   +------------------+--------------------------------------------------------------------------------------------------+
   | *mtime*          | If present returns the time of the requested copy.                                               |
   +------------------+--------------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *call put_copy(sending_pe, ens_handle, copy, vars [, mtime])*
   ::

      integer,                   intent(in)    :: sending_pe
      type(ensemble_type),       intent(inout) :: ens_handle
      integer,                   intent(in)    :: copy
      real(r8), dimension(:),    intent(in)    :: vars
      type(time_type), optional, intent(in)    :: mtime

.. container:: indent1

   Sends a state vector, in vars, from the given process to the process storing the global index copy. The time of the
   copy is also sent if mtime is present. This is generally used for operations, like IO, that require a single
   processor to do things with the entire state vector. For instance, if a single process reads in a state vector, it
   can be shipped to the storing process by this subroutine. Only the data in vars on the sending PE is processed; vars
   on all other PEs is ignored.

   +----------------+----------------------------------------------------------------------------------------------------+
   | ``sending_pe`` | This process sends the copy of the state vector.                                                   |
   +----------------+----------------------------------------------------------------------------------------------------+
   | ``ens_handle`` | Handle for ensemble.                                                                               |
   +----------------+----------------------------------------------------------------------------------------------------+
   | ``copy``       | The global index of the copy of the state vector that is to be sent.                               |
   +----------------+----------------------------------------------------------------------------------------------------+
   | ``vars``       | One dimensional array in which the requested copy of the state vector is located. Only the data in |
   |                | vars on the sending PE is processed; vars on all other PEs is ignored.                             |
   +----------------+----------------------------------------------------------------------------------------------------+
   | *mtime*        | If present send the time of the copy.                                                              |
   +----------------+----------------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *call broadcast_copy(ens_handle, copy, arraydata)*
   ::

      type(ensemble_type),    intent(in)   :: ens_handle
      integer,                intent(in)   :: copy
      real(r8), dimension(:), intent(out)  :: arraydata

.. container:: indent1

   Finds which PE has the global index copy and broadcasts that copy to all PEs. ``arraydata`` is an output on all PEs,
   even on the PE which is the owner if it is separate storage from the vars array in the ensemble handle. This is a
   collective routine, which means it must be called by all processes in the job.

   +----------------+----------------------------------------------------------------------------------------------------+
   | ``ens_handle`` | Handle for ensemble.                                                                               |
   +----------------+----------------------------------------------------------------------------------------------------+
   | ``copy``       | The global index of the copy of the state vector that is to be sent.                               |
   +----------------+----------------------------------------------------------------------------------------------------+
   | ``arraydata``  | One dimensional array into which the requested copy of the state vector will be copied on all PEs, |
   |                | including the sending PE.                                                                          |
   +----------------+----------------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *call set_ensemble_time(ens_handle, indx, mtime)*
   ::

      type(ensemble_type), intent(inout) :: ens_handle
      integer,             intent(in)    :: indx
      type(time_type),     intent(in)    :: mtime

.. container:: indent1

   Set the time of a copy to the given value. ``indx`` in this case is the local copy number for a specific task.
   get_copy_owner_index() can be called to see if you are the owning task for a given global copy number, and to get the
   local index number for that copy.

   ============== ==================================================================
   ``ens_handle`` Handle for ensemble.
   ``indx``       The local index of the copy of the state vector that is to be set.
   ``mtime``      The time to set for this copy.
   ============== ==================================================================

| 

.. container:: routine

   *call get_ensemble_time(ens_handle, indx, mtime)*
   ::

      type(ensemble_type), intent(in)   :: ens_handle
      integer,             intent(in)   :: indx
      type(time_type),     intent(out)  :: mtime

.. container:: indent1

   Get the time associated with a copy. ``indx`` in this case is the local copy number for a specific task.
   get_copy_owner_index() can be called to see if you are the owning task for a given global copy number, and to get the
   local index number for that copy.

   ============== ======================================================
   ``ens_handle`` Handle for ensemble.
   ``indx``       The local index of the copy to retrieve the time from.
   ``mtime``      The returned time value.
   ============== ======================================================

| 

.. container:: routine

   *call end_ensemble_manager(ens_handle)*
   ::

      type(ensemble_type), intent(in)  :: ens_handle

.. container:: indent1

   Frees up storage associated with an ensemble.

   ============== =======================
   ``ens_handle`` Handle for an ensemble.
   ============== =======================

| 

.. container:: routine

   *call duplicate_ens(ens1, ens2, duplicate_time)*
   ::

      type(ensemble_type), intent(in)    :: ens1
      type(ensemble_type), intent(inout) :: ens2
      logical, intent(in)                :: duplicate_time

.. container:: indent1

   Copies the contents of the vars array from ens1 into ens2. If the num_copies and num_vars are not consistent or if
   the distribution_type is not consistent, fails with an error. If duplicate_time is true, the times from ens1 are
   copied over the times of ens2. Only the vars array data is copied from the source to the destination. Transpose the
   data after duplication if you want to access the copies.

   ================== ================================================================================================
   ``ens1``           Ensemble handle of ensemble to be copies into ens2. Data from the vars array will be replicated.
   ``ens2``           Ensemble handle of ensemble into which ens1 vars data will be copied.
   ``duplicate_time`` If true, copy the times from ens1 into ens2, else leave ens2 times unchanged.
   ================== ================================================================================================

| 

.. container:: routine

   *var = get_my_num_copies(ens_handle)*
   ::

      integer                          :: get_my_num_copies
      type(ensemble_type), intent(in)  :: ens_handle

.. container:: indent1

   Returns number of copies stored by this process when storing all variables for a subset of copies. Same as num_copies
   if running with only a single process.

   ============== ======================================================================================================
   ``var``        Returns the number of copies stored by this process when storing all variables for a subset of copies.
   ``ens_handle`` Handle for an ensemble.
   ============== ======================================================================================================

| 

.. container:: routine

   *var = get_my_num_vars(ens_handle)*
   ::

      integer                         :: get_my_num_vars
      type(ensemble_type), intent(in) :: ens_handle

.. container:: indent1

   Returns number of variables stored by this process when storing all copies of a subset of variables. Same as num_vars
   if running with only a single process.

   ============== ===================================================================================================
   ``var``        Returns the number of vars stored by this process when storing all copies of a subset of variables.
   ``ens_handle`` Handle for an ensemble.
   ============== ===================================================================================================

| 

.. container:: routine

   *call get_my_copies(ens_handle, copies)*
   ::

      type(ensemble_type), intent(in) :: ens_handle
      integer, intent(out)            :: copies(:)

.. container:: indent1

   Returns a list of the global copy numbers stored on this process when storing subset of copies of all variables.

   ============== =========================================================================================
   ``ens_handle`` Handle for an ensemble.
   ``copies``     List of all copies stored by this process when storing subset of copies of all variables.
   ============== =========================================================================================

| 

.. container:: routine

   *call get_my_vars(ens_handle, vars)*
   ::

      type(ensemble_type), intent(in) :: ens_handle
      integer, intent(out)            :: vars(:)

.. container:: indent1

   Returns a list of the global variable numbers stored on this process when storing all copies of a subset of
   variables.

   ============== ==============================================================================================
   ``ens_handle`` Handle for an ensemble.
   ``vars``       List of all variables stored on this process when storing all copies of a subset of variables.
   ============== ==============================================================================================

| 

.. container:: routine

   *call get_copy_owner_index(copy_number, owner, owners_index)*
   ::

      integer, intent(in)  :: copy_number
      integer, intent(out) :: owner
      integer, intent(out) :: owners_index

.. container:: indent1

   Given the global index of a copy number, returns the PE that stores this copy when all variables of a subset of
   copies are stored and the local storage index for this copy on that process.

   ================ =============================================================================================
   ``copy_number``  Global index of a copy from an ensemble.
   ``owner``        Process Element (PE) that stores this copy when each has all variables of a subset of copies.
   ``owners_index`` Local storage index for this copy on the owning process.
   ================ =============================================================================================

| 

.. container:: routine

   *call get_var_owner_index(var_number, owner, owners_index)*
   ::

      integer, intent(in)  :: var_number
      integer, intent(out) :: owner
      integer, intent(out) :: owners_index

.. container:: indent1

   Given the global index of a variable in the vector, returns the PE that stores this variable when all copies of a
   subset of variables are stored and the local storage index for this variable on that process.

   ================ ===============================================================================================
   ``var_number``   Global index of a variable in the vector from an ensemble.
   ``owner``        Process Element (PE) that stores this variable when each has all copies of subset of variables.
   ``owners_index`` Local storage index for this variable on the owning process.
   ================ ===============================================================================================

| 

.. container:: routine

   *call all_vars_to_all_copies(ens_handle, label)*
   ::

      type(ensemble_type), intent(inout)        :: ens_handle
      character(len=*),    intent(in), optional :: label

.. container:: indent1

   Transposes data from a representation in which each PE has a subset of copies of all variables to one in which each
   has all copies of a subset of variables. In the current implementation, storage is not released so both
   representations are always available. However, one representation may be current while the other is out of date.

   Different different numbers of copies, different lengths of the vectors, different numbers of PEs and different
   implementations of the MPI parallel libraries can have very different performance characteristics. The namelist item
   ``communication_configuration`` controls one of four possible combinations of the operation order during the
   transposes. If performance is an issue the various settings on this namelist item can be explored. See the namelist
   section for more details.

   The transpose routines make both representations of the data equivalent until the next update to either the copies or
   the vars arrays, so either can be used as a data source.

   +----------------+----------------------------------------------------------------------------------------------------+
   | ``ens_handle`` | The handle of the ensemble being transposed.                                                       |
   +----------------+----------------------------------------------------------------------------------------------------+
   | ``label``      | A character string label. If present, a timestamp with this label is printed at the start and end  |
   |                | of the transpose.                                                                                  |
   +----------------+----------------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *call all_copies_to_all_vars(ens_handle, label)*
   ::

      type(ensemble_type), intent(inout)        :: ens_handle
      character(len=*),    intent(in), optional :: label

.. container:: indent1

   Transposes data from a representation in which each processor has all copies of a subset of variables to one in which
   each has a subset of copies of all variables. In the current implementation, storage is not released so both
   representations are always available. However, one representation may be current while the other is out of date.

   Different different numbers of copies, different lengths of the vectors, different numbers of PEs and different
   implementations of the MPI parallel libraries can have very different performance characteristics. The namelist item
   ``communication_configuration`` controls one of four possible combinations of the operation order during the
   transposes. If performance is an issue the various settings on this namelist item can be explored. See the namelist
   section for more details.

   The transpose routines make both representations of the data equivalent until the next update to either the copies or
   the vars arrays, so either can be used as a data source.

   +----------------+----------------------------------------------------------------------------------------------------+
   | ``ens_handle`` | The handle of the ensemble being transposed.                                                       |
   +----------------+----------------------------------------------------------------------------------------------------+
   | ``label``      | A character string label. If present, a timestamp with this label is printed at the start and end  |
   |                | of the transpose.                                                                                  |
   +----------------+----------------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *call compute_copy_mean(ens_handle, start_copy, end_copy, mean_copy)*
   ::

      type(ensemble_type), intent(inout) :: ens_handle
      integer,             intent(in)    :: start_copy
      integer,             intent(in)    :: end_copy
      integer,             intent(in)    :: mean_copy

.. container:: indent1

   Computes the mean of a contiguous subset of copies starting with global index start_copy and ending with global index
   ens_copy. Mean is written to global index mean_copy.

   When this routine is called the ensemble must have all copies of a subset of the vars. It updates the copies array
   with the mean, so after this call the copies array data is more current and the vars data is stale.

   ============== ======================================================
   ``ens_handle`` Handle for an ensemble.
   ``start_copy`` Global index of first copy in mean and sd computation.
   ``end_copy``   Global index of last copy in mean and sd computation.
   ``mean_copy``  Global index of copy into which mean is written.
   ============== ======================================================

| 

.. container:: routine

   *call compute_copy_mean_sd(ens_handle, start_copy, end_copy, mean_copy, sd_copy)*
   ::

      type(ensemble_type), intent(inout) :: ens_handle
      integer,             intent(in)    :: start_copy
      integer,             intent(in)    :: end_copy
      integer,             intent(in)    :: mean_copy
      integer,             intent(in)    :: sd_copy

.. container:: indent1

   Computes the mean and standard deviation of a contiguous subset of copies starting with global index start_copy and
   ending with global index ens_copy. Mean is written to index mean_copy and standard deviation to index sd_copy.

   When this routine is called the ensemble must have all copies of a subset of the vars. It updates the copies arrays
   with the mean and sd, so after this call the copies array data is more current and the vars data is stale.

   ============== ==============================================================
   ``ens_handle`` Handle for an ensemble.
   ``start_copy`` Global index of first copy in mean and sd computation.
   ``end_copy``   Global index of last copy in mean and sd computation.
   ``mean_copy``  Global index of copy into which mean is written.
   ``sd_copy``    Global index of copy into which standard deviation is written.
   ============== ==============================================================

| 

.. container:: routine

   *call compute_copy_mean_var(ens_handle, start_copy, end_copy, mean_copy, var_copy)*
   ::

      type(ensemble_type), intent(inout) :: ens_handle
      integer,             intent(in)  :: start_copy
      integer,             intent(in)  :: end_copy
      integer,             intent(in)  :: mean_copy
      integer,             intent(in)  :: var_copy

.. container:: indent1

   Computes the mean and variance of a contiguous subset of copies starting with global index start_copy and ending with
   global index ens_copy. Mean is written to index mean_copy and variance to index var_copy.

   When this routine is called the ensemble must have all copies of a subset of the vars. It updates the copies arrays
   with the mean and variance, so after this call the copies array data is more current and the vars data is stale.

   ============== ======================================================
   ``ens_handle`` Handle for an ensemble.
   ``start_copy`` Global index of first copy in mean and sd computation.
   ``end_copy``   Global index of last copy in mean and sd computation.
   ``mean_copy``  Global index of copy into which mean is written.
   ``var_copy``   Global index of copy into which variance is written.
   ============== ======================================================

| 

Private interfaces
------------------

== =======================
\  assign_tasks_to_pes
\  calc_tasks_on_each_node
\  create_pe_to_task_list
\  get_copy_list
\  get_max_num_copies
\  get_max_num_vars
\  get_var_list
\  round_robin
\  set_up_ens_distribution
\  simple_layout
\  sort_task_list
\  timestamp_message
== =======================

| 

.. container:: routine

   *var = get_max_num_copies(num_copies)*
   ::

      integer              :: get_max_num_copies
      integer, intent(in)  :: num_copies

.. container:: indent1

   Returns the largest number of copies that are on any pe when var complete. Depends on distribution_type with only
   option 1 currently implemented. Used to get size for creating storage to receive a list of the copies on a PE.

   ============== ============================================================================
   ``var``        Returns the largest number of copies any an individual PE when var complete.
   ``num_copies`` Total number of copies in the ensemble.
   ============== ============================================================================

| 

.. container:: routine

   *var = get_max_num_vars(num_vars)*
   ::

      integer              :: get_max_num_vars
      integer, intent(in)  :: num_vars

.. container:: indent1

   Returns the largest number of vars that are on any pe when copy complete. Depends on distribution_type with only
   option 1 currently implemented. Used to get size for creating storage to receive a list of the vars on a PE.

   ============== ===========================================================================
   ``var``        Returns the largest number of vars any an individual PE when copy complete.
   ``num_copies`` Total number of vars in an ensemble vector.
   ============== ===========================================================================

| 

.. container:: routine

   *call set_up_ens_distribution(ens_handle)*
   ::

      type(ensemble_type), intent(inout) :: ens_handle

.. container:: indent1

   Figures out how to lay out the copy complete and vars complete distributions. The distribution_type identifies
   different options. Only distribution_type 1 is implemented. This puts every Nth var or copy on a given processor
   where N is the total number of processes.

   ============== =======================
   ``ens_handle`` Handle for an ensemble.
   ============== =======================

| 

.. container:: routine

   *call get_var_list(num_vars, pe, var_list, pes_num_vars)*
   ::

      integer,   intent(in)     :: num_vars
      integer,   intent(in)     :: pe
      integer,   intent(out)    :: var_list(:)
      integer,   intent(out)    :: pes_num_vars

Returns a list of the vars stored by process pe when copy complete and the number of these vars. var_list must be
dimensioned large enough to hold all vars. Depends on distribution_type with only option 1 currently implemented.

| 

.. container:: routine

   *call get_copy_list(num_copies, pe, copy_list, pes_num_copies)*
   ::

      integer,   intent(in)     :: num_copies
      integer,   intent(in)     :: pe
      integer,   intent(out)    :: copy_list(:)
      integer,   intent(out)    :: pes_num_copies

Returns a list of the copies stored by process pe when var complete and the number of these copies. copy_list must be
dimensioned large enough to hold all copies. Depends on distribution_type with only option 1 currently implemented.

| 

.. container:: routine

   *call timestamp_message(msg [, sync] [, alltasks])*
   ::

      character(len=*), intent(in)           :: msg
      logical,          intent(in), optional :: sync
      logical,          intent(in), optional :: alltasks

.. container:: indent1

   Write current time and message to stdout and log file. If sync is present and true, sync mpi jobs before printing
   time. If alltasks is present and true, all tasks print the time. The default is only task 0 prints a timestamp.

   +------------+--------------------------------------------------------------------------------------------------------+
   | ``msg``    | character string to prepend to the time info                                                           |
   +------------+--------------------------------------------------------------------------------------------------------+
   | *sync*     | if present and true, execute an MPI_Barrier() to sync all MPI tasks before printing the time. this     |
   |            | means the time will be the value of the slowest of the tasks to reach this point.                      |
   +------------+--------------------------------------------------------------------------------------------------------+
   | *alltasks* | if present and true, have all tasks print out a timestamp. the default is for just task 0 to print.    |
   |            | the usual combination is either sync=true and alltasks=false, or sync=false and alltasks=true.         |
   +------------+--------------------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *call print_ens_handle(ens_handle, force, label)*
   ::

      type(ensemble_type),        intent(in) :: ens_handle
      logical,          optional, intent(in) :: force
      character(len=*), optional, intent(in) :: label

.. container:: indent1

   For debugging use, dump the contents of an ensemble handle derived type. If the ``debug`` namelist item is true, this
   will print in any case. If ``debug`` is false, set ``force`` to true to force printing. The optional string label can
   help provide context for the output.

   ============== =============================================================================
   ``ens_handle`` The derived type to print information about.
   ``force``      If the ``debug`` namelist item is false, set this to true to enable printing.
   ``label``      Optional string label to print to provide context for the output.
   ============== =============================================================================

| 

.. container:: routine

   *call assign_tasks_to_pes(ens_handle, nEns_members, layout_type)*
   ::

      type(ensemble_type), intent(inout)    :: ens_handle
      integer,             intent(in)       :: nEns_members
      integer,             intent(inout)    :: layout_type

.. container:: indent1

   Calulate the task layout based on the tasks per node and the total number of tasks. Allows the user to spread out the
   ensemble members as much as possible to balance memory usage between nodes. Possible options: 1. Standard task layout
   - first n tasks have the ensemble members my_pe = my_task_id() 2. Round-robin on the nodes

   ============== =======================
   ``ens_handle`` Handle for an ensemble.
   \              
   \              
   ============== =======================

| 

.. container:: routine

   *call round_robin(ens_handle)*
   ::

      type(ensemble_type), intent(inout)    :: ens_handle

.. container:: indent1

   Round-robin MPI task layout starting at the first node. Starting on the first node forces pe 0 = task 0. 

   ============== =======================
   ``ens_handle`` Handle for an ensemble.
   ============== =======================

| 

.. container:: routine

   *call create_pe_to_task_list(ens_handle)*
   ::

      type(ensemble_type), intent(inout)    :: ens_handle

.. container:: indent1

   Creates the ``ens_handle%pe_to_task_list``. ``ens_handle%task_to_pe_list`` must have been assigned first, otherwise
   this routine will just return nonsense.

   ============== =======================
   ``ens_handle`` Handle for an ensemble.
   ============== =======================

| 

.. container:: routine

   *call calc_tasks_on_each_node(nodes, last_node_task_number)*
   ::

      integer, intent(out)  :: last_node_task_number
      integer, intent(out)  :: nodes

Finds the of number nodes and how many tasks are on the last node, given the number of tasks and the tasks_per_node
(ptile). The total number of tasks is num_pes = task_count() The last node may have fewer tasks, for example, if ptile =
16 and the number of mpi tasks = 17

| 

.. container:: routine

   *call simple_layout(ens_handle, n)*
   ::

      type(ensemble_type), intent(inout) :: ens_handle
      integer,             intent(in)    :: n

.. container:: indent1

   assigns the arrays task_to_pe_list and pe_to_task list for the simple layout where my_pe = my_task_id()

   ens_handle
      Handle for an ensemble.
   n
      size

   .. container:: routine

      *call sort_task_list(i, idx, n)*
      ::

         integer, intent(in)    :: n
         integer, intent(inout) :: x(n)   ! array to be sorted
         integer, intent(out)   :: idx(n) ! index of sorted array

   sorts an array and returns the sorted array, and the index of the original array

   n
      size
   x(n)
      array to be sorted
   idx(n)
      index of sorted array

   .. container:: routine

      *call map_pe_to_task(ens_handle, p)*
      ::

         type(ensemble_type), intent(in) :: ens_handle
         integer,             intent(in) :: p

   .. container:: indent1

      Return the physical task for my_pe

      ============== =================================================
      ``ens_handle`` Handle for an ensemble.
      ``p``          The MPI task corresponding to the given PE number
      ============== =================================================

   .. container:: routine

      *call map_task_to_pe(ens_handle, t)*
      ::

         type(ensemble_type), intent(in) :: ens_handle
         integer,             intent(in) :: t

   .. container:: indent1

      Return my_pe corresponding to the physical task

      ============== =========================================================
      ``ens_handle`` Handle for an ensemble.
      ``t``          Return the PE corresponding to the given MPI task number.
      ============== =========================================================

   .. rubric:: Files
      :name: files

   -  input.nml
   -  State vector restart files, either one for all copies or one per copy.
   -  State vector output files, either one for all copies or one per copy.

   .. rubric:: References
      :name: references

   #. none

   .. rubric:: Private components
      :name: private-components

   N/A

.. |communication pattern| image:: ../../../guide/images/comm_pattern96.png
   :width: 400px
.. |communication pattern 2| image:: ../../../guide/images/comm_pattern512.png
   :width: 400px
