module mpp_mod
==============

Overview
--------

``mpp_mod``, is a set of simple calls to provide a uniform interface to different message-passing libraries. It
currently can be implemented either in the SGI/Cray native SHMEM library or in the MPI standard. Other libraries (e.g
MPI-2, Co-Array Fortran) can be incorporated as the need arises.

The data transfer between a processor and its own memory is based on ``load`` and ``store`` operations upon memory.
Shared-memory systems (including distributed shared memory systems) have a single address space and any processor can
acquire any data within the memory by ``load`` and ``store``. The situation is different for distributed parallel
systems. Specialized MPP systems such as the T3E can simulate shared-memory by direct data acquisition from remote
memory. But if the parallel code is distributed across a cluster, or across the Net, messages must be sent and received
using the protocols for long-distance communication, such as TCP/IP. This requires a \``handshaking'' between nodes of
the distributed system. One can think of the two different methods as involving ``put``\ s or ``get``\ s (e.g the SHMEM
library), or in the case of negotiated communication (e.g MPI), ``send``\ s and ``recv``\ s.

The difference between SHMEM and MPI is that SHMEM uses one-sided communication, which can have very low-latency
high-bandwidth implementations on tightly coupled systems. MPI is a standard developed for distributed computing across
loosely-coupled systems, and therefore incurs a software penalty for negotiating the communication. It is however an
open industry standard whereas SHMEM is a proprietary interface. Besides, the ``put``\ s or ``get``\ s on which it is
based cannot currently be implemented in a cluster environment (there are recent announcements from Compaq that occasion
hope).

The message-passing requirements of climate and weather codes can be reduced to a fairly simple minimal set, which is
easily implemented in any message-passing API. ``mpp_mod`` provides this API.

Features of ``mpp_mod`` include:

#. Simple, minimal API, with free access to underlying API for more complicated stuff.
#. Design toward typical use in climate/weather CFD codes.
#. Performance to be not significantly lower than any native API.

This module is used to develop higher-level calls for :doc:`./mpp_domains` and :doc:`./mpp_io`.

Parallel computing is initially daunting, but it soon becomes second nature, much the way many of us can now write
vector code without much effort. The key insight required while reading and writing parallel code is in arriving at a
mental grasp of several independent parallel execution streams through the same code (the SPMD model). Each variable you
examine may have different values for each stream, the processor ID being an obvious example. Subroutines and function
calls are particularly subtle, since it is not always obvious from looking at a call what synchronization between
execution streams it implies. An example of erroneous code would be a global barrier call (see ` <#mpp_sync>`__ below)
placed within a code block that not all PEs will execute, e.g:

::


   if( pe.EQ.0 )call mpp_sync()

Here only PE 0 reaches the barrier, where it will wait indefinitely. While this is a particularly egregious example to
illustrate the coding flaw, more subtle versions of the same are among the most common errors in parallel code.

It is therefore important to be conscious of the context of a subroutine or function call, and the implied
synchronization. There are certain calls here (e.g ``mpp_declare_pelist, mpp_init, mpp_malloc, mpp_set_stack_size``)
which must be called by all PEs. There are others which must be called by a subset of PEs (here called a ``pelist``)
which must be called by all the PEs in the ``pelist`` (e.g ``mpp_max, mpp_sum, mpp_sync``). Still others imply no
synchronization at all. I will make every effort to highlight the context of each call in the MPP modules, so that the
implicit synchronization is spelt out.

For performance it is necessary to keep synchronization as limited as the algorithm being implemented will allow. For
instance, a single message between two PEs should only imply synchronization across the PEs in question. A *global*
synchronization (or *barrier*) is likely to be slow, and is best avoided. But codes first parallelized on a Cray T3E
tend to have many global syncs, as very fast barriers were implemented there in hardware.

Another reason to use pelists is to run a single program in MPMD mode, where different PE subsets work on different
portions of the code. A typical example is to assign an ocean model and atmosphere model to different PE subsets, and
couple them concurrently instead of running them serially. The MPP module provides the notion of a *current pelist*,
which is set when a group of PEs branch off into a subset. Subsequent calls that omit the ``pelist`` optional argument
(seen below in many of the individual calls) assume that the implied synchronization is across the current pelist. The
calls ``mpp_root_pe`` and ``mpp_npes`` also return the values appropriate to the current pelist. The
``mpp_set_current_pelist`` call is provided to set the current pelist.

| 

Other modules used
------------------

.. container::

   ::

      shmem_interface
                  mpi

Public interface
----------------

.. container::

   F90 is a strictly-typed language, and the syntax pass of the compiler requires matching of type, kind and rank (TKR).
   Most calls listed here use a generic type, shown here as ``MPP_TYPE_``. This is resolved in the pre-processor stage
   to any of a variety of types. In general the MPP operations work on 4-byte and 8-byte variants of
   ``integer, real, complex, logical`` variables, of rank 0 to 5, leading to 48 specific module procedures under the
   same generic interface. Any of the variables below shown as ``MPP_TYPE_`` is treated in this way.
   ::

      use mpp_mod [, only:  mpp_max,
                            mpp_sum,
                            mpp_transmit,
                            mpp_broadcast,
                            mpp_chksum,
                            mpp_error,
                            mpp_init,
                            stdin,
                            mpp_exit,
                            mpp_pe,
                            mpp_npes,
                            mpp_declare_pelist,
                            mpp_set_current_pelist,
                            mpp_clock_set_grain,
                            mpp_sync,
                            mpp_sync_self,
                            mpp_malloc,
                            mpp_set_stack_size ]

   mpp_max:
      Reduction operations.
   mpp_sum:
      Reduction operation.
   mpp_transmit:
      Basic message-passing call.
   mpp_broadcast:
      Parallel broadcasts.
   mpp_chksum:
      Parallel checksums.
   mpp_error:
      Error handler.
   mpp_init:
      Initialize ``mpp_mod``.
   stdin:
      Standard fortran unit numbers.
   mpp_exit:
      Exit ``mpp_mod``.
   mpp_pe:
      Returns processor ID.
   mpp_npes:
      Returns processor count for current pelist.
   mpp_declare_pelist:
      Declare a pelist.
   mpp_set_current_pelist:
      Set context pelist.
   mpp_clock_set_grain:
      Set the level of granularity of timing measurements.
   mpp_sync:
      Global synchronization.
   mpp_sync_self:
      Local synchronization.
   mpp_malloc:
      Symmetric memory allocation.
   mpp_set_stack_size:
      Allocate module internal workspace.

| 

Public data
-----------

.. container::

   None.

Public routines
---------------

a. .. rubric:: Mpp_max
      :name: mpp_max

   ::

      call mpp_max ( a, pelist )

   **DESCRIPTION**
      Find the max of scalar a the PEs in pelist result is also automatically broadcast to all PEs
   **INPUT**
      +------------+--------------------------------------------------------------------------------------------------------+
      | ``a``      | ``real`` or ``integer``, of 4-byte of 8-byte kind.                                                     |
      +------------+--------------------------------------------------------------------------------------------------------+
      | ``pelist`` | If ``pelist`` is omitted, the context is assumed to be the current pelist. This call implies           |
      |            | synchronization across the PEs in ``pelist``, or the current pelist if ``pelist`` is absent.           |
      +------------+--------------------------------------------------------------------------------------------------------+

b. .. rubric:: Mpp_sum
      :name: mpp_sum

   ::

      call mpp_sum ( a, length, pelist )

   **DESCRIPTION**
      ``MPP_TYPE_`` corresponds to any 4-byte and 8-byte variant of ``integer, real, complex`` variables, of rank 0 or
      1. A contiguous block from a multi-dimensional array may be passed by its starting address and its length, as in
      ``f77``.
      Library reduction operators are not required or guaranteed to be bit-reproducible. In any case, changing the
      processor count changes the data layout, and thus very likely the order of operations. For bit-reproducible sums
      of distributed arrays, consider using the ``mpp_global_sum`` routine provided by the :doc:`./mpp_domains` module.
      The ``bit_reproducible`` flag provided in earlier versions of this routine has been removed.
      If ``pelist`` is omitted, the context is assumed to be the current pelist. This call implies synchronization
      across the PEs in ``pelist``, or the current pelist if ``pelist`` is absent.
   **INPUT**
      ``length`` ``pelist``
   **INPUT/OUTPUT**
      ``a``

c. .. rubric:: Mpp_transmit
      :name: mpp_transmit

   ::

      call mpp_transmit ( put_data, put_len, put_pe, get_data, get_len, get_pe )

   **DESCRIPTION**
      | ``MPP_TYPE_`` corresponds to any 4-byte and 8-byte variant of ``integer, real, complex, logical`` variables, of
        rank 0 or 1. A contiguous block from a multi-dimensional array may be passed by its starting address and its
        length, as in ``f77``.
      | ``mpp_transmit`` is currently implemented as asynchronous outward transmission and synchronous inward
        transmission. This follows the behaviour of ``shmem_put`` and ``shmem_get``. In MPI, it is implemented as
        ``mpi_isend`` and ``mpi_recv``. For most applications, transmissions occur in pairs, and are here accomplished
        in a single call.
      | The special PE designations ``NULL_PE``, ``ANY_PE`` and ``ALL_PES`` are provided by use association.
      | ``NULL_PE``: is used to disable one of the pair of transmissions.
      | ``ANY_PE``: is used for unspecific remote destination. (Please note that ``put_pe=ANY_PE`` has no meaning in the
        MPI context, though it is available in the SHMEM invocation. If portability is a concern, it is best avoided).
      | ``ALL_PES``: is used for broadcast operations.
      | It is recommended that ` <#mpp_broadcast>`__ be used for broadcasts.
      | The following example illustrates the use of ``NULL_PE`` and ``ALL_PES``:

      ::

             real, dimension(n) :: a
             if( pe.EQ.0 )then
                 do p = 1,npes-1
                    call mpp_transmit( a, n, p, a, n, NULL_PE )
                 end do
             else
                 call mpp_transmit( a, n, NULL_PE, a, n, 0 )
             end if
             
             call mpp_transmit( a, n, ALL_PES, a, n, 0 )

      | The do loop and the broadcast operation above are equivalent.
      | Two overloaded calls ``mpp_send`` and ``mpp_recv`` have also been provided. ``mpp_send`` calls ``mpp_transmit``
        with ``get_pe=NULL_PE``. ``mpp_recv`` calls ``mpp_transmit`` with ``put_pe=NULL_PE``. Thus the do loop above
        could be written more succinctly:

      ::

             if( pe.EQ.0 )then
                 do p = 1,npes-1
                    call mpp_send( a, n, p )
                 end do
             else
                 call mpp_recv( a, n, 0 )
             end if

d. .. rubric:: Mpp_broadcast
      :name: mpp_broadcast

   ::

      call mpp_broadcast ( data, length, from_pe, pelist )

   **DESCRIPTION**
      The ``mpp_broadcast`` call has been added because the original syntax (using ``ALL_PES`` in ``mpp_transmit``) did
      not support a broadcast across a pelist.
      ``MPP_TYPE_`` corresponds to any 4-byte and 8-byte variant of ``integer, real, complex, logical`` variables, of
      rank 0 or 1. A contiguous block from a multi-dimensional array may be passed by its starting address and its
      length, as in ``f77``.
      Global broadcasts through the ``ALL_PES`` argument to ` <#mpp_transmit>`__ are still provided for
      backward-compatibility.
      If ``pelist`` is omitted, the context is assumed to be the current pelist. ``from_pe`` must belong to the current
      pelist. This call implies synchronization across the PEs in ``pelist``, or the current pelist if ``pelist`` is
      absent.
   **INPUT**
      ``length``, ``from_pe``, ``pelist``
   **INPUT/OUTPUT**
      ``data(*)``

e. .. rubric:: Mpp_chksum
      :name: mpp_chksum

   ::

       
      mpp_chksum ( var, pelist )

   **DESCRIPTION**
      | ``mpp_chksum`` is a parallel checksum routine that returns an identical answer for the same array irrespective
        of how it has been partitioned across processors. ``LONG_KIND``\ is the ``KIND`` parameter corresponding to long
        integers (see discussion on OS-dependent preprocessor directives) defined in the header file ``os.h``.
        ``MPP_TYPE_`` corresponds to any 4-byte and 8-byte variant of ``integer, real, complex, logical`` variables, of
        rank 0 to 5.
      | Integer checksums on FP data use the F90 ``TRANSFER()`` intrinsic.
      | The `serial checksum module <http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/shared/chksum/chksum.html>`__
        is superseded by this function, and is no longer being actively maintained. This provides identical results on a
        single-processor job, and to perform serial checksums on a single processor of a parallel job, you only need to
        use the optional ``pelist`` argument.

      ::

              use mpp_mod
              integer :: pe, chksum
              real :: a(:)
              pe = mpp_pe()
              chksum = mpp_chksum( a, (/pe/) )

      | The additional functionality of ``mpp_chksum`` over serial checksums is to compute the checksum across the PEs
        in ``pelist``. The answer is guaranteed to be the same for the same distributed array irrespective of how it has
        been partitioned.
      | If ``pelist`` is omitted, the context is assumed to be the current pelist. This call implies synchronization
        across the PEs in ``pelist``, or the current pelist if ``pelist`` is absent.

   **INPUT**
      ``pelist``, ``var``

f. .. rubric:: Mpp_error
      :name: mpp_error

   ::

      call mpp_error ( errortype, routine, errormsg )

   **DESCRIPTION**
      | It is strongly recommended that all error exits pass through ``mpp_error`` to assure the program fails cleanly.
        An individual PE encountering a ``STOP`` statement, for instance, can cause the program to hang. The use of the
        ``STOP`` statement is strongly discouraged.
      | Calling mpp_error with no arguments produces an immediate error exit, i.e:

      ::

             call mpp_error
             call mpp_error(FATAL)

      | are equivalent.
      | The argument order

      ::

             call mpp_error( routine, errormsg, errortype )

      | is also provided to support legacy code. In this version of the call, none of the arguments may be omitted.
      | The behaviour of ``mpp_error`` for a ``WARNING`` can be controlled with an additional call
        ``mpp_set_warn_level``.

      ::

             call mpp_set_warn_level(ERROR)

      causes ``mpp_error`` to treat ``WARNING`` exactly like ``FATAL``.

      ::

             call mpp_set_warn_level(WARNING)

      | resets to the default behaviour described above.
      | ``mpp_error`` also has an internal error state which maintains knowledge of whether a warning has been issued.
        This can be used at startup in a subroutine that checks if the model has been properly configured. You can
        generate a series of warnings using ``mpp_error``, and then check at the end if any warnings has been issued
        using the function ``mpp_error_state()``. If the value of this is ``WARNING``, at least one warning has been
        issued, and the user can take appropriate action:

      ::

             if( ... )call mpp_error( WARNING, '...' )
             if( ... )call mpp_error( WARNING, '...' )
             if( ... )call mpp_error( WARNING, '...' )
             ...
             if( mpp_error_state().EQ.WARNING )call mpp_error( FATAL, '...' )

   **INPUT**
      ``errortype``. One of ``NOTE``, ``WARNING`` or ``FATAL`` (these definitions are acquired by use association).
      ``NOTE`` writes ``errormsg`` to ``STDOUT``. ``WARNING`` writes ``errormsg`` to ``STDERR``. ``FATAL`` writes
      ``errormsg`` to ``STDERR``, and induces a clean error exit with a call stack traceback.

g. .. rubric:: Mpp_init
      :name: mpp_init

   ::

      call mpp_init ( flags )

   **DESCRIPTION**
      Called to initialize the ``mpp_mod`` package. It is recommended that this call be the first executed line in your
      program. It sets the number of PEs assigned to this run (acquired from the command line, or through the
      environment variable ``NPES``), and associates an ID number to each PE. These can be accessed by calling
      ` <#mpp_npes>`__ and ` <#mpp_pe>`__.
   **INPUT**
      ``flags``\ <``flags`` can be set to ``MPP_VERBOSE`` to have ``mpp_mod`` keep you informed of what it's up to.
      [integer]

h. .. rubric:: Stdin
      :name: stdin

   ::

       
      stdin ()

   **DESCRIPTION**
      This function, as well as stdout(), stderr(), stdlog(), returns the current standard fortran unit numbers for
      input, output, error messages and log messages. Log messages, by convention, are written to the file
      ``logfile.out``.

i. .. rubric:: Mpp_exit
      :name: mpp_exit

   ::

      call mpp_exit ()

   **DESCRIPTION**
      Called at the end of the run, or to re-initialize ``mpp_mod``, should you require that for some odd reason.
      This call implies synchronization across all PEs.

j. .. rubric:: Mpp_pe
      :name: mpp_pe

   ::

       
      mpp_pe ()

   **DESCRIPTION**
      This returns the unique ID associated with a PE. This number runs between 0 and ``npes-1``, where ``npes`` is the
      total processor count, returned by ``mpp_npes``. For a uniprocessor application this will always return 0.

k. .. rubric:: Mpp_npes
      :name: mpp_npes

   ::

       
      mpp_npes ()

   **DESCRIPTION**
      This returns the number of PEs in the current pelist. For a uniprocessor application, this will always return 1.

l. .. rubric:: Mpp_declare_pelist
      :name: mpp_declare_pelist

   ::

      call mpp_declare_pelist ( pelist,name )

   **DESCRIPTION**
      This call is written specifically to accommodate a MPI restriction that requires a parent communicator to create a
      child communicator, In other words: a pelist cannot go off and declare a communicator, but every PE in the parent,
      including those not in pelist(:), must get together for the ``MPI_COMM_CREATE`` call. The parent is typically
      ``MPI_COMM_WORLD``, though it could also be a subset that includes all PEs in ``pelist``.
      The restriction does not apply to SMA but to have uniform code, you may as well call it.
      This call implies synchronization across the PEs in the current pelist, of which ``pelist`` is a subset.
   **INPUT**
      ``pelist``
      [integer, dimension(:)]

m. .. rubric:: Mpp_set_current_pelist
      :name: mpp_set_current_pelist

   ::

      call mpp_set_current_pelist ( pelist )

   **DESCRIPTION**
      This call sets the value of the current pelist, which is the context for all subsequent "global" calls where the
      optional ``pelist`` argument is omitted. All the PEs that are to be in the current pelist must call it.
      In MPI, this call may hang unless ``pelist`` has been previous declared using ` <#mpp_declare_pelist>`__.
      If the argument ``pelist`` is absent, the current pelist is set to the "world" pelist, of all PEs in the job.
   **INPUT**
      ``pliest``
      [integer]

n. .. rubric:: Mpp_clock_set_grain
      :name: mpp_clock_set_grain

   ::

      call mpp_clock_set_grain ( grain )

   **DESCRIPTION**
      This routine and three other routines, mpp_clock_id, mpp_clock_begin(id), and mpp_clock_end(id) may be used to
      time parallel code sections, and extract parallel statistics. Clocks are identified by names, which should be
      unique in the first 32 characters. The ``mpp_clock_id`` call initializes a clock of a given name and returns an
      integer ``id``. This ``id`` can be used by subsequent ``mpp_clock_begin`` and ``mpp_clock_end`` calls set around a
      code section to be timed. Example:

      ::

             integer :: id
             id = mpp_clock_id( 'Atmosphere' )
             call mpp_clock_begin(id)
             call atmos_model()
             call mpp_clock_end()

      | Two flags may be used to alter the behaviour of ``mpp_clock``. If the flag ``MPP_CLOCK_SYNC`` is turned on by
        ``mpp_clock_id``, the clock calls ``mpp_sync`` across all the PEs in the current pelist at the top of the timed
        code section, but allows each PE to complete the code section (and reach ``mpp_clock_end``) at different times.
        This allows us to measure load imbalance for a given code section. Statistics are written to ``stdout`` by
        ``mpp_exit``.
      | The flag ``MPP_CLOCK_DETAILED`` may be turned on by ``mpp_clock_id`` to get detailed communication profiles.
        Communication events of the types ``SEND, RECV, BROADCAST, REDUCE`` and ``WAIT`` are separately measured for
        data volume and time. Statistics are written to ``stdout`` by ``mpp_exit``, and individual PE info is also
        written to the file ``mpp_clock.out.####`` where ``####`` is the PE id given by ``mpp_pe``.
      | The flags ``MPP_CLOCK_SYNC`` and ``MPP_CLOCK_DETAILED`` are integer parameters available by use association, and
        may be summed to turn them both on.
      | While the nesting of clocks is allowed, please note that turning on the non-optional flags on inner clocks has
        certain subtle issues. Turning on ``MPP_CLOCK_SYNC`` on an inner clock may distort outer clock measurements of
        load imbalance. Turning on ``MPP_CLOCK_DETAILED`` will stop detailed measurements on its outer clock, since only
        one detailed clock may be active at one time. Also, detailed clocks only time a certain number of events per
        clock (currently 40000) to conserve memory. If this array overflows, a warning message is printed, and
        subsequent events for this clock are not timed.
      | Timings are done using the ``f90`` standard ``SYSTEM_CLOCK`` intrinsic.
      | The resolution of SYSTEM_CLOCK is often too coarse for use except across large swaths of code. On SGI systems
        this is transparently overloaded with a higher resolution clock made available in a non-portable fortran
        interface made available by ``nsclock.c``. This approach will eventually be extended to other platforms.
      | New behaviour added at the Havana release allows the user to embed profiling calls at varying levels of
        granularity all over the code, and for any particular run, set a threshold of granularity so that finer-grained
        clocks become dormant.
      | The threshold granularity is held in the private module variable ``clock_grain``. This value may be modified by
        the call ``mpp_clock_set_grain``, and affect clocks initiated by subsequent calls to ``mpp_clock_id``. The value
        of ``clock_grain`` is set to an arbitrarily large number initially.
      | Clocks initialized by ``mpp_clock_id`` can set a new optional argument ``grain`` setting their granularity
        level. Clocks check this level against the current value of ``clock_grain``, and are only triggered if they are
        *at or below ("coarser than")* the threshold. Finer-grained clocks are dormant for that run.
      | Note that subsequent changes to ``clock_grain`` do not change the status of already initiated clocks, and that
        if the optional ``grain`` argument is absent, the clock is always triggered. This guarantees backward
        compatibility.

   **INPUT**
      ``grain``
      [integer]

o. .. rubric:: Mpp_sync
      :name: mpp_sync

   ::

      call mpp_sync ( pelist )

   **DESCRIPTION**
      Synchronizes PEs at this point in the execution. If ``pelist`` is omitted all PEs are synchronized. This can be
      expensive on many systems, and should be avoided if possible. Under MPI, we do not call ``MPI_BARRIER``, as you
      might expect. This is because this call can be prohibitively slow on many systems. Instead, we perform the same
      operation as ``mpp_sync_self``, i.e all participating PEs wait for completion of all their outstanding
      non-blocking operations.
      If ``pelist`` is omitted, the context is assumed to be the current pelist. This call implies synchronization
      across the PEs in ``pelist``, or the current pelist if ``pelist`` is absent.
   **INPUT**
      ``pelist``
      [integer, dimension(:)]

p. .. rubric:: Mpp_sync_self
      :name: mpp_sync_self

   **DESCRIPTION**
      ``mpp_transmit`` is implemented as asynchronous ``put/send`` and synchronous ``get/recv``. ``mpp_sync_self``
      guarantees that outstanding asynchronous operations from the calling PE are complete. If ``pelist`` is supplied,
      ``mpp_sync_self`` checks only for outstanding puts to the PEs in ``pelist``.
      If ``pelist`` is omitted, the context is assumed to be the current pelist. This call implies synchronization
      across the PEs in ``pelist``, or the current pelist if ``pelist`` is absent.
   **INPUT**
      ``pelist`` [integer, dimension(:)]

q. .. rubric:: Mpp_malloc
      :name: mpp_malloc

   ::

      call mpp_malloc ( ptr, newlen, len )

   **DESCRIPTION**
      This routine is used on SGI systems when ``mpp_mod`` is invoked in the SHMEM library. It ensures that dynamically
      allocated memory can be used with ``shmem_get`` and ``shmem_put``. This is called *symmetric allocation* and is
      described in the ``intro_shmem`` man page. ``ptr`` is a *Cray pointer* (see the section on portability). The
      operation can be expensive (since it requires a global barrier). We therefore attempt to re-use existing
      allocation whenever possible. Therefore ``len`` and ``ptr`` must have the ``SAVE`` attribute in the calling
      routine, and retain the information about the last call to ``mpp_malloc``. Additional memory is symmetrically
      allocated if and only if ``newlen`` exceeds ``len``.
      This is never required on Cray PVP or MPP systems. While the T3E manpages do talk about symmetric allocation,
      ``mpp_mod`` is coded to remove this restriction.
      It is never required if ``mpp_mod`` is invoked in MPI.
      This call implies synchronization across all PEs.
   **INPUT**
      ``ptr``, a cray pointer, points to a dummy argument in this routine. ``newlen``, the required allocation length
      for the pointer ptr
      [integer]. ``len``, the current allocation (0 if unallocated).
      [integer].

r. .. rubric:: Mpp_set_stack_size
      :name: mpp_set_stack_size

   ::

      call mpp_set_stack_size (n)

   **DESCRIPTION**
      ``mpp_mod`` maintains a private internal array called ``mpp_stack`` for private workspace. This call sets the
      length, in words, of this array.
      The ``mpp_init`` call sets this workspace length to a default of 32768, and this call may be used if a longer
      workspace is needed.
      This call implies synchronization across all PEs.
      This workspace is symmetrically allocated, as required for efficient communication on SGI and Cray MPP systems.
      Since symmetric allocation must be performed by *all* PEs in a job, this call must also be called by all PEs,
      using the same value of ``n``. Calling ``mpp_set_stack_size`` from a subset of PEs, or with unequal argument
      ``n``, may cause the program to hang.
      If any MPP call using ``mpp_stack`` overflows the declared stack array, the program will abort with a message
      specifying the stack length that is required. Many users wonder why, if the required stack length can be computed,
      it cannot also be specified at that point. This cannot be automated because there is no way for the program to
      know if all PEs are present at that call, and with equal values of ``n``. The program must be rerun by the user
      with the correct argument to ``mpp_set_stack_size``, called at an appropriate point in the code where all PEs are
      known to be present.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``n``                                                     | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

Data sets
---------

.. container::

   None.

Error messages
--------------

.. container::

   None.

References
----------

.. container::

   None.

| 

Compiler specifics
------------------

.. container::

   Any module or program unit using ``mpp_mod`` must contain the line
   ::

          use mpp_mod

   The source file for ``mpp_mod`` is ` <ftp://ftp.gfdl.gov/pub/vb/mpp/mpp.F90>`__. Activate the preprocessor flag
   ``-Duse_libSMA`` to invoke the SHMEM library, or ``-Duse_libMPI`` to invoke the MPI library. Global translation of
   preprocessor macros is required. This required the activation of the ``-F`` flag on Cray systems and the
   ``-ftpp -macro_expand`` flags on SGI systems. On non-SGI/Cray systems, please consult the f90 manpage for the
   equivalent flag.
   On Cray PVP systems, *all* routines in a message-passing program must be compiled with ``-a taskcommon``.
   On SGI systems, it is required to use 4-byte integers and 8-byte reals, and the 64-bit ABI (``-i4 -r8 -64 -mips4``).
   It is also required on SGI systems to link the following libraries explicitly: one of ``-lmpi`` and ``-lsma``,
   depending on whether you wish to use the SHMEM or MPI implementations; and ``-lexc``). On Cray systems, all the
   required flags are default.
   On SGI, use MIPSPro f90 7.3.1.2 or higher.
   On Cray, use cf90 3.0.0.0 or higher.
   On either, use the message-passing toolkit MPT 1.2 or higher.
   The declaration ``MPI_INTEGER8`` for 8-byte integers was provided by ``mpp_mod`` because it was absent in early
   releases of the Message Passing Toolkit. It has since been included there, and the declaration in ``mpp_mod``
   commented out. This declaration may need to be reinstated if you get a compiler error from this (i.e you are using a
   superseded version of the MPT).
   By turning on the cpp flag ``-Dtest_mpp`` and compiling ``mpp_mod`` by itself, you may create a test program to
   exercise certain aspects of ``mpp_mod``, e.g
   ::

          f90 -F -Duse_libSMA -Dtest_mpp mpp.F90
          mpprun -n4 a.out

   runs a 4-PE test on a t3e.

| 

Precompiler options
-------------------

.. container::

   While the SHMEM library is currently available only on SGI/Cray systems, ``mpp_mod`` can be used on any other system
   with a standard-compliant f90 compiler and MPI library. SHMEM is now becoming available on other systems as well.
   There are some OS-dependent pre-processor directives that you might need to modify on non-SGI/Cray systems and
   compilers.
   On SGI systems, the ``f90`` standard ``SYSTEM_CLOCK`` intrinsic is overloaded with a non-portable fortran interface
   to a higher-precision clock. This is distributed with the MPP package as ``nsclock.c``. This approach will eventually
   be extended to other platforms, since the resolution of the default clock is often too coarse for our needs.

Test PROGRAM
------------

.. container::

   None.

| 

Notes
-----

.. container::

   None.

| 
