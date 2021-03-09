module mpp_domains_mod
======================

Overview
--------

``mpp_domains_mod`` is a set of simple calls for domain decomposition and domain updates on rectilinear grids. It
requires the module :doc:`./mpp`, upon which it is built.

Scalable implementations of finite-difference codes are generally based on decomposing the model domain into subdomains
that are distributed among processors. These domains will then be obliged to exchange data at their boundaries if data
dependencies are merely nearest-neighbour, or may need to acquire information from the global domain if there are
extended data dependencies, as in the spectral transform. The domain decomposition is a key operation in the development
of parallel codes.

``mpp_domains_mod`` provides a domain decomposition and domain update API for *rectilinear* grids, built on top of the
:doc:`./mpp` API for message passing. Features of ``mpp_domains_mod`` include:

Simple, minimal API, with free access to underlying API for more complicated stuff.

Design toward typical use in climate/weather CFD codes.

Domains
-------

I have assumed that domain decomposition will mainly be in 2 horizontal dimensions, which will in general be the two
fastest-varying indices. There is a separate implementation of 1D decomposition on the fastest-varying index, and 1D
decomposition on the second index, treated as a special case of 2D decomposition, is also possible. We define *domain*
as the grid associated with a *task*. We define the *compute domain* as the set of gridpoints that are computed by a
task, and the *data domain* as the set of points that are required by the task for the calculation. There can in general
be more than 1 task per PE, though often the number of domains is the same as the processor count. We define the *global
domain* as the global computational domain of the entire model (i.e, the same as the computational domain if run on a
single processor). 2D domains are defined using a derived type ``domain2D``, constructed as follows (see comments in
code for more details):

::


   type, public :: domain_axis_spec
       private
       integer :: begin, end, size, max_size
       logical :: is_global
       end type domain_axis_spec
       type, public :: domain1D
       private
       type(domain_axis_spec) :: compute, data, global, active
       logical :: mustputb, mustgetb, mustputf, mustgetf, folded
       type(domain1D), pointer, dimension(:) :: list
       integer :: pe              !PE to which this domain is assigned
       integer :: pos
       end type domain1D
   domaintypes of higher rank can be constructed from type domain1D
   typically we only need 1 and 2D, but could need higher (e.g 3D LES)
   some elements are repeated below if they are needed once per domain
       type, public :: domain2D
       private
       type(domain1D) :: x
       type(domain1D) :: y
       type(domain2D), pointer, dimension(:) :: list
       integer :: pe              !PE to which this domain is assigned
       integer :: pos
       end type domain2D
       type(domain1D), public :: NULL_DOMAIN1D
       type(domain2D), public :: NULL_DOMAIN2D

The ``domain2D`` type contains all the necessary information to define the global, compute and data domains of each
task, as well as the PE associated with the task. The PEs from which remote data may be acquired to update the data
domain are also contained in a linked list of neighbours.

Other modules used
------------------

.. container::

   ::

      mpp_mod

Public interface
----------------

.. container::

   ::

      use mpp_domains_mod [, only:  mpp_define_domains,
                                    mpp_update_domains,
                                    mpp_redistribute,
                                    mpp_global_field,
                                    mpp_global_max,
                                    mpp_global_sum,
                                    operator,
                                    mpp_get_compute_domain,
                                    mpp_get_compute_domains,
                                    mpp_get_data_domain,
                                    mpp_get_global_domain,
                                    mpp_define_layout,
                                    mpp_get_pelist,
                                    mpp_get_layout,
                                    mpp_domains_init,
                                    mpp_domains_set_stack_size,
                                    mpp_domains_exit,
                                    mpp_get_domain_components ]

   mpp_define_domains:
      Set up a domain decomposition.
   mpp_update_domains:
      Halo updates.
   mpp_redistribute:
      Reorganization of distributed global arrays.
   mpp_global_field:
      Fill in a global array from domain-decomposed arrays.
   mpp_global_max:
      Global max/min of domain-decomposed arrays.
   mpp_global_sum:
      Global sum of domain-decomposed arrays.
   operator:
      Equality/inequality operators for domaintypes.
   mpp_get_compute_domain:
      These routines retrieve the axis specifications associated with the compute domains.
   mpp_get_compute_domains:
      Retrieve the entire array of compute domain extents associated with a decomposition.
   mpp_get_data_domain:
      These routines retrieve the axis specifications associated with the data domains.
   mpp_get_global_domain:
      These routines retrieve the axis specifications associated with the global domains.
   mpp_define_layout:
      Retrieve layout associated with a domain decomposition.
   mpp_get_pelist:
      Retrieve list of PEs associated with a domain decomposition.
   mpp_get_layout:
      Retrieve layout associated with a domain decomposition.
   mpp_domains_init:
      Initialize domain decomp package.
   mpp_domains_set_stack_size:
      Set user stack size.
   mpp_domains_exit:
      Exit ``mpp_domains_mod``.
   mpp_get_domain_components:
      Retrieve 1D components of 2D decomposition.

| 

Public data
-----------

.. container::

   None.

Public routines
---------------

a. .. rubric:: Mpp_define_domains
      :name: mpp_define_domains

   ::

      call mpp_define_domains ( global_indices, ndivs, domain, & pelist, flags, halo, extent, maskmap )

   ::

      call mpp_define_domains ( global_indices, layout, domain, pelist, & xflags, yflags, xhalo, yhalo, & xextent, yextent, maskmap, name )

   **DESCRIPTION**
      There are two forms for the ``mpp_define_domains`` call. The 2D version is generally to be used but is built by
      repeated calls to the 1D version, also provided.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``global_indices``                                        | Defines the global domain.                                |
      |                                                           | [integer, dimension(2)]                                   |
      |                                                           | [integer, dimension(4)]                                   |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``ndivs``                                                 | Is the number of domain divisions required.               |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``pelist``                                                | List of PEs to which the domains are to be assigned.      |
      |                                                           | [integer, dimension(0:)]                                  |
      |                                                           | [integer, dimension(0:)]                                  |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``flags``                                                 | An optional flag to pass additional information about the |
      |                                                           | desired domain topology. Useful flags in a 1D             |
      |                                                           | decomposition include ``GLOBAL_DATA_DOMAIN`` and          |
      |                                                           | ``CYCLIC_GLOBAL_DOMAIN``. Flags are integers: multiple    |
      |                                                           | flags may be added together. The flag values are public   |
      |                                                           | parameters available by use association.                  |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``halo``                                                  | Width of the halo.                                        |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``extent``                                                | Normally ``mpp_define_domains`` attempts an even division |
      |                                                           | of the global domain across ``ndivs`` domains. The        |
      |                                                           | ``extent`` array can be used by the user to pass a custom |
      |                                                           | domain division. The ``extent`` array has ``ndivs``       |
      |                                                           | elements and holds the compute domain widths, which       |
      |                                                           | should add up to cover the global domain exactly.         |
      |                                                           | [integer, dimension(0:)]                                  |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``maskmap``                                               | Some divisions may be masked (``maskmap=.FALSE.``) to     |
      |                                                           | exclude them from the computation (e.g for ocean model    |
      |                                                           | domains that are all land). The ``maskmap`` array is      |
      |                                                           | dimensioned ``ndivs`` and contains ``.TRUE.`` values for  |
      |                                                           | any domain that must be *included* in the computation     |
      |                                                           | (default all). The ``pelist`` array length should match   |
      |                                                           | the number of domains included in the computation.        |
      |                                                           | [logical, dimension(0:)]                                  |
      |                                                           | [logical, dimension(:,:)]                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``layout``                                                | [integer, dimension(2)]                                   |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``xflags, yflags``                                        | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``xhalo, yhalo``                                          | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``xextent, yextent``                                      | [integer, dimension(0:)]                                  |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``name``                                                  | [character(len=*)]                                        |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **INPUT/OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``domain``                                                | Holds the resulting domain decomposition.                 |
      |                                                           | [type(domain1D)]                                          |
      |                                                           | [type(domain2D)]                                          |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **NOTE**
      | For example:

      ::

             call mpp_define_domains( (/1,100/), 10, domain, &
                  flags=GLOBAL_DATA_DOMAIN+CYCLIC_GLOBAL_DOMAIN, halo=2 )

      defines 10 compute domains spanning the range [1,100] of the global domain. The compute domains are
      non-overlapping blocks of 10. All the data domains are global, and with a halo of 2 span the range [-1:102]. And
      since the global domain has been declared to be cyclic, ``domain(9)%next => domain(0)`` and
      ``domain(0)%prev => domain(9)``. A field is allocated on the data domain, and computations proceed on the compute
      domain. A call to ` <#mpp_update_domains>`__ would fill in the values in the halo region:

      ::

         call mpp_get_data_domain( domain, isd, ied ) !returns -1 and 102
             call mpp_get_compute_domain( domain, is, ie ) !returns (1,10) on PE 0 ...
             allocate( a(isd:ied) )
             do i = is,ie
                a(i) = <perform computations>
             end do
             call mpp_update_domains( a, domain )

      | The call to ``mpp_update_domains`` fills in the regions outside the compute domain. Since the global domain is
        cyclic, the values at ``i=(-1,0)`` are the same as at ``i=(99,100)``; and ``i=(101,102)`` are the same as
        ``i=(1,2)``.
      | The 2D version is just an extension of this syntax to two dimensions.
      | The 2D version of the above should generally be used in codes, including 1D-decomposed ones, if there is a
        possibility of future evolution toward 2D decomposition. The arguments are similar to the 1D case, except that
        now we have optional arguments ``flags``, ``halo``, ``extent`` and ``maskmap`` along two axes.
      | ``flags`` can now take an additional possible value to fold one or more edges. This is done by using flags
        ``FOLD_WEST_EDGE``, ``FOLD_EAST_EDGE``, ``FOLD_SOUTH_EDGE`` or ``FOLD_NORTH_EDGE``. When a fold exists (e.g
        cylindrical domain), vector fields reverse sign upon crossing the fold. This parity reversal is performed only
        in the vector version of ` <#mpp_update_domains>`__. In addition, shift operations may need to be applied to
        vector fields on staggered grids, also described in the vector interface to ``mpp_update_domains``.
      | ``name`` is the name associated with the decomposition, e.g ``'Ocean model'``. If this argument is present,
        ``mpp_define_domains`` will print the domain decomposition generated to ``stdlog``.
      | Examples:

      ::


         call mpp_define_domains( (/1,100,1,100/), (/2,2/), domain, xhalo=1 )

      will create the following domain layout:

      ::

           
                        |---------|-----------|-----------|-------------|
                        |domain(1)|domain(2)  |domain(3)  |domain(4)    |
         |--------------|---------|-----------|-----------|-------------|
         |Compute domain|1,50,1,50|51,100,1,50|1,50,51,100|51,100,51,100|
         |--------------|---------|-----------|-----------|-------------|
         |Data domain   |0,51,1,50|50,101,1,50|0,51,51,100|50,101,51,100|
         |--------------|---------|-----------|-----------|-------------|

      Again, we allocate arrays on the data domain, perform computations on the compute domain, and call
      ``mpp_update_domains`` to update the halo region.

      If we wished to perfom a 1D decomposition along ``Y`` on the same global domain, we could use:

      ::

         call mpp_define_domains( (/1,100,1,100/), layout=(/4,1/), domain, xhalo=1 )

      This will create the following domain layout:

      ::


                        |----------|-----------|-----------|------------|
                        |domain(1) |domain(2)  |domain(3)  |domain(4)   |
         |--------------|----------|-----------|-----------|------------|
         |Compute domain|1,100,1,25|1,100,26,50|1,100,51,75|1,100,76,100|
         |--------------|----------|-----------|-----------|------------|
         |Data domain   |0,101,1,25|0,101,26,50|0,101,51,75|1,101,76,100|
         |--------------|----------|-----------|-----------|------------|

b. .. rubric:: Mpp_update_domains
      :name: mpp_update_domains

   ::

      call mpp_update_domains ( field, domain, flags )

   ::

      call mpp_update_domains ( fieldx, fieldy, domain, flags, gridtype )

   **DESCRIPTION**
      ``mpp_update_domains`` is used to perform a halo update of a domain-decomposed array on each PE. ``MPP_TYPE_`` can
      be of type ``complex``, ``integer``, ``logical`` or ``real``; of 4-byte or 8-byte kind; of rank up to 5. The
      vector version (with two input data fields) is only present for ``real`` types.
      For 2D domain updates, if there are halos present along both ``x`` and ``y``, we can choose to update one only, by
      specifying ``flags=XUPDATE`` or ``flags=YUPDATE``. In addition, one-sided updates can be performed by setting
      ``flags`` to any combination of ``WUPDATE``, ``EUPDATE``, ``SUPDATE`` and ``NUPDATE``, to update the west, east,
      north and south halos respectively. Any combination of halos may be used by adding the requisite flags, e.g:
      ``flags=XUPDATE+SUPDATE`` or ``flags=EUPDATE+WUPDATE+SUPDATE`` will update the east, west and south halos.
      If a call to ``mpp_update_domains`` involves at least one E-W halo and one N-S halo, the corners involved will
      also be updated, i.e, in the example above, the SE and SW corners will be updated.
      If ``flags`` is not supplied, that is equivalent to ``flags=XUPDATE+YUPDATE``.
      The vector version is passed the ``x`` and ``y`` components of a vector field in tandem, and both are updated upon
      return. They are passed together to treat parity issues on various grids. For example, on a cubic sphere
      projection, the ``x`` and ``y`` components may be interchanged when passing from an equatorial cube face to a
      polar face. For grids with folds, vector components change sign on crossing the fold.
      Special treatment at boundaries such as folds is also required for staggered grids. The following types of
      staggered grids are recognized:
      1) ``AGRID``: values are at grid centers.
      2) ``BGRID_NE``: vector fields are at the NE vertex of a grid cell, i.e: the array elements ``u(i,j)`` and
      ``v(i,j)`` are actually at (i+½,j+½) with respect to the grid centers.
      3) ``BGRID_SW``: vector fields are at the SW vertex of a grid cell, i.e: the array elements ``u(i,j)`` and
      ``v(i,j)`` are actually at (i-½,j-½) with respect to the grid centers.
      4) ``CGRID_NE``: vector fields are at the N and E faces of a grid cell, i.e: the array elements ``u(i,j)`` and
      ``v(i,j)`` are actually at (i+½,j) and (i,j+½) with respect to the grid centers.
      5) ``CGRID_SW``: vector fields are at the S and W faces of a grid cell, i.e: the array elements ``u(i,j)`` and
      ``v(i,j)`` are actually at (i-½,j) and (i,j-½) with respect to the grid centers.
      The gridtypes listed above are all available by use association as integer parameters. The scalar version of
      ``mpp_update_domains`` assumes that the values of a scalar field are always at ``AGRID`` locations, and no special
      boundary treatment is required. If vector fields are at staggered locations, the optional argument ``gridtype``
      must be appropriately set for correct treatment at boundaries.
      It is safe to apply vector field updates to the appropriate arrays irrespective of the domain topology: if the
      topology requires no special treatment of vector fields, specifying ``gridtype`` will do no harm.
      ``mpp_update_domains`` internally buffers the date being sent and received into single messages for efficiency. A
      turnable internal buffer area in memory is provided for this purpose by ``mpp_domains_mod``. The size of this
      buffer area can be set by the user by calling mpp_domains_set_stack_size in :doc:`./mpp_domains`.

c. .. rubric:: Mpp_redistribute
      :name: mpp_redistribute

   ::

      call mpp_redistribute ( domain_in, field_in, domain_out, field_out )

   **DESCRIPTION**
      ``mpp_redistribute`` is used to reorganize a distributed array. ``MPP_TYPE_`` can be of type ``integer``,
      ``complex``, or ``real``; of 4-byte or 8-byte kind; of rank up to 5.
   **INPUT**
      ============ ================================================================
      ``field_in`` ``field_in`` is dimensioned on the data domain of ``domain_in``.
      ============ ================================================================

   **OUTPUT**
      ============= ===================================================
      ``field_out`` ``field_out`` on the data domain of ``domain_out``.
      ============= ===================================================

d. .. rubric:: Mpp_global_field
      :name: mpp_global_field

   ::

      call mpp_global_field ( domain, local, global, flags )

   **DESCRIPTION**
      ``mpp_global_field`` is used to get an entire domain-decomposed array on each PE. ``MPP_TYPE_`` can be of type
      ``complex``, ``integer``, ``logical`` or ``real``; of 4-byte or 8-byte kind; of rank up to 5.
      All PEs in a domain decomposition must call ``mpp_global_field``, and each will have a complete global field at
      the end. Please note that a global array of rank 3 or higher could occupy a lot of memory.
   **INPUT**
      ========== =====================================================================================================
      ``domain`` 
      ``local``  ``local`` is dimensioned on either the compute domain or the data domain of ``domain``.
      ``flags``  ``flags`` can be given the value ``XONLY`` or ``YONLY``, to specify a globalization on one axis only.
      ========== =====================================================================================================

   **OUTPUT**
      ========== =============================================================
      ``global`` ``global`` is dimensioned on the corresponding global domain.
      ========== =============================================================

e. .. rubric:: Mpp_global_max
      :name: mpp_global_max

   ::

       
      mpp_global_max ( domain, field, locus )

   **DESCRIPTION**
      ``mpp_global_max`` is used to get the maximum value of a domain-decomposed array on each PE. ``MPP_TYPE_`` can be
      of type ``integer`` or ``real``; of 4-byte or 8-byte kind; of rank up to 5. The dimension of ``locus`` must equal
      the rank of ``field``.
      All PEs in a domain decomposition must call ``mpp_global_max``, and each will have the result upon exit.
      The function ``mpp_global_min``, with an identical syntax. is also available.
   **INPUT**
      ========== =======================================================================================
      ``domain`` 
      ``field``  ``field`` is dimensioned on either the compute domain or the data domain of ``domain``.
      ========== =======================================================================================

   **OUTPUT**
      +-----------+---------------------------------------------------------------------------------------------------------+
      | ``locus`` | ``locus``, if present, can be used to retrieve the location of the maximum (as in the ``MAXLOC``        |
      |           | intrinsic of f90).                                                                                      |
      +-----------+---------------------------------------------------------------------------------------------------------+

f. .. rubric:: Mpp_global_sum
      :name: mpp_global_sum

   ::

      call mpp_global_sum ( domain, field, flags )

   **DESCRIPTION**
      ``mpp_global_sum`` is used to get the sum of a domain-decomposed array on each PE. ``MPP_TYPE_`` can be of type
      ``integer``, ``complex``, or ``real``; of 4-byte or 8-byte kind; of rank up to 5.
   **INPUT**
      +------------+--------------------------------------------------------------------------------------------------------+
      | ``domain`` |                                                                                                        |
      +------------+--------------------------------------------------------------------------------------------------------+
      | ``field``  | ``field`` is dimensioned on either the compute domain or the data domain of ``domain``.                |
      +------------+--------------------------------------------------------------------------------------------------------+
      | ``flags``  | ``flags``, if present, must have the value ``BITWISE_EXACT_SUM``. This produces a sum that is          |
      |            | guaranteed to produce the identical result irrespective of how the domain is decomposed. This method   |
      |            | does the sum first along the ranks beyond 2, and then calls ` <#mpp_global_field>`__ to produce a      |
      |            | global 2D array which is then summed. The default method, which is considerably faster, does a local   |
      |            | sum followed by mpp_sum in :doc:`./mpp` across the domain decomposition.                               |
      +------------+--------------------------------------------------------------------------------------------------------+

   **NOTE**
      All PEs in a domain decomposition must call ``mpp_global_sum``, and each will have the result upon exit.

g. .. rubric:: Operator
      :name: operator

   **DESCRIPTION**
      | The module provides public operators to check for equality/inequality of domaintypes, e.g:

      ::

             type(domain1D) :: a, b
             type(domain2D) :: c, d
             ...
             if( a.NE.b )then
                 ...
             end if
             if( c==d )then
                 ...
             end if

      Domains are considered equal if and only if the start and end indices of each of their component global, data and
      compute domains are equal.

h. .. rubric:: Mpp_get_compute_domain
      :name: mpp_get_compute_domain

   ::

      call mpp_get_compute_domain 

   **DESCRIPTION**
      The domain is a derived type with private elements. These routines retrieve the axis specifications associated
      with the compute domains The 2D version of these is a simple extension of 1D.

i. .. rubric:: Mpp_get_compute_domains
      :name: mpp_get_compute_domains

   ::

      call mpp_get_compute_domains ( domain, xbegin, xend, xsize, & ybegin, yend, ysize )

   **DESCRIPTION**
      Retrieve the entire array of compute domain extents associated with a decomposition.
   **INPUT**
      ``domain``
   **OUTPUT**
      ``xbegin,ybegin``, ``xend,yend``, ``xsize,ysize``

j. .. rubric:: Mpp_get_data_domain
      :name: mpp_get_data_domain

   ::

      call mpp_get_data_domain 

   **DESCRIPTION**
      The domain is a derived type with private elements. These routines retrieve the axis specifications associated
      with the data domains. The 2D version of these is a simple extension of 1D.

k. .. rubric:: Mpp_get_global_domain
      :name: mpp_get_global_domain

   ::

      call mpp_get_global_domain 

   **DESCRIPTION**
      The domain is a derived type with private elements. These routines retrieve the axis specifications associated
      with the global domains. The 2D version of these is a simple extension of 1D.

l. .. rubric:: Mpp_define_layout
      :name: mpp_define_layout

   ::

      call mpp_define_layout ( global_indices, ndivs, layout )

   **DESCRIPTION**
      Given a global 2D domain and the number of divisions in the decomposition (``ndivs``: usually the PE count unless
      some domains are masked) this calls returns a 2D domain layout.
      By default, ``mpp_define_layout`` will attempt to divide the 2D index space into domains that maintain the aspect
      ratio of the global domain. If this cannot be done, the algorithm favours domains that are longer in ``x`` than
      ``y``, a preference that could improve vector performance.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``global_indices``                                        | [integer, dimension(4)]                                   |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``ndivs``                                                 | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``layout``                                                | [integer, dimension(2)]                                   |
      +-----------------------------------------------------------+-----------------------------------------------------------+

m. .. rubric:: Mpp_get_pelist
      :name: mpp_get_pelist

   **DESCRIPTION**
      The 1D version of this call returns an array of the PEs assigned to this 1D domain decomposition. In addition the
      optional argument ``pos`` may be used to retrieve the 0-based position of the domain local to the calling PE, i.e
      ``domain%list(pos)%pe`` is the local PE, as returned by mpp_pe in :doc:`./mpp`. The 2D version of this call is
      identical to 1D version.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``domain``                                                | [type(domain1D)]                                          |
      |                                                           | [type(domain2D)]                                          |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``pelist``                                                | [integer, dimension(:)]                                   |
      |                                                           | [integer, dimension(:)]                                   |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``pos``                                                   | [integer]                                                 |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

n. .. rubric:: Mpp_get_layout
      :name: mpp_get_layout

   ::

      call mpp_get_layout ( domain, layout )

   **DESCRIPTION**
      The 1D version of this call returns the number of divisions that was assigned to this decomposition axis. The 2D
      version of this call returns an array of dimension 2 holding the results on two axes.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``domain``                                                | [type(domain1D)]                                          |
      |                                                           | [type(domain2D)]                                          |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``layout``                                                | [integer]                                                 |
      |                                                           | [integer, dimension(2)]                                   |
      +-----------------------------------------------------------+-----------------------------------------------------------+

o. .. rubric:: Mpp_domains_init
      :name: mpp_domains_init

   ::

      call mpp_domains_init (flags)

   **DESCRIPTION**
      Called to initialize the ``mpp_domains_mod`` package.
      ``flags`` can be set to ``MPP_VERBOSE`` to have ``mpp_domains_mod`` keep you informed of what it's up to.
      ``MPP_DEBUG`` returns even more information for debugging.
      ``mpp_domains_init`` will call ``mpp_init``, to make sure :doc:`./mpp` is initialized. (Repeated calls to
      ``mpp_init`` do no harm, so don't worry if you already called it).
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``flags``                                                 | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

p. .. rubric:: Mpp_domains_set_stack_size
      :name: mpp_domains_set_stack_size

   ::

      call mpp_domains_set_stack_size (n)

   **DESCRIPTION**
      This sets the size of an array that is used for internal storage by ``mpp_domains``. This array is used, for
      instance, to buffer the data sent and received in halo updates.
      This call has implied global synchronization. It should be placed somewhere where all PEs can call it.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``n``                                                     | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

q. .. rubric:: Mpp_domains_exit
      :name: mpp_domains_exit

   ::

      call mpp_domains_exit ()

   **DESCRIPTION**
      Serves no particular purpose, but is provided should you require to re-initialize ``mpp_domains_mod``, for some
      odd reason.

r. .. rubric:: Mpp_get_domain_components
      :name: mpp_get_domain_components

   ::

      call mpp_get_domain_components ( domain, x, y )

   **DESCRIPTION**
      It is sometime necessary to have direct recourse to the domain1D types that compose a domain2D object. This call
      retrieves them.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``domain``                                                | [type(domain2D)]                                          |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``x,y``                                                   | [type(domain1D)]                                          |
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

   Any module or program unit using ``mpp_domains_mod`` must contain the line
   ::

           use mpp_domains_mod

   ``mpp_domains_mod`` ``use``\ s :doc:`./mpp`, and therefore is subject to the compiling and linking requirements of
   :doc:`./mpp`.

| 

Precompiler options
-------------------

.. container::

   ``mpp_domains_mod`` uses standard f90, and has no special requirements. There are some OS-dependent pre-processor
   directives that you might need to modify on non-SGI/Cray systems and compilers. The portability, as described in
   :doc:`./mpp` obviously is a constraint, since this module is built on top of it. Contact me, Balaji, SGI/GFDL, with
   questions.

| 

Loader options
--------------

The source consists of the main source file and also requires the following include files: GFDL users can check it out
of the main CVS repository as part of the CVS module. The current public tag is . External users can download the latest
package . Public access to the GFDL CVS repository will soon be made available.

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
