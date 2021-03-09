module mpp_io_mod
=================

Overview
--------

``mpp_io_mod``, is a set of simple calls for parallel I/O on distributed systems. It is geared toward the writing of
data in netCDF format. It requires the modules :doc:`./mpp_domains` and :doc:`./mpp`, upon which it is built.

In massively parallel environments, an often difficult problem is the reading and writing of data to files on disk.
MPI-IO and MPI-2 IO are moving toward providing this capability, but are currently not widely implemented. Further, it
is a rather abstruse API. ``mpp_io_mod`` is an attempt at a simple API encompassing a certain variety of the I/O tasks
that will be required. It does not attempt to be an all-encompassing standard such as MPI, however, it can be
implemented in MPI if so desired. It is equally simple to add parallel I/O capability to ``mpp_io_mod`` based on
vendor-specific APIs while providing a layer of insulation for user codes.

The ``mpp_io_mod`` parallel I/O API built on top of the :doc:`./mpp_domains` and :doc:`./mpp` API for domain
decomposition and message passing. Features of ``mpp_io_mod`` include:

#. Simple, minimal API, with free access to underlying API for more complicated stuff.
#. Self-describing files: comprehensive header information (metadata) in the file itself.
#. Strong focus on performance of parallel write: the climate models for which it is designed typically read a minimal
   amount of data (typically at the beginning of the run); but on the other hand, tend to write copious amounts of data
   during the run. An interface for reading is also supplied, but its performance has not yet been optimized.
#. Integrated netCDF capability: `netCDF <http://www.unidata.ucar.edu/packages/netcdf/>`__ is a data format widely used
   in the climate/weather modeling community. netCDF is considered the principal medium of data storage for
   ``mpp_io_mod``. But I provide a raw unformatted fortran I/O capability in case netCDF is not an option, either due to
   unavailability, inappropriateness, or poor performance.
#. May require off-line post-processing: a tool for this purpose, ``mppnccombine``, is available. GFDL users may use
   ``~hnv/pub/mppnccombine``. Outside users may obtain the source
   `here <ftp://ftp.gfdl.gov/perm/hnv/mpp/mppnccombine.c>`__. It can be compiled on any C compiler and linked with the
   netCDF library. The program is free and is covered by the `GPL license <ftp://ftp.gfdl.gov/perm/hnv/mpp/LICENSE>`__.

The internal representation of the data being written out is assumed be the default real type, which can be 4 or 8-byte.
Time data is always written as 8-bytes to avoid overflow on climatic time scales in units of seconds.

The I/O activity critical to performance in the models for which ``mpp_io_mod`` is designed is typically the writing of
large datasets on a model grid volume produced at intervals during a run. Consider a 3D grid volume, where model arrays
are stored as ``(i,j,k)``. The domain decomposition is typically along ``i`` or ``j``: thus to store data to disk as a
global volume, the distributed chunks of data have to be seen as non-contiguous. If we attempt to have all PEs write
this data into a single file, performance can be seriously compromised because of the data reordering that will be
required. Possible options are to have one PE acquire all the data and write it out, or to have all the PEs write
independent files, which are recombined offline. These three modes of operation are described in the ``mpp_io_mod``
terminology in terms of two parameters, *threading* and *fileset*, as follows:

| *Single-threaded I/O:* a single PE acquires all the data and writes it out.
| *Multi-threaded, single-fileset I/O:* many PEs write to a single file.
| *Multi-threaded, multi-fileset I/O:* many PEs write to independent files. This is also called *distributed I/O*.

| 

The middle option is the most difficult to achieve performance. The choice of one of these modes is made when a file is
opened for I/O, in ` <#mpp_open>`__.

.. _section-1:

A requirement of the design of ``mpp_io_mod`` is that the file must be entirely self-describing: comprehensive header
information describing its contents is present in the header of every file. The header information follows the model of
netCDF. Variables in the file are divided into *axes* and *fields*. An axis describes a co-ordinate variable, e.g
``x,y,z,t``. A field consists of data in the space described by the axes. An axis is described in ``mpp_io_mod`` using
the defined type ``axistype``:

::


   type, public :: axistype
     sequence
     character(len=128) :: name
     character(len=128) :: units
     character(len=256) :: longname
     character(len=8) :: cartesian
     integer :: len
     integer :: sense           !+/-1, depth or height?
     type(domain1D), pointer :: domain
     real, dimension(:), pointer :: data
     integer :: id, did
     integer :: type  ! external NetCDF type format for axis data
     integer :: natt
     type(atttype), pointer :: Att(:) ! axis attributes
   end type axistype

A field is described using the type ``fieldtype``:

::


   type, public :: fieldtype
     sequence
     character(len=128) :: name
     character(len=128) :: units
     character(len=256) :: longname
     real :: min, max, missing, fill, scale, add
     integer :: pack
     type(axistype), dimension(:), pointer :: axes
     integer, dimension(:), pointer :: size
     integer :: time_axis_index
     integer :: id
     integer :: type ! external NetCDF format for field data
     integer :: natt, ndim
     type(atttype), pointer :: Att(:) ! field metadata
   end type fieldtype

An attribute (global, field or axis) is described using the ``atttype``:

::


   type, public :: atttype
     sequence
     integer :: type, len
     character(len=128) :: name
     character(len=256)  :: catt
     real(FLOAT_KIND), pointer :: fatt(:)
   end type atttype

This default set of field attributes corresponds closely to various conventions established for netCDF files. The
``pack`` attribute of a field defines whether or not a field is to be packed on output. Allowed values of ``pack`` are
1,2,4 and 8. The value of ``pack`` is the number of variables written into 8 bytes. In typical use, we write 4-byte
reals to netCDF output; thus the default value of ``pack`` is 2. For ``pack`` = 4 or 8, packing uses a simple-minded
linear scaling scheme using the ``scale`` and ``add`` attributes. There is thus likely to be a significant loss of
dynamic range with packing. When a field is declared to be packed, the ``missing`` and ``fill`` attributes, if supplied,
are packed also.

Please note that the pack values are the same even if the default real is 4 bytes, i.e ``PACK=1`` still follows the
definition above and writes out 8 bytes.

A set of *attributes* for each variable is also available. The variable definitions and attribute information is
written/read by calling ` <#mpp_write_meta>`__ or ` <#mpp_read_meta>`__. A typical calling sequence for writing data
might be:

::


   ...
     type(domain2D), dimension(:), allocatable, target :: domain
     type(fieldtype) :: field
     type(axistype) :: x, y, z, t
   ...
     call mpp_define_domains( (/1,nx,1,ny/), domain )
     allocate( a(domain(pe)%x%data%start_index:domain(pe)%x%data%end_index, &
                 domain(pe)%y%data%start_index:domain(pe)%y%data%end_index,nz) )
   ...
     call mpp_write_meta( unit, x, 'X', 'km', 'X distance', &
         domain=domain(pe)%x, data=(/(float(i),i=1,nx)/) )
     call mpp_write_meta( unit, y, 'Y', 'km', 'Y distance', &
         domain=domain(pe)%y, data=(/(float(i),i=1,ny)/) )
     call mpp_write_meta( unit, z, 'Z', 'km', 'Z distance', &
         data=(/(float(i),i=1,nz)/) )
     call mpp_write_meta( unit, t, 'Time', 'second', 'Time' )

     call mpp_write_meta( unit, field, (/x,y,z,t/), 'a', '(m/s)', AAA', &
         missing=-1e36 )
   ...
     call mpp_write( unit, x )
     call mpp_write( unit, y )
     call mpp_write( unit, z )
   ...

In this example, ``x`` and ``y`` have been declared as distributed axes, since a domain decomposition has been
associated. ``z`` and ``t`` are undistributed axes. ``t`` is known to be a *record* axis (netCDF terminology) since we
do not allocate the ``data`` element of the ``axistype``. *Only one record axis may be associated with a file.* The call
to ` <#mpp_write_meta>`__ initializes the axes, and associates a unique variable ID with each axis. The call to
``mpp_write_meta`` with argument ``field`` declared ``field`` to be a 4D variable that is a function of ``(x,y,z,t)``,
and a unique variable ID is associated with it. A 3D field will be written at each call to ``mpp_write(field)``.

The data to any variable, including axes, is written by ``mpp_write``.

Any additional attributes of variables can be added through subsequent ``mpp_write_meta`` calls, using the variable ID
as a handle. *Global* attributes, associated with the dataset as a whole, can also be written thus. See the
` <#mpp_write_meta>`__ call syntax below for further details.

You cannot interleave calls to ``mpp_write`` and ``mpp_write_meta``: the first call to ``mpp_write`` implies that
metadata specification is complete.

A typical calling sequence for reading data might be:

::


   ...
     integer :: unit, natt, nvar, ntime
     type(domain2D), dimension(:), allocatable, target :: domain
     type(fieldtype), allocatable, dimension(:) :: fields
     type(atttype), allocatable, dimension(:) :: global_atts
     real, allocatable, dimension(:) :: times
   ...
     call mpp_define_domains( (/1,nx,1,ny/), domain )

     call mpp_read_meta(unit)
     call mpp_get_info(unit,natt,nvar,ntime)
     allocate(global_atts(natt))
     call mpp_get_atts(unit,global_atts)
     allocate(fields(nvar))
     call mpp_get_vars(unit, fields)
     allocate(times(ntime))
     call mpp_get_times(unit, times)

     allocate( a(domain(pe)%x%data%start_index:domain(pe)%x%data%end_index, &
                 domain(pe)%y%data%start_index:domain(pe)%y%data%end_index,nz) )
   ...
     do i=1, nvar
       if (fields(i)%name == 'a')  call mpp_read(unit,fields(i),domain(pe), a,
                                                 tindex)
     enddo
   ...

In this example, the data are distributed as in the previous example. The call to ` <#mpp_read_meta>`__ initializes all
of the metadata associated with the file, including global attributes, variable attributes and non-record dimension
data. The call to ``mpp_get_info`` returns the number of global attributes (``natt``), variables (``nvar``) and time
levels (``ntime``) associated with the file identified by a unique ID (``unit``). ``mpp_get_atts`` returns all global
attributes for the file in the derived type ``atttype(natt)``. ``mpp_get_vars`` returns variable types
(``fieldtype(nvar)``). Since the record dimension data are not allocated for calls to ` <#mpp_write>`__, a separate call
to ``mpp_get_times`` is required to access record dimension data. Subsequent calls to ``mpp_read`` return the field data
arrays corresponding to the fieldtype. The ``domain`` type is an optional argument. If ``domain`` is omitted, the
incoming field array should be dimensioned for the global domain, otherwise, the field data is assigned to the
computational domain of a local array.

*Multi-fileset* reads are not supported with ``mpp_read``.

Other modules used
------------------

.. container::

   ::

              mpp_mod
      mpp_domains_mod

Public interface
----------------

.. container::

   ::

      use mpp_io_mod [, only:  mpp_write_meta,
                               mpp_write,
                               mpp_read,
                               mpp_get_atts,
                               mpp_io_init,
                               mpp_io_exit,
                               mpp_open,
                               mpp_close,
                               mpp_read_meta,
                               mpp_get_info,
                               mpp_get_times,
                               mpp_flush,
                               mpp_get_ncid ]

   mpp_write_meta:
      Write metadata.
   mpp_write:
      Write to an open file.
   mpp_read:
      Read from an open file.
   mpp_get_atts:
      Get file global metdata.
   mpp_io_init:
      Initialize ``mpp_io_mod``.
   mpp_io_exit:
      Exit ``mpp_io_mod``.
   mpp_open:
      Open a file for parallel I/O.
   mpp_close:
      Close an open file.
   mpp_read_meta:
      Read metadata.
   mpp_get_info:
      Get some general information about a file.
   mpp_get_times:
      Get file time data.
   mpp_flush:
      Flush I/O buffers to disk.
   mpp_get_ncid:
      Get netCDF ID of an open file.

| 

Public data
-----------

.. container::

   None.

Public routines
---------------

a. .. rubric:: Mpp_write_meta
      :name: mpp_write_meta

   ::

      call mpp_write_meta ( unit, axis, name, units, longname, cartesian, sense, domain, data )

   ::

      call mpp_write_meta ( unit, field, axes, name, units, longname, min, max, missing, fill, scale, add, pack )

   ::

      call mpp_write_meta ( unit, id, name, rval=rval, pack=pack )

   ::

      call mpp_write_meta ( unit, id, name, ival=ival )

   ::

      call mpp_write_meta ( unit, id, name, cval=cval )

   ::

      call mpp_write_meta ( unit, name, rval=rval, pack=pack )

   ::

      call mpp_write_meta ( unit, name, ival=ival )

   ::

      call mpp_write_meta ( unit, name, cval=cval )

   **DESCRIPTION**
      This routine is used to write the metadata describing the contents of a file being written. Each file can contain
      any number of fields, which are functions of 0-3 space axes and 0-1 time axes. (Only one time axis can be defined
      per file). The basic metadata defined above for ``axistype`` and ``fieldtype`` are written in the first two forms
      of the call shown below. These calls will associate a unique variable ID with each variable (axis or field). These
      can be used to attach any other real, integer or character attribute to a variable. The last form is used to
      define a *global* real, integer or character attribute that applies to the dataset as a whole.
   **INPUT**
      ``unit``, ``name``, ``units``, ``longname``, ``cartesian``, ``sense``, ``domain``, ``data``, ``min, max``,
      ``missing``, ``fill``, ``scale``, ``add``, ``pack``, ``id``, ``cval``, ``ival``, ``rval``
   **OUTPUT**
      ``axis``, ``field``
   **NOTE**
      The first form defines a time or space axis. Metadata corresponding to the type above are written to the file on
      <unit>. A unique ID for subsequen references to this axis is returned in axis%id. If the <domain> element is
      present, this is recognized as a distributed data axis and domain decomposition information is also written if
      required (the domain decomposition info is required for multi-fileset multi-threaded I/O). If the <data> element
      is allocated, it is considered to be a space axis, otherwise it is a time axis with an unlimited dimension. Only
      one time axis is allowed per file.
      The second form defines a field. Metadata corresponding to the type above are written to the file on <unit>. A
      unique ID for subsequen references to this field is returned in field%id. At least one axis must be associated, 0D
      variables are not considered. mpp_write_meta must previously have been called on all axes associated with this
      field.
      The third form (3 - 5) defines metadata associated with a previously defined axis or field, identified to
      mpp_write_meta by its unique ID <id>. The attribute is named <name> and can take on a real, integer or character
      value. <rval> and <ival> can be scalar or 1D arrays. This need not be called for attributes already contained in
      the type.
      The last form (6 - 8) defines global metadata associated with the file as a whole. The attribute is named <name>
      and can take on a real, integer or character value. <rval> and <ival> can be scalar or 1D arrays.
      Note that ``mpp_write_meta`` is expecting axis data on the *global* domain even if it is a domain-decomposed axis.
      You cannot interleave calls to ``mpp_write`` and ``mpp_write_meta``: the first call to ``mpp_write`` implies that
      metadata specification is complete.

b. .. rubric:: Mpp_write
      :name: mpp_write

   ::

       
      mpp_write ( unit, axis )

   ::

       
      mpp_write ( unit, field, data, tstamp )

   ::

       
      mpp_write ( unit, field, domain, data, tstamp )

   **DESCRIPTION**
      ``mpp_write`` is used to write data to the file on an I/O unit using the file parameters supplied by
      ` <#mpp_open>`__. Axis and field definitions must have previously been written to the file using
      ` <#mpp_write_meta>`__. There are three forms of ``mpp_write``, one to write axis data, one to write distributed
      field data, and one to write non-distributed field data. *Distributed* data refer to arrays whose two
      fastest-varying indices are domain-decomposed. Distributed data must be 2D or 3D (in space). Non-distributed data
      can be 0-3D.
      The ``data`` argument for distributed data is expected by ``mpp_write`` to contain data specified on the *data*
      domain, and will write the data belonging to the *compute* domain, fetching or sending data as required by the
      parallel I/O mode specified in the ``mpp_open`` call. This is consistent with our definition in
      :doc:`./mpp_domains`, where all arrays are expected to be dimensioned on the data domain, and all operations
      performed on the compute domain.
      The type of the ``data`` argument must be a *default real*, which can be 4 or 8 byte.
   **INPUT**
      +------------+--------------------------------------------------------------------------------------------------------+
      | ``tstamp`` | ``tstamp`` is an optional argument. It is to be omitted if the field was defined not to be a function  |
      |            | of time. Results are unpredictable if the argument is supplied for a time- independent field, or       |
      |            | omitted for a time-dependent field. Repeated writes of a time-independent field are also not           |
      |            | recommended. One time level of one field is written per call. tstamp must be an 8-byte real, even if   |
      |            | the default real type is 4-byte.                                                                       |
      +------------+--------------------------------------------------------------------------------------------------------+

   **NOTE**
      | The type of write performed by ``mpp_write`` depends on the file characteristics on the I/O unit specified at
        the ` <#mpp_open>`__ call. Specifically, the format of the output data (e.g netCDF or IEEE), the ``threading``
        and ``fileset`` flags, etc., can be changed there, and require no changes to the ``mpp_write`` calls.
      | Packing is currently not implemented for non-netCDF files, and the ``pack`` attribute is ignored. On netCDF
        files, ``NF_DOUBLE``\ s (8-byte IEEE floating point numbers) are written for ``pack``\ =1 and ``NF_FLOAT``\ s
        for ``pack``\ =2. (``pack``\ =2 gives the customary and default behaviour). We write ``NF_SHORT``\ s (2-byte
        integers) for ``pack=4``, or ``NF_BYTE``\ s (1-byte integers) for ``pack=8``. Integer scaling is done using the
        ``scale`` and ``add`` attributes at ``pack``\ =4 or 8, satisfying the relation

      ::

             data = packed_data*scale + add

      | ``NOTE: mpp_write`` does not check to see if the scaled data in fact fits into the dynamic range implied by the
        specified packing. It is incumbent on the user to supply correct scaling attributes.
      | You cannot interleave calls to ``mpp_write`` and ``mpp_write_meta``: the first call to ``mpp_write`` implies
        that metadata specification is complete.

c. .. rubric:: Mpp_read
      :name: mpp_read

   ::

      call mpp_read ( unit, field, data, time_index )

   ::

      call mpp_read ( unit, field, domain, data, time_index )

   **DESCRIPTION**
      ``mpp_read`` is used to read data to the file on an I/O unit using the file parameters supplied by
      ` <#mpp_open>`__. There are two forms of ``mpp_read``, one to read distributed field data, and one to read
      non-distributed field data. *Distributed* data refer to arrays whose two fastest-varying indices are
      domain-decomposed. Distributed data must be 2D or 3D (in space). Non-distributed data can be 0-3D.
      The ``data`` argument for distributed data is expected by ``mpp_read`` to contain data specified on the *data*
      domain, and will read the data belonging to the *compute* domain, fetching data as required by the parallel I/O
      mode specified in the ``mpp_open`` call. This is consistent with our definition in :doc:`./mpp_domains`, where all
      arrays are expected to be dimensioned on the data domain, and all operations performed on the compute domain.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``unit``                                                  | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``field``                                                 | [type(fieldtype)]                                         |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``domain``                                                |                                                           |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``time_index``                                            | time_index is an optional argument. It is to be omitted   |
      |                                                           | if the field was defined not to be a function of time.    |
      |                                                           | Results are unpredictable if the argument is supplied for |
      |                                                           | a time- independent field, or omitted for a               |
      |                                                           | time-dependent field.                                     |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **INPUT/OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``data``                                                  | [real, dimension(:,:,:)]                                  |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **NOTE**
      The type of read performed by ``mpp_read`` depends on the file characteristics on the I/O unit specified at the
      ` <#mpp_open>`__ call. Specifically, the format of the input data (e.g netCDF or IEEE) and the ``threading``
      flags, etc., can be changed there, and require no changes to the ``mpp_read`` calls. (``fileset`` = MPP_MULTI is
      not supported by ``mpp_read``; IEEE is currently not supported).
      Packed variables are unpacked using the ``scale`` and ``add`` attributes.
      ``mpp_read_meta`` must be called prior to calling ``mpp_read.``

d. .. rubric:: Mpp_get_atts
      :name: mpp_get_atts

   ::

      call mpp_get_atts ( unit, global_atts)

   **DESCRIPTION**
      Get file global metdata.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``unit``                                                  | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``global_atts``                                           | [atttype, dimension(:)]                                   |
      +-----------------------------------------------------------+-----------------------------------------------------------+

e. .. rubric:: Mpp_io_init
      :name: mpp_io_init

   ::

      call mpp_io_init ( flags, maxunit )

   **DESCRIPTION**
      Called to initialize the ``mpp_io_mod`` package. Sets the range of valid fortran units and initializes the
      ``mpp_file`` array of ``type(filetype)``. ``mpp_io_init`` will call ``mpp_init`` and ``mpp_domains_init``, to make
      sure its parent modules have been initialized. (Repeated calls to the ``init`` routines do no harm, so don't worry
      if you already called it).
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``flags``                                                 | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``maxunit``                                               | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

f. .. rubric:: Mpp_io_exit
      :name: mpp_io_exit

   ::

      call mpp_io_exit ()

   **DESCRIPTION**
      It is recommended, though not at present required, that you call this near the end of a run. This will close all
      open files that were opened with ` <#mpp_open>`__. Files opened otherwise are not affected.

g. .. rubric:: Mpp_open
      :name: mpp_open

   ::

      call mpp_open ( unit, file, action, form, access, threading, fileset, iospec, nohdrs, recl, pelist )

   **DESCRIPTION**
      Open a file for parallel I/O.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``file``                                                  | file is the filename: REQUIRED we append .nc to filename  |
      |                                                           | if it is a netCDF file we append .<pppp> to filename if   |
      |                                                           | fileset is private (pppp is PE number)                    |
      |                                                           | [character(len=*)]                                        |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``action``                                                | action is one of MPP_RDONLY, MPP_APPEND, MPP_WRONLY or    |
      |                                                           | MPP_OVERWR.                                               |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``form``                                                  | form is one of MPP_ASCII: formatted read/write            |
      |                                                           | MPP_NATIVE: unformatted read/write with no conversion     |
      |                                                           | MPP_IEEE32: unformatted read/write with conversion to     |
      |                                                           | IEEE32 MPP_NETCDF: unformatted read/write with conversion |
      |                                                           | to netCDF                                                 |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``access``                                                | access is one of MPP_SEQUENTIAL or MPP_DIRECT (ignored    |
      |                                                           | for netCDF). RECL argument is REQUIRED for direct access  |
      |                                                           | IO.                                                       |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``threading``                                             | threading is one of MPP_SINGLE or MPP_MULTI               |
      |                                                           | single-threaded IO in a multi-PE run is done by PE0.      |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``fileset``                                               | fileset is one of MPP_MULTI and MPP_SINGLE fileset is     |
      |                                                           | only used for multi-threaded I/O if all I/O PEs in        |
      |                                                           | <pelist> use a single fileset, they write to the same     |
      |                                                           | file if all I/O PEs in <pelist> use a multi fileset, they |
      |                                                           | each write an independent file                            |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``pelist``                                                | pelist is the list of I/O PEs (currently ALL).            |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``recl``                                                  | recl is the record length in bytes.                       |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``iospec``                                                | iospec is a system hint for I/O organization, e.g         |
      |                                                           | assign(1) on SGI/Cray systems.                            |
      |                                                           | [character(len=*)]                                        |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``nohdrs``                                                | nohdrs has no effect when action=MPP_RDONLY|MPP_APPEND or |
      |                                                           | when form=MPP_NETCDF                                      |
      |                                                           | [logical]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``unit``                                                  | unit is intent(OUT): always \_returned_by\_ mpp_open().   |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **NOTE**
      | The integer parameters to be passed as flags (``MPP_RDONLY``, etc) are all made available by use association.
        The ``unit`` returned by ``mpp_open`` is guaranteed unique. For non-netCDF I/O it is a valid fortran unit number
        and fortran I/O can be directly called on the file.
      | ``MPP_WRONLY`` will guarantee that existing files named ``file`` will not be clobbered. ``MPP_OVERWR`` allows
        overwriting of files.
      | Files opened read-only by many processors will give each processor an independent pointer into the file, i.e:

      ::

               namelist / nml / ...
            ...
               call mpp_open( unit, 'input.nml', action=MPP_RDONLY )
               read(unit,nml)

      | will result in each PE independently reading the same namelist.
      | Metadata identifying the file and the version of ``mpp_io_mod`` are written to a file that is opened
        ``MPP_WRONLY`` or ``MPP_OVERWR``. If this is a multi-file set, and an additional global attribute
        ``NumFilesInSet`` is written to be used by post-processing software.
      | If ``nohdrs=.TRUE.`` all calls to write attributes will return successfully *without* performing any writes to
        the file. The default is ``.FALSE.``.
      | For netCDF files, headers are always written even if ``nohdrs=.TRUE.`` The string ``iospec`` is passed to the OS
        to characterize the I/O to be performed on the file opened on ``unit``. This is typically used for I/O
        optimization. For example, the FFIO layer on SGI/Cray systems can be used for controlling synchronicity of reads
        and writes, buffering of data between user space and disk for I/O optimization, striping across multiple disk
        partitions, automatic data conversion and the like (``man intro_ffio``). All these actions are controlled
        through the ``assign`` command. For example, to specify asynchronous caching of data going to a file open on
        ``unit``, one would do:

      ::

            call mpp_open( unit, ... iospec='-F cachea' )

      | on an SGI/Cray system, which would pass the supplied ``iospec`` to the ``assign(3F)`` system call.
      | Currently ``iospec`` performs no action on non-SGI/Cray systems. The interface is still provided, however: users
        are cordially invited to add the requisite system calls for other systems.

h. .. rubric:: Mpp_close
      :name: mpp_close

   ::

      call mpp_close ( unit, action )

   **DESCRIPTION**
      Closes the open file on ``unit``. Clears the ``type(filetype)`` object ``mpp_file(unit)`` making it available for
      reuse.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``unit``                                                  | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``action``                                                | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

i. .. rubric:: Mpp_read_meta
      :name: mpp_read_meta

   ::

      call mpp_read_meta (unit)

   **DESCRIPTION**
      This routine is used to read the metadata describing the contents of a file. Each file can contain any number of
      fields, which are functions of 0-3 space axes and 0-1 time axes. (Only one time axis can be defined per file). The
      basic metadata defined above for ``axistype`` and ``fieldtype`` are stored in ``mpp_io_mod`` and can be accessed
      outside of ``mpp_io_mod`` using calls to ``mpp_get_info``, ``mpp_get_atts``, ``mpp_get_vars`` and
      ``mpp_get_times``.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``unit``                                                  | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **NOTE**
      ``mpp_read_meta`` must be called prior to ``mpp_read``.

j. .. rubric:: Mpp_get_info
      :name: mpp_get_info

   ::

      call mpp_get_info ( unit, ndim, nvar, natt, ntime )

   **DESCRIPTION**
      Get some general information about a file.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``unit``                                                  | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``ndim``                                                  | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``nvar``                                                  | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``natt``                                                  | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``ntime``                                                 | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

k. .. rubric:: Mpp_get_times
      :name: mpp_get_times

   ::

      call mpp_get_times ( unit, time_values )

   **DESCRIPTION**
      Get file time data.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``unit``                                                  | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **INPUT/OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``time_values``                                           | [real(DOUBLE_KIND), dimension(:)]                         |
      +-----------------------------------------------------------+-----------------------------------------------------------+

l. .. rubric:: Mpp_flush
      :name: mpp_flush

   ::

      call mpp_flush (unit)

   **DESCRIPTION**
      Flushes the open file on ``unit`` to disk. Any outstanding asynchronous writes will be completed. Any buffer
      layers between the user and the disk (e.g the FFIO layer on SGI/Cray systems) will be flushed. Calling
      ``mpp_flush`` on a unit opened with the ``MPP_RDONLY`` attribute is likely to lead to erroneous behaviour.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``unit``                                                  | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

m. .. rubric:: Mpp_get_ncid
      :name: mpp_get_ncid

   ::

       
      mpp_get_ncid (unit)

   **DESCRIPTION**
      This returns the ``ncid`` associated with the open file on ``unit``. It is used in the instance that the user
      desires to perform netCDF calls upon the file that are not provided by the ``mpp_io_mod`` API itself.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``unit``                                                  | [integer]                                                 |
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

   Any module or program unit using ``mpp_io_mod`` must contain the line
   ::

           use mpp_io_mod

   If netCDF output is desired, the cpp flag ``-Duse_netCDF`` must be turned on. The loader step requires an explicit
   link to the netCDF library (typically something like ``-L/usr/local/lib -lnetcdf``, depending on the path to the
   netCDF library). `netCDF release 3 for fortran <http://www.unidata.ucar.edu/packages/netcdf/guidef>`__ is required.
   Please also consider the compiling and linking requirements of linking as described in :doc:`./mpp_domains` and
   :doc:`./mpp`, which are ``use``\ d by this module.

| 

Precompiler options
-------------------

| ``mpp_io_mod`` uses standard f90. On SGI/Cray systems, certain I/O characteristics are specified using ``assign(3F)``.
  On other systems, the user may have to provide similar capability if required.
| There are some OS-dependent pre-processor directives that you might need to modify on non-SGI/Cray systems and
  compilers.

| 

Loader options
--------------

The source consists of the main source file and also requires the following include files: (when compiled with ) GFDL
users can check it out of the main CVS repository as part of the CVS module. The current public tag is . External users
can download the latest package . Public access to the GFDL CVS repository will soon be made available.

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
