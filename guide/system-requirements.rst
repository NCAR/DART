###################
System requirements
###################

The DART software is intended to compile and run on many different
Unix/Linux operating systems with little to no change. At this point we have no
plans to port DART to Windows machines, although Windows 10
users may be interested in the free `Windows Subsystem For Linux
<https://docs.microsoft.com/en-us/windows/wsl/about>`_
which allows developers to "run a GNU/Linux environment -- including most
command-line tools, utilities, and applications -- directly on Windows,
unmodified, without the overhead of a virtual machine" (see
https://docs.microsoft.com/en-us/windows/wsl/about for more details)

.. note::

   We have tried to make the DART code as portable as possible, but we do not
   have access to all compilers on all platforms, so unfortunately we cannot
   guarantee that the code will work correctly on your particular system.
   
   We are genuinely interested in your experience building the system, so we
   welcome you to send us an email with your experiences to dart@ucar.edu.
   
   We will endeavor to incorporate your suggestions into future versions of
   this guide.

Minimally, you will need:

1.  a Fortran90 compiler,
2.  the `netCDF libraries <http://www.unidata.ucar.edu/software/netcdf/>`_
    built with the F90 interface,
3.  *perl* (just about any version),
4.  an environment that understands *csh*, *tcsh*, *sh*, and *ksh*
5.  the long-lived Unix build tool *make*
6.  and up to 1 Gb of disk space for the DART distribution.

History has shown that it is a very good idea to remove the stack and heap
limits in your run-time environment with the following terminal commands:

.. code-block:: bash

  > limit stacksize unlimited  
  > limit datasize unlimited

Additionally, the following tools have proven to be *nice* (but are not
required to run DART):

1.  `ncview <http://meteora.ucsd.edu/~pierce/ncview_home_page.html>`_: a
    great visual browser for netCDF files.
2.  `the netCDF Operators (NCO) <http://nco.sourceforge.net/>`_: tools to
    perform operations on netCDF files like concatenating, slicing, and
    dicing
3.  Some sort of MPI environment. In other words, DART does not come
    with *MPICH*, *LAM-MPI*, or *OpenMPI*, but many users of DART rely on these
    MPI distributions to run DART in a distributed-memory parallel setting. In
    order to use MPI with DART, please refer to the DART MPI introduction.
4.  If you want to use the DART diagnostic scripts, you will need a
    basic MATLAB® installation. No additional toolboxes are required, and no
    third-party toolboxes are required.
