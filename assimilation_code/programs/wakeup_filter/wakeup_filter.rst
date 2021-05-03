PROGRAM ``wakeup_filter``
=========================

.. attention::

  ``wakeup_filter`` works with versions of DART *before* Manhattan (9.x.x) and has yet to be updated. If you are interested in
  using ``wakeup_filter`` with more recent versions of DART, contact DAReS staff to assess the feasibility of an update.
  Until that time, you should consider this documentation as out-of-date.

Overview
--------

Small auxiliary program for use in the "async=4" case where the main filter program is an MPI program and the model
being run with DART is also an MPI program. The main MPI job script runs each of the model advances for the ensemble
members, and then runs this program to restart the filter program.

Modules used
------------

::

   mpi_utilities_mod

Namelist
--------

There are no namelist options for this program. It must be run as an MPI program with the same number of tasks as filter
was originally started with.

Files
-----

Named pipes (fifo) files are used to synchronize with the main MPI job run script, to ensure that the filter program and
the script do not do a "busy-wait" in which they consume CPU cycles while they are waiting for each other. The fifo
names are:

-  filter_to_model.lock
-  model_to_filter.lock
-  filter_lockNNNNN (where NNNNN is the task number with leading 0s)

References
----------

-  Anderson, J., T. Hoar, K. Raeder, H. Liu, N. Collins, R. Torn, and A. Arellano, 2009:
   The Data Assimilation Research Testbed: A Community Facility. Bull. Amer. Meteor. Soc., 90, 1283-1296.
   `DOI: 10.1175/2009BAMS2618.1 <http://dx.doi.org/10.1175%2F2009BAMS2618.1>`__
