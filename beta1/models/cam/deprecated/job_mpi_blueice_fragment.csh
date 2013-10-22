# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$

# For the IBM/blueice machine. This fragment can be dropped into the job_mpi.csh script
# to specifically schedule which mpi tasks run on which nodes. Since cores on the same
# nodes share memory, if running one mpi task per core uses too much memory, everything
# slows down drastically.  Better to run fewer tasks per node (e.g. 7 tasks on an 8 core
# node) since then each task gets more memory (the memory on these nodes is not
# partitioned by core; it is completely shared and a single task can take over all
# the memory if you let it).

   if ($parallel_cam == 'false' ) then
      if ($num_procs == 96) then
         # This environment variable tells how many processors on each node to use
         # I want 80 = 2*14 + 4*13
         echo "setenv LSB_PJL_TASK_GEOMETRY \"                  >> ${job_i}
         echo ' "{(0,1,2,3,4,5,6,7,8,9,10,11,12,13)\'           >> ${job_i}
         echo " (14,15,16,17,18,19,20,21,22,23,24,25,26,27)\"   >> ${job_i}
         echo " (28,29,30,31,32,33,34,35,36,37,38,39,40)\"      >> ${job_i}
         echo " (41,42,43,44,45,46,47,48,49,50,51,52,53)\"      >> ${job_i}
         echo " (54,55,56,57,58,59,60,61,62,63,64,65,66)\"      >> ${job_i}
         echo ' (67,68,69,70,71,72,73,74,75,76,77,78,79)}"'     >> ${job_i}
      else if ($num_procs == 48) then
         # I want 40 = 1*14 + 2*13
         echo "setenv LSB_PJL_TASK_GEOMETRY \"                  >> ${job_i}
         echo ' "{(0,1,2,3,4,5,6,7,8,9,10,11,12,13)\'           >> ${job_i}
         echo " (14,15,16,17,18,19,20,21,22,23,24,25,26)\"      >> ${job_i}
         echo ' (27,28,29,30,31,32,33,34,35,36,37,38,39)}"'     >> ${job_i}
      else
         echo "parallel_cam is false, but num_procs is not 96 or 48" >> $MASTERLOG
         exit
      endif

   endif

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

