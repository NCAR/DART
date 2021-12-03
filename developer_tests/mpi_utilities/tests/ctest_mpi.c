/* 
 * DART software - Copyright UCAR. This open source software is provided
 * by UCAR, "as is", without charge, subject to all terms of use at
 * http://www.image.ucar.edu/DAReS/DART/DART_download
 */

/* A simple MPI "c" program to test if the c interfaces for MPI work.
 *
 * DART contains no c code, but if you are having problems with either the
 * MPI or netCDF libraries and you want to diagnose whether the problem is
 * with the entire installation or with just the F90 interfaces, these
 * c programs will allow you to test the c interfaces for mpi and netCDF.
 */

#include <mpi.h>

#include <stdlib.h>
#include <stdio.h>

int main(int argc, char **argv)
{
   int i;

   MPI_Init(&argc, &argv);
   printf("returned from init\n");

   MPI_Comm_rank(MPI_COMM_WORLD, &i);
   printf("rank = %d\n", i);

   MPI_Comm_size(MPI_COMM_WORLD, &i);
   printf("task count = %d\n", i);

   MPI_Finalize();
   printf("returned from finalize\n");

   exit(0);
}
