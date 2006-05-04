/* a "c" which calls the MPI parallel communication libraries */

#include <mpi.h>

main(int argc, char **argv)
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

