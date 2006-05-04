/* a simple "c" program - no mpi calls - to test the c compiler */

main(int argc, char **argv)
{
   int i, j;

   i = 2;
   j = 3;

   printf("2 + 3 = %d\n", i+j);

   printf("'c' program ran successfully\n");

   exit(0);
}

