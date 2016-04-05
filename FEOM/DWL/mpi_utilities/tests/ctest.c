/* DART software - Copyright 2004 - 2013 UCAR. This open source software is
 * provided by UCAR, "as is", without charge, subject to all terms of use at
 * http://www.image.ucar.edu/DAReS/DART/DART_download
 *
 * $Id$
 */

/* A simple "c" program - no mpi calls, no netCDF - to test you have a      */
/* working c compiler.                                                      */
/*                                                                          */
/* DART contains no c code, but if you are having problems with either the  */
/* MPI or netCDF libraries and you want to diagnose whether the problem is  */
/* with the entire installation or with just the F90 interfaces, these      */
/* c programs will allow you to test the c interfaces for mpi and netCDF.   */

#include <stdlib.h>
#include <stdio.h>

main(int argc, char **argv)
{
   int i, j;

   i = 2;
   j = 3;

   printf("2 + 3 = %d\n", i+j);

   printf("'c' program ran successfully\n");

   exit(0);
}

/* <next few lines under version control, do not edit>
 * $URL$
 * $Id$
 * $Revision$
 * $Date$
 */
