/* Data Assimilation Research Testbed -- DART               */
/* Copyright 2004-2007, Data Assimilation Research Section  */
/* University Corporation for Atmospheric Research          */
/* Licensed under the GPL -- www.gpl.org/licenses/gpl.html  */

/* <next few lines under version control, do not edit> */
/* $URL$ */
/* $Id$ */
/* $Revision$ */
/* $Date$ */

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

