/* Data Assimilation Research Testbed -- DART               */
/* Copyright 2004-2006, Data Assimilation Research Section  */
/* University Corporation for Atmospheric Research          */
/* Licensed under the GPL -- www.gpl.org/licenses/gpl.html  */

/* <next few lines automatically updated by version control software, do not edit> */
/* $Revision$ */
/* $Date$ */
/* $Id$ */

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

