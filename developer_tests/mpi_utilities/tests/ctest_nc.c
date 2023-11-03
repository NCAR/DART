/* 
 * DART software - Copyright UCAR. This open source software is provided
 * by UCAR, "as is", without charge, subject to all terms of use at
 * http://www.image.ucar.edu/DAReS/DART/DART_download
 */

/* A simple netCDF "c" program to test that the netCDF c interfaces work.   */
/*                                                                          */
/* DART contains no c code, but if you are having problems with either the  */
/* MPI or netCDF libraries and you want to diagnose whether the problem is  */
/* with the entire installation or with just the F90 interfaces, these      */
/* c programs will allow you to test the c interfaces for mpi and netCDF.   */

#include "netcdf.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

void netcdf_error_exit(int istat);

int main(int argc, char **argv)
{

   char filename[32] = "ctestdata.nc";
   int ncfileid, istat, i;
   int test1_length = 5;
   int idata[5] = { 1, 2, 3, 4, 5 };
   int test1dimid, dataid;

   printf("program start\n");

/*  Typical sequence:                                                      */
/*  NF90_OPEN             ! create netCDF dataset: enter define mode       */
/*     NF90_def_dim       ! define dimensions: from name and length        */
/*     NF90_def_var       ! define variables: from name, type, and dims    */
/*     NF90_put_att       ! assign attribute values                        */
/*  NF90_ENDDEF           ! end definitions: leave define mode             */
/*     NF90_put_var       ! provide values for variable                    */
/*  NF90_CLOSE            ! close: save updated netCDF dataset             */

/*---------------------------------------------------------------------------*/
/*  Open/Create file                                                         */
/*---------------------------------------------------------------------------*/

   istat = nc_create(filename, NC_SHARE, &ncfileid);
   if (istat != NC_NOERR) netcdf_error_exit(istat);

   printf("successfully opened '%s'\n", filename);

/*---------------------------------------------------------------------------*/
/* Define dimension(s)                                                       */
/*---------------------------------------------------------------------------*/

   istat = nc_def_dim(ncfileid, "test1", test1_length, &test1dimid);
   if (istat != NC_NOERR) netcdf_error_exit(istat);

/*---------------------------------------------------------------------------*/
/* Write global attributes                                                   */
/*---------------------------------------------------------------------------*/

   istat = nc_put_att_text(ncfileid, NC_GLOBAL, "title", 
                           strlen("netcdf test File"), "netcdf test File");
   if (istat != NC_NOERR) netcdf_error_exit(istat);

/*---------------------------------------------------------------------------*/
/* Create variables and attributes.                                          */
/*---------------------------------------------------------------------------*/

   istat = nc_def_var(ncfileid, "data", NC_INT, 1, &test1dimid, &dataid);
   if (istat != NC_NOERR) netcdf_error_exit(istat);

   istat = nc_put_att_text(ncfileid, dataid, "long_name", 
                           strlen("test data array"), "test data array");
   if (istat != NC_NOERR) netcdf_error_exit(istat);

/*---------------------------------------------------------------------------*/
/* Leave define mode so we can fill                                          */
/*---------------------------------------------------------------------------*/
   istat = nc_enddef(ncfileid);
   if (istat != NC_NOERR) netcdf_error_exit(istat);

/*---------------------------------------------------------------------------*/
/* Fill the coordinate variables.                                            */
/*---------------------------------------------------------------------------*/

   istat = nc_put_var_int(ncfileid, dataid, idata);
   if (istat != NC_NOERR) netcdf_error_exit(istat);

/*---------------------------------------------------------------------------*/
/* Close                                                                     */
/*---------------------------------------------------------------------------*/

   istat = nc_close(ncfileid);
   if (istat != NC_NOERR) netcdf_error_exit(istat);

   printf("netcdf file successfully closed\n");

   printf("program end\n");

   exit(0);
}

void netcdf_error_exit(int istat)
{
   char error_msg[128];

   strcpy(error_msg, "netcdf error, string is: ");
   strcat(error_msg, nc_strerror(istat));

   printf("%s\n", error_msg);
   exit(-1);
}

