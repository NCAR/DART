/* 
 * DART software - Copyright UCAR. This open source software is provided
 * by UCAR, "as is", without charge, subject to all terms of use at
 * http://www.image.ucar.edu/DAReS/DART/DART_download
 *
 * DART $Id$
 */


/*
 * swap bytes in a binary/unformatted DART restart file.  
 *  (swab == the machine instruction for swap bytes)
 *
 * usage:   swabrestart < file.in > file.out
 *
 *  takes no command line arguments; always reads standard input and writes
 *  standard output
 *
 * this program changes binary data from little-endian format to big-endian
 * and back (little = intel/amd, dec alphas, pdp-11s) big = everything else).
 * this requires a specialized program which knows the data layout
 * in the file because you must swab by item size: 4 bytes at a time for
 * 4-byte ints and reals, 8 bytes at a time for 8-byte ints and reals.
 *
 * these days, the most common need for this program arises when moving 
 * binary data files between an ibm power chip platform and an intel/AMD chip 
 * platform.  (look up 'endianness' in wikipedia for more information than
 * you ever wanted to know about ordering the bits/bytes in a binary number.
 * there used to be more chip manufacturers -- sgi (mips), sun (sparc), 
 * dec (alpha) -- which made programs like this more common, but the list
 * of different cpu architectures has consolidated immensely since the 90s.)
 *
 *
 * this program now has several *compile time* options.  look below at
 * the 'user settable options' section, change the #define values, and
 * then recompile this program in order for them to take effect.
 *
 * this is a self-contained, generic C program.  all systems these days
 * come with the 'cc' compiler.  compile it with:
 *
 *   cc -O swabrestart.c -o swabrestart
 *
 * -O specifies a standard level of optimization
 * -o sets the output executable name
 * this program uses no external libraries
 * 
 * by default this program expects to be compiled and run on the target 
 * platform you want to convert the data onto!  if you need to run it on the
 * platform where the data was generated and *then* move the data to a machine
 * of the opposite byte-order, see the #defines below to change the behavior.
 *
 *
 * the fortran standard "unformatted" record format is:  
 *  4 byte integer: record byte count (*)
 *  n items of whatever type  (n = record_byte_count / item_size_in_bytes)
 *  4 byte integer: repeated record byte count
 *
 * (*) see below for comment on gfortran sometimes using 8 bytes for this
 *     record count.  note that any discussion of the gfortran compiler
 *     refers to the F90 compiler that was used to compile the DART
 *     executable either on the machine which wrote the restart files,
 *     or on the machine which is going to read the files.
 *  
 * the data inside a DART restart file is:
 *  record 1:
 *    4 byte integer:  day_number (gregorian calendar, # days since 1/1/1600)
 *    4 byte integer:  seconds  (0 <= s < 86400)
 *  record 2:
 *    real*8 data, state vector length  (or real*4 if compiled to be so)
 * 
 * unfortunately, this has gotten a lot more complicated, but here are
 * some typical examples of what you will see in a dart restart file.
 *
 * most likely scenerio:
 * 1. assume any F90 compiler except gfortran was used to compile the DART
 * executable that wrote the restart file in the first place.
 * 2. assume that this DART was compiled with the default R8 type = 8 bytes.
 * 3. assume that a single restart file was selected in the DART namelist,
 * so the restart data for all the ensemble members are in this single file.
 * 
 * the structure of this file is then:
 *
 *   4 byte integer: record byte count (will be = 8)
 *   4 byte integer: day_number (gregorian calendar, # days since 1/1/1600)
 *   4 byte integer: seconds (0 <= s < 86400)
 *   4 byte integer: record byte count (again, will be = 8)
 *   4 byte integer: record byte count (state length*8; big)
 *   N bytes real*8: state data
 *   4 byte integer: record byte count (will match prev count)
 *   <repeat of last 3 values for each ensemble in the file.>
 * 
 * another possible scenerio:
 * 1. assume any F90 compiler except gfortran was used to compile the DART
 * executable that wrote the restart file in the first place.
 * 2. assume that this DART was compiled overriding the R8 type = 4 bytes.
 * 3. assume that multiple restart files were selected in the DART namelist,
 * so this file contains only the data for a single ensemble.
 * 
 * the structure of this file is then:
 *
 *   4 byte integer: record byte count (will be = 8)
 *   4 byte integer: day_number (gregorian calendar, # days since 1/1/1600)
 *   4 byte integer: seconds (0 <= s < 86400)
 *   4 byte integer: record byte count (again, will be = 8)
 *   4 byte integer: record byte count (state length*4; big)
 *   N bytes real*4: state data
 *   4 byte integer: record byte count (will match prev count)
 * 
 * and finally, a gfortran example:
 * 1. assume a recent gfortran F90 compiler was used to compile the DART
 * executable that wrote the restart file in the first place.
 * 2. assume that this DART was compiled with the default R8 type = 8 bytes.
 * 3. assume that multiple restart files were selected in the DART namelist,
 * so this file contains only the data for a single ensemble.
 * 
 * the structure of this file is then:
 *
 *   8 byte integer: record byte count (will be = 8)
 *   4 byte integer: day number (gregorian calendar - days since 1/1/1600)
 *   4 byte integer: seconds (0 <= s < 86400)
 *   8 byte integer: record byte count (again, will be = 8)
 *   8 byte integer: record byte count (state length*8; big)
 *   N bytes real*8: state data
 *   8 byte integer: record byte count (will match prev count)
 * 
 *
 * as discussed above, the public domain gfortran compiler has changed to 
 * defaulting to use not 4 byte record sizes, but 8 bytes!  this enables 
 * larger record sizes (> 2GB) but is incompatible with every other compiler
 * on the planet at this point.  there are two solutions to the problem of
 * writing with a gfortran-compiled DART and then reading with a non-gfortran-
 * compiled DART, or vice-versa.  one is to set the gfortran flags at DART
 * compile time to set the record counts back to 4 bytes (such flags do exist).
 * the other is to look below and set the #defines to expect to read or write
 * 8 byte counts.  note that you can set the incoming and outgoing sizes
 * separately so you can translate between gfortran and non-gfortran files.
 *
 * this is not common, but if this is a 'model advance' file, then it
 * has 2 timestamp records at the start before the state vector data.
 * these files do not normally hang around, but in case you have one you
 * want to translate, there is optional code to handle a second timestamp
 * record.
 *
 * nsc 17oct06
 * nsc 28aug07 - updated comments, made R4 option easier to enable
 * nsc 29jan08 - updated (but never committed) several months ago to enable
 *               compile-time options to:
 *   - swab on source platform before moving file
 *   - read and/or write new gfortran 8 byte record counts
 *   - handle model-advance files
 *
 * ideas for additional functionality:  
 * - could add printout at end of ensemble count found in file.
 *
 * DART $Id$
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <strings.h>
#include <ctype.h>


/* start user settable options.  must recompile after changes. */

/*
 * if you run DART with R8 redefined as R4, set next line to 4 instead of 8
 * (e.g. single precision float data instead of double precision)
 * 8 byte reals is the DART default.
 */
#define R_DATA_SIZE 8     /* set to 4 if using real*4 data, 8 for real*8 */

/* 
 * if your workflow is:
 *  1. run DART and generate a binary restart file
 *  2. copy the file over to a different platform first
 *  3. run this program on the target platform after the file was moved
 *  4. run DART and read the fixed file.
 * then the existing program default value of NATIVE 0 is correct.
 *
 * instead, if you want to:
 *  1. run DART and generate a binary restart file
 *  2. run this program to switch the byte order first
 *  3. copy the file to the target platform after converting
 *  4. run DART on the target platform
 * then you must set the NATIVE flag to 1.
 *
 * (the difference is how to interpret the record byte counts around the 
 *  actual data; in one case they must be swabbed first.)
 */
#define NATIVE       0     /* set to 1 if on same platform as data generated */

/* 
 * if the data was either written by DART compiled with gfortran, and/or
 * it is going to be read in by DART compiled with gfortran, set one or both
 * of these to 8.  mixed sizes are supported to convert the rec count fields.
 * at this point in time, all other F90 compilers require a size of 4.
 * note that there are compile-time flags to force even gfortran to use
 * 4-byte record counts.  if you set that option, you can leave these at 4.
 */
#define  IN_REC_SIZE 4     /* set to 8 if gfortran was used to *write* file */
#define OUT_REC_SIZE 4     /* set to 8 if will use gfortran to *read* file  */

/*
 * this is not common, but if this is a 'model advance' file, then it has
 * a normal data timestamp record plus a 'how far the model should advance'
 * timestamp, and then the state vector.  if you set this to 1, it will
 * look for and convert 2 time records before the state vector.
 */
#define MODEL_ADV    0     /* set to 1 if 2 timestamps before state vector */

/* end user settable options */


/* my own swab routine.  libc has a built-in swab() function, 
 * but it does not work on data in-place. 
 */
void swapbytes(char *buf, int bcount);

main(int argc, char **argv) 
{
    char  in_rec[ IN_REC_SIZE];
    char out_rec[OUT_REC_SIZE];
    char fourbyte[4];
    char eightbyte[R_DATA_SIZE];
    int i, nreals, nreals2, nints, nints2, count;

    /* warn we do not process any args. */
    if (argc > 1) 
        fprintf(stderr, "warning: all arguments are ignored\n");
    
    /* record 1: */
    count = read(0, in_rec, IN_REC_SIZE);
  again:
    if (count != IN_REC_SIZE) {
        fprintf(stderr, "error: cannot read %d bytes from standard input\n",
                         IN_REC_SIZE);
        exit(-1);
    }
 
    /* number of int*4s in the first record */
    if (NATIVE != 0) {
        nints = *((int *)in_rec) / 4;
        *((int *)out_rec) = *((int *)in_rec);
        swapbytes(out_rec, OUT_REC_SIZE);
    } else {
        swapbytes(in_rec, IN_REC_SIZE);
        nints = *((int *)in_rec) / 4;
        *((int *)out_rec) = *((int *)in_rec);
    }
    write(1, out_rec, OUT_REC_SIZE); 
    fprintf(stderr, "reading %d ints\n", nints);
   
    /* simple error checking.  DART timestamps are a known size */
    if (nints != 2) {
        fprintf(stderr, 
                "file format error - first record should be 2 ints, not %d\n", 
                 nints);
        exit(-1);
    }

    /* day, secs.  DART writes these as int*4's.  */
    read(0, fourbyte, 4);
    swapbytes(fourbyte, 4);
    write(1, fourbyte, 4);

    read(0, fourbyte, 4);
    swapbytes(fourbyte, 4);
    write(1, fourbyte, 4);

    /* repeated record count */
    read(0, in_rec, IN_REC_SIZE);
    if (NATIVE != 0) {
        nints2 = *((int *)in_rec) / 4;
        *((int *)out_rec) = *((int *)in_rec);
        swapbytes(out_rec, OUT_REC_SIZE);
    } else {
        swapbytes(in_rec, IN_REC_SIZE);
        nints2 = *((int *)in_rec) / 4;
        *((int *)out_rec) = *((int *)in_rec);
    }
    write(1, out_rec, OUT_REC_SIZE);

    /* sanity check */
    if (nints != nints2) {
        fprintf(stderr, 
            "final record byte count (%d) does not match initial count (%d)\n", 
             nints2, nints);
        exit (-1);
    }

#if MODEL_ADV
    /* record 1a: */
    read(0, in_rec, IN_REC_SIZE);

    /* number of int*4s in the timestamp record */
    if (NATIVE != 0) {
        nints = *((int *)in_rec) / 4;
        *((int *)out_rec) = *((int *)in_rec);
        swapbytes(out_rec, OUT_REC_SIZE);
    } else {
        swapbytes(in_rec, IN_REC_SIZE);
        nints = *((int *)in_rec) / 4;
        *((int *)out_rec) = *((int *)in_rec);
    }
    write(1, out_rec, OUT_REC_SIZE); 
    fprintf(stderr, "reading %d ints\n", nints);

    /* simple error checking.  DART timestamps are a known size */
    if (nints != 2) {
        fprintf(stderr, 
                "file format error - second record should be 2 ints, not %d\n", 
                 nints);
        exit(-1);
    }

    /* day, secs.  DART writes these as int*4's.  */
    read(0, fourbyte, 4);
    swapbytes(fourbyte, 4);
    write(1, fourbyte, 4);

    read(0, fourbyte, 4);
    swapbytes(fourbyte, 4);
    write(1, fourbyte, 4);

    /* repeated record count */
    read(0, in_rec, IN_REC_SIZE);
    if (NATIVE != 0) {
        nints2 = *((int *)in_rec) / 4;
        *((int *)out_rec) = *((int *)in_rec);
        swapbytes(out_rec, OUT_REC_SIZE);
    } else {
        swapbytes(in_rec, IN_REC_SIZE);
        nints2 = *((int *)in_rec) / 4;
        *((int *)out_rec) = *((int *)in_rec);
    }
    write(1, out_rec, OUT_REC_SIZE);

    /* sanity check */
    if (nints != nints2) {
        fprintf(stderr, 
            "final record byte count (%d) does not match initial count (%d)\n", 
             nints2, nints);
        exit (-1);
    }
#endif

    /* record 2: */
    read(0, in_rec, IN_REC_SIZE);
    if (NATIVE != 0) {
        nreals = *((int *)in_rec) / R_DATA_SIZE;
        *((int *)out_rec) = *((int *)in_rec);
        swapbytes(out_rec, OUT_REC_SIZE);
    } else {
        swapbytes(in_rec, IN_REC_SIZE);
        nreals = *((int *)in_rec) / R_DATA_SIZE;
        *((int *)out_rec) = *((int *)in_rec);
    }
    write(1, out_rec, OUT_REC_SIZE);

    /* number of reals in the rest of the record */
    fprintf(stderr, "reading %d reals\n", nreals);

    /* actual data */
    for (i=0; i<nreals; i++) {
        read(0, eightbyte, R_DATA_SIZE);
        swapbytes(eightbyte, R_DATA_SIZE);
        write(1, eightbyte, R_DATA_SIZE);
    }

    /* final count */
    read(0, in_rec, IN_REC_SIZE);
    if (NATIVE != 0) {
        nreals2 = *((int *)in_rec) / R_DATA_SIZE;
        *((int *)out_rec) = *((int *)in_rec);
        swapbytes(out_rec, OUT_REC_SIZE);
    } else {
        swapbytes(in_rec, IN_REC_SIZE);
        nreals2 = *((int *)in_rec) / R_DATA_SIZE;
        *((int *)out_rec) = *((int *)in_rec);
    }
    write(1, out_rec, OUT_REC_SIZE);

    /* sanity check */
    if (nreals != nreals2) {
        fprintf(stderr, 
               "final data byte count (%d) does not match initial count (%d)\n",
                nreals2, nreals);
        exit (-1);
    }

    /* verify we are at the end by trying to read further.  restart files can
     * contain data for multiple ensembles, so if more data does exist in the
     * file, loop back up to the top and repeat until all data is processed.
     */
    count = read(0, in_rec, IN_REC_SIZE);
    if (count == 0) {
        fprintf(stderr, "end of data\n");
        exit (0);
    }

    /* is there more data?  if so, count already read in*/
    goto again;

    exit (0);
}


void swapbytes(char *buf, int bcount)
{
    char tmp;

    switch (bcount) {
      case 2:
        tmp = buf[0];
        buf[0] = buf[1];
        buf[1] = tmp;
        break;

      case 4:
        tmp = buf[0];
        buf[0] = buf[3];
        buf[3] = tmp;
        tmp = buf[1];
        buf[1] = buf[2];
        buf[2] = tmp;
        break;

      case 8:
        tmp = buf[0];
        buf[0] = buf[7];
        buf[7] = tmp;
        tmp = buf[1];
        buf[1] = buf[6];
        buf[6] = tmp;
        tmp = buf[2];
        buf[2] = buf[5];
        buf[5] = tmp;
        tmp = buf[3];
        buf[3] = buf[4];
        buf[4] = tmp;
        break;

      default:
        fprintf(stderr, "internal error: unsupported item size %d\n", bcount);
        exit(-1);
    }

    return;
}

/* <next few lines under version control, do not edit>
 * $URL$
 * $Revision$
 * $Date$
 */
