/*
 * swap bytes in a binary/unformatted DART restart file.
 *  (swap bytes == swab)
 *
 * the generic fortran "unformatted" record format is:  
 *  4 byte integer: record count
 *  n items of whatever type
 *  4 byte integer: repeated record count
 *  
 * the data inside a DART restart file is:
 *  record 1:
 *    4 byte integer:  day_number
 *    4 byte integer:  seconds
 *  record 2:
 *    real*8 data, state vector length
 * 
 * so all together, the bytes you will see in the file are:
 *   4 byte integer: record count and it better be 8
 *   4 byte integer: day number
 *   4 byte integer: seconds
 *   4 byte integer: record count and it better be 8
 *   4 byte integer: record count and it will be a big number (N)
 *   N bytes real*8: state data
 *   4 byte integer: record count and it better be N
 * 
 *  N better be the number of bytes in the file, minus 24, divided by 8.
 *
 * dart programs can either create individual restart files, or they can
 * concatinate the restart data for each ensemble member into the same file.
 * so at the end, this program loops to see if there is another dataset and
 * if so, processes it until it reaches the real end of file.
 *
 * this program changes data from little-endian machines to big-endian
 * and back (little = intel anything, big = everyone else except dec alphas
 * and pdp-11s).   you must swab by the item size, 4 bytes for the ints,
 * 8 bytes for the double reals (4 bytes if you have single precision reals).
 * running this program twice should produce output identical to the original 
 * input file.
 *
 * reads stdin, writes stdout.
 * takes no args.
 *
 * nsc 17oct06
 * nsc 28aug07 - updated comments, made R4 option easier to enable
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <strings.h>
#include <ctype.h>

/*
 * if you run with R8 redefined as R4, then set this next line to 4 instead of 8
 * (e.g. single precision floats instead of double precision)
 */
#define R8_SIZE 8

/* cannot call mine swab() since there is one in libc, but it does not work
 * in place. 
 */
void swapbytes(char *buf, int bcount);

main(int argc, char **argv) 
{
    char fourbyte[4];
    char eightbyte[R8_SIZE];
    int i, nreals, nints, count;

    /* warn we do not process any args. */
    if (argc > 1) 
        fprintf(stderr, "warning: all arguments are ignored\n");
    
    /* record 1: */
    count = read(0, fourbyte, 4);
again:
    if (count != 4) {
        fprintf(stderr, "error: cannot read 4 bytes from standard input\n");
        exit(-1);
    }
    swapbytes(fourbyte, 4);
    write(1, fourbyte, 4);

    /* number of int*4s in the first record */
    nints = *((int *)fourbyte) / 4;
    fprintf(stderr, "reading %d ints\n", nints);

    /* day, secs */
    read(0, fourbyte, 4);
    swapbytes(fourbyte, 4);
    write(1, fourbyte, 4);

    read(0, fourbyte, 4);
    swapbytes(fourbyte, 4);
    write(1, fourbyte, 4);

    /* repeated record count */
    read(0, fourbyte, 4);
    swapbytes(fourbyte, 4);
    write(1, fourbyte, 4);

    /* sanity check */
    if (nints != *((int *)fourbyte) / 4) {
        fprintf(stderr, "final count (%d) does not match initial count (%d)\n", 
               *((int *)fourbyte)/4, nints);
        exit (-1);
    }

    /* record 2: */
    read(0, fourbyte, 4);
    swapbytes(fourbyte, 4);
    write(1, fourbyte, 4);

    /* number of real*8s in the rest of the record */
    nreals = *((int *)fourbyte) / R8_SIZE;
    fprintf(stderr, "reading %d reals\n", nreals);

    /* actual data */
    for (i=0; i<nreals; i++) {
        read(0, eightbyte, R8_SIZE);
        swapbytes(eightbyte, R8_SIZE);
        write(1, eightbyte, R8_SIZE);
    }

    /* final count */
    read(0, fourbyte, 4);
    swapbytes(fourbyte, 4);
    write(1, fourbyte, 4);

    /* sanity check */
    if (nreals != *((int *)fourbyte) / R8_SIZE) {
        fprintf(stderr, "final count (%d) does not match initial count (%d)\n", 
               *((int *)fourbyte)/R8_SIZE, nreals);
        exit (-1);
    }

    /* verify we are at the end by trying to read further */
    count = read(0, fourbyte, 4);
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
