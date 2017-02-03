/* 
 * This code is not protected by the DART copyright agreement.
 * DART $Id$
 */


/*
  byteswap.c - byteswaps all the words of a data file.

  Written by Hans Vahlenkamp
  Geophysical Fluid Dynamics Laboratory/NOAA
  Princeton University Forrestal Campus
  Last updated: 7/1/99
*/

#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>

#define BUFSIZE 50000
#define WORDSIZE 4

void usage()
  {
   printf("Usage: byteswap infile outfile [-w n]\n\n");
   printf("  -w  Specify the size of a word as n bytes (default is 4)\n\n");
   printf("This program swaps the bytes within all the words of a data file (i.e. first\n");
   printf("byte becomes last byte, second byte becomes second-to-last byte, etc.)  By\n");
   printf("default, it assumes that the word size is the typical 4 bytes, but it could\n");
   printf("also be 2 or 8.  In order for byteswapping to be meaningful, the file size\n");
   printf("should be a multiple of the word size.  Note: the output file will be\n");
   printf("overwritten if it already exists!\n");
  }

int main(int argc, char **argv)
  {
   FILE *infid, *outfid;
   int nread, wsize=4, b, errornum=0;
   struct stat infileinfo;
   unsigned char databuf[BUFSIZE];
   unsigned char *byteptr, byteval;

   /* Check for valid command-line arguments */
   if (argc < 2)
     {
      usage(); return(1);
     }
   if (argc > 3)
     {
      if (!strcmp(argv[3],"-w"))
        {
         if (argc < 5)
           {
            usage(); return(1);
           }
         sscanf(argv[4],"%d",&wsize);
         if (wsize!=2 && wsize!=4 && wsize!=8)
           {
            fprintf(stderr,"Error: bad word size!\n"); return(1);
           }
        }
      else
        {
         usage(); return(1);
        }
     }
   if ((infid=fopen(argv[1],"rb"))==NULL)
     {
      fprintf(stderr,"Error: cannot read '%s'\n",argv[1]); return(1);
     }
   if ((outfid=fopen(argv[2],"wb"))==NULL)
     {
      fclose(infid); fprintf(stderr,"Error: cannot write to '%s'\n",argv[2]);
      return(1);
     }

   /* The file size must be a multiple of the word size */
   if (stat(argv[1],&infileinfo)!=0)
     {
      fprintf(stderr,"Error: cannot read the input file size!\n");
      fclose(infid); fclose(outfid); return(1);
     }
   if (infileinfo.st_size%wsize > 0)
     {
      fprintf(stderr,"Error: input file size is not a multiple of the word size!\n");
      fclose(infid); fclose(outfid); return(1);
     }

   /* Byteswap all the words in the file */
   switch (wsize)
     {
      case 2:
        while ((nread=fread((void *)databuf,1,BUFSIZE,infid)) > 0)
          {
           byteptr=databuf;
           for (b=0; b < nread; b+=wsize)
             {
              byteval=byteptr[0]; byteptr[0]=byteptr[1]; byteptr[1]=byteval;
              byteptr+=wsize;
             }
           if (fwrite((void *)databuf,1,nread,outfid)!=nread) break;
          }
        break;
     case 4:
        while ((nread=fread(databuf,1,BUFSIZE,infid)) > 0)
          {
           byteptr=databuf;
           for (b=0; b < nread; b+=wsize)
             {
              byteval=byteptr[0]; byteptr[0]=byteptr[3]; byteptr[3]=byteval;
              byteval=byteptr[1]; byteptr[1]=byteptr[2]; byteptr[2]=byteval;
              byteptr+=wsize;
             }
           if (fwrite(databuf,1,nread,outfid)!=nread) break;
          }
        break;
      case 8:
        while ((nread=fread(databuf,1,BUFSIZE,infid)) > 0)
          {
           byteptr=databuf;
           for (b=0; b < nread; b+=wsize)
             {
              byteval=byteptr[0]; byteptr[0]=byteptr[7]; byteptr[7]=byteval;
              byteval=byteptr[1]; byteptr[1]=byteptr[6]; byteptr[6]=byteval;
              byteval=byteptr[2]; byteptr[2]=byteptr[5]; byteptr[5]=byteval;
              byteval=byteptr[3]; byteptr[3]=byteptr[4]; byteptr[4]=byteval;
              byteptr+=wsize;
             }
           if (fwrite(databuf,1,nread,outfid)!=nread) break;
          }
        break;
     }
   if (ferror(infid))
     {
      fprintf(stderr,"Error: problem reading from the input file!\n");
      errornum=1;
     }
   if (ferror(outfid))
     {
      fprintf(stderr,"Error: problem reading from the input file!\n");
      errornum=1;
     }

   /* Done */
   fclose(infid); fclose(outfid); return(errornum);
  }

/*
 * <next few lines under version control, do not edit>
 * $URL$
 * $Id$
 * $Revision$
 * $Date$
 */
