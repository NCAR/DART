C DART software - Copyright UCAR. This open source software is provided
C by UCAR, "as is", without charge, subject to all terms of use at
C http://www.image.ucar.edu/DAReS/DART/DART_download
C
C DART $Id$

      program stat_test

C$$$  MAIN PROGRAM DOCUMENTATION BLOCK
C
C MAIN PROGRAM:  stat_test
C   PRGMMR: n.collins          ORG: ncar        DATE: 07-07-24
C
C ABSTRACT: the STAT() function returns an array of system-dependent 
C   information.  one of the entries in the array is the number of bytes
C   in a file.  the problem is that the index for this entry is different on
C   different systems.  this program makes the STAT() call and prints out
C   the array contents.  one of the numbers will match the size of a data
C   file and that will be the index required by the grabbufr program
C   (also in this directory).   if grabbufr does not work, run this program
C   and verify the STAT index is correct.
C
C PROGRAM HISTORY LOG:
C 20067-07-24  collins
C
C DART $Id$
C
C USAGE:  stat_test inputBUFRfile
C
C   INPUT FILES:
C     unit 11  - Input BUFR file.
C
C   OUTPUT FILES:
C     unit 6   - console output
C
C   SUBPROGRAMS CALLED: (LIST ALL CALLED FROM ANYWHERE IN CODES)
C     LIBRARY:
C       System   - getarg stat
C       W3LIB    - errexit 
C
C   EXIT STATES:
C     COND =   0 - SUCCESSFUL RUN
C          =   2 - Incorrect argument list
C          =   4 - Coud not allocate memory to hold Input BUFR file
C          =  99 - Could not obtain size of input BUFR file
C
C REMARKS: LIST CAVEATS, OTHER HELPFUL HINTS OR INFORMATION
C
C ATTRIBUTES:
C   LANGUAGE: Fortran 90
C   MACHINE:  IBM SP
C
C$$$
      CHARACTER(len=80) :: infile
      INTEGER(4)        :: narg,iargc,JSTAT(100)
      integer           :: i, KBYTES, rc
      integer           :: STAT

C
c liu 03/16/2005
C  GET Filename ARGUMENTS from command line.
C  IF THIS DOES NOT WORK, comment this entire section out and use the
C  hardcoded filenames below.
C
      NARG=IARGC()
      IF(NARG.NE.1) THEN
        PRINT *,'stat_test:  Incorrect usage'
        PRINT *,'Usage: stat_test inputBUFRfile'
        CALL EXIT(2)
      ENDIF

      call getarg(1,infile)
      infile = TRIM(infile)//CHAR(0)

C  If your system does not support IARGC and getarg(), comment out from
C  the previous comment to here, and comment in the following line.
C  Then link or rename your input BUFR files to match the name.

c     infile = 'prepqm.bigendian'
c
c liu 03/16/2005
C
C  Use STAT function to get size of input BUFR file
C
c  this used to be a oneliner, but the function test failed
c  on recent intel compilers.  splitting the call and test into
c  two separate lines with an explicit integer variable seems
c  to have fixed it.  also, for all the recent compilers i have
c  tested the right offset seems to be 8
      rc = STAT(infile,JSTAT)
      IF (rc.NE.0) THEN
         PRINT*,'ERROR IN FUNCTION STAT GETTING FILE INFO, RC = ',rc
         CALL EXIT(99)
      ELSE
c        Use the following print and find the index into JSTAT which is 
c        the same as the actual size of the bufr input file.  make sure
c        it matches the value of INDEXVAL in grabbufr.f
 
         do i = 1,size(JSTAT)
            if (jstat(i) /= 0) print*,'index: ',i,' = ',JSTAT(i)
         enddo

         PRINT *,''
         PRINT *,'ONE OF THESE VALUES SHOULD MATCH THE BYTE COUNT'
         PRINT *,'FROM THE ls -l COMMAND.  VERIFY THAT INDEXVAL'
         PRINT *,'IN grabbufr.f IS THE SAME'
      ENDIF
 
      stop
      end

c <next few lines under version control, do not edit>
c $URL$
c $Revision$
c $Date$
