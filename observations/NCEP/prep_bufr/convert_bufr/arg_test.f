C DART software - Copyright UCAR. This open source software is provided
C by UCAR, "as is", without charge, subject to all terms of use at
C http://www.image.ucar.edu/DAReS/DART/DART_download
C
C DART $Id$

      program arg_test
C$$$  MAIN PROGRAM DOCUMENTATION BLOCK
C
C MAIN PROGRAM:  arg_test
C   PRGMMR: collins          ORG: ncar        DATE: 07-07-27
C
C ABSTRACT: This program tests the ability to get command line arguments
C   into a Fortran program.  If it does not work, comment in the hard coded
C   filenames in the grabbufr program.
C
C PROGRAM HISTORY LOG:
C
C DART $Id$
C
C USAGE:  arg_test inputBUFRfile ouputBUFRfile
C
C   INPUT FILES:
C     unit 11  - Input BUFR file.
C
C   OUTPUT FILES:
C     unit 51  - Output BUFR file.
C
C   SUBPROGRAMS CALLED: (LIST ALL CALLED FROM ANYWHERE IN CODES)
C     UNIQUE:  
C     LIBRARY:
C       System   - getarg
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
      CHARACTER(len=80) :: infile,outfile
      INTEGER(4)        :: narg,iargc

C
C  GET Filename ARGUMENTS
C
C  if your machine does not support the iargc or getarg functions
C  comment this section out and use the hardcoded filenames below.
      NARG=IARGC()
      IF(NARG.NE.2) THEN
        PRINT *,'arg_test:  Incorrect usage'
        PRINT *,'Usage: arg_test inputBUFRfile ouputBUFRfile'
        CALL EXIT(2)
      ENDIF

      call getarg(1,infile)
      infile = TRIM(infile)//CHAR(0)
      call getarg(2,outfile)
      outfile = TRIM(outfile)//CHAR(0)

C  hardcoded filenames are the alternative to getting the filenames from
C  the command line.
c     infile = 'prepqm.bigendian'
c     outfile = 'prepqm.littleendian'
c
C
      PRINT*,'INPUT FILE: ',trim(infile)
      PRINT*,'OUTPUT FILE: ',trim(outfile)

      stop
      end 

C <next few lines under version control, do not edit>
C $URL$
C $Revision$
C $Date$
