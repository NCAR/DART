C This code is not protected by the DART copyright agreement.
C DART $Id$

      program grabbufr

C$$$  MAIN PROGRAM DOCUMENTATION BLOCK
C
C MAIN PROGRAM:  grabbufr
C   PRGMMR: Gilbert          ORG: NP11        DATE: 99-07-13
C
C ABSTRACT: This program extracts all the BUFR messages from any
C   file and writes them back out to another file.  It is being used
C   to convert the blocking structure of a BUFR file to standard unix
C   Fortran format.  Also, the program converts any BUFR edition 0 and 1
C   messages to a BUFR edition 2 message before they are written out.
C
C PROGRAM HISTORY LOG:
C 1999-07-13  Gilbert
C 1999-12-22  Gilbert  -  Made cbuf array allocatable so that there
C                         would be no hard wired size limit.
C
C DART $Id$
C
C USAGE:  grabbufr inputBUFRfile ouputBUFRfile
C
C   INPUT FILES:
C     unit 11  - Input BUFR file.
C
C   OUTPUT FILES:
C     unit 51  - Output BUFR file.
C
C   SUBPROGRAMS CALLED: (LIST ALL CALLED FROM ANYWHERE IN CODES)
C     UNIQUE:    - lenbufr
C     LIBRARY:
C       System   - getarg stat
C       W3LIB    - errexit gbyte sbyte
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
      character,allocatable :: cbuf(:)
      CHARACTER(len=80) :: infile,outfile
      character(len=4) :: bufr='BUFR',ctemp,csec0
      INTEGER(4)       narg,iargc,JSTAT(100)
      integer findbufr, i, INDEXVAL, rc
      character*1 byte(8)
      integer :: STAT
 
      data i1/11/,i2/51/,newed/2/

      call wrdlen
C
c liu 03/16/2005
C  GET Filename ARGUMENTS
C
C  if your machine does not support the iargc or getarg functions
C  comment this section out and use the hardcoded filenames below.
      NARG=IARGC()
      IF(NARG.NE.2) THEN
        PRINT *,'grabbufr:  Incorrect usage'
        PRINT *,'Usage: grabbufr inputBUFRfile ouputBUFRfile'
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
c liu 03/16/2005
C
      PRINT*,'INPUT FILE: ',infile
      PRINT*,'OUTPUT FILE: ',outfile

C
C  Use STAT function to get size of input BUFR file
C
c  the original code was:
c      IF (STAT(infile,JSTAT).NE.0) THEN
c  but this failed on the linux systems, at least with intel
c  as the compiler. but changing it so the return code was assigned to 
c  an integer variable before being tested seemed to fix the problem.
      rc = STAT(infile,JSTAT)
      IF (rc.NE.0) THEN
         PRINT*,'ERROR IN FUNCTION STAT GETTING FILE INFO, RC = ',rc
         CALL EXIT(99)
      ELSE
c        If this program has an error, or if the output files are missing data
c        run the stat_test program in this same directory and see if the
c        indexval is correct for this machine/compiler.  fix and recompile
c        if not.
c
c        The following line may need to be changed for different machines
c        and compilers (and sometimes even between 32 and 64 bit versions of
c        the same compiler).
         INDEXVAL = 8       ! seems to be right for most current compilers
c        INDEXVAL = 12      ! for old versions of pgf90 32 bytes
c        INDEXVAL = 13      ! for old versions of pgf90 64 bytes

         KBYTES = JSTAT(INDEXVAL)
         PRINT *,'NUMBER OF BYTES IN INPUT FILE = ',KBYTES
      ENDIF
C
C  Allocate array cbuf to store input file in memory.
C
      allocate(cbuf(kbytes),stat=istat)
      IF (istat.ne.0) THEN
        PRINT*,' ERROR Allocating ',kbytes,' bytes to read in file ',
     &          infile
        CALL EXIT(4)
      ENDIF
C
C  Read input BUFR file into cbuf
C
      open(i1,recl=kbytes,file=infile,access='direct')
      read(i1,rec=1) (cbuf(j),j=1,kbytes)
C
C  Open output BUFR file
C
      open(i2,file=outfile,access='sequential',form='unformatted')
 
      ibeg=1
      icnt=0
C
C  Process each BUFR message in the input file.
C
      do
C  Search for next BUFR message
        ipos=findbufr(cbuf,ibeg,kbytes)
        if (ipos.eq.0) exit
C        ibeg=ibeg+ipos-1
        ibeg=ipos
        icnt=icnt+1
C  Extract BUFR edition number
        call gbyte(cbuf(ibeg),ied,56,8)
C  Calculate length of BUFR message
        if (ied.le.1) then
          ilen=lenbufr(cbuf(ibeg))
        else
          call gbyte(cbuf(ibeg),ilen,32,24)
        endif
C  Check ending 7777 to see if we have a complete BUFR message
        iend=ibeg+ilen-1
C       ctemp=cbuf(iend-3)//cbuf(iend-2)//cbuf(iend-1)//cbuf(iend)
        CALL CHRTRNA(CTEMP,CBUF(IEND-3),4)
        if ( ctemp.eq.'7777') then
C  If BUFR message is edition 0 or 1, convert to edition 2 format
          if (ied.le.1) then
            call sbyte(ctemp,ilen+4,0,24)
            call sbyte(ctemp,newed,24,8)
            write(i2) bufr,ctemp,(cbuf(j),j=ibeg+4,iend)
          else
            write(i2) (cbuf(j),j=ibeg,iend),(byte(j),j=1,8-mod(ilen,8))
          endif
C           print *,' BUFR message ',icnt,' was copied. ',ilen,
C     &             ' bytes: from ',ibeg,' to ',iend
           ibeg=iend
        else
           print *,' Invalid BUFR message ',icnt,' at ',ibeg
           ibeg=ibeg+1
C           ibeg=ibeg+ilen
           icnt=icnt-1
        endif
      enddo
 
      print *,'grabbufr: ',icnt,' BUFR messages copied. '
 
      stop
      end
 
      integer function lenbufr(cbufr)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .                                       .
C SUBPROGRAM:    lenbufr
C   PRGMMR: Gilbert          ORG: NP11        DATE: 99-07-13
C
C ABSTRACT: Calculates the length of a given BUFR message in bytes.
C
C PROGRAM HISTORY LOG:
C 1999-07-13  Gilbert
C
C USAGE:    integer function lenbufr(cbufr)
C   INPUT ARGUMENT LIST:
C     cbufr    - Input BUFR message
C
C   RETURN VALUE:
C     lenbufr   - length of BUFR message in bytes
C
C REMARKS: LIST CAVEATS, OTHER HELPFUL HINTS OR INFORMATION
C
C ATTRIBUTES:
C   LANGUAGE: Fortran 90
C   MACHINE:  IBM SP
C
C$$$
        character*(*) cbufr
        integer ipos,itemp,isec
 
        lenbufr=4                                 !  Section 0
        ipos=32
        call gbyte(cbufr,itemp,ipos,24)           !  Section 1
        lenbufr=lenbufr+itemp                     !
 
        call gbyte(cbufr,isec,88,1)
        if (isec.eq.1) then
          ipos=lenbufr*8
          call gbyte(cbufr,itemp,ipos,24)         !  Section 2,
          lenbufr=lenbufr+itemp                   !  if included.
        endif
        ipos=lenbufr*8
        call gbyte(cbufr,itemp,ipos,24)           !  Section 3
        lenbufr=lenbufr+itemp                     !
        ipos=lenbufr*8
        call gbyte(cbufr,itemp,ipos,24)           !  Section 4
        lenbufr=lenbufr+itemp                     !
        lenbufr=lenbufr+4                         !  Section 5
 
      return
      end
 
      integer function findbufr(cbufr,ibeg,iend)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .                                       .
C SUBPROGRAM:    findbufr
C   PRGMMR: Gilbert          ORG: NP11        DATE: 99-12-22
C
C ABSTRACT: Finds the first occurence of string 'BUFR' in character
C           array cbufr starting at element ibeg and searching to
C           element iend and returns the element of the array
C           where the string begins.  If the string 'BUFR' is not found,
C           findbufr returns 0.
C
C PROGRAM HISTORY LOG:
C 1999-12-22  Gilbert
C
C USAGE:    integer function lenbufr(cbufr,ibeg,iend)
C   INPUT ARGUMENT LIST:
C     cbufr    - Input character buffer
C     ibeg     - Array element number to begin search.
C     iend     - Array element number to end search.
C
C   RETURN VALUE:
C     findbufr   - Array element number where string 'BUFR' begins
C
C REMARKS: LIST CAVEATS, OTHER HELPFUL HINTS OR INFORMATION
C
C ATTRIBUTES:
C   LANGUAGE: Fortran 90
C   MACHINE:  IBM SP
C
C$$$
      character cbufr(iend)
      character(len=4) :: bufr='BUFR',ctemp
 
      findbufr=0
      do i=ibeg,iend-3
C       ctemp=cbufr(i)//cbufr(i+1)//cbufr(i+2)//cbufr(i+3)
        CALL CHRTRNA(CTEMP,CBUFR(I),4)
        if ( ctemp .eq. bufr ) then
          findbufr=i
          return
        endif
      enddo
 
      return
      end
C-----------------------------------------------------------------------
C THIS PROGRAM WRITTEN BY.....
C             DR. ROBERT C. GAMMILL, CONSULTANT
C             NATIONAL CENTER FOR ATMOSPHERIC RESEARCH
C             MAY 1972
C
C             CHANGES FOR SiliconGraphics IRIS-4D/25
C             SiliconGraphics 3.3 FORTRAN 77
C             March 1991, RUSSELL E. JONES
C             NATIONAL WEATHER SERVICE
C
C THIS IS THE FORTRAN VERSION OF GBYTE
C
C-----------------------------------------------------------------------
C
C SUBROUTINE GBYTE (IPACKD,IUNPKD,NOFF,NBITS)
C
C PURPOSE                TO UNPACK A BYTE INTO A TARGET WORD.  THE
C                        UNPACKED BYTE IS RIGHT-JUSTIFIED IN THE
C                        TARGET WORD, AND THE REMAINDER OF THE
C                        WORD IS ZERO-FILLED.
C
C USAGE                  CALL GBYTE(IPACKD,IUNPKD,NOFF,NBITS)
C
C ARGUMENTS
C
C ON INPUT               IPACKD
C                          THE WORD OR ARRAY CONTAINING THE BYTE TO BE
C                          UNPACKED.
C
C                        IUNPKD
C                          THE WORD WHICH WILL CONTAIN THE UNPACKED
C                          BYTE.
C
C                        NOFF
C                          THE NUMBER OF BITS TO SKIP, LEFT TO RIGHT,
C                          IN 'IPACKD' IN ORDER TO LOCATE THE BYTE
C                          TO BE UNPACKED.
C
C                        NBITS
C                          NUMBER OF BITS IN THE BYTE TO BE UNPACKED.
C                          MAXIMUM OF 64 BITS ON 64 BIT MACHINE, 32
C                          BITS ON 32 BIT MACHINE.
C
C ON OUTPUT              IUNPKD
C                          CONTAINS THE REQUESTED UNPACKED BYTE.
C-----------------------------------------------------------------------
      SUBROUTINE GBYTE(IPACKD,IUNPKD,NOFF,NBITS)
 
      INTEGER    IPACKD(*)
      INTEGER    IUNPKD
      DATA IFIRST/1/
      SAVE IFIRST
 
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
 
      IF(IFIRST.EQ.1) THEN
         CALL WRDLEN
         IFIRST = 0
      ENDIF
 
      IBIT = NOFF
      CALL UPB(IUNPKD,NBITS,IPACKD,IBIT)
 
      RETURN
      END
C-----------------------------------------------------------------------
C THIS PROGRAM WRITTEN BY.....
C             DR. ROBERT C. GAMMILL, CONSULTANT
C             NATIONAL CENTER FOR ATMOSPHERIC RESEARCH
C             JULY 1972
C
C THIS IS THE FORTRAN 32 bit VERSION OF SBYTE.
C             Changes for SiliconGraphics IRIS-4D/25
C             SiliconGraphics 3.3 FORTRAN 77
C             MARCH 1991  RUSSELL E. JONES
C             NATIONAL WEATHER SERVICE
C
C-----------------------------------------------------------------------
C
C SUBROUTINE SBYTE (IPACKD,IUNPKD,NOFF,NBITS)
C
C PURPOSE                GIVEN A BYTE, RIGHT-JUSTIFIED IN A WORD, TO
C                        PACK THE BYTE INTO A TARGET WORD OR ARRAY.
C                        BITS SURROUNDING THE BYTE IN THE TARGET
C                        AREA ARE UNCHANGED.
C
C USAGE                  CALL SBYTE (IPACKD,IUNPKD,NOFF,NBITS)
C
C ARGUMENTS
C ON INPUT               IPACKD
C                          THE WORD OR ARRAY WHICH WILL CONTAIN THE
C                          PACKED BYTE.  BYTE MAY CROSS WORD BOUNDARIES.
C
C                        IUNPKD
C                          THE WORD CONTAINING THE RIGHT-JUSTIFIED BYTE
C                          TO BE PACKED.
C
C                        NOFF
C                          THE NUMBER OF BITS TO SKIP, LEFT TO RIGHT,
C                          IN 'IPACKD' IN ORDER TO LOCATE WHERE THE
C                          BYTE IS TO BE PACKED.
C
C                        NBITS
C                          NUMBER OF BITS IN THE BYTE TO BE PACKED.
C                          MAXIMUM OF 64 BITS ON 64 BIT MACHINE, 32
C                          BITS ON 32 BIT MACHINE.
C
C ON OUTPUT              IPACKD
C                          WORD OR CONSECUTIVE WORDS CONTAINING THE
C                          REQUESTED BYTE.
C-----------------------------------------------------------------------
      SUBROUTINE SBYTE(IPACKD,IUNPKD,NOFF,NBITS)
 
      INTEGER    IUNPKD
      INTEGER    IPACKD(*)
      DATA IFIRST/1/
      SAVE IFIRST
 
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
 
      IF(IFIRST.EQ.1) THEN
         CALL WRDLEN
         IFIRST = 0
      ENDIF
 
      IBIT = NOFF
      CALL PKB(IUNPKD,NBITS,IPACKD,IBIT)
 
      RETURN
      END

C <next few lines under version control, do not edit>
C $URL$
C $Id$
C $Revision$
C $Date$
