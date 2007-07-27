 
#  ------------------------------------------------------------------------
#  This script will make cwordsh.x which FORTRAN "blocks" or "unblocks" 
#  BUFR files on a number of standard computing platforms. Stictly speaking,
#  real BUFR files are "unblocked". That is, they contain a byte stream
#  containing only allowable BUFR constructs.
#
#  On some platforms it is advantagous to use the FORTRAN
#  blocked structure for I/O efficiency, and on some platforms, when
#  using FORTRAN I/O, the unblocked structure is FORTRAN UN-readable.
#
#  NOTE: The script is set up to run in the Bourne shell. If you are a
#  C-shell user, enter 'sh ./cwordsh'.
#  ------------------------------------------------------------------------

#  ------------------------------------------------------------------------
# <next few lines under version control, do not edit>
# $URL$
# $Id$
# $Revision$
# $Date$
#  ------------------------------------------------------------------------
 
set -eua
 
#  ------------------------------------------------------------------------
#  CPLAT - platform type (sgi,linux,aix,sun,hp,cray,etc.)
#  ------------------------------------------------------------------------
 
CPLAT=macosx
BUFRLIB=../lib/bufrlib.a
 
#  different platforms use different link name protocols
#  -----------------------------------------------------
 
if [ $CPLAT = sgi ]
then
   openrb=openrb_
   openwb=openwb_
   crdbfr=crdbufr_
   cwrbfr=cwrbufr_
   lenmsg=lenm_
   cc=cc; ff=f77
elif [ $CPLAT = linux ]
then
   openrb=openrb_
   openwb=openwb_
   crdbfr=crdbufr_
   cwrbfr=cwrbufr_
   lenmsg=lenm_
   cc=gcc; ff=g77
elif [ $CPLAT = intel ]
then
   openrb=openrb_
   openwb=openwb_
   crdbfr=crdbufr_
   cwrbfr=cwrbufr_
   lenmsg=lenm_
   cc=icc; ff=ifort 
elif [ $CPLAT = aix ]
then
   openrb=openrb
   openwb=openwb
   crdbfr=crdbufr
   cwrbfr=cwrbufr
   lenmsg=lenm
   cc=cc; ff=f77
elif [ $CPLAT = sun ]
then
   openrb=openrb_
   openwb=openwb_
   crdbfr=crdbufr_
   cwrbfr=cwrbufr_
   lenmsg=lenm_
   cc=cc; ff=f77
elif [ $CPLAT = hp  ]
then
   openrb=openrb
   openwb=openwb
   crdbfr=crdbufr
   cwrbfr=cwrbufr
   lenmsg=lenm
   cc=cc; ff=f77
elif [ $CPLAT = macosx  ]
then
   openrb=openrb_
   openwb=openwb_
   crdbfr=crdbufr_
   cwrbfr=cwrbufr_
   lenmsg=lenm_
   cc=gcc; ff=gfortran
elif [ $CPLAT = cray ]
then
   openrb=OPENRB
   openwb=OPENWB
   crdbfr=CRDBUFR
   cwrbfr=CWRBUFR
   lenmsg=LENM
   cc=cc; ff=f90
fi
 
#  compile the c part of the program
#  ---------------------------------
 
cat <<eof>ccwords.c; $cc -c ccwords.c
#include <stdio.h>
FILE *pb;
void $openrb (ufile) char *ufile; { pb = fopen( ufile , "rb" ); }
void $openwb (ufile) char *ufile; { pb = fopen( ufile , "wb" ); }
int  $crdbfr (bufr)
int  *bufr;
{  int  nwrd; int  nb;
   nb = sizeof(bufr);
   if((nwrd=fread(bufr,nb,8/nb,pb))!=0)
   {  nwrd = $lenmsg(bufr);
      fread(bufr+8/nb,nb,nwrd-8/nb,pb);
      return nwrd;
   }
   else
      return -1;
}
int  $cwrbfr (bufr)
int  *bufr;
{  int  nwrd; int  nb;
   nb = sizeof(bufr);
   nwrd = $lenmsg(bufr);
   fwrite(bufr,nb,nwrd,pb);
}
eof
 
#  compile the fortran part of the program
#  ---------------------------------------
 
cat <<eof>fcwords.f; $ff -c fcwords.f 
      program fcwords
      common /hrdwrd/ nb,nbitw,nrev,iord(8)
      character*80 bfile,ufile
      character*8  cword
      character*1  zbyte
      dimension    mbay(3000),iufile(20)
      equivalence  (zbyte,izero)
      equivalence  (ufile,iufile)
      integer      crdbufr,cwrbufr
      data         izero/0/

      do i=1,80
      ufile(i:i) = zbyte
      enddo
 
      read(5,'(a)') cword
      if(cword.eq.'block') then
         read(5,'(a)') ufile
         read(5,'(a)') bfile
      elseif(cword.eq.'unblk') then
         read(5,'(a)') bfile
         read(5,'(a)') ufile
      else
         print*,'cword must be block or unblk'
         call exit(8)
      endif

      do i=1,80
      if(ufile(i:i).eq.' ') ufile(i:i) = zbyte
      enddo
 
      open(8,file=bfile,form='unformatted')
      call wrdlen

      if(cword.eq.'block') then
         print*,'blocking ',ufile,' to ',bfile
         call openrb(iufile)
         do while(crdbufr(mbay).ge.0)
         write(8) (mbay(i),i=1,lenm(mbay))
         enddo
      endif
 
      if(cword.eq.'unblk') then
         print*,'unblocking ',bfile,' to ',ufile
         call openwb(iufile)
1        read(8,end=2)(mbay(i),i=1,8/nb),(mbay(i),i=1+8/nb,lenm(mbay))
         iwt = cwrbufr(mbay)
         goto 1
2        continue
      endif
 
      stop
      end
c-----------------------------------------------------------------------
      function lenm(mbay)
      common /hrdwrd/ nb,nbitw,nrev,iord(8)
      dimension mbay(*)
      lenm = (1+iupb(mbay,5,24)/8)*8/nb
      return
      end
c-----------------------------------------------------------------------
eof
 
#  link and load the executable
#  ----------------------------

CWRD=.
 
$ff -o $CWRD/cwordsh.x fcwords.o ccwords.o $BUFRLIB

#  clean up
#  --------

rm -f fcwords.[fo] ccwords.[co]
