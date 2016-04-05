MODULE byte_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

  use        types_mod, only : r4, r8

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

  public :: concat_bytes1
  public :: concat_bytes1_sign
  public :: to_positive
  public :: byte_to_word_signed
  public :: get_characteristic
  public :: to_float1
  public :: from_float1
  public :: byte_to_word
  public :: word_to_byte
  public :: get_word
  
CONTAINS

  INTEGER FUNCTION concat_bytes1(word,n,check_negative)

    ! concat_bytes1 returns the concatenated value of a set of given bytes

    INTEGER(kind=1),INTENT(in) :: word(:)
    INTEGER,INTENT(in)         :: n
    LOGICAL,INTENT(in)         :: check_negative
    INTEGER                    :: iresult,i,maxexp,work(n)

    work=word(1:n)
    iresult=0
    maxexp=n-1
    DO i=0,n-1
      IF (work(i+1).LT.0) work(i+1)=256+work(i+1)
      iresult=iresult+work(i+1)*256**(maxexp-i)
!      iresult=iresult+word(i+is)*256**(i)
    END DO
    IF (check_negative) THEN
      IF (iresult.LT.0) iresult=iresult+256
    END IF

    concat_bytes1=iresult
    RETURN

  END FUNCTION concat_bytes1

  INTEGER FUNCTION concat_bytes1_sign(word,n,check_negative)
    INTEGER(kind=1),INTENT(in) :: word(:)
    INTEGER,INTENT(in)         :: n
    LOGICAL,INTENT(in)         :: check_negative
    INTEGER                    :: iresult,i,maxexp,work(n),isign

    work=word(1:n)
    isign=1
    IF (work(1).LT.0) work(1)=256+work(1)
    IF ((work(1)/128).EQ.1) THEN
      work(1)=work(1)-128
      isign=-1
    END IF

    iresult=0
    maxexp=n-1
    DO i=0,n-1
      IF (work(i+1).LT.0) work(i+1)=256+work(i+1)
      iresult=iresult+work(i+1)*256**(maxexp-i)
!      iresult=iresult+word(i+is)*256**(i)
    END DO
    IF (check_negative) THEN
      IF (iresult.LT.0) iresult=iresult+256
    END IF

    concat_bytes1_sign=isign*iresult
    RETURN

  END FUNCTION concat_bytes1_sign

  FUNCTION to_positive(b)

    INTEGER                    :: to_positive
    INTEGER(kind=1),INTENT(in) :: b

    IF (b>0) THEN
      to_positive=b
    ELSE
      to_positive=b+256
    END IF

    RETURN
  END FUNCTION to_positive

  ! This routine converts a vector of bytes to a integer word checkin for a sign bit

  SUBROUTINE byte_to_word_signed(b,w,n)
    INTEGER(kind=1),INTENT(in) :: b(:)
    INTEGER,        INTENT(in) :: n
    INTEGER,INTENT(out)        :: w
    INTEGER                    :: ibin,work(1:n),isign

    isign=1
    work(1:n)=b(1:n)
    if (b(1)<0) work(1)=256+work(1)
    IF ((work(1)/128).EQ.1) THEN
      work(1)=work(1)-128
      isign=-1
    END IF

    w=0
    DO ibin=1,n
      w=w+work(ibin)*256**(n-ibin)
    END DO
    w=w*isign

    RETURN

  END SUBROUTINE byte_to_word_signed

  FUNCTION get_characteristic(b) result (c)
    
    IMPLICIT NONE
    
    INTEGER(kind=1)                      :: c
    INTEGER(kind=1),         INTENT(in)  :: b
    INTEGER                              :: w
    
    w=b

    IF (w<0) w=256+w

    IF (w>=128) THEN
      w=w-128
    END IF

    c=w
    RETURN
    
  END FUNCTION get_characteristic

  REAL(r8) FUNCTION to_float1(b)
    
    IMPLICIT NONE
    
    INTEGER(kind=1),         INTENT(in)  :: b(4)
    INTEGER                              :: w(4),imant,isign,ibin

    w=b

    DO ibin=1,4
      IF (w(ibin)<0) w(ibin)=256+w(ibin)
    END DO

    isign=1
    IF (w(1)>=128) THEN
      isign=-1
      w(1)=w(1)-128
    END IF

    imant=w(2)*256**2+w(3)*256**1+w(4)
    to_float1=isign*(2.**(-24))*imant*(16.**(w(1)-64))
    RETURN
    
  END FUNCTION to_float1

  FUNCTION from_float1(f,char) RESULT (b4)
    
    IMPLICIT NONE
    
    INTEGER(kind=1)     :: b4(4)

    REAL(r8),       INTENT(in) :: f
    INTEGER(kind=1),INTENT(in) :: char

    REAL(r8)            :: wf,wf0
    INTEGER             :: ibin
    INTEGER             :: ichar,isign
    INTEGER(kind=8)     :: imant,w

    if (f<0.) THEN
      isign=1
    ELSE
      isign=0
    END if

    wf=ABS(f)
    wf0=wf*(2.**24.)
!    n=int(log(wf0)/log(2.))
!    n1=4*int(n/4.+0.5)
    
!    ichar=(n1+232)/4
    ichar=char

!    imant=int((wf0/(2.**n1)*(2.**24)),8)
    imant=int((wf0/(16**(ichar-64))),8)

    b4(1)=ichar+isign*128

    w=imant

    DO ibin=4,2,-1
      b4(ibin)=MOD(w,256)
      w=w/256
    END DO

    RETURN
    
  END FUNCTION from_float1

  ! This routine converts a 2-element vector of bytes to an integer word

  INTEGER FUNCTION byte_to_word(b)
    INTEGER(kind=1),INTENT(in)  :: b(2)
    INTEGER                     :: w(2),r,n
    INTEGER                     :: ibin

    n=2
    r=0
    w=b
    DO ibin=1,n
      IF (w(ibin)<0.) THEN
        w(ibin)=256+w(ibin)
      END IF
      r=r+w(ibin)*256**(n-ibin)
    END DO
!    IF (w<0.) w=w+256**n

    byte_to_word=r
    RETURN

  END FUNCTION byte_to_word

  FUNCTION word_to_byte(d) RESULT (b2)
    INTEGER(kind=1)    :: b2(2)
    INTEGER,INTENT(in) :: d
    INTEGER            :: ibin,w

    w=d

    DO ibin=2,1,-1
      b2(ibin)=MOD(w,256)
      w=w/256
    END DO

    RETURN
  END FUNCTION word_to_byte

  INTEGER FUNCTION get_word(b)
    INTEGER(kind=1),INTENT(in) :: b(:)
    INTEGER                    :: w,n
    INTEGER                    :: ibin

    n=4
    w=0
    DO ibin=1,n
      w=w+b(ibin)*256**(n-ibin)
    END DO
    get_word=w
    RETURN

  END FUNCTION get_word

END MODULE byte_mod
