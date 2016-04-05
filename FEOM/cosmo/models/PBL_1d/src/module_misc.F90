MODULE module_misc
! miscellaneous subroutines and functions for WRF compatibility
!
! DART $Id$
!
CONTAINS

  LOGICAL FUNCTION wrf_dm_on_monitor()
    wrf_dm_on_monitor = .TRUE.
  END FUNCTION wrf_dm_on_monitor

! IBM libmassv compatibility library
! 

!#ifndef NATIVE_MASSV
      subroutine vdiv(z,x,y,n)
      real*8 x(*),y(*),z(*)
      do 10 j=1,n
      z(j)=x(j)/y(j)
   10 continue
      return
      end subroutine

      subroutine vsdiv(z,x,y,n)
      real*4 x(*),y(*),z(*)
      do 10 j=1,n
      z(j)=x(j)/y(j)
   10 continue
      return
      end subroutine

      subroutine vexp(y,x,n)
      real*8 x(*),y(*)
      do 10 j=1,n
      y(j)=exp(x(j))
   10 continue
      return
      end subroutine

      subroutine vsexp(y,x,n)
      real*4 x(*),y(*)
      do 10 j=1,n
      y(j)=exp(x(j))
   10 continue
      return
      end subroutine

      subroutine vlog(y,x,n)
      real*8 x(*),y(*)
      do 10 j=1,n
      y(j)=log(x(j))
   10 continue
      return
      end subroutine

      subroutine vslog(y,x,n)
      real*4 x(*),y(*)
      do 10 j=1,n
      y(j)=log(x(j))
   10 continue
      return
      end subroutine

      subroutine vrec(y,x,n)
      real*8 x(*),y(*)
      do 10 j=1,n
      y(j)=1.d0/x(j)
   10 continue
      return
      end subroutine

      subroutine vsrec(y,x,n)
      real*4 x(*),y(*)
      do 10 j=1,n
      y(j)=1.d0/x(j)
   10 continue
      return
      end subroutine

      subroutine vrsqrt(y,x,n)
      real*8 x(*),y(*)
      do 10 j=1,n
      y(j)=1.d0/sqrt(x(j))
   10 continue
      return
      end subroutine

      subroutine vsrsqrt(y,x,n)
      real*4 x(*),y(*)
      do 10 j=1,n
      y(j)=1.d0/sqrt(x(j))
   10 continue
      return
      end subroutine

      subroutine vsincos(x,y,z,n)
      real*8 x(*),y(*),z(*)
      do 10 j=1,n
      x(j)=sin(z(j))
      y(j)=cos(z(j))
   10 continue
      return
      end subroutine

      subroutine vssincos(x,y,z,n)
      real*4 x(*),y(*),z(*)
      do 10 j=1,n
      x(j)=sin(z(j))
      y(j)=cos(z(j))
   10 continue
      return
      end subroutine

      subroutine vsqrt(y,x,n)
      real*8 x(*),y(*)
      do 10 j=1,n
      y(j)=sqrt(x(j))
   10 continue
      return
      end subroutine

      subroutine vssqrt(y,x,n)
      real*4 x(*),y(*)
      do 10 j=1,n
      y(j)=sqrt(x(j))
   10 continue
      return
      end subroutine

      subroutine vtan(y,x,n)
      real*8 x(*),y(*)
      do 10 j=1,n
      y(j)=tan(x(j))
   10 continue
      return
      end subroutine

      subroutine vstan(y,x,n)
      real*4 x(*),y(*)
      do 10 j=1,n
      y(j)=tan(x(j))
   10 continue
      return
      end subroutine

      subroutine vatan2(z,y,x,n)
      real*8 x(*),y(*),z(*)
      do 10 j=1,n
      z(j)=atan2(y(j),x(j))
   10 continue
      return
      end subroutine

      subroutine vsatan2(z,y,x,n)
      real*4 x(*),y(*),z(*)
      do 10 j=1,n
      z(j)=atan2(y(j),x(j))
   10 continue
      return
      end subroutine

      subroutine vasin(y,x,n)
      real*8 x(*),y(*)
      do 10 j=1,n
      y(j)=asin(x(j))
   10 continue
      return
      end subroutine

      subroutine vsin(y,x,n)
      real*8 x(*),y(*)
      do 10 j=1,n
      y(j)=sin(x(j))
   10 continue
      return
      end subroutine

      subroutine vssin(y,x,n)
      real*4 x(*),y(*)
      do 10 j=1,n
      y(j)=sin(x(j))
   10 continue
      return
      end subroutine

      subroutine vacos(y,x,n)
      real*8 x(*),y(*)
      do 10 j=1,n
      y(j)=acos(x(j))
   10 continue
      return
      end subroutine

      subroutine vcos(y,x,n)
      real*8 x(*),y(*)
      do 10 j=1,n
      y(j)=cos(x(j))
   10 continue
      return
      end subroutine

      subroutine vscos(y,x,n)
      real*4 x(*),y(*)
      do 10 j=1,n
      y(j)=cos(x(j))
   10 continue
      return
      end subroutine

      subroutine vcosisin(y,x,n)
      complex*16 y(*)
      real*8 x(*)
      do 10 j=1,n
      y(j)=dcmplx(cos(x(j)),sin(x(j)))
   10 continue
      return
      end subroutine

      subroutine vscosisin(y,x,n)
      complex*8 y(*)
      real*4 x(*)
      do 10 j=1,n
      y(j)= cmplx(cos(x(j)),sin(x(j)))
   10 continue
      return
      end subroutine

      subroutine vdint(y,x,n)
      real*8 x(*),y(*)
      do 10 j=1,n
!     y(j)=dint(x(j))
      y(j)=int(x(j))
   10 continue
      return
      end subroutine

      subroutine vdnint(y,x,n)
      real*8 x(*),y(*)
      do 10 j=1,n
!     y(j)=dnint(x(j))
      y(j)=nint(x(j))
   10 continue
      return
      end subroutine

      subroutine vlog10(y,x,n)
      real*8 x(*),y(*)
      do 10 j=1,n
      y(j)=log10(x(j))
   10 continue
      return
      end subroutine

!      subroutine vlog1p(y,x,n)
!      real*8 x(*),y(*)
!      interface
!        real*8 function log1p(%val(x))
!          real*8 x
!        end function log1p subroutine
!      end interface subroutine
!      do 10 j=1,n
!      y(j)=log1p(x(j))
!   10 continue
!      return
!      end subroutine

      subroutine vcosh(y,x,n)
      real*8 x(*),y(*)
      do 10 j=1,n
      y(j)=cosh(x(j))
   10 continue
      return
      end subroutine

      subroutine vsinh(y,x,n)
      real*8 x(*),y(*)
      do 10 j=1,n
      y(j)=sinh(x(j))
   10 continue
      return
      end subroutine

      subroutine vtanh(y,x,n)
      real*8 x(*),y(*)
      do 10 j=1,n
      y(j)=tanh(x(j))
   10 continue
      return
      end subroutine

!      subroutine vexpm1(y,x,n)
!      real*8 x(*),y(*)
!      interface
!        real*8 function expm1(%val(x))
!          real*8 x
!        end function expm1 subroutine
!      end interface  subroutine
!      do 10 j=1,n
!      y(j)=expm1(x(j))
!   10 continue
!      return
!      end subroutine


      subroutine vsasin(y,x,n)
      real*4 x(*),y(*)
      do 10 j=1,n
      y(j)=asin(x(j))
   10 continue
      return
      end subroutine

      subroutine vsacos(y,x,n)
      real*4 x(*),y(*)
      do 10 j=1,n
      y(j)=acos(x(j))
   10 continue
      return
      end subroutine

      subroutine vscosh(y,x,n)
      real*4 x(*),y(*)
      do 10 j=1,n
      y(j)=cosh(x(j))
   10 continue
      return
      end subroutine

!      subroutine vsexpm1(y,x,n)
!      real*4 x(*),y(*)
!      interface
!        real*8 function expm1(%val(x))
!          real*8 x
!        end function expm1 subroutine
!      end interface subroutine
!      do 10 j=1,n
!      y(j)=expm1(real(x(j),8))
!   10 continue
!      return
!      end subroutine

      subroutine vslog10(y,x,n)
      real*4 x(*),y(*)
      do 10 j=1,n
      y(j)=log10(x(j))
   10 continue
      return
      end subroutine

!      subroutine vslog1p(y,x,n)
!      real*4 x(*),y(*)
!      interface
!        real*8 function log1p(%val(x))
!          real*8 x
!        end function log1p subroutine
!      end interface subroutine
!      do 10 j=1,n
!      y(j)=log1p(real(x(j),8))
!   10 continue
!      return
!      end subroutine


      subroutine vssinh(y,x,n)
      real*4 x(*),y(*)
      do 10 j=1,n
      y(j)=sinh(x(j))
   10 continue
      return
      end subroutine

      subroutine vstanh(y,x,n)
      real*4 x(*),y(*)
      do 10 j=1,n
      y(j)=tanh(x(j))
   10 continue
      return
      end subroutine
!#endif subroutine

      subroutine vspow(z,y,x,n)
      real*4 x(*),y(*),z(*)
      do 10 j=1,n
      z(j)=y(j)**x(j)
   10 continue
      return
      end subroutine

      subroutine vpow(z,y,x,n)
      real*8 x(*),y(*),z(*)
      do 10 j=1,n
      z(j)=y(j)**x(j)
   10 continue
      return
      end subroutine

END MODULE module_misc

