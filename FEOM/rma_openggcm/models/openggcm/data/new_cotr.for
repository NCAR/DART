!This module contains (pure) f90 code in fixed-format with modules.  
!There is no fppn necessary to compile it. For the sake of future 
!projects, it might be nice to keep it that way :-) -- MLG
! 
!As such, it needs to be
!compiled BEFORE any other code which uses it.  It has no dependencies,
!so there is no harm in compiling it first.

!doc -
!doc - These modules contains stuff originally in mhd-cotr.for
!doc -
!doc - Original AUTHOR:
!doc -           Joachim Raeder, IGPP/UCLA, June 1994
!doc - The new modules are updated to minimize global data
!doc - and to use some of the neater features of fortran90
!doc -
!doc - Updating AUTHOR:
!doc -           Matthew L Gilson, SSC/UNH, January 2013


      module date_lib

      integer,private,parameter,dimension(12,2) :: mdays = RESHAPE(
     *     (/0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334,
     *       0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335/),
     *     (/12,2/))


      contains
cc----------------------------------------------------------------
      subroutine epoch1966(dsecs,iy,mo,id,ih,mi,sec,iimod)
c----------------------------------------------------------------
cdoc-
cdoc-FUNCTION:
cdoc-  converts seconds since JAN 1, 1966 0UT to year,month,day,hour,
cdoc-  minute and second or vice versa
cdoc-
cdoc-PARAMETERS:
cdoc-  dsec (in/out): seconds of epoch 1966
cdoc-  iy,mo,id,ih,mi,se (in):  time in UT
cdoc-  iimod (in,integer): flag
cdoc-                     iimod=0: convert iy,... to dsec
cdoc-                     iimod!=0: convert dsec to iy,...
cdoc-
cdoc-NOTE:
cdoc-  works only for times after JAN 1, 1966, i.e. dsec positive
cdoc-
      implicit none

      integer iy,mo,id,ih,mi,iimod
      real*8 dsecs
      real sec

      integer i,ij,jj
      real*8 ss,ds,dso
      integer,parameter,dimension(0:11) ::  
     *     rdays=( /31,28,31,30,31,30,31,31,30,31,30,31/ )
      integer,parameter,dimension(0:11) ::  
     *     jdays = ( /31,29,31,30,31,30,31,31,30,31,30,31/)
c
      if(iimod.eq.0) then
      jj=iy
      if(jj.lt.100)jj=jj+1900
      ss=0.0
      do i=1966,jj-1
         if(mod(i,4).eq.0) ss=ss+366.0d0*24.0d0*3600.0d0
         if(mod(i,4).ne.0) ss=ss+365.0d0*24.0d0*3600.0d0
      end do
      do i=1,mo-1
         if( mod(jj,4).eq.0 ) ss=ss+dble(jdays(i-1))*24.0d0*3600.0d0
         if( mod(jj,4).ne.0 ) ss=ss+dble(rdays(i-1))*24.0d0*3600.0d0
      end do
      ss=ss+dble(id-1)*3600.0d0*24.0d0
      ss=ss+dble(ih)*3600.0d0
      ss=ss+dble(mi)*60.0d0
      ss=ss+dble(sec)
      dsecs=ss
c
      else
      ss=dsecs
      ds=0.0d0
      dso=0.0d0
      jj=1966
      do i=0,100
         if(mod(jj,4).eq.0) then
            ds=ds+(366.0d0*24.0d0*3600.0d0)
         else
            ds=ds+(365.0d0*24.0d0*3600.0d0)
         endif
         if(ds.gt.ss) exit
         dso=ds
         jj=jj+1
      end do

      ss=ss-dso
      iy=jj-1900
      ij=0
      if(mod(jj,4).eq.0)ij=1
      ds=0.0d0
      dso=0.0d0
      jj=0
      do i=0,11
         if(ij.eq.0) then
            ds=ds+dble(rdays(i))*24.0d0*3600.0d0
         else
            ds=ds+dble(jdays(i))*24.0d0*3600.0d0
         endif
         if(ds.gt.ss) exit
         dso=ds
         jj=jj+1
      end do

      ss=ss-dso
      mo=jj+1
      ds=ss/(24.0d0*3600.0d0)
      id=int(ds)
      ss=ss-dble(id)*24.0d0*3600.0d0
      id=id+1
      ds=ss/3600.0d0
      ih=int(ds)
      ss=ss-dble(ih)*3600.0d0
      ds=ss/60.0d0
      mi=int(ds)
      ss=ss-dble(mi)*60.0d0
      sec=ss
      endif
c
      return
      end subroutine epoch1966
      
c---------------------------------------------------------------
      real*8 function djul(iy,mo,id,ih,mi,se)
c---------------------------------------------------------------
cdoc-
cdoc-FUNCTION:
cdoc-  converts time to fractional julian days
cdoc-
cdoc-PARAMETERS:
cdoc-  iy,mo,id,ih,mi,se (in):  time in UT
cdoc-  djul (out,real*8):   fractional julian day
cdoc-
cdoc-NOTE:
cdoc-  expected to work for any time after 4500 BC
cdoc-
      implicit none

      integer iy,mo,id,ih,mi
      real*4 se

      integer jy
      integer i1,i2,i3,i4,ier
      real*8 d1,zero,one,xjd,uthour
      parameter(zero=0.0d0,one=1.0d0)
      jy=iy
      if(jy.lt.100)jy=jy+1900
      ier=0
      if(jy.lt.1901) ier=1
      if(jy.gt.2049) ier=7
      if(mo.lt.1.or.mo.gt.12) ier=2
      if(id.lt.1.or.id.gt.31) ier=3
      if(ih.lt.0.or.ih.gt.23) ier=4
      if(mi.lt.0.or.mi.gt.63) ier=5
!      if(se.lt.0.0.or.se.gt.60.0) ier=6
      if(ier.ne.0) then
         write(0,*)'error djul, ier= ',ier
!      stop
      endif
      uthour=dble(float(ih))+(dble(float(mi))/60.0d0)+
     *     dble(se)/3600.0d0
      i1=367*jy
      i2=(mo+9)/12
      i2=(i2+jy)*7
      i2=i2/4
      i3=(275*mo)/9
      i4=100*jy+mo
      d1=one
      if((dble(i4)-190002.5d0).lt.zero)d1=(-d1)
      xjd=dble(i1) - dble(i2) + dble(i3) + dble(id) + 1721013.5d0 
     *     + uthour/24.0d0 - 0.5d0*d1 + 0.5
      djul=xjd
      return
      end function djul

c----------------------------------------------------
      integer function julday(iy,mo,id)
c----------------------------------------------------
cdoc-
cdoc-FUNCTION:
cdoc-  convert date from day - month system to julian day system
cdoc-
cdoc-PARAMETERS:
cdoc-  iy (in,integer):  year (nn or nnnn notation)
cdoc-  mo (in,integer): month
cdoc-  id (in,integer): day of month
cdoc-  julday (out,integer): julian day of year
cdoc-
      implicit none

      integer iy,mo,id
      integer leap,jmo

      jmo=max0(1,min0(12,mo))
      leap=1
      if(mod(iy,4).eq.0)leap=2
      julday=mdays(jmo,leap)+id
      return
      end function julday


c----------------------------------------------------
      subroutine daymon(iy,jd,mo,id)
c----------------------------------------------------
cdoc-
cdoc-FUNCTION:
cdoc-  convert date from julian day system to day - month system
cdoc-
cdoc-PARAMETERS:
cdoc-  iy (in,integer):  year (nn or nnnn notation)
cdoc-  jd (in,integer):  julian day (1<=jd<=366), january 1 <> jd=1
cdoc-  mo (out,integer): month
cdoc-  id (out,integer): day of month
cdoc-
      implicit none

      integer iy,jd,mo,id
      integer leap,i,k

      leap=1
      if(mod(iy,4).eq.0)leap=2
      k=1
      do i=1,12
         if(jd.le.mdays(i,leap)) exit
         k=i
      end do

      mo=k
      id=jd-mdays(k,leap)
      return
      end subroutine daymon

      end module date_lib

! ================================================================

      module new_cotr

      real*8,parameter,private :: rad=17.45329252d-3
      real*8,parameter,private :: deg=57.29577951d0


      type transform
      integer :: iscalled
      integer :: icalled = 1
      real*8,dimension(3,3,2) :: tmpmat
      real*8,dimension(3,3) :: cmat
      real*8,dimension(3,3,-7:7) :: pmat
      real*8 :: t0
      integer :: year,month,day,hour,minute
      real :: seconds
      character(len=4) :: ccfr,ccto
      end type transform
      
      interface cotr_set
        module procedure cotr_set_1  !tranform from explicit times
        module procedure cotr_set_r8 !transform from r*8 time
      end interface cotr_set

      interface diptil
        module procedure diptil_1
        module procedure diptil_tran
      end interface diptil

      contains

      function transform_2_r8(tran)
      use date_lib
      implicit none
      type(transform),intent(IN) :: tran
      real*8 :: transform_2_r8
      call epoch1966(transform_2_r8,
     *     tran%year,tran%month,tran%day,
     *     tran%hour,tran%minute,tran%seconds,0)
      return
      end function transform_2_r8
      

!----------------------------------------------------
      subroutine xyzdeg(x,y,z,r,p,t)
!----------------------------------------------------
!doc-
!doc-FUNCTION:
!doc-  convert from cartesian to polar coordinates, angles in degrees
!doc-
!doc-PARAMETERS:
!doc-  x,y,z (in,real): cartesian coordinates
!doc-  r,p,t (out,real): polar coordinates, 
!doc-                    -180.0<=p<=180.0, 0<=t<=180.0
!doc-                    i.e. colatitude system
!doc-
      implicit none
      real x,y,z,r,p,t
      r=sqrt(x*x+y*y+z*z)
      p=0.0
      if(abs(x)+abs(y).gt.0.0)p=deg*atan2(y,x)
      t=0.0
      if(r.gt.0.0)t=deg*acos(z/r)
      return
      end subroutine xyzdeg

!----------------------------------------------------
      subroutine radxyz(r,p,t,x,y,z)
!----------------------------------------------------
!doc-
!doc-FUNCTION:
!doc-  convert from polar coordinates (angles in radian) to cartesian
!doc-
!doc-CALL SEQUENCE:
!doc-
!doc-PARAMETERS:
!doc-  r,p,t (in,real): polar coordinates, -pi<=p<=pi, 0<=t<=pi
!doc-                    i.e. colatitude system
!doc-  x,y,z (out,real): cartesian coordinates
!doc-
      implicit none
      
      real r,p,t,x,y,z

      real st

      st=sin(t)
      x=r*cos(p)*st
      y=r*sin(p)*st
      z=r*cos(t)
      return
      end subroutine radxyz

!----------------------------------------------------
      subroutine degxyz(r,p,t,x,y,z)
!----------------------------------------------------
!doc-
!doc-FUNCTION:
!doc-  convert from polar coordinates (angles in degrees) to cartesian
!doc-
!doc-PARAMETERS:
!doc-  r,p,t (in,real): polar coordinates,
!doc-                    -180.0<=p<=180.0, 0<=t<=180.0
!doc-                    i.e. colatitude system
!doc-  x,y,z (out,real): cartesian coordinates
!doc-
      implicit none

      real r,p,t,x,y,z
      real pp,tt,st

      pp=p*rad
      tt=t*rad
      st=sin(tt)
      x=r*cos(pp)*st
      y=r*sin(pp)*st
      z=r*cos(tt)
      return
      end subroutine degxyz


!----------------------------------------------------
      subroutine cotr(tran,cfr,cto,x1,y1,z1,x2,y2,z2)
!----------------------------------------------------
!     doc-  
!     doc-FUNCTION:
!     doc-  transform vector (x1,y1,z1) in coordinate system 'cfr'
!     doc-  to vector (x2,y2,z2) in coordinate system 'cto'
!     doc-  coordinate systems 'cfr' and 'cto' can be any of the following:
!     doc-  'gei'  :  geocentric equatorial inertial
!     doc-  'geo'  :  geographic
!     doc-  'gse'  :  geocentric solar ecliptic
!     doc-  'gsm'  :  geocentric solar magnetospheric
!     doc-  'sm '  :  solar magnetic
!     doc-  'mag'  :  geomagnetic
!     doc-  'mhd'  :  global mhd simulation system, like 'gse' with
!     doc-            x and y axes mirrored
!     doc-
!     doc-PARAMETERS:
!     doc-  cfr (in,character*3)   'from' coordinate system, for ex. 'gse'
!     doc-  cto (in,character*3)   'to' coordinate system, for ex. 'gsm'
!     doc-  x1,y1,z1 (in,real)     input vector
!     doc-  x2,y2,z2 (out,real)    transformed vector
!     doc-
      implicit none
      
      type(transform) :: tran
      
      character*3 cfr,cto
      real*4 x1,y1,z1,x2,y2,z2    
      integer i,j,k,l
      integer ifr,ito
      integer kmat
      integer imat(7,7,7)
      real*8 xx1,xx2,yy1,yy2,zz1,zz2
      character*3 clfr,clto
      
      data (imat(l,1,1),l=1,7)/ 0, 0, 0, 0, 0, 0, 0/
      data (imat(l,2,1),l=1,7)/-1, 0, 0, 0, 0, 0, 0/
      data (imat(l,3,1),l=1,7)/-2, 0, 0, 0, 0, 0, 0/
      data (imat(l,4,1),l=1,7)/-3,-2, 0, 0, 0, 0, 0/
      data (imat(l,5,1),l=1,7)/-4,-3,-2, 0, 0, 0, 0/
      data (imat(l,6,1),l=1,7)/-5,-1, 0, 0, 0, 0, 0/
      data (imat(l,7,1),l=1,7)/ 6,-2, 0, 0, 0, 0, 0/
      data (imat(l,1,2),l=1,7)/ 1, 0, 0, 0, 0, 0, 0/
      data (imat(l,2,2),l=1,7)/ 0, 0, 0, 0, 0, 0, 0/
      data (imat(l,3,2),l=1,7)/-2, 1, 0, 0, 0, 0, 0/
      data (imat(l,4,2),l=1,7)/-3,-2, 1, 0, 0, 0, 0/
      data (imat(l,5,2),l=1,7)/-4,-3,-2, 1, 0, 0, 0/
      data (imat(l,6,2),l=1,7)/-5, 0, 0, 0, 0, 0, 0/
      data (imat(l,7,2),l=1,7)/ 6,-2, 1, 0, 0, 0, 0/
      data (imat(l,1,3),l=1,7)/ 2, 0, 0, 0, 0, 0, 0/
      data (imat(l,2,3),l=1,7)/-1, 2, 0, 0, 0, 0, 0/
      data (imat(l,3,3),l=1,7)/ 0, 0, 0, 0, 0, 0, 0/
      data (imat(l,4,3),l=1,7)/-3, 0, 0, 0, 0, 0, 0/
      data (imat(l,5,3),l=1,7)/-4,-3, 0, 0, 0, 0, 0/
      data (imat(l,6,3),l=1,7)/-5,-1, 2, 0, 0, 0, 0/
      data (imat(l,7,3),l=1,7)/ 6, 0, 0, 0, 0, 0, 0/
      data (imat(l,1,4),l=1,7)/ 2, 3, 0, 0, 0, 0, 0/
      data (imat(l,2,4),l=1,7)/-1, 2, 3, 0, 0, 0, 0/
      data (imat(l,3,4),l=1,7)/ 3, 0, 0, 0, 0, 0, 0/
      data (imat(l,4,4),l=1,7)/ 0, 0, 0, 0, 0, 0, 0/
      data (imat(l,5,4),l=1,7)/-4, 0, 0, 0, 0, 0, 0/
      data (imat(l,6,4),l=1,7)/-5,-1, 2, 3, 0, 0, 0/
      data (imat(l,7,4),l=1,7)/ 6, 3, 0, 0, 0, 0, 0/
      data (imat(l,1,5),l=1,7)/ 2, 3, 4, 0, 0, 0, 0/
      data (imat(l,2,5),l=1,7)/-1, 2, 3, 4, 0, 0, 0/
      data (imat(l,3,5),l=1,7)/ 3, 4, 0, 0, 0, 0, 0/
      data (imat(l,4,5),l=1,7)/ 4, 0, 0, 0, 0, 0, 0/
      data (imat(l,5,5),l=1,7)/ 0, 0, 0, 0, 0, 0, 0/
      data (imat(l,6,5),l=1,7)/-5,-1, 2, 3, 4, 0, 0/
      data (imat(l,7,5),l=1,7)/ 6, 3, 4, 0, 0, 0, 0/
      data (imat(l,1,6),l=1,7)/ 1, 5, 0, 0, 0, 0, 0/
      data (imat(l,2,6),l=1,7)/ 5, 0, 0, 0, 0, 0, 0/
      data (imat(l,3,6),l=1,7)/-2, 1, 5, 0, 0, 0, 0/
      data (imat(l,4,6),l=1,7)/-3,-2, 1, 5, 0, 0, 0/
      data (imat(l,5,6),l=1,7)/-4,-3,-2, 1, 5, 0, 0/
      data (imat(l,6,6),l=1,7)/ 0, 0, 0, 0, 0, 0, 0/
      data (imat(l,7,6),l=1,7)/ 6,-2, 1, 5, 0, 0, 0/
      data (imat(l,1,7),l=1,7)/ 2,-6, 0, 0, 0, 0, 0/
      data (imat(l,2,7),l=1,7)/-1, 2, 6, 0, 0, 0, 0/
      data (imat(l,3,7),l=1,7)/-6, 0, 0, 0, 0, 0, 0/
      data (imat(l,4,7),l=1,7)/-3,-6, 0, 0, 0, 0, 0/
      data (imat(l,5,7),l=1,7)/-4,-3,-6, 0, 0, 0, 0/
      data (imat(l,6,7),l=1,7)/-5,-1, 2,-6, 0, 0, 0/
      data (imat(l,7,7),l=1,7)/ 0, 0, 0, 0, 0, 0, 0/
!     
      if(tran%icalled.eq.0.or.tran%iscalled.ne.0) then
         tran%icalled=1
         tran%ccfr='XXX'
         tran%ccto='xyz'
      endif
      tran%iscalled=0
!     
      xx1=dble(x1)
      yy1=dble(y1)
      zz1=dble(z1)
      clfr=cfr
      clto=cto
!     
      if(clfr.ne.tran%ccfr(1:3).or.clto.ne.tran%ccto(1:3)) then
         ifr=0
         if(clfr.eq.'GEI'.or.clfr.eq.'gei')ifr=1
         if(clfr.eq.'GEO'.or.clfr.eq.'geo')ifr=2
         if(clfr.eq.'GSE'.or.clfr.eq.'gse')ifr=3
         if(clfr.eq.'GSM'.or.clfr.eq.'gsm')ifr=4
         if(clfr(1:2).eq.'SM'.or.clfr(1:2).eq.'sm')ifr=5
         if(clfr.eq.'MAG'.or.clfr.eq.'mag')ifr=6
         if(clfr.eq.'MHD'.or.clfr.eq.'mhd')ifr=7
         ito=0
         if(clto.eq.'GEI'.or.clto.eq.'gei')ito=1
         if(clto.eq.'GEO'.or.clto.eq.'geo')ito=2
         if(clto.eq.'GSE'.or.clto.eq.'gse')ito=3
         if(clto.eq.'GSM'.or.clto.eq.'gsm')ito=4
         if(clto(1:2).eq.'SM'.or.clto(1:2).eq.'sm')ito=5
         if(clto.eq.'MAG'.or.clto.eq.'mag')ito=6
         if(clto.eq.'MHD'.or.clto.eq.'mhd')ito=7
         if(ifr.lt.1.or.ito.lt.1) then
            write(0,*)' error cotr: ifr,ito ',ifr,ito
            write(0,*)'             clfr=',clfr
            write(0,*)'             clto=',clto
            stop
         endif
!     
         tran%ccfr=clfr
         tran%ccto=clto
         do i=1,3
            do j=1,3
               tran%cmat(i,j)=0.0d0
               if(i.eq.j) tran%cmat(i,j)=1.0d0
            end do
         end do
         
         do k=1,7
            kmat=imat(k,ifr,ito)
            if(kmat.eq.0) goto 100
            do i=1,3
               do j=1,3
                  tran%tmpmat(i,j,1)=tran%pmat(i,1,kmat)*tran%cmat(1,j) 
     *            +tran%pmat(i,2,kmat)*tran%cmat(2,j)                   
     *            +tran%pmat(i,3,kmat)*tran%cmat(3,j)
               end do
            end do
            do i=1,3
               do j=1,3
                  tran%cmat(i,j)=tran%tmpmat(i,j,1)
               end do
            end do
         end do
      endif
!
 100  continue
      xx2=tran%cmat(1,1)*xx1+tran%cmat(1,2)*yy1+tran%cmat(1,3)*zz1
      yy2=tran%cmat(2,1)*xx1+tran%cmat(2,2)*yy1+tran%cmat(2,3)*zz1
      zz2=tran%cmat(3,1)*xx1+tran%cmat(3,2)*yy1+tran%cmat(3,3)*zz1
      x2=sngl(xx2)
      y2=sngl(yy2)
      z2=sngl(zz2)
!
      return
      end subroutine cotr

      
!----------------------------------------------------
      subroutine cotr_set_1(iy,mo,id,ih,mi,se,tran)
!----------------------------------------------------
!doc-  
!doc-FUNCTION:
!doc-  computes coordinate transformation matrices
!doc-  for cotr() for a given time
!doc-
!doc-PARAMETERS:
!doc-  iy (in,integer):  year (nn or nnnn notation) (1900<iy<2000)
!doc-  mo (in,integer):  month             (1<=mo<=12)
!doc-  id (in,integer):  day of month      (1<=id<=31)
!doc-  ih (in,integer):  hour of day       (0<=ih<24)
!doc-  mi (in,integer):  minute of hour    (0<=mi<60)
!doc-  se (in,real):     seconds of minute (0.0<=se<60.0)
!doc-  tran(out,transform): Transform to be passed to cotr
!doc-
!doc-NOTES:
!doc-  1.  dipole coordinates are computed from IGRF g10,g11,h11
!doc-      coefficients between 1945 and 1990 with linear time
!doc-      interpolation.  Before 1945 coefficients for 1945 are
!doc-      used and after 1990 coefficients for 1990 are used.
!doc-  2.  internal arithmetic is handled in double precision
!doc-
!doc-  3.  for iy=1967 mo=1 id=1 ih=0 mi=0 se=any all transformations
!doc-      are set to identity
!doc-
!doc-
!doc-REFERENCES:
!doc-
!doc-        Almanac for Computers 1991
!doc-        Nautical Almanac Office
!doc-        US Naval Observatory
!doc-        Washington, DC
!doc-
!doc-        M A Hapgood
!doc-        Space physics coordinate transformations: A user guide
!doc-        Planet. Space Sci., 40, 711, 1992
!doc-
      implicit none

      type(transform) :: tran
      
      integer iy,mo,id,ih,mi
      real*4 se
      
      integer ier
      integer i,j,k
      integer i1,i2,i3,i4
      integer jy, mjd
      real*8 zero, one, two, twopi
      real*8 eps
      real*8 d1
      real*8 psi,pih
      real*8 t, tmp, tmpco, tmpsi
      real*8 gg11, gg10, hh11
      real*8 gmsth, gmstd, gang
      real*8 uthour
      real*8 w1
      real*8 xmjd, xjd0, xjd
      real*8 xxp, xmu, xxl, xm, xy, xl, xls
      real*8 qx,qy,qz,qxg,qyg,qzg
      real*8 tmpx,tmpy,tmpz,tmp1,tmp2
      
      integer mdays(12,2)
      data mdays / 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 
     *334, 0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335 /
      real*4 tig,g10,g11,h11  !These are arrays.  
      integer nyrsIGRF
      parameter (nyrsIGRF=17)
      dimension tig(nyrsIGRF)  !year
      dimension g10(nyrsIGRF),g11(nyrsIGRF),h11(nyrsIGRF) !IGRF coefficients
      data tig( 1) /1000.00/
      data g10( 1),g11( 1),h11( 1)/ -30594.00, -2285.00, 5810.00/
      data tig( 2) /1945.00/
      data g10( 2),g11( 2),h11( 2)/ -30594.00, -2285.00, 5810.00/
      data tig( 3) /1950.00/
      data g10( 3),g11( 3),h11( 3)/ -30554.00, -2250.00, 5815.00/
      data tig( 4) /1955.00/
      data g10( 4),g11( 4),h11( 4)/ -30500.00, -2215.00, 5820.00/
      data tig( 5) /1960.00/
      data g10( 5),g11( 5),h11( 5)/ -30421.00, -2169.00, 5791.00/
      data tig( 6) /1965.00/
      data g10( 6),g11( 6),h11( 6)/ -30334.00, -2119.00, 5776.00/
      data tig( 7) /1970.00/
      data g10( 7),g11( 7),h11( 7)/ -30220.00, -2068.00, 5737.00/
      data tig( 8) /1975.00/
      data g10( 8),g11( 8),h11( 8)/ -30100.00, -2013.00, 5675.00/
      data tig( 9) /1980.00/
      data g10( 9),g11( 9),h11( 9)/ -29992.00, -1956.00, 5604.00/
      data tig(10) /1985.00/
      data g10(10),g11(10),h11(10)/ -29873.00, -1905.00, 5500.00/
      data tig(11) /1990.00/
      data g10(11),g11(11),h11(11)/ -29775.00, -1848.00, 5406.00/
      data tig(12) /1995.00/
      data g10(12),g11(12),h11(12)/ -29692.00, -1784.20, 5306.00/
      data tig(13) /2000.00/
      data g10(13),g11(13),h11(13)/ -29619.40, -1728.20, 5186.10/
      data tig(14) /2005.00/
      data g10(14),g11(14),h11(14)/ -29554.63, -1669.05, 5077.99/
      data tig(15) /2010.00/
      data g10(15),g11(15),h11(15)/ -29496.50, -1585.90, 4945.10/
! extrapolate: 2010 - 2015
      data tig(16) /2015.00/
      data g10(16),g11(16),h11(16)/ -29438.37, -1502.73, 4812.21/
! flat after 2015
      data tig(17) /9999.00/
      data g10(17),g11(17),h11(17)/ -29438.37, -1502.73, 4812.21/
      parameter(pih=1.570796372d0,twopi=6.283185307179586476925287d0)
      parameter(zero=0.0d0,one=1.0d0,two=2.0d0)
!
!..... check input consistency
!

      tran%year = iy
      tran%month = mo
      tran%day = id
      tran%hour = ih
      tran%minute = mi
      tran%seconds = se

      jy=iy
      if(jy.lt.200)jy=jy+1900
      ier=0
      if(jy.lt.1901) ier=1
      if(jy.gt.2016) ier=7
      if(mo.lt.1.or.mo.gt.12) ier=2
      if(id.lt.1.or.id.gt.31) ier=3
      if(ih.lt.0.or.ih.gt.23) ier=4
      if(mi.lt.0.or.mi.gt.63) ier=5
      if(se.lt.0.0.or.se.gt.60.0) ier=6
      if(ier.ne.0) then
         write(0,*)'error cotr_set, ier= ',ier
         write(0,*)jy,mo,id,ih,mi,se
         stop
      endif
!
      tran%iscalled=1
!
!.....  universal time of day in hours (UT1)
!
      uthour=dble(float(ih))+(dble(float(mi))/60.0d0)+dble(se)/3600.0d0
!
!..... full julian day
!
      i1=367*jy
      i2=(mo+9)/12
      i2=(i2+jy)*7
      i2=i2/4
      i3=(275*mo)/9
      i4=100*jy+mo
      d1=one
      if((dble(i4)-190002.5d0).lt.zero)d1=(-d1)
      xjd=dble(i1)-dble(i2)+dble(i3)+dble(id)+
     *     1721013.5d0+uthour/24.0d0-0.5d0*d1+0.5
!
!..... julian day at 0 UT
!
      xjd0=xjd-(uthour/24.0d0)
!
!.....  modified julian day
!
      xmjd=xjd-2400000.5d0
      mjd=int(xmjd)
!
!.....  set to identity transformation if iy=1967,....
!       except for the MHD transformation which still flips x- y- axes
!
      if(jy.eq.1967.and.mo.eq.1.and.id.eq.1.and.ih.eq.0.and.mi.eq.0)then
         do k=-7,7
            do j=1,3
               do i=1,3
                  tran%pmat(i,j,k)=0.0d0
                  if(i.eq.j)tran%pmat(i,j,k)=1.0d0
               end do
            end do
         end do
         do i=1,3
            do j=1,3
               tran%pmat(j,i,6)=zero
            end do
         end do
         tran%pmat(1,1,6)=(-1.0d0)
         tran%pmat(2,2,6)=(-1.0d0)
         tran%pmat(3,3,6)=(1.0d0)
         do i=1,3
            do j=1,3
               tran%pmat(j,i,-6)=tran%pmat(i,j,6)
            end do
         end do
         
         return
      endif
!
!
!..... solar coordinates
!
      t=(xjd-2451545.0d0)/36525.0d0
      tran%t0=(xjd0-2451545.0d0)/36525.0d0
      xm=357.528d0+35999.050d0*t
      xl=280.460d0+36000.772d0*t
      xm=dmod(xm+360.0d3,360.0d0)
      xl=dmod(xl+360.0d3,360.0d0)
      tmp1=rad*xm
      tmp2=2.0d0*rad*xm
      xls=xl+((1.915d0-0.0048*tran%t0)*dsin(tmp1))+(0.020d0*dsin(tmp2))
      eps=23.439d0-0.013d0*t
      eps=dmod(eps+360.0d3,360.0d0)
!
!..... greenwich mean siderial time
!
      gmsth=6.69737456d0 +2400.051336d0*tran%t0 +0.0000258622d0*
     *     tran%t0*tran%t0 +1.002737909d0*uthour
      gmsth=dmod(gmsth+24.0d3,24.0d0)
      gmstd=15.0d0*gmsth
      gang=rad*gmstd
!
      eps=rad*eps
      xls=rad*xls
!
!.... gei to geo
!
      do i=1,3
         do j=1,3
            tran%pmat(j,i,1)=zero
         end do
      end do
      tmpsi=dsin(gang)
      tmpco=dcos(gang)
      tran%pmat(1,1,1)=tmpco
      tran%pmat(2,2,1)=tmpco
      tran%pmat(3,3,1)=one
      tran%pmat(2,1,1)=-tmpsi
      tran%pmat(1,2,1)=tmpsi
      do i=1,3
         do j=1,3
            tran%pmat(j,i,-1)=tran%pmat(i,j,1)
         end do
      end do
!
!.... gei to gse
!
      do i=1,3
         do j=1,3
            tran%tmpmat(j,i,1)=zero
         end do
      end do
      tmpsi=dsin(xls)
      tmpco=dcos(xls)
      tran%tmpmat(1,1,1)=tmpco
      tran%tmpmat(2,2,1)=tmpco
      tran%tmpmat(3,3,1)=one
      tran%tmpmat(2,1,1)=-tmpsi
      tran%tmpmat(1,2,1)=tmpsi
!
      do i=1,3
         do j=1,3
            tran%tmpmat(j,i,2)=zero
         end do
      end do
      tmpsi=dsin(eps)
      tmpco=dcos(eps)
      tran%tmpmat(1,1,2)=one
      tran%tmpmat(2,2,2)=tmpco
      tran%tmpmat(3,3,2)=tmpco
      tran%tmpmat(3,2,2)=-tmpsi
      tran%tmpmat(2,3,2)=tmpsi
!
      do i=1,3
         do j=1,3
            tran%pmat(i,j,2)=tran%tmpmat(i,1,1)*tran%tmpmat(1,j,2) 
     *           +tran%tmpmat(i,2,1)*tran%tmpmat(2,j,2)  
     *           +tran%tmpmat(i,3,1)*tran%tmpmat(3,j,2)
         end do
      end do
      
      do i=1,3
         do j=1,3
            tran%pmat(j,i,-2)=tran%pmat(i,j,2)
         end do
      end do
!
!
! .... gse to gsm
!
!     CHECKME-CHECKME-CHECKME
!     These were uninitialized
      gg10 = 0.d0
      gg11 = 0.d0
      hh11 = 0.d0
!
      xy=100.0d0*(20.0d0+t) 
      do k=1,nyrsIGRF-1
         if(xy.ge.tig(k).and.xy.le.tig(k+1)) then
            tmp=(xy-tig(k))/(tig(k+1)-tig(k))
            gg10=(one-tmp)*g10(k)+tmp*g10(k+1)
            gg11=(one-tmp)*g11(k)+tmp*g11(k+1)
            hh11=(one-tmp)*h11(k)+tmp*h11(k+1)
            exit
         endif
      end do
      
      tmp1=hh11/gg11
      xxl=datan(tmp1)
      tmp1=(gg11*dcos(xxl)+hh11*dsin(xxl))/gg10
      xxp=pih-dasin( tmp1 )
!
      qxg=dcos(xxp)*dcos(xxl)
      qyg=dcos(xxp)*dsin(xxl)
      qzg=dsin(xxp)
      tmpx=tran%pmat(1,1,-1)*qxg+tran%pmat(1,2,-1)*qyg + 
     *     tran%pmat(1,3,-1)*qzg
      tmpy=tran%pmat(2,1,-1)*qxg+tran%pmat(2,2,-1)*qyg + 
     *     tran%pmat(2,3,-1)*qzg
      tmpz=tran%pmat(3,1,-1)*qxg+tran%pmat(3,2,-1)*qyg + 
     *     tran%pmat(3,3,-1)*qzg
      qx=tran%pmat(1,1,2)*tmpx+tran%pmat(1,2,2)*tmpy + 
     *     tran%pmat(1,3,2)*tmpz
      qy=tran%pmat(2,1,2)*tmpx+tran%pmat(2,2,2)*tmpy + 
     *     tran%pmat(2,3,2)*tmpz
      qz=tran%pmat(3,1,2)*tmpx+tran%pmat(3,2,2)*tmpy + 
     *     tran%pmat(3,3,2)*tmpz
      psi=datan2(qy,qz)
      do i=1,3
         do j=1,3
            tran%pmat(j,i,-3)=zero
         end do
      end do
      tmpsi=dsin(psi)
      tmpco=dcos(psi)
      tran%pmat(1,1,-3)=one
      tran%pmat(2,2,-3)=tmpco
      tran%pmat(3,3,-3)=tmpco
      tran%pmat(3,2,-3)=-tmpsi
      tran%pmat(2,3,-3)=tmpsi
      do i=1,3
         do j=1,3
            tran%pmat(j,i,3)=tran%pmat(i,j,-3)
         end do
      end do
!
!.... gsm to sm
!
      xmu=datan(qx/sqrt(qy*qy+qz*qz))
      do i=1,3
         do j=1,3
            tran%pmat(j,i,-4)=zero
         end do
      end do
      tmpsi=dsin(xmu)
      tmpco=dcos(xmu)
      tran%pmat(1,1,-4)=tmpco
      tran%pmat(2,2,-4)=one
      tran%pmat(3,3,-4)=tmpco
      tran%pmat(3,1,-4)=-tmpsi
      tran%pmat(1,3,-4)=tmpsi
      do i=1,3
         do j=1,3
            tran%pmat(j,i,4)=tran%pmat(i,j,-4)
         end do
      end do
!
!.... geo to mag
!
      w1=xxp-pih
      do i=1,3
         do j=1,3
            tran%tmpmat(j,i,1)=zero
         end do
      end do
      do i=1,3
         do j=1,3
            tran%tmpmat(j,i,2)=zero
         end do
      end do
      tmpsi=dsin(w1)
      tmpco=dcos(w1)
      tran%tmpmat(1,1,1)=tmpco
      tran%tmpmat(2,2,1)=one
      tran%tmpmat(3,3,1)=tmpco
      tran%tmpmat(3,1,1)=-tmpsi
      tran%tmpmat(1,3,1)=tmpsi
      tmpsi=dsin(xxl)
      tmpco=dcos(xxl)
      tran%tmpmat(1,1,2)=tmpco
      tran%tmpmat(2,2,2)=tmpco
      tran%tmpmat(3,3,2)=one
      tran%tmpmat(2,1,2)=-tmpsi
      tran%tmpmat(1,2,2)=tmpsi
      do i=1,3
         do j=1,3
            tran%pmat(i,j,5)=tran%tmpmat(i,1,1)*tran%tmpmat(1,j,2) 
     *           +tran%tmpmat(i,2,1)*tran%tmpmat(2,j,2) 
     *           +tran%tmpmat(i,3,1)*tran%tmpmat(3,j,2)
         end do
      end do
      do i=1,3
         do j=1,3
            tran%pmat(j,i,-5)=tran%pmat(i,j,5)
         end do
      end do
!
!.....  gse to mhd
!
      do i=1,3
         do j=1,3
            tran%pmat(j,i,6)=zero
         end do
      end do
      tran%pmat(1,1,6)=(-1.0d0)
      tran%pmat(2,2,6)=(-1.0d0)
      tran%pmat(3,3,6)=(1.0d0)
      do i=1,3
         do j=1,3
            tran%pmat(j,i,-6)=tran%pmat(i,j,6)
         end do
      end do
!      
      return
      end subroutine cotr_set_1

!     Convenience functions:
!----------------------------------------------------
      subroutine cotr_set_r8(rtime,tran)
!----------------------------------------------------
      use date_lib
      implicit none
      type(transform),intent(out) :: tran
      real*8,intent(in) :: rtime

      integer iy,mo,id,ih,mi
      real sec
      call epoch1966(rtime,iy,mo,id,ih,mi,sec,1)
      if(iy.lt.200)then
         iy = iy+1900
      endif
      call cotr_set_1(iy,mo,id,ih,mi,sec,tran)
      end subroutine cotr_set_r8

c---------------------------------------------------------------
      subroutine mhdmlt(iy,mo,id,ih,mi,se,pmhd,tmhd,xmlt,colat)
c---------------------------------------------------------------
cdoc-  
cdoc-FUNCTION:
cdoc-  convert MHD logitude (-180,180) and colatitude (0,180)
cdoc-  to magnetic local time (0-24) and colatitude
cdoc-
cdoc-PARAMETERS:
cdoc-  iy,mo,id,ih,mi,se (in):  time in UT
cdoc-  pmhd (in,real):  deg longitude in MHD system (-180.0,180.0)
cdoc-  tmhd (in,real):  deg colatitude in MHD system (0.0,180.0)
cdoc-  xmlt (out,real): hours magnetic local time (0.0-24.0)
cdoc-  colat (out,real): deg magnetic colatitude (0.0,180.0)
cdoc-
      implicit none

      type(transform) :: tran
      integer iy,mo,id,ih,mi
      real se,pmhd,tmhd,xmlt,colat

      real x1,y1,z1,x2,y2,z2
      real rr,pp,tt

      call cotr_set(iy,mo,id,ih,mi,se,tran)
      call degxyz(1.0,pmhd,tmhd,x1,y1,z1)
      call cotr(tran,'mhd','sm ',x1,y1,z1,x2,y2,z2)
      call xyzdeg(x2,y2,z2,rr,pp,tt)
      colat=tt
      xmlt=(pp+180.0)/15.0
      return
      end subroutine mhdmlt

c---------------------------------------------------------
      subroutine diptil_1(iy,mo,id,ih,mi,se,degsun,degy)
c---------------------------------------------------------
cdoc-  
cdoc-FUNCTION:
cdoc-  compute dipole tilt angles for a given time
cdoc-
cdoc-PARAMETERS:
cdoc-  iy,mo,id,ih,mi,se (in):  time in UT
cdoc-  degsun (out,real): deg tilt in GSE x-z plane
cdoc-  degy (out,real):   deg tilt in GSE y-z plane
cdoc-

      implicit none

      type(transform) :: tran
      integer iy,mo,id,ih,mi
      real se,degsun,degy

      call cotr_set(iy,mo,id,ih,mi,se,tran)
      call diptil_tran(tran,degsun,degy)
      return
      end subroutine diptil_1

      subroutine diptil_tran(tran,degsun,degy)
      implicit none
      real,intent(out) :: degsun,degy
      type(transform),intent(in) :: tran
      real :: xsm,ysm,zsm,xgse,ygse,zgse

      xsm = 0.0; ysm = 0.0; zsm = 1.0
      call cotr(tran,'sm ','gse',xsm,ysm,zsm,xgse,ygse,zgse)
      degsun=deg*atan2(xgse,zgse)
      degy=deg*atan2(ygse,zgse)
      end subroutine diptil_tran
      
      end module new_cotr
