!------------------------------
! MODULE:       coamps_intrinsic_mod
! AUTHOR:       T. R. Whitcomb
!               Naval Research Laboratory
! MODIFIED:     07 July 2008
! DART VERSION: Jamaica
!
! Collect the routines taken directly from the COAMPS utility
! package that have been translated to Fortran 90 and patched
! to use real(kind=r8) instead of whatever was there before. 
! This module is automatically generated.
!------------------------------
module coamps_intrinsic_mod
  use types_mod, only : r8

  implicit none

  !------------------------------
  ! BEGIN PUBLIC INTERFACE
  !------------------------------

  ! All routines are accessible
  public
  
  !------------------------------
  ! END PUBLIC INTERFACE
  !------------------------------

  !------------------------------
  ! BEGIN EXTERNAL INTERFACES
  !------------------------------
  ! [none]
  !------------------------------
  ! END EXTERNAL INTERFACES
  !------------------------------
  !------------------------------
  ! BEGIN TYPES AND CONSTANTS 
  !------------------------------
      real(kind=r8), parameter :: cp   = 1004.64
      real(kind=r8), parameter :: g    = 9.80616
      real(kind=r8), parameter :: rgas = 287.04
      real(kind=r8), parameter :: cv   = cp-rgas
      real(kind=r8), parameter :: rocp = rgas/cp
      real(kind=r8), parameter :: cpor = cp/rgas
      real(kind=r8), parameter :: p00  = 1.e5
      real(kind=r8), parameter :: pi   = 3.141592741012573
      real(kind=r8), parameter :: c27=2.5e6/1004.
      real(kind=r8), parameter :: eps=0.622
      real(kind=r8), parameter :: e0=6.11
      real(kind=r8), parameter :: xl=2.5e6
      real(kind=r8), parameter :: gamma=6.5e-3
      real(kind=r8), parameter :: heatlv=597.3*4186.0
      real(kind=r8), parameter :: heatls=677.0*4186.0
      real(kind=r8), parameter :: hlvocp=heatlv/cp
      real(kind=r8), parameter :: hlsocp=heatls/cp
  !------------------------------
  ! END TYPES AND CONSTANTS 
  !------------------------------

  !------------------------------
  ! BEGIN MODULE VARIABLES
  !------------------------------
  ! [none]
  !------------------------------
  ! END MODULE VARIABLES
  !------------------------------

contains
subroutine ij2ll(igrid,reflat,reflon,iref,jref,stdlt1,stdlt2&
                      ,stdlon,delx,dely,grdi,grdj,npts&
                      ,grdlat,grdlon)

! rcs keywords: $RCSfile: ij2ll.f,v $ 
!               $Revision: 1.1 $ $Date: 1997/12/05 20:19:58 $c
!***********************************************************************
!
      implicit none
!
!***********************************************************************
!           parameters:
!***********************************************************************
!
      integer igrid
      integer iref
      integer jref
      integer npts
!
      real(kind=r8) delx
      real(kind=r8) dely
      real(kind=r8) grdi   (npts)
      real(kind=r8) grdj   (npts)
      real(kind=r8) grdlat (npts)
      real(kind=r8) grdlon (npts)
      real(kind=r8) reflat
      real(kind=r8) reflon
      real(kind=r8) stdlon
      real(kind=r8) stdlt1
      real(kind=r8) stdlt2
!
!***********************************************************************
!          local variables and dynamic storage:
!***********************************************************************
!
      integer i
      integer ihem
!
      real(kind=r8) angle
      real(kind=r8) cn1
      real(kind=r8) cn2
      real(kind=r8) cn3
      real(kind=r8) cn4
      real(kind=r8) con1
      real(kind=r8) con2
      real(kind=r8) d2r
      real(kind=r8) deg
      real(kind=r8) gcon
      real(kind=r8) ogcon
      real(kind=r8) omega4
      real(kind=r8) onedeg
      real(kind=r8) pi
      real(kind=r8) pi2
      real(kind=r8) pi4
      real(kind=r8) r2d
      real(kind=r8) radius
      real(kind=r8) rih
      real(kind=r8) rr
      real(kind=r8) x
      real(kind=r8) xih
      real(kind=r8) xx
      real(kind=r8) y
      real(kind=r8) yih
      real(kind=r8) yy
!
!***********************************************************************
!
!          subroutine: ij2ll
!
!          purpose: to compute latitude and longitude of specified
!                   i- and j-points on a grid.  all latitudes in
!                   this routine start with -90.0 at the south pole
!                   and increase northward to +90.0 at the north
!                   pole.  the longitudes start with 0.0 at the
!                   greenwich meridian and increase to the east,
!                   so that 90.0 refers to 90.0e, 180.0 is the
!                   international dateline and 270.0 is 90.0w.
!
!          input variables:
!
!            igrid:  type of grid projection:
!
!                       = 1, mercator projection
!                       = 2, lambert conformal projection
!                       = 3, polar stereographic projection
!                       = 4, cartesian coordinates
!                       = 5, spherical projection
!            reflat: latitude at reference point (iref,jref)
!
!            reflon: longitude at reference point (iref,jref)
!
!            iref:   i-coordinate value of reference point
!
!            jref:   j-coordinate value of reference point
!
!            stdlt1: standard latitude of grid
!
!            stdlt2: second standard latitude of grid (only required
!                    if igrid = 2, lambert conformal)
!
!            stdlon: standard longitude of grid (longitude that
!                     points to the north)
!
!            delx:   grid spacing of grid in x-direction
!                    for igrid = 1,2,3 or 4, delx must be in meters
!                    for igrid = 5, delx must be in degrees
!
!            dely:   grid spacing (in meters) of grid in y-direction
!                    for igrid = 1,2,3 or 4, delx must be in meters
!                    for igrid = 5, dely must be in degrees
!
!            grdi:   i-coordinate(s) that this routine will generate
!                    information for
!
!            grdj:   j-coordinate(s) that this routine will generate
!                    information for
!
!            npts:   number of points to find information for
!
!          output variables:
!
!            grdlat: latitude of point (grdi,grdj)
!
!            grdlon: longitude of point (grdi,grdj)
!
!********************************************************************
!
!********************************************************************
!          make sure igrid is an acceptable value
!********************************************************************
!
      if (igrid.lt.1.or.igrid.gt.5) then
        print 800
        print 805,igrid
        print 800
        return
      endif
!
!********************************************************************
!          local constants
!********************************************************************
!
      pi=4.0*atan(1.0)
      pi2=pi/2.0
      pi4=pi/4.0
      d2r=pi/180.0
      r2d=180.0/pi
      radius=6371229.0
      omega4=4.0*pi/86400.0
!
!********************************************************************
!          mercator projection (igrid=1)
!********************************************************************
!
      if (igrid.eq.1) then
        gcon=0.0
        deg=abs(stdlt1)*d2r
        con1=cos(deg)
        con2=radius*con1
        deg=reflat*0.5*d2r
        rih=con2*dlog(tan(pi4+deg))
        do i=1,npts
          rr=rih +(grdj(i)-jref)*dely
          grdlat(i)=(2.0*atan(exp(rr/con2))-pi2)*r2d
          grdlon(i)=reflon +(grdi(i)-iref)*r2d*delx/con2
          if (grdlon(i).gt.360.0) grdlon(i)=grdlon(i)-360.0
          if (grdlon(i).lt.  0.0) grdlon(i)=grdlon(i)+360.0
        enddo
        return
!
!********************************************************************
!          lambert conformal (igrid=2) or
!          polar stereographic (igrid=3)
!********************************************************************
!
      else if (igrid.eq.2.or.igrid.eq.3) then
        if (igrid.eq.2) then
          if (stdlt1.eq.stdlt2) then
            gcon=sin(abs(stdlt1)*d2r)
          else
            gcon=(log(sin((90.0-abs(stdlt1))*d2r))&
                 -log(sin((90.0-abs(stdlt2))*d2r)))&
                /(log(tan((90.0-abs(stdlt1))*0.5*d2r))&
                 -log(tan((90.0-abs(stdlt2))*0.5*d2r)))
          endif
        else
          gcon=1.0
        endif
        ogcon=1.0/gcon
        ihem=nint(abs(stdlt1)/stdlt1)
        deg=(90.0-abs(stdlt1))*d2r
        cn1=sin(deg)
        cn2=radius*cn1*ogcon
        deg=deg*0.5
        cn3=tan(deg)
        deg=(90.0-abs(reflat))*0.5*d2r
        cn4=tan(deg)
        rih=cn2*(cn4/cn3)**gcon
        deg=(reflon-stdlon)*d2r*gcon
        xih= rih*sin(deg)
        yih=-rih*cos(deg)*ihem
        do i=1,npts
          x=xih+(grdi(i)-iref)*delx
          y=yih+(grdj(i)-jref)*dely
          rr=sqrt(x*x+y*y)
          grdlat(i)=r2d*(pi2-2.0*atan(cn3*(rr/cn2)**ogcon))*ihem
          xx= x
          yy=-y*ihem
          if (yy.eq.0.0) then
            if (xx.le.0.0) then
              angle=-90.0
            else if (xx.gt.0.0) then
              angle=90.0
            endif
          else
            angle=atan2(xx,yy)*r2d
          endif
          grdlon(i)=stdlon+angle*ogcon
          if (grdlon(i).gt.360.0) grdlon(i)=grdlon(i)-360.0
          if (grdlon(i).lt.  0.0) grdlon(i)=grdlon(i)+360.0
        enddo
        return
!
!********************************************************************
!          analytic grid (igrid=4)
!********************************************************************
!
      else if (igrid.eq.4) then
        onedeg=radius*2.0*pi/360.0
        cn2=delx/onedeg
        do i=1,npts
          grdlat(i)=reflat+(grdj(i)-jref)*cn2
          grdlon(i)=reflon+(grdi(i)-iref)*cn2
          if (grdlon(i).gt.360.0) grdlon(i)=grdlon(i)-360.0
          if (grdlon(i).lt.  0.0) grdlon(i)=grdlon(i)+360.0
        enddo
        return
!
!********************************************************************
!          spherical grid (igrid=5)
!********************************************************************
!
      else if (igrid.eq.5) then
        do i=1,npts
          grdlat(i)= (grdj(i)-float(jref))*dely+reflat
          grdlon(i)= (grdi(i)-float(iref))*delx+reflon
          if (grdlon(i).gt.360.0) grdlon(i)=grdlon(i)-360.0
          if (grdlon(i).lt.  0.0) grdlon(i)=grdlon(i)+360.0
        enddo
        return
      endif
!
!********************************************************************
!          format statements
!********************************************************************
!
  800 format(/,' ',72('-'),/)
  805 format(/,' ERROR from subroutine ij2ll:',/&
              ,'   igrid must be one of the following values:',//&
              ,'            1: mercator projection',/&
              ,'            2: lambert conformal projection',/&
              ,'            3: polar steographic projection',/&
              ,'            4: cartesian coordinates',/&
              ,'            5: spherical projection',//&
              ,' Your entry was:',i6,', Correct and try again',/)
!
!********************************************************************
!
      return
end subroutine ij2ll

subroutine ll2ij(igrid,reflat,reflon,iref,jref,stdlt1,stdlt2&
                      ,stdlon,delx,dely,grdi,grdj,npts&
                      ,grdlat,grdlon)

! rcs keywords: $RCSfile: ll2ij.f,v $ 
!               $Revision: 1.2 $ $Date: 1998/11/24 23:38:08 $c
!********************************************************************
!********************************************************************
!          SUBROUTINE: ll2ij
!
!          PURPOSE: To compute i- and j-coordinates of a specified
!                   grid given the latitude and longitude points.
!                   All latitudes in this routine start
!                   with -90.0 at the south pole and increase
!                   northward to +90.0 at the north pole.  The
!                   longitudes start with 0.0 at the Greenwich
!                   meridian and increase to the east, so that
!                   90.0 refers to 90.0E, 180.0 is the inter-
!                   national dateline and 270.0 is 90.0W.
!
!          INPUT VARIABLES:
!
!            igrid:  type of grid projection:
!
!                       = 1, mercator projection
!                       = 2, lambert conformal projection
!                       = 3, polar stereographic projection
!                       = 4, cartesian coordinates
!                       = 5, spherical projection
!
!            reflat: latitude at reference point (iref,jref)
!
!            reflon: longitude at reference point (iref,jref)
!
!            iref:   i-coordinate value of reference point
!
!            jref:   j-coordinate value of reference point
!
!            stdlt1: standard latitude of grid
!
!            stdlt2: second standard latitude of grid (only required
!                    if igrid = 2, lambert conformal)
!
!            stdlon: standard longitude of grid (longitude that
!                     points to the north)
!
!            delx:   grid spacing of grid in x-direction
!                    for igrid = 1,2,3 or 4, and 5 delx must be in meters
!
!            dely:   grid spacing (in meters) of grid in y-direction
!                    for igrid = 1,2,3 or 4, and 5 delx must be in meters
!
!            grdlat: latitude of point (grdi,grdj)
!
!            grdlon: longitude of point (grdi,grdj)
!
!            npts:   number of points to find information for
!
!          OUTPUT VARIABLES:
!
!            grdi:   i-coordinate(s) that this routine will generate
!                    information for
!
!            grdj:   j-coordinate(s) that this routine will generate
!                    information for
!
!********************************************************************
!********************************************************************
!
      implicit none
!
!********************************************************************
!
      integer i
      integer igrid
      integer ihem
      integer iref
      integer jref
      integer npts
!
      real(kind=r8) alnfix
      real(kind=r8) alon
      real(kind=r8) check
      real(kind=r8) cn1
      real(kind=r8) cn2
      real(kind=r8) cn3
      real(kind=r8) cn4
      real(kind=r8) cnx
      real(kind=r8) cny
      real(kind=r8) con1
      real(kind=r8) con2
      real(kind=r8) d2r
      real(kind=r8) deg
      real(kind=r8) delx
      real(kind=r8) dely
      real(kind=r8) gcon
      real(kind=r8) grdi  (npts)
      real(kind=r8) grdj  (npts)
      real(kind=r8) grdlat(npts)
      real(kind=r8) grdlon(npts)
      real(kind=r8) ogcon
      real(kind=r8) omega4
      real(kind=r8) onedeg
      real(kind=r8) pi
      real(kind=r8) pi2
      real(kind=r8) pi4
      real(kind=r8) r2d
      real(kind=r8) radius
      real(kind=r8) reflat
      real(kind=r8) reflon
      real(kind=r8) rih
      real(kind=r8) rrih
      real(kind=r8) stdlon
      real(kind=r8) stdlt1
      real(kind=r8) stdlt2
      real(kind=r8) x
      real(kind=r8) xih
      real(kind=r8) y
      real(kind=r8) yih
!
!********************************************************************
!          make sure igrid is an acceptable value
!********************************************************************
!
      if (igrid.lt.1.or.igrid.gt.5) then
        print 800
        print 805,igrid
        print 800
        return
      endif
!
!********************************************************************
!          local constants
!********************************************************************
!
      pi=4.0*atan(1.0)
      pi2=pi/2.0
      pi4=pi/4.0
      d2r=pi/180.0
      r2d=180.0/pi
      radius=6371229.0
      omega4=4.0*pi/86400.0
      onedeg=radius*2.0*pi/360.0
!
!********************************************************************
!          mercator projection (igrid=1)
!********************************************************************
!
      if (igrid.eq.1) then
        deg=abs(stdlt1)*d2r
        con1=cos(deg)
        con2=radius*con1
        deg=reflat*0.5*d2r
        rih=con2*dlog(tan(pi4+deg))
        do i=1,npts
          alon=grdlon(i)+180.0-reflon
          if (alon.lt.  0.0) alon=alon+360.0
          if (alon.gt.360.0) alon=alon-360.0
          grdi(i)=iref+(alon-180.0)*con2/(r2d*delx)
          deg=grdlat(i)*d2r+pi2
          deg=deg*0.5
          grdj(i)=jref+(con2*dlog(tan(deg))-rih)/dely
        enddo
        return
!
!********************************************************************
!          lambert conformal (igrid=2) or
!          polar stereographic (igrid=3)
!********************************************************************
!
      else if (igrid.eq.2.or.igrid.eq.3) then
        if (igrid.eq.2) then
          if (stdlt1.eq.stdlt2) then
            gcon=sin(abs(stdlt1)*d2r)
          else
            gcon=(log(sin((90.0-abs(stdlt1))*d2r))&
                 -log(sin((90.0-abs(stdlt2))*d2r)))&
                /(log(tan((90.0-abs(stdlt1))*0.5*d2r))&
                 -log(tan((90.0-abs(stdlt2))*0.5*d2r)))
          endif
        else
          gcon=1.0
        endif
        ogcon=1.0/gcon
        ihem=nint(abs(stdlt1)/stdlt1)
        deg=(90.0-abs(stdlt1))*d2r
        cn1=sin(deg)
        cn2=radius*cn1*ogcon
        deg=deg*0.5
        cn3=tan(deg)
        deg=(90.0-abs(reflat))*0.5*d2r
        cn4=tan(deg)
        rih=cn2*(cn4/cn3)**gcon
        deg=(reflon-stdlon)*d2r*gcon
        xih= rih*sin(deg)
        yih=-rih*cos(deg)*ihem
        do i=1,npts
          deg=(90.0-grdlat(i)*ihem)*0.5*d2r
          cn4=tan(deg)
          rrih=cn2*(cn4/cn3)**gcon
          check=180.0-stdlon
          alnfix=stdlon+check
          alon=grdlon(i)+check
          if (alon.lt.  0.0) alon=alon+360.0
          if (alon.gt.360.0) alon=alon-360.0
          deg=(alon-alnfix)*gcon*d2r
          x= rrih*sin(deg)
          y=-rrih*cos(deg)*ihem
          grdi(i)=iref+(x-xih)/delx
          grdj(i)=jref+(y-yih)/dely
        enddo
        return
!
!********************************************************************
!          analytic grid (igrid=4)
!********************************************************************
!
      else if (igrid.eq.4) then
        cnx=delx/onedeg
        cny=dely/onedeg
        do i=1,npts
          grdi(i)=iref+(grdlon(i)-reflon)/cnx
          grdj(i)=jref+(grdlat(i)-reflat)/cny
        enddo
        return
!
!********************************************************************
!          spherical grid (igrid=5)
!********************************************************************
!
      else if (igrid.eq.5) then
        cnx=delx/onedeg
        cny=dely/onedeg
        do i=1,npts
          grdi(i)=iref+(grdlon(i)-reflon)/cnx
          grdj(i)=jref+(grdlat(i)-reflat)/cny
        enddo
        return
      endif
!
!********************************************************************
!          format statements
!********************************************************************
!
  800 format(/,' ',72('-'),/)
  805 format(/,' ERROR from subroutine ll2ij:',/&
              ,'   igrid must be one of the following values:',//&
              ,'            1: mercator projection',/&
              ,'            2: lambert conformal projection',/&
              ,'            3: polar steographic projection',/&
              ,'            4: cartesian coordinates',/&
              ,'            5: spherical projection',//&
              ,' Your entry was:',i6,', Correct and try again',/)
!
!********************************************************************
!
      return
end subroutine ll2ij

subroutine s2pint(din,dout,zin,zout,kin,kout,len,missing,value)

! rcs keywords: $RCSfile: s2pint.f,v $ 
!               $Revision: 1.1 $ $Date: 1997/12/05 20:20:02 $
!
      implicit none
!
!***********************************************************************
!           parameters:
!***********************************************************************
!
      integer kin
      integer kout
      integer len
      integer missing
!
      real(kind=r8) din (len,kin)
      real(kind=r8) dout (len,kout)
      real(kind=r8) value
      real(kind=r8) zin (len,kin)
      real(kind=r8) zout (kout)
!
!***********************************************************************
!          local variables and dynamic storage:
!***********************************************************************
!
      integer i
      integer ki
      integer ko
!
      real(kind=r8) temp
!
!***********************************************************************
!          end of definitions
!***********************************************************************
!***********************************************************************
!          interpolate from sigma levels to p-levels
!          set up interpolation loops
!***********************************************************************
!
      do ko=1,kout
        do i=1,len
          if (zout(ko).le.zin(i,1)) dout(i,ko)=din(i,1)
!          if (zout(ko).ge.zin(i,kin)) dout(i,ko)=din(i,kin)
!***********************************************************************
!         assign under ground (missing) values
!***********************************************************************
          if (zout(ko).ge.zin(i,kin)) then
            if(missing.eq.1) then
              dout(i,ko)=din(i,kin)
            else
              dout(i,ko)=value
            endif
          endif

        enddo
      enddo
      do ki=2,kin
        do ko=1,kout
          do i=1,len
            if ((zout(ko).gt.zin(i,ki-1)).and.(zout(ko).le.zin(i,ki)))&
             then
              temp=(dlog(zout(ko))-dlog(zin(i,ki-1)))&
                  /(dlog(zin(i,ki))-dlog(zin(i,ki-1)))
              dout(i,ko)=din(i,ki-1)+temp*(din(i,ki)-din(i,ki-1))
            endif
          enddo
        enddo
      enddo
!
!***********************************************************************
!
      return
end subroutine s2pint

subroutine sfcp(th,p,thbm,exbw,dsigw,sigmwa,zsfc,m,n,kk,ps)

      implicit none
      
      integer i
      integer m
      integer n
      integer mn
      integer kk
      integer kkp1

      real(kind=r8) th(m,n,kk)
      real(kind=r8) p(m,n,kk)
      real(kind=r8) thbm(m,n,kk)
      real(kind=r8) dsigw(kk+1)
      real(kind=r8) aoz
      real(kind=r8) exbw(m,n,kk+1)
      real(kind=r8) ppsfc
      real(kind=r8) ps(m,n)
      real(kind=r8) sigmwa(kk+1)
      real(kind=r8) zsfc(m,n)


      kkp1=kk+1
      mn=m*n
!
!  compute surface pressure (mb)
!
      do i=1,mn
        aoz=sigmwa(1)/(sigmwa(1)-zsfc(i,1))
        ppsfc=p(i,1,kk)&
             -dsigw(kkp1)*g*(th(i,1,kk)-thbm(i,1,kk))&
            /(aoz*cp*thbm(i,1,kk)*thbm(i,1,kk))
        ps(i,1)=(exbw(i,1,kkp1)+ppsfc)**(cp/rgas)*1000.0
        if(i.eq.1) then
           print *,'aoz,ppsfc,ps= ',aoz,ppsfc,ps(1,1)
        endif
      enddo

      print *,'ps(10,10)= ',ps(10,10)
      return
end subroutine sfcp

subroutine z2zint(din,dout,ifin,zin,zout,zsfc,kin,kout,len&
        ,missing,value)

! rcs keywords: $RCSfile: z2zint.f,v $ 
!               $Revision: 1.1 $ $Date$
!
      implicit none

!
!***********************************************************************
!           parameters:
!***********************************************************************
!
      integer kin
      integer kout
      integer len
      integer ifin
!
      real(kind=r8) din    (len,kin)
      real(kind=r8) dout   (len,kout)
      real(kind=r8) zin    (len,kin)
!      real(kind=r8) zout   (kout)
      real(kind=r8) zout   (maxlevels)
      real(kind=r8) zsfc   (len)
      real(kind=r8) zagl   (len,kout)
!
!***********************************************************************
!          local variables and dynamic storage:
!***********************************************************************
!
      integer i
      integer ki
      integer ko
      integer missing
!
      real(kind=r8) temp
      real(kind=r8) value
!
!***********************************************************************
!         add in terrain height for agl calculation
!***********************************************************************
!
      if (ifin.eq.6) then
        do i=1,len
          do ko=1,kout
            zagl(i,ko)=zout(ko)+zsfc(i)
          enddo
        enddo
      else
        do i=1,len
          do ko=1,kout
            zagl(i,ko)=zout(ko)
          enddo
        enddo
      endif
!
!***********************************************************************
!          interpolate from sigma-z levels to constant z levels
!***********************************************************************
!
      do i=1,len
        do ko=1,kout
          if (zagl(i,ko).ge.zin(i,1)) dout(i,ko)=din(i,1)
!          if (zout(ko).le.zin(i,kin)) dout(i,ko)=din(i,kin)
!***********************************************************************
!         assign under ground (missing) values
!***********************************************************************
          if (zagl(i,ko).le.zin(i,kin)) then
            if(missing.eq.1) then
              dout(i,ko)=din(i,kin)
            else
              dout(i,ko)=value
            endif
          endif
!
        enddo
      enddo
!
      do ki=2,kin
        do ko=1,kout
          do i=1,len
            if ((zagl(i,ko).lt.zin(i,ki-1)).and.&
             (zagl(i,ko).ge.zin(i,ki))) then
              temp=(zagl(i,ko)-zin(i,ki-1))/(zin(i,ki)-zin(i,ki-1))
              dout(i,ko)=din(i,ki-1)+temp*(din(i,ki)-din(i,ki-1))
            endif
          enddo
        enddo
      enddo
!
!***********************************************************************
!
      return
end subroutine z2zint

subroutine utom(u,m,n,kk)

!
! rcs keywords: $RCSfile: utom.f,v $ 
!               $Revision: 1.1 $ $Date: 1997/12/05 20:20:05 $
!

!
!***********************************************************************
!           subroutine to destagger u-component
!***********************************************************************
!
      integer m
      integer n
      integer kk
      integer m1
      integer i
      integer j
      integer k

      real(kind=r8) u(m,n,kk)
      real(kind=r8) work(m,n,kk)
!
      m1=m-1
!
      do k=1,kk
        do j=1,n
          do i=2,m1
            work(i,j,k)=(u(i-1,j,k)+u(i,j,k))*0.5
          enddo
          work(1,j,k)=u( 1,j,k)
          work(m,j,k)=u(m1,j,k)
        enddo
        do i=1,m
          do j=1,n
            u(i,j,k)=work(i,j,k)
          enddo
        enddo
      enddo 
!
!***********************************************************************
!
      return
end subroutine utom

subroutine vtom(v,m,n,kk)

!
! rcs keywords: $RCSfile: vtom.f,v $ 
!               $Revision: 1.1 $ $Date: 1997/12/05 20:20:07 $
!

!
!***********************************************************************
!           subroutine to destagger v-component
!***********************************************************************
!
      implicit none
      integer i
      integer j
      integer k
      integer m
      integer n
      integer kk
      integer n1
      real(kind=r8) v(m,n,kk)      
      real(kind=r8) work(m,n,kk)
!
      n1=n-1
!
      do k=1,kk
        do j=2,n1
          do i=1,m
            work(i,j,k)=(v(i,j-1,k)+v(i,j,k))*0.5
          enddo
        enddo
        do i=1,m
          work(i,1,k)=v(i, 1,k)
          work(i,n,k)=v(i,n1,k)
        enddo
        do i=1,m
          do j=1,n
            v(i,j,k)=work(i,j,k)
          enddo
        enddo 
      enddo
!
!***********************************************************************
!
      return
end subroutine vtom

subroutine uvg2uv (u, v, kk, mn, rot,utru,vtru)
!
!****************************************************************** 
! arguments:
!
!        u:   u-component of winds in grid coordinate
!
!        v:   v-component of winds in grid coordinate 
!
!        mn:  number of points in the x and y plane
!
!        kk:  number of points in the z (or k) direction
!
!        utru: u-component of winds in world (lat/long) coordinate
!
!        vtru: v-component of winds in world (lat/long) coordinate
!
!  convert grid u and v to real(kind=r8) u and v
!  assuming grid u = real(kind=r8) u and grid v = real(kind=r8) v along
!  the standard longitude and rot is the rotation array
!
!  Author: Pedro T.H. Tsai   30 Nov 1998
!
!********************************************************************
!
        implicit none
        integer mn, i, k, kk
        real(kind=r8) u(mn,kk), v(mn,kk), rot(mn), utru(mn,kk), vtru(mn,kk)
        real(kind=r8) r2d, d2r, ff, dd, ndd, uu, vv
  
        d2r = atan(1.0)/45.0
        r2d = 1.0/d2r
        do 10 k = 1, kk
          do 10 i = 1, mn
!  calc rotation angle between grid lon and standard lon
            uu   = u(i,k)
            vv   = v(i,k)
            ff   = sqrt((uu*uu) + (vv*vv))
            if (uu .eq. 0.0) uu = 1.0E-6
            dd   = 270.0 - (atan2(vv,uu)*r2d)
            dd   = dmod(dd, 360.0)
            ndd  = dd - rot(i) + 360.0
            ndd  = dmod(ndd, 360.0)
            utru(i,k) = -sin(ndd*d2r)*ff
            vtru(i,k) = -cos(ndd*d2r)*ff
10    continue
  
      return
end subroutine uvg2uv

subroutine pstd(m,zout,pout)

! SCCS IDENTIFICATION:  %W% %G%
!
      implicit none
!
!***********************************************************************
!           parameters:
!***********************************************************************
!
      integer m
!
      real(kind=r8) pout (m)
      real(kind=r8) zout (m)
!
!***********************************************************************
!          local variables and dynamic storage:
!***********************************************************************
!
      integer i
      integer ki
!
      real(kind=r8) dz
      real(kind=r8) pin (223)
      real(kind=r8) zin (223)
!
!***********************************************************************
!          end of definitions
!***********************************************************************
!***********************************************************************
!          table lookup for standard atmosphere pressures as a function
!          of height.  table has resolution of 10mbs below 100mb, 1mb
!          from 100mbs to 1 mb.  interpolation is linear in z.
!***********************************************************************
!th
!***********************************************************************
!***********************************************************************
!          new: above 1 mb, add levels to 0.1 mb 
!          every 0.05 mb using 1976 standard atmosphere calculator
!          ref:  www.digitaldutch.com/atmoscalc
!
!          total of 223 levels
!***********************************************************************
!***********************************************************************
!th
!
      data zin/&
        636870. ,609463. ,589448.&
       ,573612. ,560442. ,549165. ,539271. ,530456.&
       ,522493. ,515237. ,508569. ,502381. ,496625.&
       ,491271. ,486250. ,481543. ,477101. ,472904.&
       ,469048. ,416281. ,386712. ,366358. ,350923.&
       ,338509. ,328134. ,319230. ,311444. ,304532.&
       ,298311. ,292647. ,287449. ,282648. ,278187.&
       ,274022. ,270117. ,266442. ,262971. ,259683.&
       ,256561. ,253588. ,250750. ,248037. ,245438.&
       ,242944. ,240546. ,238239. ,236014. ,233867.&
       ,231793. ,229786. ,227843. ,225960. ,224133.&
       ,222359. ,220634. ,218957. ,217325. ,215736.&
       ,214186. ,212676. ,211201. ,209762. ,208356.&
       ,206982. ,205638. ,204323. ,203036. ,201776.&
       ,200542. ,199332. ,198146. ,196983. ,195842.&
       ,194721. ,193620. ,192539. ,191476. ,190431.&
       ,189403. ,188391. ,187396. ,186417. ,185453.&
       ,184503. ,183568. ,182647. ,181739. ,180844.&
       ,179962. ,179092. ,178235. ,177389. ,176554.&
       ,175730. ,174917. ,174115. ,173323. ,172540.&
       ,171768. ,171005. ,170251. ,169506. ,168770.&
       ,168043. ,167324. ,166613. ,165911. ,165216.&
       ,164529. ,163849. ,163177. ,162512. ,161854.&
       ,161202. ,160558. ,159920. ,159289. ,158664.&
       ,152737. ,147326. ,142348. ,137739. ,133449.&
       ,129436. ,125665. ,122111. ,118749. ,115559.&
       ,112525. ,109632. ,106866. ,104200. ,101623.&
       , 99128. , 96710. , 94363. , 92083. , 89866.&
       , 87707. , 85605. , 83555. , 81555. , 79601.&
       , 77693. , 75827. , 74001. , 72214. , 70464.&
       , 68748. , 67066. , 65417. , 63798. , 62209.&
       , 60647. , 59114. , 57606. , 56124. , 54666.&
       , 53231. , 51819. , 50429. , 49060. , 47711.&
       , 46382. , 45072. , 43780. , 42507. , 41251.&
       , 40011. , 38788. , 37581. , 36389. , 35212.&
       , 34050. , 32902. , 31768. , 30647. , 29539.&
       , 28444. , 27362. , 26291. , 25233. , 24186.&
       , 23150. , 22125. , 21111. , 20107. , 19113.&
       , 18130. , 17156. , 16192. , 15237. , 14292.&
       , 13355. , 12427. , 11508. , 10597. ,  9695.&
       ,  8800. ,  7914. ,  7035. ,  6164. ,  5300.&
       ,  4443. ,  3594. ,  2752. ,  1917. ,  1088.&
       ,   266. ,  -549. , -1358. , -2160. , -2957.&
       , -3747. , -4531. , -5309. , -6081. , -6848.&
       , -7609. , -8364. , -9114. , -9859. ,-10598./
!
      data pin/&
           0.10,    0.15,    0.20,&
           0.25,    0.30,    0.35,    0.40,    0.45,&
           0.50,    0.55,    0.60,    0.65,    0.70,&
           0.75,    0.80,    0.85,    0.90,    0.95,&
             1.,      2.,      3.,      4.,      5.,&
             6.,      7.,      8.,      9.,     10.,&
            11.,     12.,     13.,     14.,     15.,&
            16.,     17.,     18.,     19.,     20.,&
            21.,     22.,     23.,     24.,     25.,&
            26.,     27.,     28.,     29.,     30.,&
            31.,     32.,     33.,     34.,     35.,&
            36.,     37.,     38.,     39.,     40.,&
            41.,     42.,     43.,     44.,     45.,&
            46.,     47.,     48.,     49.,     50.,&
            51.,     52.,     53.,     54.,     55.,&
            56.,     57.,     58.,     59.,     60.,&
            61.,     62.,     63.,     64.,     65.,&
            66.,     67.,     68.,     69.,     70.,&
            71.,     72.,     73.,     74.,     75.,&
            76.,     77.,     78.,     79.,     80.,&
            81.,     82.,     83.,     84.,     85.,&
            86.,     87.,     88.,     89.,     90.,&
            91.,     92.,     93.,     94.,     95.,&
            96.,     97.,     98.,     99.,    100.,&
           110.,    120.,    130.,    140.,    150.,&
           160.,    170.,    180.,    190.,    200.,&
           210.,    220.,    230.,    240.,    250.,&
           260.,    270.,    280.,    290.,    300.,&
           310.,    320.,    330.,    340.,    350.,&
           360.,    370.,    380.,    390.,    400.,&
           410.,    420.,    430.,    440.,    450.,&
           460.,    470.,    480.,    490.,    500.,&
           510.,    520.,    530.,    540.,    550.,&
           560.,    570.,    580.,    590.,    600.,&
           610.,    620.,    630.,    640.,    650.,&
           660.,    670.,    680.,    690.,    700.,&
           710.,    720.,    730.,    740.,    750.,&
           760.,    770.,    780.,    790.,    800.,&
           810.,    820.,    830.,    840.,    850.,&
           860.,    870.,    880.,    890.,    900.,&
           910.,    920.,    930.,    940.,    950.,&
           960.,    970.,    980.,    990.,   1000.,&
          1010.,   1020.,   1030.,   1040.,   1050.,&
          1060.,   1070.,   1080.,   1090.,   1100.,&
          1110.,   1120.,   1130.,   1140.,   1150./
!
!***********************************************************************
!           for every entry in the height table
!***********************************************************************
!
      do ki=2,223
!
!***********************************************************************
!           for every value of the input height array
!     search the height table for spanning entries
!***********************************************************************
!
        do i=1,m
!
!***********************************************************************
!           if spanning entries found then calculate delta height adjust
!***********************************************************************
!
          if (zout(i).lt.zin(ki-1).and.zout(i).ge.zin(ki)) then
            dz=(zin(ki-1)-zout(i))/(zin(ki-1)-zin(ki))
!
!***********************************************************************
!           and calculate output pressure
!***********************************************************************
!
            pout(i)=exp(log(pin(ki-1))&
                   +(log(pin(ki))-log(pin(ki-1)))*dz)
          endif
        enddo
      enddo
!
!***********************************************************************
!
      return
end subroutine pstd

subroutine tstd(m,pin,ttin)

! SCCS IDENTIFICATION:  %W% %G%
!
      implicit none
!
!***********************************************************************
!           parameters:
!***********************************************************************
!
      integer m
!
      real(kind=r8) pin    (m)
      real(kind=r8) ttin   (m)
!
!***********************************************************************
!          local variables and dynamic storage:
!***********************************************************************
!
      integer i
      integer ii     (m)
!
      real(kind=r8) hndrd
      real(kind=r8) one
      real(kind=r8) pw     (m,2)
      real(kind=r8) tenth
      real(kind=r8) tt     (230)
      real(kind=r8) tt1    (70)
      real(kind=r8) tt2    (40)
      real(kind=r8) tt3    (70)
      real(kind=r8) tt4    (40)
      real(kind=r8) tt5    (10)
!
!***********************************************************************
!          lookup table for standard atmosphere temperature as a
!          function of pressure.  table has resolution of 10 mbs below
!          100mb, 1mb from 100 to 1 mb. interpolation is linear in log p
!th
!          added data from 1 mb to 0.1 mb every 0.1 mb based on
!          www.digitaldutch.com/atmoscalc
!  
!          these are in tt5 (from 0.1 mb to 1.0 mb)
!          10 values
!
!  also modified existing temperatures from 10 mb to 1 mb
!  because they differed from new ref above
!  these are the old ones
!
!      data tt3/
!     1  272.0700 ,258.3600 ,249.7100 ,243.3700 ,238.6400
!     2 ,235.7700 ,233.3300 ,231.2300 ,229.3700 ,227.7100
! 
!***********************************************************************
!
      equivalence (tt(111),tt3)
      equivalence (tt(181),tt4)
      equivalence (tt,tt1)
      equivalence (tt(71),tt2)
      equivalence (tt(221),tt5)
!
!***********************************************************************
!          local constants
!***********************************************************************
!
      data tt1/&
        227.7088 ,223.1353 ,220.5026 ,218.6536 ,217.2301&
       ,216.6500 ,216.6500 ,216.6500 ,216.6500 ,216.6500&
       ,216.6500 ,216.6500 ,216.6500 ,216.6500 ,216.6500&
       ,216.6500 ,216.6500 ,216.6500 ,216.6500 ,216.6500&
       ,216.6500 ,216.6500 ,217.3175 ,219.0844 ,220.7927&
       ,222.4464 ,224.0495 ,225.6052 ,227.1165 ,228.5862&
       ,230.0167 ,231.4104 ,232.7692 ,234.0951 ,235.3898&
       ,236.6548 ,237.8917 ,239.1019 ,240.2865 ,241.4467&
       ,242.5838 ,243.6985 ,244.7920 ,245.8651 ,246.9186&
       ,247.9533 ,248.9700 ,249.9693 ,250.9519 ,251.9184&
       ,252.8693 ,253.8053 ,254.7268 ,255.6343 ,256.5284&
       ,257.4093 ,258.2776 ,259.1337 ,259.9779 ,260.8106&
       ,261.6321 ,262.4428 ,263.2430 ,264.0329 ,264.8129&
       ,265.5833 ,266.3442 ,267.0961 ,267.8390 ,268.5733/
      data tt2/&
        269.2991 ,270.0166 ,270.7262 ,271.4279 ,272.1220&
       ,272.8087 ,273.4880 ,274.1603 ,274.8256 ,275.4841&
       ,276.1360 ,276.7814 ,277.4205 ,278.0533 ,278.6801&
       ,279.3010 ,279.9160 ,280.5253 ,281.1291 ,281.7274&
       ,282.3203 ,282.9080 ,283.4905 ,284.0679 ,284.6405&
       ,285.2081 ,285.7710 ,286.3292 ,286.8828 ,287.4319&
       ,287.9766 ,288.5170 ,289.0530 ,289.5849 ,290.1126&
       ,290.6363 ,291.1560 ,291.6718 ,292.1837 ,292.6918/
      data tt3/&
        270.6500 ,257.9200 ,249.4700 ,243.6600 ,239.2500&
       ,235.7100 ,232.7400 ,230.2000 ,228.4100 ,227.7100&
       ,227.0744 ,226.4968 ,225.9667 ,225.4771 ,225.0222&
       ,224.5975 ,224.1993 ,223.8245 ,223.4705 ,223.1353&
       ,222.8168 ,222.5136 ,222.2243 ,221.9476 ,221.6826&
       ,221.4282 ,221.1837 ,220.9484 ,220.7216 ,220.5026&
       ,220.2911 ,220.0865 ,219.8883 ,219.6963 ,219.5099&
       ,219.3290 ,219.1532 ,218.9822 ,218.8157 ,218.6536&
       ,218.4956 ,218.3416 ,218.1912 ,218.0444 ,217.9011&
       ,217.7609 ,217.6239 ,217.4898 ,217.3586 ,217.2301&
       ,217.1042 ,216.9808 ,216.8599 ,216.7413 ,216.6500&
       ,216.6500 ,216.6500 ,216.6500 ,216.6500 ,216.6500&
       ,216.6500 ,216.6500 ,216.6500 ,216.6500 ,216.6500&
       ,216.6500 ,216.6500 ,216.6500 ,216.6500 ,216.6500/
      data tt4/&
        216.6500 ,216.6500 ,216.6500 ,216.6500 ,216.6500&
       ,216.6500 ,216.6500 ,216.6500 ,216.6500 ,216.6500&
       ,216.6500 ,216.6500 ,216.6500 ,216.6500 ,216.6500&
       ,216.6500 ,216.6500 ,216.6500 ,216.6500 ,216.6500&
       ,216.6500 ,216.6500 ,216.6500 ,216.6500 ,216.6500&
       ,216.6500 ,216.6500 ,216.6500 ,216.6500 ,216.6500&
       ,216.6500 ,216.6500 ,216.6500 ,216.6500 ,216.6500&
       ,216.6500 ,216.6500 ,216.6500 ,216.6500 ,216.6500/
      data tt5/   231.6000 ,245.1400 ,253.4200 ,259.4700&
       ,264.2600 ,268.2300 ,270.6500 ,270.6500 ,270.6500&
       ,270.6500/   

      data one/1.0/
      data hndrd/100.0/
      data tenth/0.1/
!
!***********************************************************************
!          initialize real(kind=r8) lookup index for pressures >= 100 mb
!***********************************************************************
!
      do i=1,m
        ttin(i)=pin(i)*tenth-one
      enddo
!
!***********************************************************************
!          initialize real(kind=r8) lookup index for pressures < 100 mb
!***********************************************************************
!
      do i=1,m
        if (pin(i).lt.hndrd) ttin(i)=pin(i)-one
      enddo
!
!***********************************************************************
!          initialize real(kind=r8) lookup index for pressures < 1 mb
!***********************************************************************
!
      do i=1,m
        if (pin(i).lt.one) ttin(i)=10.*pin(i)-one
      enddo
!
!***********************************************************************
!          convert real(kind=r8) lookup index to integer and calculate 
!      pressure weighting factor
!***********************************************************************
!
      do i=1,m
        ii(i)=int(ttin(i))
        pw(i,2)=ttin(i)+one
        ttin(i)=float(ii(i)+1)
        pw(i,2)=pw(i,2)/ttin(i)
        ttin(i)=(ttin(i)+one)/ttin(i)
        pw(i,2)=log(pw(i,2))/log(ttin(i))
      enddo
!
      do i=1,m
        if (pin(i).lt.hndrd.and.pin(i).ge.one) ii(i)=ii(i)+110
      enddo
      do i=1,m
        if (pin(i).lt.one) ii(i)=ii(i)+220
      enddo
!
!***********************************************************************
!          lookup temperature corresponding to integer lookup index
!***********************************************************************
!
      do i=1,m
        ttin(i)=tt(ii(i)+1)
        pw(i,1)=tt(ii(i)+2)
      enddo
!
!***********************************************************************
!          weight the temperature looked up and return
!***********************************************************************
!
      do i=1,m
        ttin(i)=ttin(i)+pw(i,2)*(pw(i,1)-ttin(i))
      enddo
!
!***********************************************************************
!
      return
end subroutine tstd



subroutine vor(u,v,m,n,kk,delx,dely,hxm,hym,fm,sigmma,sigmwa&
                     ,dsigw,zsfc,iflag,rvor)
!
!     compute relative rvor
!
      implicit none

      integer i
      integer j
      integer k
      integer m
      integer n
      integer kk
      integer iflag

      real(kind=r8) delx
      real(kind=r8) dely
      real(kind=r8) scale
      real(kind=r8) dvdx
      real(kind=r8) dudy

      real(kind=r8) hxu(m,n)
      real(kind=r8) hyu(m,n)
      real(kind=r8) hxv(m,n)
      real(kind=r8) hyv(m,n)

      real(kind=r8) hxm(m,n)
      real(kind=r8) hym(m,n)
      real(kind=r8) fm(m,n)

      real(kind=r8) u(m,n,kk)
      real(kind=r8) v(m,n,kk)
      real(kind=r8) rvor(m,n,kk)
      real(kind=r8) aoz(m,n)
      real(kind=r8) sigmwa(kk+1)
      real(kind=r8) sigmma(kk)
      real(kind=r8) dsigw(kk+1)
      real(kind=r8) zsfc(m,n)
      real(kind=r8) zinm(m,n,kk)
      real(kind=r8) dudzu
      real(kind=r8) dudzl
      real(kind=r8) dudz1
      real(kind=r8) dudz2
      real(kind=r8) dudz
      real(kind=r8) dvdzu
      real(kind=r8) dvdzl
      real(kind=r8) dvdz1
      real(kind=r8) dvdz2
      real(kind=r8) dvdz
      real(kind=r8) dzdx
      real(kind=r8) dzdy

      call hm2uv(hxm,hym,hxu,hyu,hxv,hyv,m,n)

! zinm is the height of the sigmma levels.

      do i=1,m
        do j=1,n
         aoz(i,j)=sigmwa(1)/(sigmwa(1)-zsfc(i,j))
        enddo
      enddo

      do k=1,kk
        do i=1,m
          do j=1,n
             zinm(i,j,k)=sigmma(k)*(sigmwa(1)-zsfc(i,j))&
                /sigmwa(1)+zsfc(i,j)
          enddo
        enddo
      enddo

      do k=1,kk
        do i=2,m-1
          do j=2,n-1
!
! compute du/dz
!
            if(k.eq.1) then
              dudzu=(u(i,j,1)-u(i,j,2))*aoz(i,j)/dsigw(2)
              dudzl=(u(i,j,1)-u(i,j,2))*aoz(i,j)/dsigw(2)
            elseif(k.eq.kk) then
              dudzu=(u(i,j,k-1)-u(i,j,k))*aoz(i,j)/dsigw(k)
              dudzl=(u(i,j,k-1)-u(i,j,k))*aoz(i,j)/dsigw(k)
            else
              dudzu=(u(i,j,k-1)-u(i,j,k  ))*aoz(i,j)/dsigw(k  )
              dudzl=(u(i,j,k  )-u(i,j,k+1))*aoz(i,j)/dsigw(k+1)
            endif
            dudz1=(dudzu+dudzl)/2.0
            if(k.eq.1) then
              dudzu=(u(i-1,j,1)-u(i-1,j,2))*aoz(i-1,j)/dsigw(2)
              dudzl=(u(i-1,j,1)-u(i-1,j,2))*aoz(i-1,j)/dsigw(2)
            elseif(k.eq.kk) then
              dudzu=(u(i-1,j,k-1)-u(i-1,j,k))*aoz(i-1,j)/dsigw(k)
              dudzl=(u(i-1,j,k-1)-u(i-1,j,k))*aoz(i-1,j)/dsigw(k)
            else
              dudzu=(u(i-1,j,k-1)-u(i-1,j,k  ))*aoz(i-1,j)/dsigw(k  )
              dudzl=(u(i-1,j,k  )-u(i-1,j,k+1))*aoz(i-1,j)/dsigw(k+1)
            endif
            dudz2=(dudzu+dudzl)/2.0
            dudz=(dudz1+dudz2)/2.0
!
! compute dv/dz
!
            if(k.eq.1) then
              dvdzu=(v(i,j,1)-v(i,j,2))*aoz(i,j)/dsigw(2)
              dvdzl=(v(i,j,1)-v(i,j,2))*aoz(i,j)/dsigw(2)
            elseif(k.eq.kk) then
              dvdzu=(v(i,j,k-1)-v(i,j,k))*aoz(i,j)/dsigw(k)
              dvdzl=(v(i,j,k-1)-v(i,j,k))*aoz(i,j)/dsigw(k)
            else
              dvdzu=(v(i,j,k-1)-v(i,j,k  ))*aoz(i,j)/dsigw(k  )
              dvdzl=(v(i,j,k  )-v(i,j,k+1))*aoz(i,j)/dsigw(k+1)
            endif
            dvdz1=(dvdzu+dvdzl)/2.0
            if(k.eq.1) then
              dvdzu=(v(i,j-1,1)-v(i,j-1,2))*aoz(i,j-1)/dsigw(2)
              dvdzl=(v(i,j-1,1)-v(i,j-1,2))*aoz(i,j-1)/dsigw(2)
            elseif(k.eq.kk) then
              dvdzu=(v(i,j-1,k-1)-v(i,j-1,k))*aoz(i,j-1)/dsigw(k)
              dvdzl=(v(i,j-1,k-1)-v(i,j-1,k))*aoz(i,j-1)/dsigw(k)
            else
              dvdzu=(v(i,j-1,k-1)-v(i,j-1,k  ))*aoz(i,j-1)/dsigw(k  )
              dvdzl=(v(i,j-1,k  )-v(i,j-1,k+1))*aoz(i,j-1)/dsigw(k+1)
            endif
            dvdz2=(dvdzu+dvdzl)/2.0
            dvdz=(dvdz1+dvdz2)/2.0

             dzdx=(zinm(i+1,j,k)*hxm(i+1,j)-zinm(i-1,j,k)*hxm(i-1,j))&
                  /delx
             dzdy=(zinm(i,j+1,k)*hym(i,j+1)-zinm(i,j-1,k)*hym(i,j-1))&
                  /dely
             dvdx=(v(i+1,j,k)*hyv(i+1,j)-v(i,j,k)*hyv(i,j))/delx&
                  -dzdx*dvdz
             dudy=(u(i,j+1,k)*hxu(i,j+1)-u(i,j,k)*hxu(i,j))/dely&
                  -dzdy*dudz
             rvor(i,j,k)=(dvdx-dudy)/(hxm(i,j)*hym(i,j))
          enddo
        enddo

        do j=2,n-1
           rvor(1,j,k)=rvor(2,j,k)*2.0-rvor(3,j,k)
           rvor(m,j,k)=rvor(m-1,j,k)*2.0-rvor(m-2,j,k)
        enddo

        do i=2,m-1
           rvor(i,1,k)=rvor(i,2,k)*2.0-rvor(i,3,k)
           rvor(i,n,k)=rvor(i,n-1,k)*2.0-rvor(i,n-2,k)
        enddo

        rvor(1,1,k)=(rvor(1,2,k)+rvor(2,1,k))*0.5
        rvor(m,1,k)=(rvor(m-1,1,k)+rvor(m,2,k))*0.5
        rvor(1,n,k)=(rvor(2,n,k)+rvor(1,n-1,k))*0.5
        rvor(m,n,k)=(rvor(m-1,n,k)+rvor(m,n-1,k))*0.5

!
!  average to the mass point
!

        do j=2,n-1
           do i=2,m-1
              rvor(i,j,k)=(rvor(i,j,k)+rvor(i-1,j,k)&
              +rvor(i,j-1,k)+rvor(i-1,j-1,k))*0.25
           enddo
        enddo
!
!***********************************************************************
!  iflag: 1 - relative vorticity
!         2 - absolute vorticity
!***********************************************************************
!
      enddo

      if(iflag.eq.2) then
         do k=1,kk
          do i=1,m
            do j=1,n
              rvor(i,j,k)=rvor(i,j,k)+fm(i,j)
            enddo
          enddo
         enddo
      endif
!
      return
end subroutine vor



subroutine filt9(z1,m,n,nfpass,isdd)

! rcs keywords: $RCSfile: filt9.f,v $ 
!               $Revision: 1.1 $ $Date: 1997/12/05 20:19:57 $c
      implicit none
!
!***********************************************************************
!           parameters:
!***********************************************************************
!
      integer isd
      integer isdd
      integer m
      integer n
      integer nfpass
!
      real(kind=r8) z1     (m,n)
!
!***********************************************************************
!          local variables and dynamic storage:
!***********************************************************************
!
      integer i
      integer j
      integer lenx
      integer leny
      integer m1
      integer nf
      integer np
!
      real(kind=r8) s1
      real(kind=r8) z2     (m,n)
!
!***********************************************************************
!        9 point smoother/desmoother
!
!        input variables:
!
!          z1     : first word address of array to be filtered
!          m      : first dimension of array z1
!          n      : second dimension of array z1
!          nfpass : number of passes through filter
!          isd=1  : 9 point smoother
!          isd=2  : 9 point smoother-desmoother
!***********************************************************************
!
      if (nfpass.le.0) return
!
!***********************************************************************
!        local variables
!***********************************************************************
!
      isd=isdd
      if (isd.le.1) isd=1
      if (isd.ge.2) isd=2
!
      m1=m-1
      lenx=m*n-1
      leny=m*(n-2)
!
!***********************************************************************
!        start loop for number of passes
!***********************************************************************
!
      do np=1,nfpass
!
!***********************************************************************
!        filter the field
!***********************************************************************
!
        do nf=1,isd
          if (nf.eq.1) then
            s1=0.25
          else
            s1=-0.25
          endif
!
!***********************************************************************
!        filter in the x-direction
!***********************************************************************
!
          do i=2,lenx
            z2(i,1)=z1(i,1)+(z1(i-1,1)+z1(i+1,1)-2.0*z1(i,1))*s1
          enddo
!
          do j=1,n
            do i=2,m1
              z1(i,j)=z2(i,j)
            enddo
          enddo
!
!***********************************************************************
!        filter in the y-direction
!***********************************************************************
!
          do i=1,leny
            z2(i,2)=z1(i,2)+(z1(i,1)+z1(i,3)-2.0*z1(i,2))*s1
          enddo
!
          do i=1,leny
            z1(i,2)=z2(i,2)
          enddo
!
        enddo
!
!***********************************************************************
!        end of loop over number of passes
!***********************************************************************
!
      enddo
!
!***********************************************************************
!
      return
end subroutine filt9

  ! define_mean_exner
  ! ---------------------
  ! Defines the mean exner based on the standard 
  ! atmosphere on sigma surfaces in a 3D domain.
  !  PARAMETERS
  !   OUT  exner   Mean exner function 
  !   IN   mn      Horizontal dimension
  !   IN   kk      Vertical levels
  !   IN   ztop    Domain top
  !   IN   zsfc    Topography height
  !   IN   sigma   Sigma levels
  subroutine define_mean_exner(exner, mn, kk, ztop, zsfc, sigma)
    real(kind=r8), dimension(:,:), intent(out) :: exner
    integer,                       intent(in)  :: mn
    integer,                       intent(in)  :: kk
    real(kind=r8),                 intent(in)  :: ztop
    real(kind=r8), dimension(:),   intent(in)  :: zsfc
    real(kind=r8), dimension(:),   intent(in)  :: sigma

    real(kind=r8) :: zout(mn)
    real(kind=r8) :: pres(mn)
    integer       :: k

    do k=1,kk
      zout(:)=( sigma(k)*(ztop-zsfc(:))/ztop + zsfc(:) )*g
      call pstd(mn, zout, pres)
      exner(:,k)=( (pres(:)*100.d0 )/p00)**rocp
    enddo
  end subroutine define_mean_exner

  ! define_mean_theta
  ! ---------------------
  ! Defines the mean theta based on the standard 
  ! atmosphere on sigma surfaces in a 3D domain.
  !  PARAMETERS
  !   OUT  theta   Mean theta function 
  !   IN   mn      Horizontal dimension
  !   IN   kk      Vertical levels
  !   IN   ztop    Domain top
  !   IN   zsfc    Topography height
  !   IN   sigma   Sigma levels
  subroutine define_mean_theta(theta, mn, kk, ztop, zsfc, sigma)
    real(kind=r8), dimension(:,:), intent(out) :: theta
    integer,                       intent(in)  :: mn
    integer,                       intent(in)  :: kk
    real(kind=r8),                 intent(in)  :: ztop
    real(kind=r8), dimension(:),   intent(in)  :: zsfc
    real(kind=r8), dimension(:),   intent(in)  :: sigma

    real(kind=r8) :: zout(mn)
    real(kind=r8) :: pres(mn)
    integer       :: k

    do k=1,kk
      zout(:)=( sigma(k)*(ztop-zsfc(:))/ztop + zsfc(:) )*g
      call pstd(mn, zout, pres)
      call tstd(mn, pres, zout)
      theta(:,k)=zout(:)*(p00/(pres(:)*100.d0))**rocp
    enddo
  end subroutine define_mean_theta

  ! compute_altimeter
  ! ---------------------
  ! Computes the altimeter of a pressure value at
  ! a gemotric height
  !  PARAMETERS
  !   IN   pres    pressure value  (Pa)
  !   IN   zhgt    pressure height (m)
  function compute_altimeter(pres, zhgt) result(altimeter)
    real(kind=r8), intent(in) :: pres
    real(kind=r8), intent(in) :: zhgt
    real(kind=r8)             :: altimeter

    real(kind=r8), parameter :: k1 = 0.190284_r8
    real(kind=r8), parameter :: k2 = 8.4228807E-5_r8

    real(kind=r8), parameter  :: hPa_to_Pa= 100.0_r8

    altimeter = hPa_to_Pa * ((pres/hPa_to_Pa - 0.3_r8) ** k1 + k2 * zhgt) ** (1.0_r8 / k1)

    return
  end function compute_altimeter

end module coamps_intrinsic_mod
