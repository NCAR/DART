!
! inside function is by Alexander F. Shchepetkin
! provided by ali.aydogdu@cmcc.it
!
module utilities

! integer precision:
integer, parameter :: i2 = SELECTED_INT_KIND(3)
integer, parameter :: i4 = SELECTED_INT_KIND(8)
integer, parameter :: i8 = SELECTED_INT_KIND(13)

! real precision:
! TO RUN WITH REDUCED PRECISION REALS (and use correspondingly less memory)
! comment OUT the r8 definition below and use the second one:
integer, parameter :: r4 = SELECTED_REAL_KIND(6,30)
integer, parameter :: r8 = SELECTED_REAL_KIND(12)   ! real r8

public             :: skip_read_line, &  ! skip a number of lines while reading
                      to_radian,      &  ! degree to radian
                      earthdistance,  &  ! compute earth distance of two locations
                      inside             ! check if inside polygon by Alexander F. Shchepetkin

contains

 subroutine skip_read_line(ir,num)
   integer,    intent(in) :: num
   integer,    intent(in) :: ir
   SKIP: do i=1,num
       read(ir,*)
       end do SKIP
 end subroutine skip_read_line

 function to_radian(degree) result(rad)

!degrees to radians
   real(r8),intent(in) :: degree
   real, parameter :: deg_to_rad = atan(1.0)/45 ! exploit intrinsic atan to generate pi/180 runtime constant
   real :: rad

   rad = degree*deg_to_rad

 end function to_radian

!compute distance on earth
 real function earthdistance(zon1,zon2,mer1,mer2)

 use o_param, only : rad

 implicit none

  real(r8),intent(in) :: zon1, zon2, mer1, mer2
  real(r8)            :: lat1, lat2, dlat, dlon
  real(r8)            :: a, c
  !  real,parameter :: radius = 6372.8
  real,parameter :: radius = 6369.161

  lat1 = to_radian(zon1)
  lat2 = to_radian(zon2)
  dlat = to_radian(zon1-zon2)
  dlon = to_radian(mer1-mer2)
  a = (sin(dlat/2))**2 + cos(lat1)*cos(lat2)*(sin(dlon/2))**2
  c = 2*asin(sqrt(a))
  earthdistance = radius*c*1000 ! in meters
  !      distno=sin(dlat/2)**2+cos(lat1)*cos(lat2)*sin(dlon/2)**2
  !      earthdistance= ( acos(distno) / rad ) * 60 * 1.609344 * 1.1515
 return
 end function earthdistance


 function inside (Xo,Yo,Xb,Yb,Nb)
!
!svn $Id: inside.F 435 2010-01-02 15:55:10Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  Given the vectors Xb and Yb of size Nb, defining the coordinates    !
!  of a closed polygon,  this function find if the point (Xo,Yo) is    !
!  inside the polygon.  If the point  (Xo,Yo)  falls exactly on the    !
!  boundary of the polygon, it still considered inside.                !
!                                                                      !
!  This algorithm does not rely on the setting of  Xb(Nb)=Xb(1) and    !
!  Yb(Nb)=Yb(1).  Instead, it assumes that the last closing segment    !
!  is (Xb(Nb),Yb(Nb)) --> (Xb(1),Yb(1)).                               !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Reid, C., 1969: A long way from Euclid. Oceanography EMR,         !
!      page 174.                                                       !
!                                                                      !
!  Algorithm:                                                          !
!                                                                      !
!  The decision whether the point is  inside or outside the polygon    !
!  is done by counting the number of crossings from the ray (Xo,Yo)    !
!  to (Xo,-infinity), hereafter called meridian, by the boundary of    !
!  the polygon.  In this counting procedure,  a crossing is counted    !
!  as +2 if the crossing happens from "left to right" or -2 if from    !
!  "right to left". If the counting adds up to zero, then the point    !
!  is outside.  Otherwise,  it is either inside or on the boundary.    !
!                                                                      !
!  This routine is a modified version of the Reid (1969) algorithm,    !
!  where all crossings were counted as positive and the decision is    !
!  made  based on  whether the  number of crossings is even or odd.    !
!  This new algorithm may produce different results  in cases where    !
!  Xo accidentally coinsides with one of the (Xb(k),k=1:Nb) points.    !
!  In this case, the crossing is counted here as +1 or -1 depending    !
!  of the sign of (Xb(k+1)-Xb(k)).  Crossings  are  not  counted if    !
!  Xo=Xb(k)=Xb(k+1).  Therefore, if Xo=Xb(k0) and Yo>Yb(k0), and if    !
!  Xb(k0-1) < Xb(k0) < Xb(k0+1),  the crossing is counted twice but    !
!  with weight +1 (for segments with k=k0-1 and k=k0). Similarly if    !
!  Xb(k0-1) > Xb(k0) > Xb(k0+1), the crossing is counted twice with    !
!  weight -1 each time.  If,  on the other hand,  the meridian only    !
!  touches the boundary, that is, for example, Xb(k0-1) < Xb(k0)=Xo    !
!  and Xb(k0+1) < Xb(k0)=Xo, then the crossing is counted as +1 for    !
!  segment k=k0-1 and -1 for segment k=k0, resulting in no crossing.   !
!                                                                      !
!  Note 1: (Explanation of the logical condition)                      !
!                                                                      !
!  Suppose  that there exist two points  (x1,y1)=(Xb(k),Yb(k))  and    !
!  (x2,y2)=(Xb(k+1),Yb(k+1)),  such that,  either (x1 < Xo < x2) or    !
!  (x1 > Xo > x2).  Therefore, meridian x=Xo intersects the segment    !
!  (x1,y1) -> (x2,x2) and the ordinate of the point of intersection    !
!  is:                                                                 !
!                                                                      !
!                 y1*(x2-Xo) + y2*(Xo-x1)                              !
!             y = -----------------------                              !
!                          x2-x1                                       !
!                                                                      !
!  The mathematical statement that point  (Xo,Yo)  either coinsides    !
!  with the point of intersection or lies to the north (Yo>=y) from    !
!  it is, therefore, equivalent to the statement:                      !
!                                                                      !
!         Yo*(x2-x1) >= y1*(x2-Xo) + y2*(Xo-x1),   if   x2-x1 > 0      !
!  or                                                                  !
!         Yo*(x2-x1) <= y1*(x2-Xo) + y2*(Xo-x1),   if   x2-x1 < 0      !
!                                                                      !
!  which, after noting that  Yo*(x2-x1) = Yo*(x2-Xo + Xo-x1) may be    !
!  rewritten as:                                                       !
!                                                                      !
!        (Yo-y1)*(x2-Xo) + (Yo-y2)*(Xo-x1) >= 0,   if   x2-x1 > 0      !
!  or                                                                  !
!        (Yo-y1)*(x2-Xo) + (Yo-y2)*(Xo-x1) <= 0,   if   x2-x1 < 0      !
!                                                                      !
!  and both versions can be merged into  essentially  the condition    !
!  that (Yo-y1)*(x2-Xo)+(Yo-y2)*(Xo-x1) has the same sign as x2-x1.    !
!  That is, the product of these two must be positive or zero.         !
!                                                                      !
!=======================================================================
!
      implicit none
      integer Nstep
      parameter (Nstep=128)
!
      logical inside
      integer Nb, crossings, i, inc, k, kk, nc
      integer index(Nstep)
      real(r8) Xo, Yo, dx1, dx2, dxy
      real(r8)  Xb(Nb+1), Yb(Nb+1)
!
!-----------------------------------------------------------------------
!  Find intersections.
!-----------------------------------------------------------------------
!
!  Set crossings counter and close the contour of the polygon.
!
      crossings=0
      Xb(Nb+1)=Xb(1)
      Yb(Nb+1)=Yb(1)
!
!  The search is optimized.  First select the indices of segments
!  where Xb(k) is different from Xb(k+1) and Xo falls between them.
!  Then, further investigate these segments in a separate loop.
!  Doing it in two stages takes less time because the first loop is
!  pipelined.
!
      DO kk=0,Nb-1,Nstep
        nc=0
        DO k=kk+1,MIN(kk+Nstep,Nb)
          IF (((Xb(k+1)-Xo)*(Xo-Xb(k)).ge.0.0).and.                     &
     &        (Xb(k).ne.Xb(k+1))) THEN
            nc=nc+1
            index(nc)=k
          END IF
        END DO
        DO i=1,nc
          k=index(i)
          IF (Xb(k).ne.Xb(k+1)) THEN
            dx1=Xo-Xb(k)
            dx2=Xb(k+1)-Xo
            dxy=dx2*(Yo-Yb(k))-dx1*(Yb(k+1)-Yo)
            inc=0
            IF ((Xb(k).eq.Xo).and.(Yb(k).eq.Yo)) THEN
              crossings=1
              GO TO 10
            ELSE IF (((dx1.eq.0.0).and.(Yo.ge.Yb(k  ))).or.             &
     &              ((dx2.eq.0.0).and.(Yo.ge.Yb(k+1)))) THEN
              inc=1
            ELSE IF ((dx1*dx2.gt.0.0).and.                              &
     &              ((Xb(k+1)-Xb(k))*dxy.ge.0.0)) THEN
              inc=2
            END IF
            IF (Xb(k+1).gt.Xb(k)) THEN
              crossings=crossings+inc
            ELSE
              crossings=crossings-inc
            END IF
          END IF
        END DO
      END DO
!
!  Determine if point (Xo,Yo) is inside of closed polygon.
!
  10  IF (crossings.eq.0) THEN
        inside=.false.
      ELSE
        inside=.true.
      END IF
      RETURN
 end function inside

end module utilities
