      subroutine trajbc(x,y,z,nxext,ny,nzp2,p2,p,ph,dxrad,dyrad,ymin)
c
c--------------------------------------------------------------------
c  trajectories for those departure points
c  whose latitudinal angle is outside the domain
c--------------------------------------------------------------------
c
      implicit none

      integer :: nxext, ny, nzp2
      real, dimension(nxext,ny,nzp2) :: x, y, z
      real :: p2, p, ph, dxrad, dyrad, ymin, temp0, temp1, temp2, tempx
      integer :: i, j, k
c
c  conversion from grid units to longitudinal and latitudinal angles

      do i=2,nxext-1
         do j=1,ny
            do k=1,nzp2
               x(i,j,k) = (x(i,j,k)-1.)*dxrad
               y(i,j,k) = ymin + dyrad*( y(i,j,k) - 1. )
            end do
         end do
      end do

c map x to an angle between 0 and 2 pi
      do i=2,nxext-1
         do j=1,ny
            do k=1,nzp2
               temp0 = x(i,j,k)
               temp0 = temp0 - p2*ifix(temp0/p2)
               if(temp0.lt.0.) temp0 = temp0 + p2
               x(i,j,k) = temp0
            end do
         end do
      end do

c  temp1 is rotated about the north pole by pi 
c  temp2 reflects the latitudinal angle about pi/2 (eg. 91 to 89 degrees)
      do i=2,nxext-1
         do j=1,ny
            do k=1,nzp2
               temp0 = y(i,j,k)
               if(abs(temp0).gt.ph)then
                  tempx = x(i,j,k)
                  if (tempx.lt.p) temp1 = tempx + p
                  if (tempx.ge.p) temp1 = tempx - p
                  temp2 = sign(p,temp0) - temp0
                  x(i,j,k)=temp1
                  y(i,j,k)=temp2
               endif
            end do
         end do
      end do

c  convert angles back to grid units.
      do i=2,nxext-1
         do j=1,ny
            do k=1,nzp2
               x(i,j,k) = x(i,j,k)/dxrad+1.
               y(i,j,k) = (y(i,j,k)-ymin)/dyrad+1.
            end do
         end do
      end do

c vertical b.c. -  maps departure points outside the upper and lower
c levels to the top and bottom levels
      temp0 = float(nzp2)+0.49
      do i=2,nxext-1
         do j=1,ny
            do k=1,nzp2
               temp1 = amin1( z(i,j,k), temp0 )
               z(i,j,k) = amax1( temp1, 0.51 )
            end do
         end do
      end do

      call trxbc(x,nxext,ny,nzp2)
      call trxbc(y,nxext,ny,nzp2)
      call trxbc(z,nxext,ny,nzp2)

      end


      subroutine traject(delt,nbiter)
c
c-----------------------------------------------------------------------
c   set nbiter to 1 or 2 (higher is more accurate and costlier)
c-----------------------------------------------------------------------
c

      use params
      use dynam

      implicit none

      integer, parameter :: nxext=nx+2

      real :: delt
      integer :: nbiter
      real, dimension(nxext,ny,nzp2) :: u1, u2, u3
      real, dimension(nxext,ny,nzp2) :: vl1, vl2, vl3
      real, dimension(nxext,ny,nzp2) :: xm, ym, zm
      real, dimension(nxext,ny,nzp2) :: x0, y0, z0
      common / dep / x0, y0, z0

      real :: pi2, pih, re, deltn, temp1, temp2, temp3
      integer :: i, i1, j, k, iter, it, nuscal, nuvect
c
c----------------------------------------------------------
c   dxrad : angular size of a longitudinal interval, (set in msetfix)
c   dyrad : angular size of a latitude interval.
c----------------------------------------------------------
c
      pi2 = 2.*pi
      pih = .5*pi
c
c-----------------------------------------------------------
c   re is radius of the earth in meters.
c-----------------------------------------------------------
c
      re = 6.37e6
      deltn = delt/nbiter
c
c-----------------------------------------------------------
c  from the winds u,v,w the actual distance these winds travel
c    in a timestep deltn is calculated.
c  u1 : u1 = (u*deltn)/(re*cosphi*dxrad)
c    distance, in longitudinal grid intervals, covered by
c    an air parcel in a deltn timestep.
c  u2 : u2 = (-v*deltn)/(re*dyrad) 
c    fraction of a latitudinal interval traveled by a
c    air parcel with meridional wind v.
c  u3 : number of vertical grids covered by the wind in
c     a timestep.
c-----------------------------------------------------------

      temp2 = deltn*gcy
      temp3 = deltn*gcz
      do j=1,ny
         temp1 = deltn*gcx(j)
         do i=2,nxext-1
            i1 = i-1
            do k=1,nzp2
               u1(i,j,k) = uusm(k,i1,j)*temp1
               u2(i,j,k) = vvsm(k,i1,j)*temp2
               u3(i,j,k) = wwsm(k,i1,j)*temp3
            end do
         end do
      end do

      call trxbc (u1,nxext,ny,nzp2)
      call trxbc (u2,nxext,ny,nzp2)
      call trxbc (u3,nxext,ny,nzp2)
      do k=1,nzp2
         do j=1,ny
            do i=2,nxext-1
               x0(i,j,k)=i
               y0(i,j,k)=j
               z0(i,j,k)=k
            end do
         end do
      end do

      do iter=1,nbiter
         do i=2,nxext-1
            do j=1,ny
               do k=1,nzp2
                  xm(i,j,k)=x0(i,j,k)-0.5*u1(i,j,k)
                  ym(i,j,k)=y0(i,j,k)-0.5*u2(i,j,k)
                  zm(i,j,k)=z0(i,j,k)-0.5*u3(i,j,k)
               end do
            end do
         end do

         call trajbc(xm,ym,zm,nxext,ny,nzp2,pi2,pi,pih,dxrad,dyrad,ymin)
         
         do it=1,2

c  resets vl1,vl2,vl3 equal to the wind distances at the grid points
            do i=1,nxext
               do j=1,ny
                  do k=1,nzp2
                     vl1(i,j,k)= u1(i,j,k)
                     vl2(i,j,k)= u2(i,j,k)
                     vl3(i,j,k)= u3(i,j,k)
                  end do
               end do
            end do

c   for interpolation use trord1, trord2 or trord3 depending on the
c   order desired (higher order is more accurate, slower)
c   the vl1,vl2,vl3 go into the interpolation being the wind distances 
c   defined at the grid points, but come out being the winds defined at 
c   the midpoints xm,ym,zm

            nuscal=1
            nuvect=-1
            call trle3(vl1,xm,ym,zm,nuvect,nuvect,1,1)
            call trle3(vl2,xm,ym,zm,nuvect,nuvect,1,1)
            call trle3(vl3,xm,ym,zm,nuscal,nuscal,1,1)
            do i=2,nxext-1
               do j=1,ny
                  do k=1,nzp2
                     xm(i,j,k)=x0(i,j,k)-0.5*vl1(i,j,k)
                     ym(i,j,k)=y0(i,j,k)-0.5*vl2(i,j,k)
                     zm(i,j,k)=z0(i,j,k)-0.5*vl3(i,j,k)
                  end do
               end do
            end do

            call trajbc(xm,ym,zm,nxext,ny,nzp2,pi2,pi,pih,dxrad,
     $                  dyrad,ymin)
         end do

         do i=2,nxext-1
            do j=1,ny
               do k=1,nzp2
                  x0(i,j,k)=2.*xm(i,j,k)-x0(i,j,k)
                  y0(i,j,k)=2.*ym(i,j,k)-y0(i,j,k)
                  z0(i,j,k)=2.*zm(i,j,k)-z0(i,j,k)
               end do
            end do
         end do

         call trajbc(x0,y0,z0,nxext,ny,nzp2,pi2,pi,pih,dxrad,dyrad,ymin)

      end do

      end


      subroutine transport
c
c----------------------------------------------------------------------
c     calculates transport of species
c     semi-lagrangian transport scheme (slt)
c     by smolarkiewicz, ncar
c
c     WARNING -- species with index = 0 are NOT TRANSPORTED
c
c...  constituents (transport and chemistry)
c      1  n2o       2  ch4        3  h2o       4  o1d
c      5  hno3      6  n2o5       7  h         8  oh
c      9  co       10  hcl       11  clono2   12  hocl
c     13  h2o2     14  ho2       15  ho2no2   16  h2
c     17  ch2o     18  o         19  o3       20  cl
c     21  clo      22  n         23  no       24  no2
c     25  no3 
c----------------------------------------------------------------------

      use params
      use chem
      use dynam

      implicit none

      integer, parameter :: nxext=nx+2
      real :: dttrns 
      integer :: nbiter, nc, i, i1, j, k, nuvect, nuscal
      real, dimension(nxext,ny,nzp2) :: x0, y0, z0
      common / dep / x0, y0, z0

      real, dimension(nxext,ny,nzp2) :: xx 
      integer, dimension(nbcon) :: iop, index

      data index /1, 1, 1, 0,    ! n2o    ch4     h2o      o1d 
     $            1, 1, 1, 1,    ! hno3   n2o5    h        oh
     $            1, 1, 1, 1,    ! co     hcl     clono2   hocl
     $            1, 1, 1, 1,    ! h2o2   ho2     ho2no2   h2
     $            1, 1, 1, 1,    ! ch2o   o       o3       cl
     $            1, 1, 1, 1,    ! clo    n       no       no2
     $            1, 1 /         ! no3    o2

      data iop   /1, 1, 1, 0,    ! n2o    ch4     h2o      o1d 
     $            1, 1, 1, 1,    ! hno3   n2o5    h        oh
     $            1, 1, 1, 1,    ! co     hcl     clono2   hocl
     $            1, 1, 1, 1,    ! h2o2   ho2     ho2no2   h2
     $            1, 1, 1, 1,    ! ch2o   o       o3       cl
     $            1, 1, 1, 1,    ! clo    n       no       no2
     $            1, 1 /         ! no3    o2

c----------------------------------------------------------------------

      dttrns = ntrans * deltat
!.....................................................................
c  Extend dynamics arrays in vertical
      do j=1,ny
         do i=1,nx
            uusm(1,i,j) = ubc(i,j)
            vvsm(1,i,j) = vbc(i,j)
            wwsm(1,i,j) = ww(1,i,j)
            do k=2,nz+1
               uusm(k,i,j) = un1(k-1,i,j)
               vvsm(k,i,j) = vn1(k-1,i,j)
               wwsm(k,i,j) = ww(k-1,i,j)
            end do
            uusm(nzp2,i,j) = uubc(j)
            vvsm(nzp2,i,j) = 0.
            wwsm(nzp2,i,j) = 0.
         end do
      end do
!.....................................................................

c  get departure points for all species.
      nbiter = 1
      call traject(dttrns,nbiter)

      nuvect=-1
      nuscal=1


      do nc = 1,nbcon
         if(index(nc).ne.0)then
            do i=2,nxext-1
               i1 = i-1
               do j=1,ny
                  do k=1,nz
                     xx(i,j,k+1) = qn1(k,i1,j,nc)
                  end do
                  xx(i,j,1) = qlbc(i1,j,nc)
                  xx(i,j,nzp2) = qubc(i1,j,nc)
               end do
            end do

            call trxbc(xx,nxext,ny,nzp2)
            call trle3(xx,x0,y0,z0,nuscal,nuscal,1,iop(nc))

            do i=2,nxext-1
               i1 = i-1
               do j=1,ny
                  do k=1,nz
                     qn1(k,i1,j,nc) = xx(i,j,k+1)
                  end do
               end do
            end do
         end if
      end do

      end


      subroutine trle3(xx,xd,yd,zd,nu1,nu2,ior,nonos)
c
c----------------------------------------------------------------------------
c  trle3 interpolates to find the field values at the points xd,yd,zd.
c
c  nonos=1 is the "non-oscillatory option", i.e., the monotonicity of
c   the solution is guaranteed and the required 
c   computer time is increased by a factor of about 3
c   recommended:  nonos = 1 for short-lived species (e.g. nox)
c                 nonos = 0 for the long-lived species (e.g. ch4 and cfcs).
c
c  ior indicates roughly to what order the expansion goes. or how many grid
c   distances away from the nearest grid point you go to use in the 
c   interpolation, ior = 1 only uses the 14 grid points immediately 
c   surrounding the nearest grid point.
c--------------------------------------------------------------------------
c
c

      use params

      implicit none

      integer :: ior, nonos
      integer, parameter :: io=3
      integer, parameter :: iom=-io+1
      integer, parameter :: nxext=nx+2 
      integer, parameter :: inx=nxext+io, iny=ny+io, inz=nzp2+io
      integer, parameter :: nxyzext=nxext*ny*nzp2 
      integer, parameter :: nxyext = nxext*ny
      real, dimension(nxyzext) :: xf, xd1, yd1, zd1
      integer, dimension(nxext) :: ip
      real, dimension(nxext,ny,nzp2) :: xx, xd, yd, zd
      real, dimension(nxyzext,-io:io) :: z, y
      real :: ep = 1.e-30
      real :: mx, mn
      integer, dimension(nxyzext) :: ig0, jg0, kg0
      real, dimension(iom:inx,iom:iny,iom:inz) :: x
ccc      real, dimension(-2:37,-2:39,-2:43) :: x

      real :: y1, y2, a, xi

      real :: tr2
      real :: u, ym1, y0, yp1, f0, f1, fl1, fl0, w
      real :: zzm, zz0, zzp, ov, un, yy0, yyp, yym
      real :: pp, pn, pp0, pp1, pn0, pn1
      integer :: i, ii, ii0, ii1, is, is1, isz, isp, ism, itemp
      integer :: j, jior, jg0j
      integer :: k, kior, kg0k
      integer :: nu1, nu2, nior1

c statement functions

      tr2(y1,y2,a) = a*.5*(y1+y2-a*(y2-y1))

c     tr4(ym1,y0,yp1,yp2,a)=
c    1    a/12.*(7.*(yp1+y0)-(yp2+ym1)
c    2   -a*((15.*(yp1-y0)-(yp2-ym1))/2.+a*((yp1+y0)
c    3   -(yp2+ym1)-a/2.*(3.*(yp1-y0)-(yp2-ym1)))))
c
c      tr6(ym2,ym1,y0,yp1,yp2,yp3,a) =
c    1   -a*((-ym2+8.*(ym1+yp2)-37.*(y0+yp1)-yp3)/60.
c    2   +a*((2.*(yp3-ym2)+25.*(ym1-yp2)-245.*(y0-yp1))/360.
c    3   +a*((ym2-7.*(ym1+yp2)+6.*(y0+yp1)+yp3)/48.
c    4   +a*((ym2-11.*(ym1-yp2)+28.*(y0-yp1)-yp3)/144.
c    5   +a*((-ym2+3.*(ym1+yp2)-2.*(y0+yp1)-yp3)/240.
c    6   +a*((-ym2+5.*(ym1-yp2)-10.*(y0-yp1)+yp3)/720.))))))

      pp(xi) = amax1(0.,xi)
      pn(xi) = amin1(0.,xi)

      nior1 = -ior + 1
      jior = ny+ior
      kior = nzp2+ior
      
      do i=1,nxext
         ip(i) = mod(i+nxhlf-1,nx) + 1
      end do

      do k=1,nzp2
         ii0 = (k-1)*nxyext
         do i=1,nxext
            ii1 = ii0+i
            do j=1,ny
               ii = ii1 + (j-1)*nxext
               xd1(ii) = xd(i,j,k)
               yd1(ii) = yd(i,j,k)
               zd1(ii) = zd(i,j,k)
               xf(ii) = xx(i,j,k) 
            end do
         end do
      end do

      do k=1,nxyzext
         ig0(k) = nint(xd1(k))
         jg0(k) = nint(yd1(k))
         kg0(k) = nint(zd1(k))
      end do

      do k=1,nzp2      
         ii0 = (k-1)*nxyext
         do i=1,nxext
            ii1 = ii0 + i
            do j=1,ny
               x(i,j,k) = xf (ii1+(j-1)*nxext) 
            end do
         end do
      end do

      do is=1,ior
         is1 = 1-is
         isz = nzp2+is
         do j=1,ny
            do i=1,nxext
               x(i,j,is1) = x(i,j,1)
               x(i,j,isz) = x(i,j,nzp2)
            end do
         end do
      end do

      do is=nior1,0
         is1 = 1-is
         do k=nior1,kior
            do i=1,nxext
               x(i,is,k) = x(ip(i),is1,k) * nu1
            end do
         end do
      end do

      do is=1,ior
         isp = ny + is
         ism = nyp1 - is
         do k=nior1,kior
            do i=1,nxext
               x(i,isp,k) = x(ip(i),ism,k) * nu2
            end do
         end do
      end do

      do k=nior1,kior
         do j=nior1,jior
            do is=nior1,0
               x(is,j,k) = x(nx+is,j,k)
            end do
            do is=1,ior
               x(nxext+is,j,k) = x(is+2,j,k)
            end do
         end do
      end do
c
c  end of grid extension
c  here starts residual advection
c
      do k=-ior,ior
         do j=-ior,ior

            if(nonos.eq.1) then
               do ii=1,nxyzext
                  itemp = ig0(ii)
                  jg0j = jg0(ii) + j
                  kg0k = kg0(ii) + k
                  u = itemp - xd1(ii)
                  ym1 = x(itemp-1,jg0j, kg0k)
                  y0 =  x(itemp,  jg0j, kg0k)
                  yp1 = x(itemp+1,jg0j, kg0k)
                  f0 = tr2(ym1, y0 ,u)
                  f1 = tr2(y0  ,yp1,u)
                  if(u.lt.0.)then
                     fl0 = y0*u
                     fl1 = yp1*u
                  else
                     fl0 = ym1*u
                     fl1 = y0*u
                  endif
                  w = y0 -(fl1-fl0) 
                  mx = amax1(ym1,y0 ,yp1,w)
                  mn = amin1(ym1,y0 ,yp1,w)
                  f0 = f0-fl0
                  f1 = f1-fl1
                  pp0 = pp(f0)
                  pp1 = pp(f1)
                  pn0 = pn(f0)
                  pn1 = pn(f1)
                  ov = (mx-w)/(-pn1+pp0+ep)
                  un = (w-mn)/( pp1-pn0+ep)
                  ov = amin1(1.,ov)
                  un = amin1(1.,un)
                  f0 = pp0*ov + pn0*un
                  f1 = pp1*un + pn1*ov
                  z(ii,j) = w-(f1-f0)
               end do 
            else
               do ii=1,nxyzext
                  itemp = ig0(ii)
                  jg0j = jg0(ii) + j
                  kg0k = kg0(ii) + k
                  u = itemp - xd1(ii)
                  ym1 = x(itemp-1,jg0j,kg0k)
                  y0 =  x(itemp,jg0j,kg0k)
                  yp1 = x(itemp+1,jg0j,kg0k)
                  f0 = tr2(ym1, y0 ,u)
                  f1 = tr2(y0  ,yp1,u)
                  z(ii,j) = y0 - (f1-f0) 
               end do
            endif
         end do

         if(nonos.eq.1) then
            do ii=1,nxyzext
               u = jg0(ii)-yd1(ii)
               zzm = z(ii,-1)
               zz0 = z(ii,0)
               zzp = z(ii,1)
               f0 = tr2(zzm,zz0,u)
               f1 = tr2(zz0,zzp,u)
               if(u.lt.0.)then
                  fl0 = zz0*u
                  fl1 = zzp*u
               else
                  fl0 = zzm*u
                  fl1 = zz0*u
               endif
               w = zz0 - (fl1-fl0) 
               mx = amax1(zzm,zz0,zzp,w)
               mn = amin1(zzm,zz0,zzp,w)
               f0 = f0-fl0
               f1 = f1-fl1
               pp0 = pp(f0)
               pp1 = pp(f1)
               pn0 = pn(f0)
               pn1 = pn(f1)
               ov = (mx-w)/(-pn1+pp0+ep)
               un = (w-mn)/( pp1-pn0+ep)
               ov = amin1(1.,ov)
               un = amin1(1.,un)
               f0 = pp0*ov + pn0*un
               f1 = pp1*un + pn1*ov
               y(ii,k) = w - (f1-f0) 
            end do
         else
            do ii=1,nxyzext
               zz0 = z(ii,0)
               u = jg0(ii)-yd1(ii)
               f0 = tr2(z(ii,-1),zz0,u)
               f1 = tr2(zz0,z(ii,1),u)
               y(ii,k) = zz0-(f1-f0) 
            end do
         endif
      end do

      if(nonos.eq.1) then
         do ii=1,nxyzext
            u = kg0(ii)-zd1(ii)
            yym = y(ii,-1)
            yy0 = y(ii,0)
            yyp = y(ii,1)
            f0 = tr2(yym,yy0,u)
            f1 = tr2(yy0,yyp,u)
            if(u.lt.0.)then
               fl0 = yy0*u
               fl1 = yyp*u
            else
               fl0 = yym*u
               fl1 = yy0*u
            endif
            w = yy0 - (fl1-fl0) 
            mx = amax1(yym,yy0,yyp,w)
            mn = amin1(yym,yy0,yyp,w)
            f0 = f0-fl0
            f1 = f1-fl1
            pp0 = pp(f0)
            pp1 = pp(f1)
            pn0 = pn(f0)
            pn1 = pn(f1)
            ov = (mx-w)/(-pn1+pp0+ep)
            un = (w-mn)/( pp1-pn0+ep)
            ov = amin1(1.,ov)
            un = amin1(1.,un)
            f0 = pp0*ov + pn0*un
            f1 = pp1*un + pn1*ov
            xf(ii) = w-(f1-f0) 
         end do
      else
         do ii=1,nxyzext
            u = kg0(ii)-zd1(ii)
            yy0 = y(ii,0)
            f0 = tr2(y(ii,-1),yy0,u)
            f1 = tr2(yy0,y(ii,1),u)
            xf(ii) = yy0 - (f1-f0) 
         end do
      endif

      do k=1,nzp2
         ii0 = (k-1) * nxyext
         do i=1,nxext
            ii1 = ii0 + i
            do j=1,ny
               xx(i,j,k) = xf((j-1)*nxext+ii1) 
            end do
         end do
      end do
      end   


      subroutine trxbc(x,nxext,ny,nzp2)
c
c   periodic longitudinal boundary conditions for transport
c

      implicit none

      integer :: nxext, ny, nzp2
      real, dimension(nxext,ny,nzp2) :: x
      integer :: i, j, k

      do k=1,nzp2
         do j=1,ny
            x(1,j,k)     = x(nxext-1,j,k)
            x(nxext,j,k) = x(2,j,k)
         end do
      end do
      end
