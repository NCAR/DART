      subroutine euvheat(z,heuv,ij)
c
c----------called from abcr1 for every height, z
c Fills heuv array with real heuv heating for migrating tide in thermosphere
c jhackney 8/98
c
      implicit none
c
      integer i,m,n,iz,mx,ny,ij,f107
      real euvlat(36),reuv(36,36),ieuv(36,36),euvalt(36)
      real sigmat,dtheta, RtoDEG
      real deuv(36,36,3),heuv(91)
      real colat,z,heuv1d,heuvdx,heuvdy,heuvdxx,heuvdxy
      real heuvdyy
      real clatt(91),xlat(91),snclat(91),csclat(91),tnclat(91)
      real ctnclat(91),sin2i(91),gmlat(91),dip(91)
c
      common/euvi/mx,ny,iz,f107
      common/euvr/euvlat,euvalt,reuv,ieuv,deuv,sigmat
c
      common/latvect/dtheta,clatt,xlat,snclat,csclat,
     +  tnclat,ctnclat,sin2i,gmlat,dip
c
      RtoDEG=57.2958
c
      do i=1,ij
         colat=clatt(i)*RtoDEG
c
         call surfd(colat,z,heuv1d,heuvdx,heuvdy,heuvdxx,heuvdxy,
     +        heuvdyy,mx,ny,euvlat,euvalt,reuv,iz,deuv,sigmat)
c
         heuv(i)=heuv1d
c
      end do
      return
      end
