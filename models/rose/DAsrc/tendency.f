      subroutine tendency
c     ------
c  compute time tendency for u, v, T
c     ------

      use params
      use dynam
      use phys

      implicit none

      integer :: i, j, k, icall, klo, ip1, im1, jp1, jm1
      real :: alphar, xdiv, dxhm1, dycos, deltfi
      real :: upx, umx, upy, umy, vpx, vmx, vpy, vmy
      real :: tpx, tmx, tpy, tmy, vpyc, vmyc
      real, dimension(nz,nx,ny) :: fu, fv, ft, ptemp, vertadv
      real, dimension(nz) :: fnpole, fspole, scoef, umid, vmid, tmid
      real :: sum, ysum, divx
      save ysum, divx
      data icall/0/, scoef/nz*0./
      data alphar/0.02/

      if(icall.eq.0)then
         icall = 1
         xdiv=1./float(nx)
c  sponge coefficient
         klo = 30
         do k=1,klo-1
            scoef(k) = 0.
         end do
         do k=klo,nz
            scoef(k) = 2. * (1. + tanh((zkm(k)-105.)/10.))*1.e-5
         end do

         ysum = 0.
         do j=1,ny
            ysum = ysum + cosfi(j)
         end do
         divx = 1./float(nx)
      end if

      do k=1,nz
         fnpole(k) = 0.
         fspole(k) = 0.
         do i=1,nx
            fnpole(k) = fnpole(k) + fin1(k,i,ny)*xdiv
            fspole(k) = fspole(k) + fin1(k,i,1)*xdiv
         end do
      end do

c  convert to potential temperature
      do j=1,ny
         do i=1,nx
            do k=1,nz
               ptemp(k,i,j) = (tn1(k,i,j)+tref(k))*thconv(k)
            end do
         end do
      end do

      do j=1,ny
         jp1=j+1
         if(j.eq.ny)jp1=ny
         jm1=j-1
         if(j.eq.1)jm1=1
         dxhm1 = 1./dx(j)
         dycos = 1./(cosfi(j)*dy)
         do i=1,nx
            ip1=i+1
            im1=i-1
            if(i.eq.1) im1=nx
            if(i.eq.nx)ip1=1

c  fields at midpoints
            do k=2,nz
               umid(k) = (un1(k,i,j)+un1(k-1,i,j))*.5
               vmid(k) = (vn1(k,i,j)+vn1(k-1,i,j))*.5
               tmid(k) = (ptemp(k,i,j)+ptemp(k-1,i,j))*.5
            end do
            umid(1) = (un1(1,i,j)+ubc(i,j))*.5
            vmid(1) = (vn1(1,i,j)+vbc(i,j))*.5
            tmid(1) = (ptemp(1,i,j)+thbc(i,j))*.5

c  tendency terms involving vertical derivatives and/or boundaries
c  thermal equation in potential temperature (converted to normal
c  temperature below for advection)
            do k=1,nz-1
               fu(k,i,j) = - (umid(k+1)*wn1(k+1,i,j)
     +                        - umid(k)*wn1(k,i,j))/dz/rou(k)
               fv(k,i,j) = - (vmid(k+1)*wn1(k+1,i,j)
     +                        - vmid(k)*wn1(k,i,j))/dz/rou(k)
               ft(k,i,j) = - (tmid(k+1)*wn1(k+1,i,j)
     +                        - tmid(k)*wn1(k,i,j))/dz/rou(k)
               vertadv(k,i,j) = ft(k,i,j) / thconv(k)
            end do

            fu(nz,i,j) = umid(nz)*wn1(nz,i,j)/dz/rou(nz)
            fv(nz,i,j) = vmid(nz)*wn1(nz,i,j)/dz/rou(nz)
            ft(nz,i,j) = tmid(nz)*wn1(nz,i,j)/dz/rou(nz)
            vertadv(nz,i,j) = ft(nz,i,j) / thconv(nz)

c  compute and add in all other tendency terms
            do k=1,nz
               upx = (un1(k,ip1,j) + un1(k,i,j))  *.5
               umx = (un1(k,i,j)   + un1(k,im1,j))*.5
               upy = (un1(k,i,jp1) + un1(k,i,j))  *.5
               umy = (un1(k,i,j)   + un1(k,i,jm1))*.5
               vpx = (vn1(k,ip1,j) + vn1(k,i,j))  *.5
               vmx = (vn1(k,i,j)   + vn1(k,im1,j))*.5
               vpy = (vn1(k,i,jp1) + vn1(k,i,j))  *.5
               vmy = (vn1(k,i,j)   + vn1(k,i,jm1))*.5
               tpx = (tn1(k,ip1,j) + tn1(k,i,j))  *.5
               tmx = (tn1(k,i,j)   + tn1(k,im1,j))*.5
               tpy = (tn1(k,i,jp1) + tn1(k,i,j))  *.5
               tmy = (tn1(k,i,j)   + tn1(k,i,jm1))*.5
               vpyc = vpy * cosf2(j+1)
               vmyc = vmy * cosf2(j)
c  zonal wind
               fu(k,i,j) = fu(k,i,j)
     +                   - (upx*upx - umx*umx)*dxhm1
     +                   - (upy*vpyc - umy*vmyc)*dycos
     +                   + vn1(k,i,j)*cor(j) 
     +                   + un1(k,i,j)*vn1(k,i,j)*tgfia(j)
     +                   - (fin1(k,ip1,j) - fin1(k,im1,j))*.5*dxhm1
     +                   + fxh(k,i,j)
c     +                   + cgwfx(k,i,j)
     +                   + falph(k,i,j)
c  meridional wind
               if(j.gt.1.and.j.lt.ny) deltfi = (fin1(k,i,j+1)
     +                                         -fin1(k,i,j-1))*.5
               if(j.eq.1) deltfi = (fin1(k,i,2)
     +                             +fin1(k,i,1))*.5 - fspole(k)
               if(j.eq.ny) deltfi = fnpole(k) - (fin1(k,i,ny-1)
     +                             +fin1(k,i,ny))*.5
               fv(k,i,j) = fv(k,i,j)
     +                   - (vpx*upx - vmx*umx)*dxhm1
     +                   - (vpy*vpyc - vmy*vmyc)*dycos
     +                   - un1(k,i,j)*cor(j) 
     +                   - un1(k,i,j)*un1(k,i,j)*tgfia(j)
     +                   - deltfi/dy
     +                   + fyh(k,i,j)
c     +                   + cgwfy(k,i,j)
c  temperature
               ft(k,i,j) = ft(k,i,j) / thconv(k)
     +                   - (tpx*upx - tmx*umx)*dxhm1
     +                   - tref(k) * (upx - umx)*dxhm1
     +                   - (tpy*vpyc - tmy*vmyc)*dycos
     +                   - tref(k) * (vpyc - vmyc)*dycos
     +                   + qir(k+nztrop,i,j)
     +                   + solht(k,i,j)
     +                   + ch_heat(k,i,j)
     +                   + q_resid(k,i,j)
            end do
         end do
      end do

c  Robert time-stepping
      do j=1,ny
         do i=1,nx
            do k=1,nz
c..............................................................
c  no sponge layer
c  with Lindzen gwd
c               un2(k,i,j) = (un0(k,i,j) + dtleap * (fu(k,i,j)
c     +                                  - fcgr(k,i,j)))
c     +                      /(1. - fgr(k,i,j) * dtleap)
c               vn2(k,i,j) = vn0(k,i,j) + dtleap * fv(k,i,j)
c               tn2(k,i,j) = tn0(k,i,j) + dtleap * ft(k,i,j)
c               fsponge(k,i,j) = 0.
c..............................................................
c  with sponge layer
c  no Lindzen gwd
c               un2(k,i,j) = (un0(k,i,j) + dtleap * (fu(k,i,j)
c     +                                  + scoef(k) * uclimo(k,j))) 
c     +                      /(1. + scoef(k) * dtleap)
c               vn2(k,i,j) = (vn0(k,i,j) + dtleap * (fv(k,i,j)))
c     +                      /(1. + scoef(k)*dtleap)
c               tn2(k,i,j) = tn0(k,i,j) + dtleap * ft(k,i,j)
c               fsponge(k,i,j) = -scoef(k) * (un2(k,i,j) - uclimo(k,j))
c..............................................................
c  with sponge layer
c  with Lindzen gwd
               un2(k,i,j) = (un0(k,i,j) + dtleap * (fu(k,i,j)
     +                                  + scoef(k) * uclimo(k,j)
     +                                  - fcgr(k,i,j)))
     +            /(1. + (scoef(k) - fgr(k,i,j)) * dtleap)
               vn2(k,i,j) = (vn0(k,i,j) + dtleap * (fv(k,i,j)))
     +                      /(1. + scoef(k)*dtleap)
               tn2(k,i,j) = tn0(k,i,j) + dtleap * ft(k,i,j)
               fsponge(k,i,j) = -scoef(k) * (un2(k,i,j) - uclimo(k,j))
c..............................................................

               un1(k,i,j) = un1(k,i,j) 
     +               + alphar * (un2(k,i,j)-2.*un1(k,i,j)+un0(k,i,j))
               vn1(k,i,j) = vn1(k,i,j) 
     +               + alphar * (vn2(k,i,j)-2.*vn1(k,i,j)+vn0(k,i,j))
               tn1(k,i,j) = tn1(k,i,j) 
     +               + alphar * (tn2(k,i,j)-2.*tn1(k,i,j)+tn0(k,i,j))
            end do
         end do
      end do

! temperature tendency due to vertical advection
      do k=1,nz
         dtadv(k) = 0.
         do j=1,ny
            sum = 0.
            do i=1,nx
               sum = sum + vertadv(k,i,j)
            end do
            dtadv(k) = dtadv(k) + sum*cosfi(j)
         end do
         dtadv(k) = dtadv(k)*divx / ysum
      end do

      return
      end
