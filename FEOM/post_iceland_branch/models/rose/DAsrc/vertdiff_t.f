      subroutine vertdiff_t
c-----------------------------------------------------------------------
c
c driver routine to compute vertical diffusion of momentum
c and potential temperature.
c
c free atmosphere diffusivities are computed first; then
c passed to parameterization svdiff
c
c upper and lower boundary information included in diffusion calculation
c
c global mean tendency due to diffusion saved and used in hadjust
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c     !   Basic grid point resolution parameters
c
      use params
      use dynam
      use phys

      implicit none

c-----------------------------------------------------------------------
c
c
      integer :: i, j, k
      real up1(nx,ny,nz),       ! u-wind after vertical diffusion
     $     vp1(nx,ny,nz),       ! v-wind after vertical diffusion
     $     thp(nx,ny,nz),       ! potential temperature after vertical dif
     $     dzztint(nz),         ! molecular diffusion at midpoints
     $     dzzgint(nx,ny,nz),   ! g.w. diffusion at midpoints
     $     tupperbc(nx,ny),     ! T upper boundary
     $     uupperbc(nx,ny),     ! u upper boundary
     $     vupperbc(nx,ny)      ! v upper boundary
c
c local workspace 
c
      real thm(nx,ny,nz),       ! potential temperature
     $     umm(nx,ny,nz),       ! zonal wind
     $     vmm(nx,ny,nz)        ! meridional wind
      real cah(nx,ny,nz),       ! -upper diagonal for heat
     $     cam(nx,ny,nz),       ! -upper diagonal for momentum
     $     cch(nx,ny,nz),       ! -lower diagonal for heat
     $     ccm(nx,ny,nz),       ! -lower diagonal for momentum
     $     dvdz2(nx,ny),        ! (du/dz)**2 + (dv/dz)**2
     $     fstab(nx,ny),        ! stable f(ri)
     $     funst(nx,ny),        ! unstable f(ri)
     $     xkvg(nx,ny,nzp1),    ! free atmosphere xkv including g.w.
     $     xkvt(nx,ny,nzp1),    ! xkv including molecular and g.w.
     $     rinub(nx,ny),        ! richardson number = (g/theta)(dtheta/dz)/
c                               !                     (du/dz**2+dv/dz**2)
     $     sstab(nx,ny),        ! static stability = g/th  * dth/dz
     $     tmp2,                ! temporary storage
     $     prndd,               ! inverse Prandtl number for temperature
     $     factr,
     $     rhoratio

c
      real xkvn                 ! neutral Xkv
      real tnew(nz)             ! global mean T after diffusion
      real ysum                 ! global mean area corrector

      real divx, sum
c
c-------------------------------------------------------------------------

      ysum = 0.
      do j=1,ny
         ysum = ysum + cosfi(j)
      end do
      divx = 1./float(nx)

      prndd = .3       ! inverse Prandtl number for dynamics

c set potential temperatures and determine diffusion at midpoints
c
      do j=1,ny
         do i=1,nx
            do k=1,nz
               thm(i,j,k) = (tn2(k,i,j)+tref(k)) * thconv(k)
               umm(i,j,k) = un2(k,i,j)
               vmm(i,j,k) = vn2(k,i,j)
            end do
            do k=2,nz
               dzzgint(i,j,k) = fdzzh(k,i,j)
            end do
            tupperbc(i,j) = tubc(j)
            uupperbc(i,j) = uubc(j)
            vupperbc(i,j) = 0.
         end do
      end do

      do k=2,nz
         dzztint(k) = dtzz(k)
      end do
c
c set the vertical diffusion coefficient above/below the diffusion leve
c
      do j=1,ny
         do i=1,nx
            xkvg(i,j,1) = 0.0
            xkvt(i,j,1) = 0.0
            xkvg(i,j,nzp1) = 0.0
            xkvt(i,j,nzp1) = 0.0
         end do
      end do

c compute the free atmosphere vertical diffusion coefficients
c xkvh = xkvq = xkvm.
c
      do k=1,nz-1
         do j=1,ny
            do i=1,nx
c
c vertical shear squared, min value of (delta v)**2 prevents
c zero shear.
c
               dvdz2(i,j) = (umm(i,j,k+1)-umm(i,j,k))**2 +
     $                      (vmm(i,j,k+1)-vmm(i,j,k))**2
               dvdz2(i,j) = amax1(dvdz2(i,j),1.e-10)
               dvdz2(i,j) = dvdz2(i,j)/dz**2
c
c static stability (use virtual potential temperature)
c
               sstab(i,j) = gravit*2.0/(thm(i,j,k) + thm(i,j,k+1))
     $                      *(thm(i,j,k+1) - thm(i,j,k))/dz
c
c richardson number, stable and unstable modifying functions
c
               rinub(i,j) = sstab(i,j)/dvdz2(i,j)
               fstab(i,j) = 1.0 / (1.0 + 10.0*rinub(i,j) * 
     $                      (1.0+8.0*rinub(i,j)))
               funst(i,j) = amax1(1. - 18.*rinub(i,j),0.)
c
c select the appropriate function of the richardson number
c
               if (rinub(i,j) .lt. 0) fstab(i,j) = sqrt(funst(i,j))
c
c neutral diffusion coefficient
c compute mixing length (z), where z is the interface height estimated
c with an 8 km scale height.
c
               xkvn      = xml2(k) * sqrt(dvdz2(i,j))
c
c full diffusion coefficient (modified by f(ri)),
c note index k for xkvg refers to 1/2 level above the standard grid
c
               xkvg(i,j,k+1) = amax1(zkmin(k+1),  xkvn*fstab(i,j))
               xkvt(i,j,k+1) = amax1(zkmin(k+1) 
     $                              + prndd * dzzgint(i,j,k+1)
     $                              + dzztint(k+1), xkvn*fstab(i,j))

            end do
         end do
      end do

c determine superdiagonal (ca(k)) and subdiagonal (cc(k)) coefficients
c of the tridiagonal diffusion matrix. the diagonal elements are a
c combination of ca and cc; they are not explicitly provided to the solv

      factr = deltat / (dz**2)
      rhoratio = exp(-dz/2./h)
      do k=1,nz-1
         do j=1,ny
            do i=1,nx
               cch(i,j,k+1) = xkvt(i,j,k+1) * factr / rhoratio
               cah(i,j,k  ) = xkvt(i,j,k+1) * factr * rhoratio
               ccm(i,j,k+1) = xkvg(i,j,k+1) * factr / rhoratio
               cam(i,j,k  ) = xkvg(i,j,k+1) * factr * rhoratio
            end do
         end do
      end do
c
c the last element of upper and lower diagonals
c
      do j=1,ny
         do i=1,nx
            cah(i,j,nz) = cah(i,j,nz-1)
            cam(i,j,nz) = cam(i,j,nz-1)
            cch(i,j,1) = cch(i,j,2)
            ccm(i,j,1) = ccm(i,j,2)
         end do
      end do

c
c diffuse momentum
c
      call svdiffbc(umm,uupperbc,ubc,cam,ccm,up1)
      call svdiffbc(vmm,vupperbc,vbc,cam,ccm,vp1)
c
c diffuse potential temperature
c
      call svdiffbc(thm,tupperbc,thbc,cah,cch,thp)

c  save input global mean temperature

      do k=1,nz
         dtglob(k) = 0.
         do j=1,ny
            sum = 0.
            do i=1,nx
               sum = sum + tn2(k,i,j)
            end do
            dtglob(k) = dtglob(k) + sum*cosfi(j)
         end do
         dtglob(k) = dtglob(k)*divx / ysum
      end do

c  move back to original grids and convert potential 
c  temperatures back to temperatures

      do j=1,ny
         do i=1,nx
            do k=1,nz
               un2(k,i,j) = up1(i,j,k)
               vn2(k,i,j) = vp1(i,j,k)
               tn2(k,i,j) = thp(i,j,k)/thconv(k) - tref(k)
            end do
         end do
      end do

c  calculate and save global mean temperature tendency

      do k=1,nz
         tnew(k) = 0.
         do j=1,ny
            sum = 0.
            do i=1,nx
               sum = sum + tn2(k,i,j)
            end do
            tnew(k) = tnew(k) + sum*cosfi(j)
         end do
         dtglob(k) = (tnew(k)*divx / ysum - dtglob(k)) / deltat
      end do
c
c done
c
      return
      end
