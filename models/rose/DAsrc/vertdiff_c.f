      subroutine vertdiff_c
c-----------------------------------------------------------------------
c
c driver routine to compute vertical diffusion of chemicals
c upper and lower boundary information included in diffusion calculation
c Effective vertical velocity due to molecular diffusion is included
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      use params
      use chem
      use dynam
      use phys

      implicit none

c-----------------------------------------------------------------------
c
c
      integer :: i, j, k, nc
      real dzztint(nz),         ! molecular diffusion at midpoints
     $     dzzgint(nz),         ! g.w. diffusion at midpoints
     $     qq(nz),              ! temporary chem array
     $     qp1(nz),             ! chem after vertical diffusion
     $     qupperbc,            ! q upper boundary
     $     qlowerbc             ! q lower boundary
c
c local workspace 
c
      real ca(nz),              ! -upper diagonal
     $     cb(nz),              ! -diagonal
     $     cc(nz),              ! -lower diagonal
     $     xkvt(nzp1),          ! xkv including molecular and g.w.
     $     tmp2,                ! temporary storage
     $     prndc,               ! inverse Prandtl number for chemistry
     $     factr_d,
     $     factr_a,
     $     df(0:nz),
     $     adv(0:nz),
     $     rhoratio,
     $     rhoratio2

      integer, dimension(nbcon) :: index
      data index /1, 1, 1, 0,    ! n2o    ch4     h2o      o1d
     $            1, 1, 1, 1,    ! hno3   n2o5    h        oh
     $            1, 1, 1, 1,    ! co     hcl     clono2   hocl
     $            1, 1, 1, 1,    ! h2o2   ho2     ho2no2   h2
     $            1, 1, 1, 1,    ! ch2o   o       o3       cl
     $            1, 1, 1, 1,    ! clo    n       no       no2
     $            1, 0  /        ! no3    o2
c
c-----------------------------------------------------------------------
      prndc = .3       ! inverse Prandtl number for chemistry

      factr_d = deltat / (dz**2)
      factr_a = deltat / (2.*dz)
      rhoratio = exp(-dz/2./h)
      rhoratio2 = exp(-dz/h)

c set the vertical diffusion coefficient
c Includes both Hines and ccm gw terms and molecular diffusion 
      do nc=1,nbcon

         if(index(nc).ne.0)then

            do j=1,ny
               do i=1,nx
                  xkvt(1) = zkmin(1) 
     $                    + prndc * fdzzh(1,i,j)
                  xkvt(nzp1) = zkmin(nz) 
     $                       + prndc * fdzzh(nz,i,j) 
     $                       + dmolch(nz,nc)
                  do k=1,nz-1
                     xkvt(k+1) = zkmin(k+1) 
     $                         + prndc * fdzzh(k+1,i,j) 
     $                         + dmolch(k+1,nc)
                  end do

c indexing for terms:
c    for diffusive terms k -> 1/2 gridpoint above grid of same index
c    for advective terms k -> normal grid position

                  df(0)  = xkvt(1) * factr_d / rhoratio
                  adv(0) = 0.
                  do k=1,nz-1
                     df(k)  = xkvt(k+1) * factr_d
                     adv(k) = wdif(k,nc)  * factr_a
                  end do
                  df(nz)  = df(nz-1)
                  adv(nz) = wdif(nz,nc) * factr_a

c determine superdiagonal (ca), diagonal (bb) and subdiagonal (cc)
c coefficients of the tridiagonal diffusion matrix. 
                  cc(1) = df(0)
                  do k=1,nz-1
                     ca(k) = df(k)*rhoratio - adv(k+1)*rhoratio2
                     cb(k) = 1. + df(k) * rhoratio
     $                          + df(k-1) / rhoratio
                     cc(k) = df(k-1)/rhoratio + adv(k-1)/rhoratio2 
                  end do
                  ca(nz) = df(nz)*rhoratio - adv(nz)*rhoratio2
                  cb(nz) = 1. + df(nz) * rhoratio 
     $                        + df(nz-1)/rhoratio 
                  cc(nz) = df(nz-1)/rhoratio + adv(nz-1)/rhoratio2 
c
c set chemical fields to pass to svdiffbc
c
                  do k=1,nz
                     qq(k) = qn1(k,i,j,nc)
                  end do
                  qupperbc = qubc(i,j,nc)
                  qlowerbc = qlbc(i,j,nc)
c
c diffuse tracers
c
                  call svdiffch(qq,qupperbc,qlowerbc,ca,cb,cc,qp1)

c  back to original grids
                  do k=1,nz
                     qn1(k,i,j,nc) = qp1(k)
                  end do
               end do
            end do
         end if
      end do
c
c done
c
      return
      end
