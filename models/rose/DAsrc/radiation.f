      subroutine radiation (doy, gmt_frac, tfull)
!------------------------------------------------------------------------
!  calculate radiative heating and cooling from progams supplied by X. Zhu
!------------------------------------------------------------------------

      use params
      use dynam
      use phys
      use chem, only : hnm

      implicit none

      integer, intent(in) :: doy
      real, intent(in) :: gmt_frac
      real, dimension(nz,nx,ny), intent(in) :: tfull

      real, dimension(nztrop) :: ptrop
      real, dimension(nzz,nx,ny) :: tzu, den_m3
      real :: aa
      integer icall, i, j, k, kk
      logical :: matrx = .true.

      save icall, ptrop
      data icall /0/

      icall = icall + 1
      if(icall.le.1) then
!... tropospheric pressure
         do k=1,nztrop
            ptrop(k) = 1013. * exp(-(k-1)*2.5/7.)
         end do
      end if

!-----------------------------------------------------------
! density & temperature on extended altitude grid nzz 
      do j=1,ny
         do k=1,nztrop
            aa = ptrop(k)/(tnmc(k,j)*1.38e-19) * 1.e6
            do i=1,nx
               tzu(k,i,j) = tnmc(k,j)
               den_m3(k,i,j) = aa
            end do
         end do
         do k=1,nz
            kk = k + nztrop
            do i=1,nx
               tzu(kk,i,j) = tfull(k,i,j)
               den_m3(kk,i,j) = hnm(k,i,j)*1.e6
            end do
         end do
      end do

      print *, 'cool_ir'
      call cool_ir (den_m3, tzu, co2mr, cl_o3zz, cl_h2o,
     $                   tcltop, qir, matrx)

      print *, 'heat_uv'        ! solar heating
      call heat_uv (tfull)

      print *, 'chemheat'       ! chemical heating
      call chemheat

      print *, 'hadjust_fix'    ! ensure zero net global mean heating
      call hadjust_fix

      end
