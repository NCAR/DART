      subroutine svdiffbc(qm1,qupper,qlower,ca,cc,qp1)
c-----------------------------------------------------------------------
c     ! solve vertical diffusion equation
c     ! procedure for solution of the implicit equation follows
c     ! Richtmyer and Morton (1967, pp 198-199).
c     ! This version includes upper and lower boundary information (assumes
c     !  that diffusion rate is constant beyond boundary).
c-----------------------------------------------------------------------
c
c     !   basic grid point resolution parameters
c
      use params

      implicit none

c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c
c     ! input variables
c
      real qm1(nx,ny,nz),   ! initial field
     $     qupper(nx,ny),   ! field at upper boundary
     $     qlower(nx,ny),   ! field at lower boundary
     $     ca(nx,ny,nz),    ! -upper diagonal coeff. of tri-diagonal matrix
     $     cc(nx,ny,nz)     ! -lower diagonal coeff. of tri-diagonal matrix
c
c     ! output variables
c
      real qp1(nx,ny,nz)    ! final field
c
c     ! local workspace
c
      integer :: i, j, k
      real :: tmp0, e0, f0, tmp1 
      real ze(nx,ny,nz),   ! terms appearing in soln of tridiagonal system
     $    zfq(nx,ny,nz),   ! terms appearing in soln of tridiagonal system
     $    tmp              ! temporary workspace
c
c     ! calculate e(k) and fq(k).
c     ! these are terms required in solution of tridiagonal matrix
c     ! defined by implicit diffusion eqn.
c
      do j=1,ny
         do i=1,nx
            tmp0 = 1. + ca(i,j,1)
            e0 = ca(i,j,1)/tmp0
            f0 = qlower(i,j)/tmp0
            tmp1 = 1. + ca(i,j,1) + cc(i,j,1) - cc(i,j,1) * e0
            ze(i,j,1) = ca(i,j,1)/tmp1
            zfq(i,j,1) = (qm1(i,j,1) + cc(i,j,1)*f0) / tmp1
         end do
      end do

      do k=2,nz
         do j=1,ny
            do i=1,nx
               tmp = 1. + ca(i,j,k) + cc(i,j,k) - cc(i,j,k)*ze(i,j,k-1)
               ze(i,j,k) = ca(i,j,k)/tmp
               zfq(i,j,k) = (qm1(i,j,k) + cc(i,j,k)*zfq(i,j,k-1))/tmp
            end do
         end do
      end do
c
c     ! perform back substitution
c
      do j=1,ny
         do i=1,nx
            qp1(i,j,nz) = zfq(i,j,nz) + ze(i,j,nz)*qupper(i,j)
         end do
      end do
c
      do k=nz-1,1,-1
         do j=1,ny
            do i=1,nx
               qp1(i,j,k) = zfq(i,j,k) + ze(i,j,k)*qp1(i,j,k+1)
            end do
         end do
      end do

      return
      end

