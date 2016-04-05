      subroutine svdiffch (qm1,qupper,qlower,ca,cb,cc,qp1)
c-----------------------------------------------------------------------
c     ! solve vertical diffusion equation
c     ! procedure for solution of the implicit equation follows
c     ! Richtmyer and Morton (1967, pp 198-199).
c     ! This version includes effective vertical velocity from 
c     !  molecular diffusion.
c     !  Flux upper boundary condition
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
      real qm1(nz),   ! initial field
     $     qupper,    ! field at upper boundary
     $     qlower,    ! field at lower boundary
     $     ca(nz),    ! -upper diagonal coeff. of tri-diagonal matrix
     $     cb(nz),    ! -diagonal coeff. of tri-diagonal matrix
     $     cc(nz)     ! -lower diagonal coeff. of tri-diagonal matrix
c
c     ! output variables
c
      real qp1(nz)    ! final field
c
c     ! local workspace
c
      integer :: i, j, k
      real :: tmp0, e0, f0, tmp1
      real ze(nz),   ! terms appearing in soln of tridiagonal system
     $    zfq(nz),   ! terms appearing in soln of tridiagonal system
     $    tmp        ! temporary workspace
c
c     ! calculate e(k) and fq(k).
c     ! these are terms required in solution of tridiagonal matrix
c     ! defined by implicit diffusion eqn.
c
            tmp0 = 1. + ca(1)
            e0 = ca(1)/tmp0
            f0 = qlower/tmp0
            tmp1 = cb(1) - cc(1) * e0
            ze(1) = ca(1)/tmp1
            zfq(1) = (qm1(1) + cc(1)*f0) / tmp1

            do k=2,nz
               tmp = cb(k) - cc(k)*ze(k-1)
               ze(k) = ca(k)/tmp
               zfq(k) = (qm1(k) + cc(k)*zfq(k-1))/tmp
            end do
c
c     ! perform back substitution
c
            qp1(nz) = zfq(nz) + ze(nz)*qupper
c
            do k=nz-1,1,-1
               qp1(k) = zfq(k) + ze(k)*qp1(k+1)
            end do

      return
      end

