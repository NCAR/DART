!----------------------------------------------------------------------
      module dynam 
!----------------------------------------------------------------------
             
      use params

      implicit none

      save

!... variables for program control, dynamical variables, and filtering

      integer :: ndiabat, ntrans, ninterp

      integer, dimension (nz) :: nifx,nify

      real :: deltat, dtleap, h, defi, pi, dy, dz, dzrh2, 
     +     treflb, trefub, dxrad, dyrad, dzz, gcy, gcz, ymin

      real, dimension (nz) :: rou, row, z, thconv, tref, ssref, 
     +                        zkm, pmb, tex, tey

c      real, dimension (ny) :: cor, tgfia, sinfi, cosfi, dx, tubc,
c     +                        uubc, tnzp1, unzp1, fi0, u0, t0, v0, 
c     +                        pwav1c, pwav1s, pwav2c, pwav2s,
c     +                        phideg
      real, dimension (ny) :: cor, tgfia, sinfi, cosfi, dx, tubc,
     +                        uubc, tnzp1, unzp1, phideg, gcx

      real, dimension (nyp1) :: cosf2

      real, dimension (nx,ny) :: fi_lbc, t_lbc, u_lbc, v_lbc

      real, dimension (nx,ny) :: tbc, fibc, ubc, vbc, thbc

      real, dimension (nz,nx,ny) :: un1, vn1, tn1, un0, vn0, tn0,
     +                              un2, vn2, tn2, fin1, wn1, ww

      real, dimension (nzp2,nx,ny) :: uusm, vvsm, wwsm

      real, dimension(ny) :: dau, dbu, dav, dbv, dat, dbt, daz, dbz,
     $                       sau, sbu, sav, sbv, sat, sbt, saz, sbz

c  i/o file names
      character*128 :: namf10,namf12,namf13,namf14,namf16,namf17,
     +                 namf18,namf19,namf20,namf21,namf22,namf23,
     +                 namf24,namf25,namf26,namf27,namf28,namf29,
     +                 namf30,namf40,namf41,namf52,namf53

      complex, dimension (nxhlf) :: wfft

      end module dynam 


