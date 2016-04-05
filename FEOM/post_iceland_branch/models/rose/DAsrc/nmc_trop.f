c-------------------------------------------------------------------
      subroutine nmc_lbc (dayofyr, yearno)
c-------------------------------------------------------------------
c  interpolate dynamical tropospheric and lower boundary data in time
c  7-day averages from Sep 25, 2001  to Dec 30, 2002
c  yearno is 0 for 1999, etc.
c-------------------------------------------------------------------

      use params
      use dynam
      use phys
      use chem
      use trop

      implicit none

      integer, parameter :: ntimes=52
      integer, intent(in) :: dayofyr, yearno

      integer :: i, j, k, kk, m, mm, m1, m2, ispan, mid, nn
      integer :: iday, icall, nbar
      integer :: kc, ks, icount, mdref
      integer :: nmc_doy(ntimes), nmc_yr(ntimes), nmc_m1, nmc_m2
      real :: g, f1, f2, phibar
      real, dimension(nztrop,nx,ny,ntimes) :: gphin, tin, uin, vin
      integer :: inp_dates(2,ntimes)
      real, dimension(ny) :: u_zonal, v_zonal, t_zonal, p_zonal 
      real, parameter :: gcp = 9.8e-03

      real :: uu_zm, vv_zm, tt_zm, zz_zm, f_reduce

      save icall, nmc_doy, nmc_yr

      data icall/0/, g/9.8/

c-------------------------------------------------------------------
c  read all input files once and save for interpolation
c-------------------------------------------------------------------
      icall = icall+1
      if(icall.le.1)then

         call trop_open( namf20 )
         print *, 'trop_open complete ', namf20

         call trop_retrieve('phi', inp_dates, gphin )
         print *, 'gph read'
         call trop_retrieve('u  ', inp_dates, uin )
         print *, 'u read'
         call trop_retrieve('v  ', inp_dates, vin )
         print *, 'v read'
         call trop_retrieve('t  ', inp_dates, tin )
         print *, 'T read'

         call trop_close()
         do m=1,ntimes
            nmc_doy(m) = inp_dates(2,m)
            nmc_yr(m)  = inp_dates(1,m)
            print *, nmc_doy(m), nmc_yr(m)
         end do

c reference (global annual mean) temperature
         nbar = nx * ny * ntimes
         do k=1,nztrop
            troptref(k) = 0.
            do m=1,ntimes
               do j=1,ny
                  do i=1,nx
                     troptref(k) = troptref(k) + tin(k,i,j,m)
                  end do
               end do
            end do
            troptref(k) = troptref(k) / float(nbar)
         end do
c referrence static stability
         do k=2,nztrop-1
            tropsref(k) = (troptref(k+1) - troptref(k-1))
     $                    / (2. * 7.e3) + gcp
         end do
         tropsref(1) = tropsref(2)
         tropsref(nztrop) = tropsref(nztrop-1)
c referrence geopotential
         phibar = 0.
         do m=1,ntimes
            do j=1,ny
               do i=1,nx
                  phibar = phibar + gphin(nztrop,i,j,m)
               end do
            end do
         end do
         phibar = phibar / float(nbar)
      end if

c______________________________________________________________________
c______________________________________________________________________

c  Now interpolate all data to day of year
c  Also, convert temperature to perturbation temperature

c______________________________________________________________________
c______________________________________________________________________

      ispan = 7

c  interpolation parameters
      iday = dayofyr + 365 * yearno
      icount = iday - nmc_doy(1)
      if(mod(icount,ispan).eq.0) then
         m1 = icount/ispan + 1
         m2 = m1
         f1 = 1.
      else
         m1 = icount/ispan + 1
         mdref = ispan*(m1-1)
         m2 = m1 + 1
         f1 = 1. - (icount-mdref)/float(ispan)
         if(m1.lt.1.or.m2.lt.1) then
            m1 = 1
            m2 = 1
            f1 = 1.
         end if
         if(m2.gt.ntimes.or.m2.gt.ntimes) then
            m1 = ntimes
            m2 = ntimes
            f1 = 0.
         end if
      end if

      nmc_m1 = nmc_doy(m1) + 365 * (nmc_yr(m1)-2002)
      nmc_m2 = nmc_doy(m2) + 365 * (nmc_yr(m2)-2002)

      if((iday.lt.nmc_m1).or.(iday.gt.nmc_m2)) then
         print *, '*** ERROR in subr. nmc_lbc ***'
         print *, 'm1, m2, iday, dayofyr, yearno', 
     $             m1, m2, iday, dayofyr, yearno 
         print *, 'nmc_m1, nmc_doy(m1), nmc_yr(m1)', 
     $             nmc_m1, nmc_doy(m1), nmc_yr(m1)
         print *, 'nmc_m2, nmc_doy(m2), nmc_yr(m2)', 
     $             nmc_m2, nmc_doy(m2), nmc_yr(m2)
         stop
      end if

      f2 = 1.-f1

c______________________________________________________________________
c______________________________________________________________________

c  fields returned to model

c______________________________________________________________________
c  dynamical lower boundary

!  perturbation T
!  perturbation phi
      do j=1,ny
         do i=1,nx
            u_lbc(i,j) = f1*uin(nztrop,i,j,m1) + f2*uin(nztrop,i,j,m2)
            v_lbc(i,j) = f1*vin(nztrop,i,j,m1) + f2*vin(nztrop,i,j,m2)
            t_lbc(i,j) = f1*tin(nztrop,i,j,m1) + f2*tin(nztrop,i,j,m2)
     $                 - treflb 
            fi_lbc(i,j) =  f1*gphin(nztrop,i,j,m1)
     $                    +f2*gphin(nztrop,i,j,m2) - phibar
         end do
      end do


!______________________________________________________________________
!  tropospheric fields - zonally averaged
      do j=1,ny
         do k=1,nztrop
            tnmc(k,j) = 0.
            unmc(k,j) = 0.
            do i=1,nx
!  full T for radiation
               tnmc(k,j) = tnmc(k,j) + (f1*tin(k,i,j,m1) 
     $                               +  f2*tin(k,i,j,m2))/float(nx)
               unmc(k,j) = unmc(k,j) + (f1*uin(k,i,j,m1)  
     $                               +  f2*uin(k,i,j,m2))/float(nx)
            end do
!  perturbation T
            tpnmc(k,j) = tnmc(k,j) - troptref(k)
         end do
      end do

      end
