      subroutine chembc_intrp (dayofyr, gmt_frac)
c  interpolate chemical upper and lower boundary fields in time
c  NOTE: for now, b.c. on O2 taken from initial condidtions

      use params
      use chem
      use dynam

      use utilities_mod, only : open_file, close_file

      implicit none

      integer, intent(in) :: dayofyr
      real,    intent(in) :: gmt_frac

      integer, parameter :: nmonths = 12
      integer :: iunit
      integer :: icall, iday, m1, m2, mm, nfirst(nmonths), i, j
      real :: ddays(nmonths), f1, f2, smult, waterlbc(ny)

c  monthly means from SOCRATES
      real,dimension(ny,nmonths) :: ub_n2o, ub_ch4, ub_h2o, ub_hno3,
     $                              ub_n2o5, ub_co, ub_hcl, ub_clono2,
     $                              ub_hocl, ub_h2o2, ub_ho2no2,ub_ch2o,
     $                              ub_o1d, ub_oh, ub_cl, ub_o, ub_o3, 
     $                              ub_ho2, ub_no2, ub_no, ub_n, ub_clo,
     $                              ub_no3, ub_h
      real,dimension(ny,nmonths) :: lb_n2o, lb_ch4, lb_h2o, lb_hno3,
     $                              lb_n2o5, lb_co, lb_hcl, lb_clono2,
     $                              lb_hocl, lb_h2o2, lb_ho2no2,lb_ch2o,
     $                              lb_o1d, lb_oh, lb_cl, lb_o, lb_o3, 
     $                              lb_ho2, lb_no2, lb_no, lb_n, lb_clo,
     $                              lb_no3, lb_h
c   interpolated daily values
      real, dimension (ny) :: ubc_n2o, ubc_ch4, ubc_h2o, ubc_hno3,
     $                        ubc_n2o5, ubc_co, ubc_hcl, ubc_clono2,
     $                        ubc_hocl, ubc_h2o2, ubc_ho2no2, ubc_ch2o, 
     $                        ubc_o1d, ubc_oh, ubc_cl, ubc_o, ubc_o3, 
     $                        ubc_ho2, ubc_no2, ubc_no, ubc_n, ubc_clo,
     $                        ubc_no3, ubc_h
      real, dimension (ny) :: lbc_n2o, lbc_ch4, lbc_h2o, lbc_hno3,
     $                        lbc_n2o5, lbc_co, lbc_hcl, lbc_clono2,
     $                        lbc_hocl, lbc_h2o2, lbc_ho2no2, lbc_ch2o, 
     $                        lbc_o1d, lbc_oh, lbc_cl, lbc_o, lbc_o3, 
     $                        lbc_ho2, lbc_no2, lbc_no, lbc_n, lbc_clo,
     $                        lbc_no3, lbc_h

      save ub_n2o, ub_ch4, ub_h2o, ub_hno3, ub_n2o5, ub_co, ub_hcl, 
     $     ub_clono2, ub_hocl, ub_h2o2, ub_ho2no2, ub_ch2o, ub_o1d,
     $     ub_oh, ub_cl, ub_o, ub_o3, ub_ho2, ub_no2, ub_no, ub_n, 
     $     ub_clo, ub_no3, ub_h, lb_n2o, lb_ch4, lb_h2o, lb_hno3,
     $     lb_n2o5, lb_co, lb_hcl, lb_clono2, lb_hocl, lb_h2o2, 
     $     lb_ho2no2, lb_ch2o, lb_o1d, lb_oh, lb_cl, lb_o, lb_o3, 
     $     lb_ho2, lb_no2, lb_no, lb_n, lb_clo, lb_no3, lb_h

      data nfirst/15,46,74,105,135,166,196,226,258,288,319,349/
      data ddays/31.,28.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31./
      data icall/0/

c-------------------------------------------------------------------
c  read all SOCRATES input files once and save for interpolation
c-------------------------------------------------------------------
      icall = icall+1

      if(icall.le.1)then

c  read zonal mean boundary data
c  lower chemical

         iunit = open_file(namf27,form='formatted',action='read')

c        open (unit = 27,
c    $         file = namf27,
c    $         form = 'formatted',
c    $         status = 'old')

         read(iunit,'(1x,6e13.5)') lb_n2o, lb_ch4, lb_h2o, lb_hno3,
     $                             lb_n2o5, lb_co, lb_hcl, lb_clono2,
     $                             lb_hocl, lb_h2o2, lb_ho2no2, lb_ch2o,
     $                             lb_o1d, lb_oh, lb_cl, lb_o, lb_o3, 
     $                             lb_ho2, lb_no2, lb_no, lb_n, lb_clo,
     $                             lb_no3, lb_h

         call close_file(iunit)

c  upper chemical

         iunit = open_file(namf26,form='formatted',action='read')

c        open (unit = 26,
c    $         file = namf26,
c    $         form = 'formatted',
c    $         status = 'old')

         read(iunit,'(1x,6e13.5)') ub_n2o, ub_ch4, ub_h2o, ub_hno3,
     $                             ub_n2o5, ub_co, ub_hcl, ub_clono2,
     $                             ub_hocl, ub_h2o2, ub_ho2no2, ub_ch2o,
     $                             ub_o1d, ub_oh, ub_cl, ub_o, ub_o3, 
     $                             ub_ho2, ub_no2, ub_no, ub_n, ub_clo,
     $                             ub_no3, ub_h

         call close_file(iunit)

c  water bc for dry tropopause
         do j=1,ny
cc            if(abs(phideg(j)).le.30.) waterlbc(j) = 1.4 e-6
cc            if(abs(phideg(j)).ge.60.) waterlbc(j) = 3. e-6
cc            if((abs(phideg(j)).gt.30.).and.(abs(phideg(j)).lt.60.)) 
cc     +         waterlbc(j) = 1.4e-6 + 1.6e-6 
cc     +            * (1. - (cos((abs(phideg(j))-30.)/30.*pi/2.))**2)
            waterlbc(j) = 2.2 e-6
         end do

      end if
c______________________________________________________________________
c______________________________________________________________________

c  Once per day, interpolate SOCRATES data to day of year

c______________________________________________________________________
c______________________________________________________________________

      if(icall.le.1.or.gmt_frac.eq.0.)then
c  interpolation parameters
         iday = dayofyr
         if(iday.le.nfirst(1))then
            m1 = nmonths
            m2 = 1
            f1 = float(15-iday)/ddays(1)
         end if
         if(iday.gt.nfirst(nmonths))then
            m1 = nmonths
            m2 = 1
            f1 = float(380-iday)/ddays(nmonths)
         end if
         if(iday.gt.nfirst(1).and.iday.le.nfirst(nmonths))then
            do mm=1,11
               if(iday.gt.nfirst(mm).and.iday.le.nfirst(mm+1))then
                  m1 = mm
                  m2 = mm+1
                  f1 = float(nfirst(m2)-iday)/ddays(m1)
               end if
            end do
         end if
         f2 = 1.-f1

         do j=1,ny
            ubc_n2o(j)    = f1*ub_n2o(j,m1)    + f2*ub_n2o(j,m2)
            ubc_ch4(j)    = f1*ub_ch4(j,m1)    + f2*ub_ch4(j,m2)
            ubc_h2o(j)    = f1*ub_h2o(j,m1)    + f2*ub_h2o(j,m2)
            ubc_hno3(j)   = f1*ub_hno3(j,m1)   + f2*ub_hno3(j,m2)
            ubc_n2o5(j)   = f1*ub_n2o5(j,m1)   + f2*ub_n2o5(j,m2)
            ubc_co(j)     = f1*ub_co(j,m1)     + f2*ub_co(j,m2)
            ubc_hcl(j)    = f1*ub_hcl(j,m1)    + f2*ub_hcl(j,m2)
            ubc_clono2(j) = f1*ub_clono2(j,m1) + f2*ub_clono2(j,m2)
            ubc_hocl(j)   = f1*ub_hocl(j,m1)   + f2*ub_hocl(j,m2)
            ubc_h2o2(j)   = f1*ub_h2o2(j,m1)   + f2*ub_h2o2(j,m2)
            ubc_ho2no2(j) = f1*ub_ho2no2(j,m1) + f2*ub_ho2no2(j,m2)
            ubc_ch2o(j)   = f1*ub_ch2o(j,m1)   + f2*ub_ch2o(j,m2)
            ubc_o1d(j)    = f1*ub_o1d(j,m1)    + f2*ub_o1d(j,m2)
            ubc_oh(j)     = f1*ub_oh(j,m1)     + f2*ub_oh(j,m2)
            ubc_cl(j)     = f1*ub_cl(j,m1)     + f2*ub_cl(j,m2)
            ubc_o(j)      = f1*ub_o(j,m1)      + f2*ub_o(j,m2)
            ubc_o3(j)     = f1*ub_o3(j,m1)     + f2*ub_o3(j,m2)
            ubc_ho2(j)    = f1*ub_ho2(j,m1)    + f2*ub_ho2(j,m2)
            ubc_no2(j)    = f1*ub_no2(j,m1)    + f2*ub_no2(j,m2)
            ubc_no(j)     = 7.e-6
            ubc_n(j)      = f1*ub_n(j,m1)      + f2*ub_n(j,m2)
            ubc_clo(j)    = f1*ub_clo(j,m1)    + f2*ub_clo(j,m2)
            ubc_no3(j)    = f1*ub_no3(j,m1)    + f2*ub_no3(j,m2)
            ubc_h(j)      = f1*ub_h(j,m1)      + f2*ub_h(j,m2)

            lbc_n2o(j)    = f1*lb_n2o(j,m1)    + f2*lb_n2o(j,m2)
            lbc_ch4(j)    = f1*lb_ch4(j,m1)    + f2*lb_ch4(j,m2)
            lbc_h2o(j)    = waterlbc(j)
            lbc_hno3(j)   = f1*lb_hno3(j,m1)   + f2*lb_hno3(j,m2)
            lbc_n2o5(j)   = f1*lb_n2o5(j,m1)   + f2*lb_n2o5(j,m2)
            lbc_co(j)     = f1*lb_co(j,m1)     + f2*lb_co(j,m2)
            lbc_hcl(j)    = f1*lb_hcl(j,m1)    + f2*lb_hcl(j,m2)
            lbc_clono2(j) = f1*lb_clono2(j,m1) + f2*lb_clono2(j,m2)
            lbc_hocl(j)   = f1*lb_hocl(j,m1)   + f2*lb_hocl(j,m2)
            lbc_h2o2(j)   = f1*lb_h2o2(j,m1)   + f2*lb_h2o2(j,m2)
            lbc_ho2no2(j) = f1*lb_ho2no2(j,m1) + f2*lb_ho2no2(j,m2)
            lbc_ch2o(j)   = f1*lb_ch2o(j,m1)   + f2*lb_ch2o(j,m2)
            lbc_o1d(j)    = f1*lb_o1d(j,m1)    + f2*lb_o1d(j,m2)
            lbc_oh(j)     = f1*lb_oh(j,m1)     + f2*lb_oh(j,m2)
            lbc_cl(j)     = f1*lb_cl(j,m1)     + f2*lb_cl(j,m2)
            lbc_o(j)      = f1*lb_o(j,m1)      + f2*lb_o(j,m2)
            lbc_o3(j)     = f1*lb_o3(j,m1)     + f2*lb_o3(j,m2)
            lbc_ho2(j)    = f1*lb_ho2(j,m1)    + f2*lb_ho2(j,m2)
            lbc_no2(j)    = f1*lb_no2(j,m1)    + f2*lb_no2(j,m2)
            lbc_no(j)     = f1*lb_no(j,m1)     + f2*lb_no(j,m2)
            lbc_n(j)      = f1*lb_n(j,m1)      + f2*lb_n(j,m2)
            lbc_clo(j)    = f1*lb_clo(j,m1)    + f2*lb_clo(j,m2)
            lbc_no3(j)    = f1*lb_no3(j,m1)    + f2*lb_no3(j,m2)
            lbc_h(j)      = f1*lb_h(j,m1)      + f2*lb_h(j,m2)
         end do
      end if

c______________________________________________________________________
c______________________________________________________________________
c
c______________________________________________________________________
c 
c  Sort into model fields (arrays qubc, qlbc)
c______________________________________________________________________
c   Impose zenith angle control for several short-lived species
c______________________________________________________________________

      do j=1,ny
         do i=1,nx
            if(sc2d(i,j).gt.50.) then
               smult = 0.
            else
               smult = 2.
            end if
            qubc(i,j,1) = ubc_n2o(j)
            qubc(i,j,2) = ubc_ch4(j)
            qubc(i,j,3) = ubc_h2o(j)
            qubc(i,j,4) = 0.0 
            qubc(i,j,5) = ubc_hno3(j)
            qubc(i,j,6) = ubc_n2o5(j)
            qubc(i,j,7) = ubc_h(j)
            qubc(i,j,8) = ubc_oh(j)
            qubc(i,j,9) = ubc_co(j)
            qubc(i,j,10) = ubc_hcl(j)
            qubc(i,j,11) = ubc_clono2(j)
            qubc(i,j,12) = ubc_hocl(j)
            qubc(i,j,13) = ubc_h2o2(j)
            qubc(i,j,14) = ubc_ho2(j)
            qubc(i,j,15) = ubc_ho2no2(j)
            qubc(i,j,16) = 5.e-7
            qubc(i,j,17) = ubc_ch2o(j)
            qubc(i,j,18) = ubc_o(j)
            qubc(i,j,19) = ubc_o3(j)
            qubc(i,j,20) = ubc_cl(j)
            qubc(i,j,21) = ubc_clo(j)
            qubc(i,j,22) = ubc_n(j)
            qubc(i,j,23) = ubc_no(j)
            qubc(i,j,24) = ubc_no2(j)
            qubc(i,j,25) = ubc_no3(j)
            qubc(i,j,26) = q_o2(nz,j) 

            qlbc(i,j,1) = lbc_n2o(j)
            qlbc(i,j,2) = lbc_ch4(j)
            qlbc(i,j,3) = lbc_h2o(j)
            qlbc(i,j,4) = 0.0 
            qlbc(i,j,5) = lbc_hno3(j)
            qlbc(i,j,6) = lbc_n2o5(j)
            qlbc(i,j,7) = lbc_h(j) * smult
            qlbc(i,j,8) = lbc_oh(j) * smult
            qlbc(i,j,9) = lbc_co(j)
            qlbc(i,j,10) = lbc_hcl(j)
            qlbc(i,j,11) = lbc_clono2(j)
            qlbc(i,j,12) = lbc_hocl(j)
            qlbc(i,j,13) = lbc_h2o2(j)
            qlbc(i,j,14) = lbc_ho2(j)  * smult
            qlbc(i,j,15) = lbc_ho2no2(j)
            qlbc(i,j,16) = 5.e-7
            qlbc(i,j,17) = lbc_ch2o(j)
            qlbc(i,j,18) = lbc_o(j) * smult
            qlbc(i,j,19) = lbc_o3(j)
            qlbc(i,j,20) = lbc_cl(j) * smult
            qlbc(i,j,21) = lbc_clo(j)
            qlbc(i,j,22) = lbc_n(j) * smult
            qlbc(i,j,23) = lbc_no(j) * smult
            qlbc(i,j,24) = lbc_no2(j)
            qlbc(i,j,25) = lbc_no3(j)
            qlbc(i,j,26) = q_o2(1,j) 
         end do
      end do
c______________________________________________________________________
c______________________________________________________________________

      end
