c----------------------------------------------------------------------
      subroutine chem_gs( gs_convrg, qnew )
c----------------------------------------------------------------------
c
c     Semi-implicit solver using Gauss-Seidel iteration
c
c----------------------------------------------------------------------
c                                                                     c
c     A semi-implicit solver is used to calculate concentrations      c
c     of long-lived species                                           c
c                                                                     c
c      c(n+1) - c(n) = dt * (P - L*c(n+1))                            c
c                                                                     c
c     Gauss-Seidel interation improves the accuracy of the solution   c
c                                                                     c
c----------------------------------------------------------------------
c
c    Species solved for in this routine:
c
c    N2O, CH4, H2O, HNO3, N2O5, CO, HCl, 
c    H2O2, HO2NO2, H2, CH2O, NO3, O2
c
c----------------------------------------------------------------------

      use chem0d

      implicit none

      logical :: gs_convrg
      real, dimension(nblong) :: qnew 

c  local variables 

      logical, dimension(nblong) :: convrg
      real :: pr, pe
      real :: absdif

      character (LEN=6), dimension(nblong) :: gs_species = 
     $  (/ 'CH4   ', 'H2O   ', 'H2    ', 'N2O   ', 'CO    ', 'HCl   ', 
     $     'HNO3  ', 'H2O2  ', 'N2O5  ', 'HO2NO2', 'CH2O  ', 'NO3   ', 
     $     'O2    ' /)

      integer :: nn, iter

c----------------------------------------------------------------------

      convrg(:) = .false. 

c----------------------------------------------------------------------

c... Here starts the iteration loop (niter is set in subr. msetvar)

      do iter = 1,niter

c----------------------------------------------------------------------
c
c            chemistry of long-lived species  
c            Solve for longest-lived first: based on the chemical 
c            lifetimes around 50 km, from Brasseur & Solomon.
c
c----------------------------------------------------------------------
c
c              o2
c
            pr = hk_1 * c_o * c_o              ! production via Ox rxns
     $         + ( 2 * hk_3 * c_o
     $           + tj0(3) + tj0(4)
     $           + a_2 * c_h ) * c_o3          ! production via HOx rxns
     $         + ( a_5 * c_o 
     $           + a_6 * c_o3 ) * c_oh
     $         + ( 2 * a_6b * c_o3
     $           + a_7 * c_o
     $           + a_17 * c_oh
     $           + a_23b * c_h
     $           + a_27 * c_ho2 ) * c_ho2

            pe = tj0(1)
     $         + tj0(2)
     $         + hk_2 * c_o
     $         + a_1 * c_h

            absdif = o2 
            o2 = (c_prev(26) + dtchem * pr)
     $                      / (1. + dtchem * pe)
            o2 = amax1(o2, dnmin)

            absdif = abs(o2 - absdif)
            convrg(13) = (absdif .le. o2*convlim)

c----------------------------------------------------------------------
c
c              ch4
c
            pe = tj0(5) 
     $         + c_2 * c_oh
     $         + c_1 * c_o1d 
     $         + d_5 * c_cl

            absdif = c_ch4
            c_ch4 = c_prev(2) / (1. + dtchem * pe)
            c_ch4 = amax1(c_ch4, dnmin)

            absdif = abs(c_ch4 - absdif)
            convrg(1) = (absdif .le. c_ch4*convlim)

c----------------------------------------------------------------------
c
c              h2o
c
            pr = (a_17 * c_oh 
     $            + a_23c * c_h) * c_ho2
     $         + (a_19 * c_h2
     $            + a_82 * c_oh
     $            + a_30 * c_h2o2
     $            + b_27 * c_hno3
     $            + b_28 * c_ho2no2
     $            + c_2 * c_ch4
     $            + c_8 * c_ch2o
     $            + d_11 * c_hcl
     $            + d_34 * c_hocl) * c_oh
     $            + .05 * tj0(5) * c_ch4
            pe = tj0(31) 
     $         + tj0(32) 
     $         + tj0(33) 
     $         + a_1et * c_o1d

            absdif = c_h2o
            c_h2o = (c_prev(3) + dtchem * pr)
     $                      / (1. + dtchem * pe)
            c_h2o = amax1(c_h2o, dnmin)

            absdif = abs(c_h2o - absdif)
            convrg(2) = (absdif .le. c_h2o*convlim)

c----------------------------------------------------------------------
c
c             h2 
c
	    pr = a_23b * c_h * c_ho2
     $         + 1.44 * tj0(5) * c_ch4
     $         + tj0(25) * c_ch2o
     $         + tj0(32) * c_h2o

            pe = a_3et * c_o1d + a_19 * c_oh
     $         + d_6 * c_cl

            absdif = c_h2
            c_h2 = (c_prev(16) + dtchem * pr)
     $                      / (1. + dtchem * pe)
            c_h2 = amax1(c_h2, dnmin)

            absdif = abs(c_h2 - absdif)
            convrg(3) =  (absdif .le. c_h2 * convlim)

c----------------------------------------------------------------------
c
c              n2o
c
            pr = hk_21 * c_o1d * n2
            pe = tj0(19)
     $         + (b_38 + b_39) * c_o1d

            absdif = c_n2o
            c_n2o = (c_prev(1) + dtchem * pr)
     $                      / (1. + dtchem * pe)
            c_n2o = amax1(c_n2o, dnmin)

            absdif = abs(c_n2o - absdif)
            convrg(4) = (absdif .le. c_n2o*convlim)

c----------------------------------------------------------------------
c
c              co
c
            pr = ( c_8 * c_oh
     $           + c_9 * c_o 
     $           + d_10 * c_cl
     $           + tj0(24)
     $           + tj0(25) )  *  c_ch2o    
     $         + tj0(18) * c_co2
     $         + .38 * tj0(5) * c_ch4
            pe = a_36 * c_oh

            absdif = c_co
            c_co = (c_prev(9) + dtchem * pr)
     $                      / (1. + dtchem * pe)
            c_co = amax1(c_co, dnmin)

            absdif = abs(c_co - absdif)
            convrg(5) = (absdif .le. c_co*convlim)

c----------------------------------------------------------------------
c
c              hcl
c
            pr = ( d_5 * c_ch4
     $           + d_7 * c_ho2
     $           + d_6 * c_h2
     $           + d_10 * c_ch2o
     $           + d_84 * c_h2o2) * c_cl
     $         + d_46 * c_clo * c_oh
            pe = tj0(16)
     $         + d_11 * c_oh
     $         + d_75 * c_o1d
     $         + d_85 * c_o

            absdif = c_hcl
            c_hcl = (c_prev(10) + dtchem * pr)
     $                      / (1. + dtchem * pe)
            c_hcl = amax1(c_hcl, dnmin)

            absdif = abs(c_hcl - absdif)
            convrg(6) = (absdif .le. c_hcl*convlim)
 
c----------------------------------------------------------------------
c
c              hno3
c
            pr = b_22 * c_oh * c_no2
     $         + b_73b * c_ho2 * c_no3
            pe = tj0(7) 
     $         + b_27 * c_oh

            absdif = c_hno3
            c_hno3 = (c_prev(5) + dtchem * pr)
     $                      / (1. + dtchem * pe)
            c_hno3 = amax1(c_hno3, dnmin)

            absdif = abs(c_hno3 - absdif)
            convrg(7) = (absdif .le. c_hno3*convlim)
 
c----------------------------------------------------------------------
c
c              h2o2
c
            pr = a_27 * c_ho2 * c_ho2
     $         + a_83 * c_oh * c_oh
            pe = tj0(13)
     $         + a_30 * c_oh
     $         + a_81 * c_o
     $         + d_84 * c_cl

            absdif = c_h2o2
            c_h2o2 = (c_prev(13) + dtchem * pr)
     $                      / (1. + dtchem * pe)
            c_h2o2 = amax1(c_h2o2, dnmin)

            absdif = abs(c_h2o2 - absdif)
            convrg(8) = (absdif .le. c_h2o2*convlim)
 
c----------------------------------------------------------------------
c
c              n2o5
c
            pr = b_12 * c_no3 * c_no2

            pe = tj0(11)
     $         + tj0(12)
     $         + b_32

            absdif = c_n2o5
            c_n2o5 = (c_prev(6) + dtchem * pr)
     $                      / (1. + dtchem * pe)
            c_n2o5 = amax1(c_n2o5, dnmin)

            absdif = abs(c_n2o5 - absdif)
            convrg(9) = (absdif .le. c_n2o5*convlim)
 
c----------------------------------------------------------------------
c
c              ho2no2
c
            pr = b_23 * c_ho2 * c_no2
            pe = tj0(9)
     $         + tj0(10)
     $         + b_28 * c_oh
     $         + b_24

            absdif = c_ho2no2
            c_ho2no2 = (c_prev(15) + dtchem * pr)
     $                      / (1. + dtchem * pe)
            c_ho2no2 = amax1(c_ho2no2, dnmin)

            absdif = abs(c_ho2no2 - absdif)
            convrg(10) = (absdif .le. c_ho2no2*convlim)
 
c----------------------------------------------------------------------
c
c              ch2o
c
            pr = ( c_1 * c_o1d
     $           + c_2 * c_oh
     $           + d_5 * c_cl
     $           + .18 * tj0(5) ) * c_ch4
            pe = tj0(24)
     $         + tj0(25)
     $         + c_8 * c_oh
     $         + c_9 * c_o
     $         + d_10 * c_cl

            absdif = c_ch2o
            c_ch2o = (c_prev(17) + dtchem * pr)
     $                      / (1. + dtchem * pe)
            c_ch2o = amax1(c_ch2o, dnmin)

            absdif = abs(c_ch2o - absdif)
            convrg(11) = (absdif .le. c_ch2o*convlim)
c
c----------------------------------------------------------------------
c
c             no3 
c
            pr = b_9 * c_o3 * c_no2
     $         + b_32 * c_n2o5
     $         + b_27 * c_oh * c_hno3 
     $         + b_81 * c_no2 * c_o 
     $         + d_32 * c_o * c_clono2
     $         + d_87 * c_oh * c_clono2
     $         + (tj0(11) + tj0(12)) * c_n2o5
     $         + tj0(10) * c_ho2no2
     $         + tj0(27) * c_clono2
            pe = b_12 * c_no2
     $         + b_71 * c_o 
     $         + b_72 * c_oh 
     $         + b_84 * c_no
     $         + b_73a * c_ho2 
     $         + b_73b * c_ho2
     $         + d_73 * c_cl
     $         + tj0(29)
     $         + tj0(30)

            absdif = c_no3
            c_no3 = (c_prev(25) + dtchem * pr)
     $                      / (1. + dtchem * pe)
            c_no3 = amax1(c_no3, dnmin)	

            absdif = abs(c_no3 - absdif)
            convrg(12) = (absdif .le. c_no3 * convlim)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     if convergence has been met, exit iteration loop
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         gs_convrg = all(convrg)
         if (gs_convrg) go to 850 

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     end of the iteration loop
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      end do

      if ( .not. gs_convrg ) then
        print * , 'GS failed to converge after ', 
     $              niter, ' iterations'
	do nn = 1, nblong
	  if ( .not. convrg(nn) ) print *, gs_species(nn)
	enddo
      endif

 850  continue

      qnew(1)  = c_ch4
      qnew(2)  = c_h2o
      qnew(3)  = c_h2
      qnew(4)  = c_n2o
      qnew(5)  = c_co
      qnew(6)  = c_hcl
      qnew(7)  = c_hno3
      qnew(8)  = c_h2o2
      qnew(9)  = c_n2o5
      qnew(10) = c_ho2no2
      qnew(11) = c_ch2o
      qnew(12) = c_no3
      qnew(13) = o2

      end
