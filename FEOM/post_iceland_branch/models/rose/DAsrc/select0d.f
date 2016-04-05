
c----------------------------------------------------------------------
      subroutine select0d( k, i, j)
c----------------------------------------------------------------------

      use params, only : nphot, nztrop
      use chem
      use chem0d
      use phys,   only : co2mr
      use time_gcm

      implicit none

c... common parameters and arrays

      integer, intent(in) :: i, j, k
      integer ::  nn

      logical :: firstcall

      data firstcall /.true./

c... temperature independent rate constants 

      if (firstcall)then
         firstcall = .false.
         a_1et = a1et  
         a_3et = a3et  
         a_23a = a23a
         a_23b = a23b
         a_23c = a23c
         a_23  = a23

         b_7a = b7a
         b_38 = b38
         b_39 = b39
         b_71 = b71
         b_72 = b72
         b_73a = b73a
         b_73b = b73b

         c_1 = c1
         c_8 = c8

         d_35 = d35
         d_73 = d73
         d_75 = d75

         hk_7 = hk7

	 nrstep = dtchem

      end if

c... temperature dependent rate constants 

      a_1 = a1(k,i,j)
      a_2 = a2(k,i,j)
      a_5 = a5(k,i,j)
      a_6 = a6(k,i,j)
      a_6b = a6b(k,i,j)
      a_7 = a7(k,i,j)
      a_17 = a17(k,i,j)
      a_19 = a19(k,i,j)
      a_26 = a26(k,i,j)
      a_27 = a27(k,i,j)
      a_30 = a30(k,i,j)
      a_36 = a36(k,i,j)
      a_81 = a81(k,i,j)
      a_82 = a82(k,i,j)
      a_83 = a83(k,i,j)

      b_3 = b3(k,i,j)
      b_4 = b4(k,i,j)
      b_6 = b6(k,i,j)
      b_7 = b7(k,i,j)
      b_9 = b9(k,i,j)
      b_12 = b12(k,i,j)
      b_22 = b22(k,i,j)
      b_23 = b23(k,i,j)
      b_24 = b24(k,i,j)
      b_27 = b27(k,i,j)
      b_28 = b28(k,i,j)
      b_32 = b32(k,i,j)
      b_81 = b81(k,i,j)
      b_82 = b82(k,i,j)
      b_84 = b84(k,i,j)

      c_2 = 0.0
      c_9 = c9(k,i,j)

      d_2 = d2(k,i,j)
      d_3 = d3(k,i,j)
      d_4 = d4(k,i,j)
      d_5 = d5(k,i,j)
      d_6 = d6(k,i,j)
      d_7 = d7(k,i,j)
      d_8 = d8(k,i,j)

      d_10 = d10(k,i,j)
      d_11 = d11(k,i,j)
      d_31 = d31(k,i,j)
      d_32 = d32(k,i,j)
      d_33 = d33(k,i,j)
      d_34 = d34(k,i,j)
      d_46 = d46(k,i,j)
      d_47 = d47(k,i,j)

      d_82 = 0.0
      d_83 = d83(k,i,j)
      d_84 = d84(k,i,j)
      d_85 = d85(k,i,j)
      d_87 = d87(k,i,j)

      hk_1 = hk1(k,i,j)
      hk_2 = hk2(k,i,j)
      hk_3 = hk3(k,i,j)
      hk_4 = hk4(k,i,j)
      hk_5 = hk5(k,i,j)
      hk_21 = hk21(k,i,j)

c... photolysis rates

      do nn = 1, nphot
        tj0(nn) = tj(k,i,j,nn)
      enddo

c... concentrations

      den      = hnm(k,i,j) 

      c_n2o    = qn1(k,i,j,1) * den
      c_ch4    = qn1(k,i,j,2) * den
      c_h2o    = qn1(k,i,j,3) * den
      c_o1d    = qn1(k,i,j,4) * den
      c_hno3   = qn1(k,i,j,5) * den
      c_n2o5   = qn1(k,i,j,6) * den
      c_h      = qn1(k,i,j,7) * den
      c_oh     = qn1(k,i,j,8) * den
      c_co     = qn1(k,i,j,9) * den
      c_hcl    = qn1(k,i,j,10) * den
      c_clono2 = qn1(k,i,j,11) * den
      c_hocl   = qn1(k,i,j,12) * den
      c_h2o2   = qn1(k,i,j,13) * den
      c_ho2    = qn1(k,i,j,14) * den
      c_ho2no2 = qn1(k,i,j,15) * den
      c_h2     = qn1(k,i,j,16) * den
      c_ch2o   = qn1(k,i,j,17) * den
      c_o      = qn1(k,i,j,18) * den
      c_o3     = qn1(k,i,j,19) * den
      c_cl     = qn1(k,i,j,20) * den
      c_clo    = qn1(k,i,j,21) * den
      c_n      = qn1(k,i,j,22) * den
      c_no     = qn1(k,i,j,23) * den
      c_no2    = qn1(k,i,j,24) * den
      c_no3    = qn1(k,i,j,25) * den

      o2       = qn1(k,i,j,26) * den
      n2       = ( 1. - qn1(k,i,j,26) - qn1(k,i,j,18)) * den

      c_co2    = co2mr(k+nztrop) * den   

      c_prev   = qn1(k,i,j,:) * den

      if(sc2d(i,j).le.50.) then
         c_n2d = n2d(k) * den
      else
         c_n2d = 0.
      end if

      dnmin = 1e-6 

      end subroutine select0d

