
c----------------------------------------------------------------------
      subroutine update0d( k, i, j)
c----------------------------------------------------------------------

      use chem
      use chem0d

      implicit none

c... common parameters and arrays

      integer, intent(in) :: i, j, k

      qn1(k,i,j,1)  = c_n2o    / den
      qn1(k,i,j,2)  = c_ch4    / den
      qn1(k,i,j,3)  = c_h2o    / den
      qn1(k,i,j,4)  = c_o1d    / den
      qn1(k,i,j,5)  = c_hno3   / den
      qn1(k,i,j,6)  = c_n2o5   / den
      qn1(k,i,j,7)  = c_h      / den
      qn1(k,i,j,8)  = c_oh     / den
      qn1(k,i,j,9)  = c_co     / den
      qn1(k,i,j,10) = c_hcl    / den
      qn1(k,i,j,11) = c_clono2 / den
      qn1(k,i,j,12) = c_hocl   / den
      qn1(k,i,j,13) = c_h2o2   / den
      qn1(k,i,j,14) = c_ho2    / den
      qn1(k,i,j,15) = c_ho2no2 / den
      qn1(k,i,j,16) = c_h2     / den
      qn1(k,i,j,17) = c_ch2o   / den
      qn1(k,i,j,18) = c_o      / den
      qn1(k,i,j,19) = c_o3     / den
      qn1(k,i,j,20) = c_cl     / den
      qn1(k,i,j,21) = c_clo    / den
      qn1(k,i,j,22) = c_n      / den
      qn1(k,i,j,23) = c_no     / den
      qn1(k,i,j,24) = c_no2    / den
      qn1(k,i,j,25) = c_no3    / den
      qn1(k,i,j,26) = o2       / den 


      end subroutine update0d

