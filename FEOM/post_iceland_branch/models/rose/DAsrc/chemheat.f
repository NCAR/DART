      subroutine chemheat

c     -------------------------------------------------------------------
c                          Chemical Heating Model
c                             Amanda Szymczak
c                              June 13, 2000
c
c      The purpose of this subroutine is to calculate the chemical heating of 
c      the atmosphere.
c      
c      Modified by AKS to calculate HOx heating terms directly (not assuming
c          photochimical equilibrium)
c     -------------------------------------------------------------------

      use params
      use chem
      use phys
     
      implicit none

c     local arrays
      real er1, er2, er3, er4, er5, er6, er7, const, SpD
      integer i, j, k
      
c      rate coff: hk1 hk2 hk3 a2 a5 a7 a1  - found in rates.f
c      chemicals: qn1(19) qn1(18) qn1(7) qn1(8) qn1(14)
c                  o3       o       h      oh    ho2

c....................................................................
c  constants
c      hcp = 1004.       ! J/(Kelvin*kg) heat capacity of air
c      Jev = 1.6021e-19  ! Joules per ev
c      dd = 4.81069e-26  ! (.02897 kg/mole air)/(6.022e23 molecules/mole)
      const = 1.6021e-19 / 4.81069e-26 / 1004. ! convert units to K/sec

      er1 = 5.11 * const  !   O   + O  + M -> O2  + M
      er2 = 1.05 * const  !   O   + O2 + M -> O3  + M
      er3 = 4.06 * const  !   O   + O3     -> O2  + O2
      er4 = 3.34 * const  !   H   + O3     -> OH  + O2
      er5 = 0.72 * const  !   OH  + O      -> H   + O2
      er6 = 2.39 * const  !   HO2 + O      -> OH  + O2
      er7 = 2.0  * const  !   H   + O2 + M -> HO2 + M
 
      E2 = 0.
      ch_heat = 0.

      do k=18,nz
         do j=1,ny
            do i=1,nx
            
c     units of E and ch_heat are Kelvin/sec
               E1(k,i,j) = hnm(k,i,j) * er1 * hk1(k,i,j)
     $                     * qn1(k,i,j,18) * qn1(k,i,j,18)
               if(k.ge.22) then
                  E2(k,i,j) = hnm(k,i,j) * er2 * hk2(k,i,j)
     $                        * qn1(k,i,j,18) * q_o2(k,j)
               end if
               E3(k,i,j) = hnm(k,i,j) * er3 * hk3(k,i,j)
     $                     * qn1(k,i,j,18) * qn1(k,i,j,19)
               E4(k,i,j) = 0.6 * hnm(k,i,j) * er4 * a2(k,i,j)
     $                     * qn1(k,i,j,7) * qn1(k,i,j,19)
               E5(k,i,j) = hnm(k,i,j) * er5 * a5(k,i,j)
     $                     * qn1(k,i,j,8) * qn1(k,i,j,18)
               E6(k,i,j) = hnm(k,i,j) * er6 * a7(k,i,j)
     $                     * qn1(k,i,j,14) * qn1(k,i,j,18)
               E7(k,i,j) = hnm(k,i,j) * er7 * a1(k,i,j)
     $                     * qn1(k,i,j,7) * q_o2(k,j)

               ch_heat(k,i,j) = E1(k,i,j) + E2(k,i,j)
     $                        + E3(k,i,j) + E4(k,i,j)
     $                        + E5(k,i,j) + E6(k,i,j)
     $                        + E7(k,i,j)
    
            end do
         end do
      end do
      
      end










