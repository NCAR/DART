c----------------------------------------------------------------------
      subroutine chem_nr( nr_convrg, q1 )
c----------------------------------------------------------------------
c
c     Newton-Raphson iteration 
c
c----------------------------------------------------------------------
c
c   Begin solving the fully-implicit backward Euler system of
c   non-linear equations using Newton-Raphson iteration. Done
c   for fast reacting members of the Ox, HOx and ClOx families.
c
c   Define G = q(n+1) - q(n) - nrstep * f(q(n+1))
c
c   solve system G=0 for q(n+1) by Newton-Raphson iteration:
c
c   q(m+1) = q(m) - y  where  y solves G(m) = J(m) * y
c
c----------------------------------------------------------------------
c
c   Species solved for in this routine:
c 
c   O, O3, H, OH, HO2, Cl, ClO, N, NO, NO2
c
c----------------------------------------------------------------------

      use chem0d

      implicit none

      logical :: nr_convrg
      real, dimension (nbshort) :: q1

c... local variables
      real, dimension (nbshort, nbshort) :: jacob
      real, dimension (nbshort) :: q0, qprev, f1
      real, dimension (nbshort) :: g1
      integer :: iter, info, nn

      real, dimension (nbshort, nbshort) :: a1
      real, dimension (nbshort) :: b1, x1

      integer, parameter :: itermax = 32 
      logical, dimension (nbshort) :: convrg
      
      character (LEN=6), dimension(nbshort) :: nr_species = 
     $  (/ 'O     ', 'O3    ', 'H     ', 'OH    ', 'HO2    ', 
     $     'Cl    ', 'ClO   ', 'N     ', 'NO    ', 'NO2    ',
     $     'HOCl  ', 'ClONO2' /)

      convrg = .false.

      q0(1)  = c_prev(18)  ! o(3p)
      q0(2)  = c_prev(19)  ! o3
      q0(3)  = c_prev(7)   ! h
      q0(4)  = c_prev(8)   ! oh
      q0(5)  = c_prev(14)  ! ho2
      q0(6)  = c_prev(20)  ! cl 
      q0(7)  = c_prev(21)  ! clo 
      q0(8)  = c_prev(22)  ! n 
      q0(9)  = c_prev(23)  ! no 
      q0(10) = c_prev(24)  ! no2 
      q0(11) = c_prev(12)  ! hocl
      q0(12) = c_prev(11)  ! clono2

      q1(1)  = c_o 
      q1(2)  = c_o3 
      q1(3)  = c_h
      q1(4)  = c_oh
      q1(5)  = c_ho2
      q1(6)  = c_cl
      q1(7)  = c_clo
      q1(8)  = c_n
      q1(9)  = c_no
      q1(10) = c_no2
      q1(11) = c_hocl
      q1(12) = c_clono2

      call jfixed()

      do iter = 1, itermax  ! Newton-Raphson iteration

	qprev = q1

	call eval_f ( q1, f1 )        ! calc. chemical forcing

	g1 = q1 - q0 - f1
	      
	call eval_j ( q1, jacob )     ! calc. Jacobian
        a1 = jacob
	b1 = g1

c... solve linear equation G(m) = J(m) * y

        call NR_SOLVE( jacob, g1, nbshort, INFO )

        if (info.ne.0) then
          print 40,info
 40       format('NR_SOLVE error',i3, ' halting')
          do nn = 1, nbshort
	    print *, nr_species(nn), q1(nn)
	  enddo
          stop
        end if

c... update q(m+1) 
        q1 = q1 - g1

	where( q1 < dnmin )           ! limit iterate
	  q1 = dnmin 
	endwhere
              
	convrg = ( abs(q1 - qprev) .le. qprev * convlim )
	nr_convrg = all( convrg )
	if ( nr_convrg ) exit

      enddo

      if ( .not. nr_convrg ) then
         print * , 'NR failed to converge after ', 
     $              itermax, ' iterations'
         do nn = 1, nbshort
	   if ( .not. convrg(nn) ) print *, nr_species(nn)
	 enddo
      endif

      c_o   =  q1(1)  
      c_o3  =  q1(2)
      c_h   =  q1(3)
      c_oh  =  q1(4)
      c_ho2 =  q1(5)
      c_cl  =  q1(6)
      c_clo =  q1(7)
      c_n   =  q1(8)
      c_no  =  q1(9)
      c_no2 =  q1(10)
      c_hocl   = q1(11)
      c_clono2 = q1(12)

      return

      end subroutine chem_nr

c----------------------------------------------------------------------
      subroutine jfixed
c----------------------------------------------------------------------
c
c   calculate constants in the Jacobian
c
c----------------------------------------------------------------------

      use chem0d

      cj11 = a_81  * c_h2o2
     $     + hk_2 * o2
     $     + b_71 * c_no3
     $     + d_85 * c_hcl
     $     + c_9 * c_ch2o

      cj12 = tj0(3) + tj0(4)

      cj18 = b_7 * o2

      cj19 = tj0(28) 

      cj110 = tj0(6) 

      cj21 = hk_2 * o2 

      cj22 = tj0(3) + tj0(4)
      
      cj33 = a_1 * o2 

      cj34 = a_36 * c_co 
     $     + a_19 * c_h2
      
      cj36 = d_6 * c_h2

      cj41 = a_81 * c_h2o2
     $     + c_9 * c_ch2o
     $     + d_85 * c_hcl 

      cj44 = a_19 * c_h2
     $     + a_30 * c_h2o2
     $     + b_72 * c_no3
     $     + b_27 * c_hno3
     $     + b_28 * c_ho2no2
     $     + a_36 * c_co
     $     + c_8 * c_ch2o
     $     + d_11 * c_hcl
      
      cj45 = b_73a * c_no3

      cj46 = d_5 * c_ch4   

      cj411 = tj0(8)

      cj51 = c_9 * c_ch2o
     $     + a_81 * c_h2o2

      cj53 = a_1 * o2 

      cj54 = a_30 * c_h2o2
     $     + b_72 * c_no3
     $     + c_8 * c_ch2o 

      cj55 = ( b_73a + b_73b ) * c_no3

      cj56 = d_10 * c_ch2o + d_84 * c_h2o2 

      cj61 = d_85 * c_hcl

      cj64 = d_11 * c_hcl

      cj66 = d_5 * c_ch4 
     $     + d_10 * c_ch2o
     $     + d_6 * c_h2
     $     + d_84 * c_h2o2 
     $     + d_73 * c_no3 

      cj611 = tj0(8)

      cj612 = tj0(27)

      cj76 = d_73 * c_no3

      cj712 = tj0(26)

      cj89 = tj0(28)

      cj98 = b_7 * o2 

      cj99 = b_84 * c_no3 + tj0(28)

      cj910 = tj0(6)

      cj101 = b_71 * c_no3 

      cj104 = b_72 * c_no3 + b_28 * c_ho2no2

      cj105 = b_73a * c_no3

      cj106 = d_73 * c_no3

      cj109 = 2 * b_84 * c_no3

      cj1010 = b_12 * c_no3 + tj0(6)

      cj1012 = tj0(26)

      cj1111 = tj0(8)

      cj1212 = tj0(26) + tj0(27)

      end subroutine jfixed 

c----------------------------------------------------------------------
      subroutine eval_f ( q, f )
c----------------------------------------------------------------------
c
c     calc. chemical forcing
c
c----------------------------------------------------------------------
      
      use chem0d

      implicit none

      real, dimension (nbshort), intent(in) :: q
      real, dimension (nbshort), intent(out) :: f

c... local variables

      real :: pr, pe 

      c_o = q(1)
      c_o3 = q(2)
      c_h = q(3)
      c_oh = q(4)
      c_ho2 = q(5)
      c_cl = q(6)
      c_clo = q(7)
      c_n = q(8)
      c_no = q(9)
      c_no2 = q(10)
      c_hocl = q(11)
      c_clono2 = q(12)

c... update O(1D)
      c_o1d = ( tj0(4) * c_o3
     $        + tj0(2) * o2
     $        + tj0(19) * c_n2o
     $        + tj0(32) * c_h2o )
     $        / ( hk_4 * n2 + hk_5 * o2 + hk_7 * c_o3)

c... f(o)

      pr = a_82 * c_oh * c_oh
     $   + a_23c * c_h * c_ho2
     $   + b_7 * c_n * o2
     $   + b_7a * c_n2d * o2
     $   + b_6 * c_n * c_no 
     $   + 2 * (tj0(1) + tj0(2)) * o2
     $   + (tj0(3) + tj0(4)) * c_o3 
     $   + tj0(28) * c_no 
     $   + tj0(6) * c_no2 
     $   + tj0(29) * c_no3
     $   + tj0(19) * c_n2o
     $   + .18 * tj0(5) * c_ch4
     $   + (tj0(32) + tj0(33)) * c_h2o 
     $   + tj0(12) * c_n2o5 
     $   + tj0(18) * c_co2

      pe = a_5 * c_oh
     $   + a_7 * c_ho2
     $   + a_81 * c_h2o2 
     $   + b_3 * c_no2 
     $   + hk_2 * o2
     $   + hk_3 * c_o3 
     $   + 2 * hk_1 * c_o
     $   + d_3 * c_clo
     $   + b_82 * c_no
     $   + b_81 * c_no2
     $   + b_71 * c_no3
     $   + c_9 * c_ch2o
     $   + d_32 * c_clono2
     $   + d_35 * c_hocl
     $   + d_85 * c_hcl

      f(1) = pr - pe * c_o

c... f(o3)

      pr = hk_2 * c_o * o2
c
      pe = a_2 * c_h
     $   + a_6 * c_oh
     $   + a_6b * c_ho2
     $   + hk_3 * c_o
     $   + d_2 * c_cl 
     $   + b_4 * c_no
     $   + b_9 * c_no2
     $   + tj0(3) + tj0(4)

      f(2) = pr - pe * c_o3

c... f(h)

      pr = ( a_5 * c_o 
     $   + a_36 * c_co) *c_oh
     $   + ( a_3et * c_o1d  
     $   + a_19 * c_oh 
     $   + d_6 * c_cl ) * c_h2 
     $   + (tj0(31) + 2 * tj0(33)) * c_h2o 
     $   + tj0(16) * c_hcl
     $   + tj0(24) * c_ch2o

      pe = a_1 * o2
     $   + a_2 * c_o3
     $	 + ( a_23a + a_23b + a_23c ) * c_ho2

      f(3) = pr - pe * c_h

c... f(oh)

      pr = ( 2.*a_1et * c_h2o
     $     + 2.*c_1 * c_ch4
     $     + a_3et * c_h2
     $     + d_75 * c_hcl ) * c_o1d
     $   + ( a_6b * c_o3
     $     + a_7 * c_o
     $     + 2.*a_23a * c_h
     $     + a_26 * c_no
     $     + b_73a * c_no3
     $     + d_83  * c_cl ) * c_ho2
     $   + ( a_81  * c_h2o2
     $     + c_9   * c_ch2o 
     $     + d_35  * c_hocl
     $     + d_85  * c_hcl ) * c_o
     $   + a_2 * c_h * c_o3
     $   + d_5 * c_ch4 * c_cl
     $   + tj0(31) * c_h2o
     $   + tj0(7) * c_hno3
     $   + tj0(8) * c_hocl
     $   + tj0(10) * c_ho2no2
     $   + 2. * tj0(13) * c_h2o2
     $   + .33 * tj0(5) * c_ch4

      pe = (2.*a_82 + 2.*a_83) * c_oh
     $   + a_17 * c_ho2
     $   + a_5 * c_o
     $   + a_6 * c_o3
     $   + a_19 * c_h2
     $   + a_30 * c_h2o2
     $   + a_36 * c_co
     $   + b_22 * c_no2
     $   + b_27 * c_hno3
     $   + b_28 * c_ho2no2
     $   + b_72 * c_no3
     $   + c_8 * c_ch2o
     $   + ( d_8 + d_46 ) * c_clo
     $   + d_11 * c_hcl
     $   + d_34 * c_hocl
     $   + d_87 * c_clono2

      f(4) = pr - pe * c_oh

c... f(ho2)

      pr = a_1 * o2 * c_h
     $   + ( a_6 * c_o3
     $     + a_30 * c_h2o2
     $     + b_72 * c_no3
     $     + c_8 * c_ch2o
     $     + d_8 * c_clo) * c_oh
     $   + ( a_81 * c_h2o2
     $     + c_9 * c_ch2o ) * c_o
     $   + ( d_10 * c_ch2o 
     $     + d_84 * c_h2o2 )* c_cl
     $   + ( b_24
     $     + tj0(9) ) * c_ho2no2 
     $   + tj0(24) * c_ch2o
     $   + .33 * tj0(5) * c_ch4

      pe = 2. * a_27 * c_ho2
     $   + a_6b * c_o3
     $   + a_7 * c_o
     $   + a_17 * c_oh
     $   + a_23 * c_h
     $   + a_26 * c_no
     $   + b_23 * c_no2
     $   + b_73a * c_no3
     $   + b_73b * c_no3
     $   + ( d_83 + d_7 ) * c_cl
     $   + d_33 * c_clo

      f(5) = pr - pe * c_ho2

c... f(cl)

      pr = ( d_3 * c_o
     $     + d_8 * c_oh
     $     + d_4 * c_no ) * c_clo 
     $   + ( d_11 * c_oh
     $     + d_75 * c_o1d
     $     + d_85 * c_o ) * c_hcl 
     $   + tj0(16) * c_hcl
     $   + tj0(8) * c_hocl
     $   + tj0(27) * c_clono2

      pe = d_5 * c_ch4 
     $   + d_10 * c_ch2o
     $   + d_2 * c_o3
     $   + d_6 * c_h2
     $   + ( d_7 + d_83 ) * c_ho2
     $   + d_84 * c_h2o2 
     $   + d_73 * c_no3 
     $   + d_82 * c_clono2

      f(6) = pr - pe * c_cl

c... f(clo)

      pr = ( d_2 * c_o3 
     $     + d_73 * c_no3
     $     + d_83 * c_ho2 ) * c_cl 
     $   + ( d_32 * c_clono2
     $     + d_35 * c_hocl ) * c_o
     $   + d_34 * c_oh * c_hocl
     $   + tj0(26) * c_clono2

      pe = d_3 * c_o
     $   + ( d_8 + d_46 ) * c_oh
     $   + d_33 * c_ho2
     $   + d_4 * c_no
     $   + d_31 * c_no2

      f(7) = pr - pe * c_clo

c... f(n) 

      pr = tj0(28) * c_no 

      pe = b_7 * o2 + b_6 * c_no

      f(8) = pr - pe * c_n

c... f(no) 

      pr = b_3 * c_o * c_no2
     $   + b_7 * c_n * o2 
     $   + b_7a * c_n2d * o2 
     $   + 2 * b_39 * c_o1d * c_n2o
     $   + tj0(6) * c_no2
     $   + tj0(30) * c_no3 
     $   + tj0(12) * c_n2o5 

      pe = b_4 * c_o3
     $   + b_82 * c_o
     $   + a_26 * c_ho2
     $   + b_6 * c_n 
     $   + b_84 * c_no3
     $   + d_4 * c_clo
     $   + tj0(28)

      f(9) = pr - pe * c_no

c... f(no2)

      pr = ( b_4 * c_o3 
     $     + a_26 * c_ho2
     $     + 2 * b_84 * c_no3
     $     + d_4 * c_clo
     $     + b_82 * c_o ) * c_no
     $   + ( b_71 * c_o
     $     + b_72 * c_oh
     $     + d_73 * c_cl
     $     + b_73a * c_ho2 ) * c_no3 
     $   + b_32 * c_n2o5  
     $   + ( b_28 * c_oh
     $     + b_24) * c_ho2no2
     $   + tj0(9) * c_ho2no2 
     $   + tj0(11) * c_n2o5 
     $   + tj0(7) * c_hno3
     $   + tj0(26) * c_clono2
     $   + tj0(29) * c_no3 

      pe = ( b_3 + b_81 ) * c_o
     $   + b_9 * c_o3
     $   + b_12 * c_no3
     $   + b_22 * c_oh
     $   + b_23 * c_ho2 
     $   + d_31 * c_clo
     $   + tj0(6)

      f(10) = pr - pe * c_no2

c... f(hocl)

      pr =   d_33 * c_ho2 * c_clo 
     $     + d_87 * c_oh * c_clono2

      pe = d_34 * c_oh
     $   + d_35 * c_o
     $   + tj0(8)

      f(11) = pr - pe * c_hocl

c... f(clono2)

      pr =   d_31 * c_no2 * c_clo 

      pe = d_32 * c_o
     $   + d_87 * c_oh
     $   + tj0(26)
     $   + tj0(27)

      f(12) = pr - pe * c_clono2

!..........................

      f = f * nrstep

      end subroutine eval_f

c----------------------------------------------------------------------
      subroutine eval_j ( q1, jacob )
c----------------------------------------------------------------------
c
c     calc. Jacobian
c
c----------------------------------------------------------------------
      
      use chem0d 

      implicit none

      real, dimension (nbshort), intent(in) :: q1
      real, dimension (nbshort,nbshort), intent(out) :: jacob

      integer :: nn

      jacob = 0

c... df(o)/dqi

      jacob(1,1) =  - (4 * hk_1 * q1(1) + hk_3 * q1(2)
     $          + a_5 * q1(4) + a_7 * q1(5)
     $          + (b_3 + b_81) * q1(10)
     $          + b_82 * q1(9) 
     $          + d_3 * q1(7) 
     $          + d_32 * q1(12)
     $          + d_35 * q1(11)
     $          + cj11 ) 

      jacob(1,2) = cj12 - hk_3 * q1(1)

      jacob(1,3) = a_23c * q1(5)

      jacob(1,4) = 2 * a_82 * q1(4) - a_5 * q1(1) 

      jacob(1,5) = a_23c * q1(3) - a_7 * q1(1)

      jacob(1,7) = - d_3 * q1(1) 

      jacob(1,8) = cj18 + b_6 * q1(9)

      jacob(1,9) = cj19 + b_6 * q1(8) - b_82 * q1(1) 

      jacob(1,10) = cj110 - (b_3 + b_81) * q1(1) 

      jacob(1,11) = - d_35 * q1(1)

      jacob(1,12) = - d_32 * q1(1)

c... df(o3)/dqi

      jacob(2,1) = - hk_3 * q1(2) + cj21 

      jacob(2,2) = - (hk_3 * q1(1) + a_6 * q1(4)
     $          + a_2 * q1(3) + a_6b * q1(5)  
     $          + d_2 * q1(6) +  b_4 * q1(9) 
     $          + b_9 * q1(10) + cj22 )

      jacob(2,3) = - a_2 * q1(2) 

      jacob(2,4) = - a_6 * q1(2) 

      jacob(2,5) = - a_6b * q1(2) 

      jacob(2,6) = - d_2 * q1(2)

      jacob(2,9) = - b_4 * q1(2)

      jacob(2,10) = - b_9 * q1(2)

c... df(h)/dqi

      jacob(3,1) = a_5 * q1(4) 

      jacob(3,2) = - a_2 * q1(3) 

      jacob(3,3) = - (a_2 * q1(2)  
     $           + a_23 * q1(5)
     $           + cj33)

      jacob(3,4) = a_5 * q1(1) + cj34

      jacob(3,5) = - a_23 * q1(3) 

      jacob(3,6) = cj36

c... df(oh)/dqi

      jacob(4,1) = - a_5 * q1(4) 
     $           + a_7 * q1(5) 
     $           + d_35 * q1(11)
     $           + cj41

      jacob(4,2) = - a_6 * q1(4)
     $           + a_2 * q1(3) 
     $           + a_6b * q1(5) 

      jacob(4,3) = a_2 * q1(2)
     $           + 2 * a_23a * q1(5) 

      jacob(4,4) = - ( 4 * ( a_82 + a_83 ) * q1(4)
     $           + a_5 * q1(1) + a_6 * q1(2) 
     $           + (d_8 + d_46) * q1(7)
     $           + b_22 * q1(10) 
     $           + a_17 * q1(5)
     $           + d_34 * q1(11)
     $           + d_87 * q1(12)
     $           + cj44 )

      jacob(4,5) = -a_17 * q1(4)
     $           + a_6b * q1(2)
     $           + a_7 * q1(1)
     $           + 2 * a_23a * q1(3)
     $           + d_83 * q1(6) 
     $           + a_26 * q1(9) + cj45

      jacob(4,6) = d_83 * q1(5) + cj46

      jacob(4,7) = - ( d_8 + d_46 ) * q1(4)

      jacob(4,9) = a_26 * q1(5)

      jacob(4,10) = - b_22 * q1(4)  

      jacob(4,11) = - d_34 * q1(4)
     $              + cj411 + d_35 * q1(1)

      jacob(4,12) = - d_87 * q1(4)

c... df(ho2)/dqi

      jacob(5,1) = cj51 - a_7 * q1(5)

      jacob(5,2) = a_6 * q1(4)
     $           - a_6b * q1(5)

      jacob(5,3) = -a_23 * q1(5) + cj53
      
      jacob(5,4) = a_6 * q1(2) - a_17 * q1(5)
     $           + d_8 * q1(7) + cj54

      jacob(5,5) = - ( a_23 * q1(3) + a_17 * q1(4)
     $           + a_7 * q1(1) + a_6b * q1(2)
     $           + ( d_83 + d_7 ) * q1(6)
     $           + d_33 * q1(7) 
     $           + 4 * a_27 * q1(5) + a_26 * q1(9) 
     $           + b_23 * q1(10) + cj55 )

      jacob(5,6) = cj56 - ( d_83 + d_7 ) * q1(5)

      jacob(5,7) =  d_8 * q1(4) - d_33 * q1(5)

      jacob(5,9) =  - a_26 * q1(5)

      jacob(5,10) = - b_23 * q1(5) 

c... df(cl)/dqi

      jacob(6,1) = d_3 * q1(7) + cj61
      
      jacob(6,2) = - d_2 * q1(6)

      jacob(6,4) = d_8 * q1(7) + cj64

      jacob(6,5) = - ( d_7 + d_83 ) * q1(6)

      jacob(6,6) = - ( ( d_7 + d_83 ) * q1(5) + d_2 * q1(2) 
     $                 + d_82 * q1(12)
     $                 + cj66 )
      
      jacob(6,7) = d_3 * q1(1) + d_8 * q1(4) + d_4 * q1(9) 

      jacob(6,9) =  d_4  * q1(7)

      jacob(6,11) = cj611 

      jacob(6,12) = cj612 

c... df(clo)/dqi

      jacob(7,1) = - d_3 * q1(7) + d_32 * q1(12) + d_35 * q1(11) 
      
      jacob(7,2) = d_2 * q1(6)

      jacob(7,4) = d_34 * q1(11) - ( d_8 + d_46 ) * q1(7)

      jacob(7,5) = d_83 * q1(6) - d_33 * q1(7)

      jacob(7,6) = d_2 * q1(2) + d_83 * q1(5) + cj76
      
      jacob(7,7) = - ( d_3 * q1(1) + ( d_8 + d_46 ) * q1(4) 
     $               + d_33 * q1(5)  
     $               + d_4 * q1(9) + d_31 * q1(10) )

      jacob(7,9) =  - d_4  * q1(7)

      jacob(7,10) =  - d_31  * q1(7)

      jacob(7,11) = d_34 * q1(4) + d_35 * q1(1) 

      jacob(7,12) = d_32 * q1(1) + cj712

c... df(n)/dqi

      jacob(8,8) = - ( b_7 * o2 + b_6 * q1(9) )

      jacob(8,9) = cj89 - b_6 * q1(8) 

c... df(no)/dqi

      jacob(9,1) = b_3 * q1(10) - b_82 * q1(9) 

      jacob(9,2) = - b_4 * q1(9) 

      jacob(9,5) = - a_26 * q1(9) 

      jacob(9,7) = - d_4 * q1(9) 

      jacob(9,8) = cj98 - b_6 * q1(9) 

      jacob(9,9) = - ( b_82 * q1(1) + b_4 * q1(2) 
     $   + a_26 * q1(5) + d_4 * q1(7) + b_6 * q1(8)
     $   + cj99 ) 

      jacob(9,10) = cj910 + b_3 * q1(1)  

c... df(no2)/dqi

      jacob(10,1) = b_82 * q1(9) + cj101 - ( b_3 + b_81 ) * q1(10) 

      jacob(10,2) = b_4 * q1(9) - b_9 * q1(10) 

      jacob(10,4) = cj104 - b_22 * q1(10) 

      jacob(10,5) = a_26 * q1(9) + cj105 - b_23 * q1(10) 

      jacob(10,6) = cj106 

      jacob(10,7) = d_4 * q1(9) - d_31 * q1(10) 

      jacob(10,9) = b_4 * q1(2) + a_26 * q1(5)
     $            + b_82 * q1(1) + d_4 * q1(7) + cj109

      jacob(10,10) = - ( (b_3 + b_81) * q1(1) + b_9 * q1(2)
     $                 + b_22 * q1(4) + b_23 * q1(5) + 
     $                 + d_31 * q1(7) + cj1010 )   

      jacob(10,12) =   cj1012

c... df(hocl)/dqi

      jacob(11,1) = - d_35 * q1(11)

      jacob(11,4) = d_87 * q1(12)
     $            - d_34 * q1(11)

      jacob(11,5) = d_33 * q1(7)

      jacob(11,7) = d_33 * q1(5)  

      jacob(11,11) = - (d_34 * q1(4) + d_35 * q1(1) + cj1111)

      jacob(11,12) = d_87 * q1(4)

c... df(clono2)/dqi

      jacob(12,1) = - d_32 * q1(12)

      jacob(12,4) = - d_87 * q1(12)

      jacob(12,7) = d_31 * q1(10)  

      jacob(12,10) = d_31 * q1(7)

      jacob(12,12) = - (d_32 * q1(1) + d_87 * q1(4) + cj1212)

c... calculate Jacobian Jij = I - dt * dfi/dqj

      jacob = - nrstep * jacob
       
      do nn = 1,nbshort
        jacob(nn,nn) = 1 + jacob(nn,nn)
      enddo

      end subroutine eval_j

