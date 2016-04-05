c----------------------------------------------------------------------
      module chem0d 
c----------------------------------------------------------------------
       
      use params, only : nphot, nbcon

      implicit none

      save

      integer, parameter :: nblong = 13
      integer, parameter :: nbshort = 12

      real :: dtchem       ! chemical time step

c... one element densities
      real :: den, c_n2o, c_ch4, c_h2o, c_noy, c_hno3, 
     $        c_n2o5, c_clx, c_co, c_hcl,
     $        c_h2o2, c_ho2no2, c_h2, c_ch2o,
     $        c_no3, c_o1d, o2, n2, c_n2d  

      real :: c_h, c_oh, c_ho2, c_o, c_o3, c_cl, c_clo, 
     $        c_n, c_no, c_no2, c_hocl, c_clono2

      real :: c_co2

      real :: c_prev(nbcon)

      real, dimension (nphot) :: tj0

c... rate constants
      real :: cj11, cj12, cj18, cj19, cj110,
     $        cj21, cj22, 
     $        cj33, cj34, cj36,
     $        cj41, cj44, cj45, cj46, cj411,
     $        cj51, cj53, cj54, cj55, cj56

      real :: cj61, cj64, cj66, cj611, cj612,
     $        cj74, cj76, cj712, 
     $        cj89, cj98, cj99, cj910,
     $        cj101, cj104, cj105, cj106, cj109, cj1010, cj1012, 
     $        cj1111,
     $        cj1212

      real :: nrstep 

c... one element of 3-d rate arrays

      real :: a_1, a_2, a_5, a_6, a_6b, a_7, a_17, a_19, a_26,
     $        a_27, a_30, a_36, a_81, a_82, a_83

      real :: b_3, b_4, b_6, b_7, b_9, b_12, b_22, b_23, b_24,
     $        b_27, b_28, b_32, b_81, b_82, b_84

      real :: c_2, c_9, d_2, d_3, d_4, d_5, d_6, d_7, d_8,
     $        d_10, d_11, d_31, d_32, d_33, d_34, d_35, d_36,
     $        d_37, d_46, d_47, d_48, d_60, d_61, d_62, d_63,
     $        d_82, d_83, d_84, d_85, d_87

      real :: hk_1, hk_2, hk_3, hk_4, hk_5, hk_21

c... temperature independent rate constant 

      real :: a_1et, a_3et, a_23a, a_23b, a_23c, a_23

      real :: b_7a, b_38, b_39, b_71, b_72, b_73a, b_73b

      real :: c_1, c_1et, c_8, d_90, d_91, d_71, d_72,
     $        d_73, d_74, d_75 

      real :: hk_7

c... convergence parameters
      real :: convlim      
      integer :: niter

      real :: dnmin

      logical :: particle

      contains

c----------------------------------------------------------------------
      subroutine chemsolv( err_flag )
c----------------------------------------------------------------------

      implicit none

      logical :: err_flag  ! error flag for failed convergence

      logical :: nr_convrg, gs_convrg

      real, dimension (nbshort) :: nr_iter, nr_dif
      real, dimension (nblong) :: gs_iter, gs_dif

      err_flag = .false.

      convlim = .1

      call chem_nr( nr_convrg, nr_dif )
      call chem_gs( gs_convrg, gs_dif )

      convlim = 1.e-3

      call chem_nr( nr_convrg, nr_iter )
      call chem_gs( gs_convrg, gs_iter )

      nr_dif = abs( nr_iter - nr_dif )
      gs_dif = abs( gs_iter - gs_dif )

      if ((all( nr_dif .le. nr_iter*convlim*10. )) .and. 
     $    (all( gs_dif .le. gs_iter*convlim*10. ))) then 
        err_flag = (.not. nr_convrg ) .or. (.not. gs_convrg )
        return
      endif

      nr_dif = nr_iter
      gs_dif = gs_iter

      call chem_nr( nr_convrg, nr_iter )
      call chem_gs( gs_convrg, gs_iter )

      nr_dif = abs( nr_iter - nr_dif )
      gs_dif = abs( gs_iter - gs_dif )

      if ((all( nr_dif .le. nr_iter*convlim*10. )) .and. 
     $    (all( gs_dif .le. gs_iter*convlim*10. ))) then 
        err_flag = (.not. nr_convrg ) .or. (.not. gs_convrg )
        return
      endif

      print *, 'chemsolv failed to converge'
      err_flag = .true.

      return

      end subroutine chemsolv

c----------------------------------------------------------------------
      end module chem0d 
c----------------------------------------------------------------------

