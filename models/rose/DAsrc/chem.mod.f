      
!----------------------------------------------------------------------
      module chem
!----------------------------------------------------------------------
!                                                                     c
!  This module calculates the concentrations of short-lived           c
!  species (assumed to be in photochemical equilibrium) and           c
!  long-lived species.                                                c
!                                                                     c
!  PSC's and heterogeneous reactions are not included                 c
!                                                                     c
!----------------------------------------------------------------------
!
!... constituent mixing ratios (qn)
!
!      1  n2o       2  ch4        3  h2o       4  o1d
!      5  hno3      6  n2o5       7  h         8  oh
!      9  co       10  hcl       11  clono2   12  hocl
!     13  h2o2     14  ho2       15  ho2no2   16  h2
!     17  ch2o     18  o         19  o3       20  cl
!     21  clo      22  n         23  no       24  no2
!     25  no3      26  o2
!
!----------------------------------------------------------------------
!
!... photodissociation rates (tj) 
!
!      1 o2,           2 o2->o1d       3 o3          4 o3->o1d
!      5 ch4           6 no2           7 hno3        8 hocl 
!      9 ho2no2->ho2  10 ho2no2->oh   11 n2o5->no2  12 n2o5->no
!     13 h2o2         14 oclo         15 cl2o2      16 hcl 
!     17 cl2          18 co2          19 n2o        20 brono2
!     21 brcl         22 hobr         23 hno2       24 ch2o->h+ho2+co
!     25 ch2o->h2     26 clono2->clo  27 clono2->cl 28 no  
!     29 no3->no2     30 no3->no      31 h2o->h+oh  32 h2o->h2+o(1d) 
!     33 h2o->2h+o
!
!... note -> branching is already applied
!
!----------------------------------------------------------------------
           
      use params

      implicit none

      save

      integer :: nrates         !  frequency of rate constant update

!... constituent mixing ratios
      real, dimension (nz,nx,ny,nbcon) :: qn1
      real, dimension (nz,ny) :: q_o2, q_n2
      real, dimension (nz,nx,ny) :: hnm, o3t, o2t

!... boundary conditions
      real, dimension (nx,ny,nbcon) :: qlbc, qubc

!... rate coefficients
      real :: a1et, a3et, a23a, a23b, a23c, a23, a23p
      real :: b7a, b38, b39, b71, b72, b73a, b73b
      real :: c1, c8, d35, d73, d74, d75, hk7

      real, dimension (nz,nx,ny) :: a1, a2, a5, a6, a6b, 
     $      a7, a17, a19, a26, a27, a30, a36, a81, a82, a83

      real, dimension (nz,nx,ny) :: b3, b4, b6, b7, b9, 
     $      b12, b22, b23, b24, b27, b28, b32, b81, b82, b84

      real, dimension (nz,nx,ny) :: c2, c9, d0, d2, d3, d4, d5,
     $      d6, d7, d8, d10, d11, d31, d32, d33, d34

      real, dimension (nz,nx,ny) :: d46, d47, d83, d84, d85, d87,
     $      hk1, hk2, hk3, hk4, hk5, hk21

!... photolysis
      real, dimension (nz,nx,ny,nphot) :: tj
      real, dimension (nx,ny) :: sc2d

      contains

!----------------------------------------------------------------------
      subroutine chem3d
!----------------------------------------------------------------------
!
!  Main chemistry driver. Uses Newton-Raphson iteration for non-linear
!  constituents. Gauss-Seidel iteration for remaining constituents
!
!----------------------------------------------------------------------

      use chem0d, only : chemsolv
      use dynam, only : zkm

      implicit none

!...local variables
      logical :: err
      integer :: i, j, k

!... loop over all model grid-points
      do k=1,nz
        do j=1,ny
          do i=1,nx
	    
	    call select0d( k, i, j)

	    call chemsolv(err)              ! solve chemistry at ea. grid-point

	    if (err) then
              print *, (i-1)*360./nx, j, zkm(k)
	    endif

	    call update0d( k, i, j) 

          enddo
        enddo
      enddo
   
      return

!----------------------------------------------------------------------
      end subroutine chem3d
!----------------------------------------------------------------------

!----------------------------------------------------------------------
      end module chem
!----------------------------------------------------------------------
