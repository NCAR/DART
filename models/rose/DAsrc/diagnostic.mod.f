c-------------------------------------------------------------------
      module diagnostic
c-------------------------------------------------------------------

      use params
      use ncdf

      implicit none

      save

      private
      public :: diag_init, diag_update, diag_close 

c... parameters
      integer, parameter ::  ndiag = 6

c... netCDF id
      integer ::  diag_id

c... variable ids
      integer :: nstp_id
      integer :: time_id
      integer :: mtime_id

      integer, dimension(ndiag) ::  dvar_ids

c... timing variables
      integer :: diag_ctr 

      contains
      
c-------------------------------------------------------------------
      subroutine diag_init( nc_name ) 
c-------------------------------------------------------------------

      implicit none     

      character(len=*), intent(in) :: nc_name

      character (LEN=6), dimension(ndiag) :: var_names = 
     1  (/ 'v5577 ', 'vo200 ', 'e_v_oh', 
     2     'v_oh  ', 'chm_ht', 'sol_ht'/)
      
      character (LEN=9), dimension(ndiag) :: var_units = 
     1  (/ 'ph/cm3/s ', 'ph/cm3/s ', 'erg/cm3/s', 
     2     'ph/cm3/s ', 'K/day    ', 'K/day    '/)

      diag_ctr = 0

      call nc_init( nc_name, var_names, var_units, dvar_ids, ndiag,
     $              diag_id, nstp_id, time_id, mtime_id )

      call nc_add_global( diag_id )

      end subroutine diag_init

c-------------------------------------------------------------------
      subroutine diag_update( iyear, doy, utsec )
c-------------------------------------------------------------------

      use params
      use chem
      use phys
      use dynam
      use airglow, only: eoh83, eo200, e5577, oh_ver 

      implicit none
      
      integer, intent(in) :: iyear
      integer, intent(in) :: doy
      integer, intent(in) :: utsec

      real, dimension(nz, nx, ny, ndiag) :: vars 
      integer, dimension(3) :: mtime
      
c... airglow variables
      real, dimension(nz,nx,ny) :: c_o, c_o2, c_n2, c_h, c_o3
      real, dimension(nz,nx,ny) :: tn
      real, dimension(nz,nx,ny) :: v5577, vo200, e_v_oh, v_oh 

c... local variables
      integer :: i, j, k

      print *, 'writing diagnostic variables: '
      print *, 'diag_update', iyear, doy, utsec, nstep 

c... increment counter
      diag_ctr = diag_ctr + 1

      mtime(1) = iyear
      mtime(2) = doy
      mtime(3) = utsec

      c_o = qn1(:,:,:,18) * hnm(:,:,:)
      c_h = qn1(:,:,:,7) * hnm(:,:,:)
      c_o2 = qn1(:,:,:,26) * hnm(:,:,:)
      c_o3 = qn1(:,:,:,19) * hnm(:,:,:)
      c_n2 = hnm(:,:,:) - c_o2(:,:,:) - c_o2(:,:,:)

c... total temperature
      do k=1,nz
         do i=1,nx
            do j=1,ny
               tn(k,i,j) = tn1(k,i,j) + tref(k)
            end do
         end do
      end do

      call e5577( tn, c_o, c_o2, c_n2, v5577 )

      call eo200( tn, c_o, c_o2, c_n2, vo200 )

!!    call eoh83( tn, c_o, c_o2, c_n2, voh83 )

      call oh_ver( tn, c_o, c_o2, c_n2, c_o3, c_h, v_oh, e_v_oh )

      vars(:,:,:,1) = v5577 
      vars(:,:,:,2) = vo200
      vars(:,:,:,3) = e_v_oh 
      vars(:,:,:,4) = v_oh 

      vars(:,:,:,5) = ch_heat * 86400.
      vars(:,:,:,6) = solht * 86400.

      call nc_update( diag_id, diag_ctr, nstp_id, mtime_id, time_id,
     $                mtime, ndiag, dvar_ids, vars) 

      end subroutine diag_update
       
c----------------------------------------------------------------------
      subroutine diag_close 
c----------------------------------------------------------------------

      implicit none
      include 'netcdf.inc'

      integer :: iret

      iret = nf_close(diag_id)
      call check_err(iret)

      end subroutine diag_close

c-------------------------------------------------------------------
      end module diagnostic
c-------------------------------------------------------------------
