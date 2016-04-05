
c----------------------------------------------------------------------
      module output 
c----------------------------------------------------------------------

      use params
      use ncdf

      implicit none

      private
      public :: main_init, main_update, main_close

c... netCDF id
      integer ::  main_id

c... variable ids
      integer :: nstp_id
      integer :: time_id
      integer :: mtime_id

      integer, dimension(nvars) ::  mvar_ids

c... timing variables
      integer :: main_ctr 

      contains
       
c----------------------------------------------------------------------
      subroutine main_init( nc_name )
c----------------------------------------------------------------------

      implicit none     

      character(len=*), intent(in) :: nc_name

      character (LEN=6), dimension(nvars) :: var_names = 
     1  (/ 'un    ', 'vn    ', 'wn    ', 'tn    ', 'N2O   ', 'CH4   ',
     2     'H2O   ', 'O_1D  ', 'HNO3  ', 'N2O5  ', 'H     ', 'OH    ',
     3     'CO    ', 'HCl   ', 'ClONO2', 'HOCl  ', 'H2O2  ', 'HO2   ',
     4     'HO2NO2', 'H2    ', 'CH2O  ', 'O     ', 'O3    ', 'Cl    ',
     5     'ClO   ', 'N     ', 'NO    ', 'NO2   ', 'NO3   ', 'O2    ' /)
      
      character (LEN=6), dimension(nvars) :: var_units = 
     1  (/ 'm/s   ', 'm/s   ', 'm/s   ', 'deg k ', 'cm3-qn', 'cm3-qn',
     2     'cm3-qn', 'cm3-qn', 'cm3-qn', 'cm3-qn', 'cm3-qn', 'cm3-qn',
     3     'cm3-qn', 'cm3-qn', 'cm3-qn', 'cm3-qn', 'cm3-qn', 'cm3-qn',
     4     'cm3-qn', 'cm3-qn', 'cm3-qn', 'cm3-qn', 'cm3-qn', 'cm3-qn',
     5     'cm3-qn', 'cm3-qn', 'cm3-qn', 'cm3-qn', 'cm3-qn', 'cm3-qn' /)

      main_ctr = 0

      call nc_init( nc_name, var_names, var_units, mvar_ids, nvars, 
     $              main_id, nstp_id, time_id, mtime_id )

      call nc_add_global( main_id )

      end subroutine main_init        

c----------------------------------------------------------------------
      subroutine main_update( vars, iyear, doy, utsec )
c----------------------------------------------------------------------

      implicit none

      real, dimension(nz, nx, ny, nvars), intent(in) :: vars 
      integer, intent(in) :: iyear
      integer, intent(in) :: doy
      integer, intent(in) :: utsec
      
      integer, dimension(3) :: mtime

c... increment counter
      main_ctr = main_ctr + 1
 
      print *, 'main_update', iyear, doy, utsec, nstep 

      mtime(1) = iyear
      mtime(2) = doy
      mtime(3) = utsec

      call nc_update( main_id, main_ctr, nstp_id, mtime_id, time_id, 
     $                mtime, nvars, mvar_ids, vars ) 

      end subroutine main_update
       
c----------------------------------------------------------------------
      subroutine main_close 
c----------------------------------------------------------------------

      implicit none
      include 'netcdf.inc'

      integer :: iret

      iret = nf_close(main_id)
      call check_err(iret)

      end subroutine main_close

c----------------------------------------------------------------------
      end module output 
c----------------------------------------------------------------------
