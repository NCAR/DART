
      program surfdata_xform
!------------------------------------------------------------------------------------------
!	... Reorder megan dataset coordinates and data
!           for compatability with WRF
!------------------------------------------------------------------------------------------

      implicit none

      include "netcdf.inc"

!------------------------------------------------------------------------------------------
!	... Local variables
!------------------------------------------------------------------------------------------
      integer :: astat
      integer :: i, j, n
      integer :: dimid, varid
      integer :: ncid
      integer :: nshift
      integer :: nlon_srf
      integer :: nlat_srf
      integer :: ntime_srf
      real, allocatable :: srf_lons(:)
      real, allocatable :: srf_data(:,:,:)
      character(len=32) :: var_name
      character(len=80) :: base_name
      character(len=80) :: filename
      character(len=80) :: err_string

!-----------------------------------------------------------------------
!  	... get filename
!-----------------------------------------------------------------------
      write(*,*) 'Enter filename'
      read(*,*) base_name
!-----------------------------------------------------------------------
!  	... Open "raw" surfdata netcdf file
!-----------------------------------------------------------------------
      filename = trim( base_name ) // '.nc'
      err_string = 'surfdata_xform: Failed to open ' // trim( filename )
      call handle_ncerr( nf_open( trim(filename), nf_write, ncid ), err_string )
      write(*,*) 'Opened file ',trim(filename)
!-----------------------------------------------------------------------
!     	... get surface data dimensions
!-----------------------------------------------------------------------
      err_string = 'Failed to get lon dimension id'
      call handle_ncerr( nf_inq_dimid( ncid, 'lon', dimid ), err_string )
      err_string = 'Failed to get lon dimension'
      call handle_ncerr( nf_inq_dimlen( ncid, dimid, nlon_srf ), err_string )
      err_string = 'Failed to get lat dimension id'
      call handle_ncerr( nf_inq_dimid( ncid, 'lat', dimid ), err_string )
      err_string = 'Failed to get lat dimension'
      call handle_ncerr( nf_inq_dimlen( ncid, dimid, nlat_srf ), err_string )
      err_string = 'Failed to get time dimension id'
      call handle_ncerr( nf_inq_dimid( ncid, 'time', dimid ), err_string )
      err_string = 'Failed to get lat dimension'
      call handle_ncerr( nf_inq_dimlen( ncid, dimid, ntime_srf ), err_string )
!-----------------------------------------------------------------------
!     	... allocate arrays
!-----------------------------------------------------------------------
      allocate( srf_lons(nlon_srf),stat=astat )
      if( astat /= 0 ) then
         write(*,*) 'surfdata_xform: failed to allocate megan_lats; error = ',astat
         stop 'Allocation error'
      end if
      allocate( srf_data(nlon_srf,nlat_srf,ntime_srf),stat=astat )
      if( astat /= 0 ) then
         write(*,*) 'surfdata_xform: failed to allocate megan_var; error = ',astat
         stop 'Allocation error'
      end if
!---------------------------------------------------------------------
!   	... read srf longitudes
!---------------------------------------------------------------------
      err_string = 'Failed to get lon variable id'
      call handle_ncerr( nf_inq_varid( ncid, 'lon', varid ), err_string )
      err_string = 'Failed to read lon variable'
      call handle_ncerr( nf_get_var_real( ncid, varid, srf_lons ), err_string )
!---------------------------------------------------------------------
!   	... reorder longitudes
!---------------------------------------------------------------------
      do nshift = 1,nlon_srf
         if( srf_lons(nshift) < 0. ) then
            exit
         end if
      end do
      nshift = nshift - 1
      srf_lons(:) = cshift( srf_lons(:), nshift )

!---------------------------------------------------------------------
!   	... write longitudes
!---------------------------------------------------------------------
      err_string = 'Failed to get lon variable id'
      call handle_ncerr( nf_inq_varid( ncid, 'lon', varid ), err_string )
      err_string = 'Failed to write lon variable'
      call handle_ncerr( nf_put_var_real( ncid, varid, srf_lons ), err_string )

!---------------------------------------------------------------------
!   	... read surface data
!---------------------------------------------------------------------
      var_name   = trim( base_name ) // '_AVE'
      err_string = 'Failed to get ' // trim( var_name ) // ' variable id'
      call handle_ncerr( nf_inq_varid( ncid, trim( var_name ), varid ), err_string )
      err_string = 'Failed to read ' // trim( var_name ) // ' data variable'
      call handle_ncerr( nf_get_var_real( ncid, varid, srf_data ), err_string )

!---------------------------------------------------------------------
!   	... reorder data
!---------------------------------------------------------------------
      do n = 1,ntime_srf
         do j = 1,nlat_srf
            srf_data(:,j,n) = cshift( srf_data(:,j,n), nshift )
         end do
      end do

!---------------------------------------------------------------------
!   	... write surface data
!---------------------------------------------------------------------
      err_string = 'Failed to get ' // trim( var_name ) // ' variable id'
      call handle_ncerr( nf_inq_varid( ncid, trim( var_name ), varid ), err_string )
      err_string = 'Failed to write ' // trim( var_name ) // ' data variable'
      call handle_ncerr( nf_put_var_real( ncid, varid, srf_data ), err_string )
!-----------------------------------------------------------------------
!     	... Close file
!-----------------------------------------------------------------------
      err_string = 'surfdata_xform: Failed to close ' // trim( filename )
      call handle_ncerr( nf_close( ncid ), err_string )

      write(*,*) 'surfdata_xform: Successfully converted surfdata file'

      end program surfdata_xform

      subroutine handle_ncerr( ret, mes )
!-----------------------------------------------------------------------
! 	... check netcdf library function return code.  if error detected 
!           issue error message then abort.
!-----------------------------------------------------------------------

      implicit none

      include "netcdf.inc"

!-----------------------------------------------------------------------
! 	... dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in) :: ret            ! return code from netcdf library routine
      character(len=*), intent(in) :: mes   ! message to be printed if error detected

      if( ret /= nf_noerr ) then
         write(*,*) mes
         write(*,*) nf_strerror( ret )
         stop
      end if

      end subroutine handle_ncerr
