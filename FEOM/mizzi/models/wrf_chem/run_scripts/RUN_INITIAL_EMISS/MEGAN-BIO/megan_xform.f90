
      program megan_xform
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
      integer :: i, ii, ij, j, jj, n, nn
      integer :: dimid, varid
      integer :: ncid
      integer :: nlon_megan
      integer :: nlat_megan
      real    :: xx
      real, allocatable :: megan_lats(:)
      real, allocatable :: megan_var(:,:)
      real, allocatable :: wrk(:)
      character(len=80) :: filename
      character(len=80) :: err_string
      character(len=64) :: varname

!-----------------------------------------------------------------------
!  	... Open "raw" megan netcdf file
!-----------------------------------------------------------------------
      filename = ' '
      write(*,*) 'Enter filename'
      read(*,*) filename
      err_string = 'megan_xform: Failed to open ' // trim( filename )
      call handle_ncerr( nf_open( trim(filename), nf_write, ncid ), err_string )
!-----------------------------------------------------------------------
!     	... get megan dimensions
!-----------------------------------------------------------------------
      err_string = 'Failed to get lon dimension id'
      call handle_ncerr( nf_inq_dimid( ncid, 'lon', dimid ), err_string )
      err_string = 'Failed to get lon dimension'
      call handle_ncerr( nf_inq_dimlen( ncid, dimid, nlon_megan ), err_string )
      err_string = 'Failed to get lat dimension id'
      call handle_ncerr( nf_inq_dimid( ncid, 'lat', dimid ), err_string )
      err_string = 'Failed to get lat dimension'
      call handle_ncerr( nf_inq_dimlen( ncid, dimid, nlat_megan ), err_string )
!-----------------------------------------------------------------------
!     	... allocate arrays
!-----------------------------------------------------------------------
      allocate( wrk(nlon_megan),stat=astat )
      if( astat /= 0 ) then
         write(*,*) 'megan_xform: failed to allocate wrk; error = ',astat
         stop 'Allocation error'
      end if
      allocate( megan_lats(nlat_megan),stat=astat )
      if( astat /= 0 ) then
         write(*,*) 'megan_xform: failed to allocate megan_lats; error = ',astat
         stop 'Allocation error'
      end if
      allocate( megan_var(nlon_megan,nlat_megan),stat=astat )
      if( astat /= 0 ) then
         write(*,*) 'megan_xform: failed to allocate megan_var; error = ',astat
         stop 'Allocation error'
      end if
!---------------------------------------------------------------------
!   	... read megan latitude variable
!---------------------------------------------------------------------
      err_string = 'Failed to get lat variable id'
      call handle_ncerr( nf_inq_varid( ncid, 'lat', varid ), err_string )
      err_string = 'Failed to read lat variable'
      call handle_ncerr( nf_get_var_real( ncid, varid, megan_lats ), err_string )
!---------------------------------------------------------------------
!   	... reorder latitudes
!---------------------------------------------------------------------
      do j = 1,nlat_megan/2
         jj = nlat_megan - j + 1
         xx = megan_lats(jj)
         megan_lats(jj) =  megan_lats(j)
         megan_lats(j)  =  xx
      end do

      write(*,*) ' '
      write(*,*) 'megan_xform: transformed latitudes'
      write(*,'(1p5g15.7)') megan_lats(1:5)
      write(*,'(1p5g15.7)') megan_lats(nlat_megan/2-2:nlat_megan/2+2)
      write(*,'(1p5g15.7)') megan_lats(nlat_megan-4:nlat_megan)
      write(*,*) ' '

!---------------------------------------------------------------------
!   	... write megan latitude variable
!---------------------------------------------------------------------
      err_string = 'Failed to get lat variable id'
      call handle_ncerr( nf_inq_varid( ncid, 'lat', varid ), err_string )
      err_string = 'Failed to write lat variable'
      call handle_ncerr( nf_put_var_real( ncid, varid, megan_lats ), err_string )

!---------------------------------------------------------------------
!   	... read megan data
!---------------------------------------------------------------------
      varname = ' '
      write(*,*) 'Enter variable name'
      read(*,*) varname
      err_string = 'Failed to get ' // trim(varname) // ' variable id'
      call handle_ncerr( nf_inq_varid( ncid, trim(varname), varid ), err_string )
      err_string = 'Failed to read ' // trim(varname) // ' variable'
      call handle_ncerr( nf_get_var_real( ncid, varid, megan_var ), err_string )
!---------------------------------------------------------------------
!   	... reorder data
!---------------------------------------------------------------------
      do j = 1,nlat_megan/2
         jj = nlat_megan - j + 1
         wrk(:) = megan_var(:,jj)
         megan_var(:,jj) =  megan_var(:,j)
         megan_var(:,j)  =  wrk(:)
      end do
!---------------------------------------------------------------------
!   	... write megan data
!---------------------------------------------------------------------
      err_string = 'Failed to get ' // trim(varname) // ' variable id'
      call handle_ncerr( nf_inq_varid( ncid, trim(varname), varid ), err_string )
      err_string = 'Failed to write ' // trim(varname) // ' variable'
      call handle_ncerr( nf_put_var_real( ncid, varid, megan_var ), err_string )
!-----------------------------------------------------------------------
!     	... Close file
!-----------------------------------------------------------------------
      err_string = 'megan_xform: Failed to close ' // trim( filename )
      call handle_ncerr( nf_close( ncid ), err_string )

      write(*,*) 'megan_xform: Successfully converted megan file'

      end program megan_xform

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
