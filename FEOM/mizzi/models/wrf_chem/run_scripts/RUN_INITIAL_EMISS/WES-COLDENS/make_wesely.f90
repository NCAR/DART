
      program make_wes
!------------------------------------------------------------------------------------------
!	... read mozart netcdf wesely drydep file and interpolate
!           to specific wrf-chem horizontal grid
!------------------------------------------------------------------------------------------

      use mo_calendar, only : caldayr

      implicit none

      include "netcdf.inc"

!------------------------------------------------------------------------------------------
!	... Local variables
!------------------------------------------------------------------------------------------
      integer :: domain
      integer :: domains = 1
      integer :: astat, istat
      integer :: west_east, west_east_id
      integer :: south_north, south_north_id
      integer :: nlats, npfts, ntimes
      integer :: wrf_ntimes, times_id
      integer :: i, j, n
      integer :: inp_ncid
      integer :: wrf_ncid
      integer :: dimid
      integer :: varid
      integer :: jy(0:1)
      integer, allocatable :: seasonal_pft(:,:,:,:)
      
      real :: by(0:1)
      real, allocatable :: wrf_lats(:,:,:)
      real, allocatable :: wrf_lons(:,:,:)
      real, allocatable :: usgs_seasonal_pft(:,:,:)
      real, allocatable :: lats(:)

      character(len=132) :: filespec
      character(len=80) :: pft_dir = '.'
      character(len=80) :: wrf_dir = '.'
      character(len=80) :: pft_flnm = 'season_wes_usgs.nc'
      character(len=80) :: filename
      character(len=80) :: err_msg
      character(len=80) :: attribute
      character(len=3)  :: numc

      namelist /control/ pft_dir, wrf_dir, pft_flnm, domains

!-----------------------------------------------------------------
!     	... read control variables
!-----------------------------------------------------------------
      read(*,nml=control,iostat=istat)
      if( istat /= 0 ) then
        write(*,*) 'make_wes: failed to read namelist; error = ',istat
        stop
      end if
!-----------------------------------------------------------------
!     	... check namelist inputs
!-----------------------------------------------------------------
      if( domains < 1 ) then
        write(*,*) 'make_wes: domains must be >= 1'
        stop
      end if

!-----------------------------------------------------------------------
!  	... loop over domains
!-----------------------------------------------------------------------
domain_loop : &
      do domain = 1,domains
        write(numc,'(i3)') 100+domain
!-----------------------------------------------------------------------
!  	... open wrfinput file
!-----------------------------------------------------------------------
        filename = 'wrfinput_d' // numc(2:3)
        filespec = trim(wrf_dir) // '/' // trim(filename)
        err_msg = 'make_wes: Failed to open ' // trim(filespec)
        call handle_ncerr( nf_open( trim(filespec), nf_nowrite, wrf_ncid ), trim(err_msg) )
!------------------------------------------------------------------------------
!       ... get wrf dimensions
!------------------------------------------------------------------------------
        err_msg = 'make_wes: failed to get wrf lat id'
        call handle_ncerr( nf_inq_dimid( wrf_ncid, 'south_north', south_north_id ), trim(err_msg) )
        err_msg = 'make_wes: failed to get wrf lat dimension'
        call handle_ncerr( nf_inq_dimlen( wrf_ncid, south_north_id, south_north ), trim(err_msg) )
        err_msg = 'make_wes: failed to get wrf lon id'
        call handle_ncerr( nf_inq_dimid( wrf_ncid, 'west_east', west_east_id ), trim(err_msg) )
        err_msg = 'make_wes: failed to get wrf lon dimension'
        call handle_ncerr( nf_inq_dimlen( wrf_ncid, west_east_id, west_east ), trim(err_msg) )
        err_msg = 'make_wes: failed to get wrf Time id'
        call handle_ncerr( nf_inq_dimid( wrf_ncid, 'Time', dimid ), trim(err_msg) )
        err_msg = 'make_wes: failed to get wrf lon dimension'
        call handle_ncerr( nf_inq_dimlen( wrf_ncid, dimid, wrf_ntimes ), trim(err_msg) )
!------------------------------------------------------------------------------
!       ... allocate wrf arrays
!------------------------------------------------------------------------------
        allocate( wrf_lats(west_east,south_north,wrf_ntimes),stat=astat )
        if( astat /= nf_noerr) then 
           write(*,*) 'make_wes: failed to allocate wrf_lats ; error = ',astat
           stop
        end if
        allocate( wrf_lons(west_east,south_north,wrf_ntimes),stat=astat )
        if( astat /= nf_noerr) then 
           write(*,*) 'make_wes: failed to allocate wrf_lons ; error = ',astat
           stop
        end if
!------------------------------------------------------------------------------
!       ... read wrf latitudes
!------------------------------------------------------------------------------
        err_msg = 'make_wes: failed to get wrf lat variable id'
        call handle_ncerr( nf_inq_varid( wrf_ncid, 'XLAT', varid ), trim(err_msg) )
        err_msg = 'make_wes: failed to read wrf lat variable'
        call handle_ncerr( nf_get_var_real( wrf_ncid, varid, wrf_lats ), trim(err_msg) )
!------------------------------------------------------------------------------
!       ... read wrf longitudes
!------------------------------------------------------------------------------
        err_msg = 'make_wes: failed to get wrf lon variable id'
        call handle_ncerr( nf_inq_varid( wrf_ncid, 'XLONG', varid ), trim(err_msg) )
        err_msg = 'make_wes: failed to read wrf lon variable'
        call handle_ncerr( nf_get_var_real( wrf_ncid, varid, wrf_lons ), trim(err_msg) )
!-----------------------------------------------------------------------
!     	... close wrf template file
!-----------------------------------------------------------------------
        call handle_ncerr( nf_close( wrf_ncid ), 'make_wes: Failed to close netcdf file ' )

initial_domain : &
        if( domain == 1 ) then
!-----------------------------------------------------------------------
!  	... open and read the dry dep seasonal file
!-----------------------------------------------------------------------
          err_msg = 'make_wes: Failed to open ' // trim(pft_dir) // '/' // trim(pft_flnm)
          call handle_ncerr( nf_open( trim(pft_dir)//'/'//trim(pft_flnm), nf_noclobber, inp_ncid ), trim(err_msg) )
!------------------------------------------------------------------------------
!       ... get input dimensions
!------------------------------------------------------------------------------
          err_msg = 'make_wes: failed to get lat dimension id'
          call handle_ncerr( nf_inq_dimid( inp_ncid, 'lat', dimid ), trim(err_msg) )
          err_msg = 'make_wes: failed to get lat dimension'
          call handle_ncerr( nf_inq_dimlen( inp_ncid, dimid, nlats ), trim(err_msg) )
          err_msg = 'make_wes: failed to get time dimension id'
          call handle_ncerr( nf_inq_dimid( inp_ncid, 'time', dimid ), trim(err_msg) )
          err_msg = 'make_wes: failed to get time dimension'
          call handle_ncerr( nf_inq_dimlen( inp_ncid, dimid, ntimes ), trim(err_msg) )
          err_msg = 'make_wes: failed to get pft dimension id'
          call handle_ncerr( nf_inq_dimid( inp_ncid, 'pft', dimid ), trim(err_msg) )
          err_msg = 'make_wes: failed to get pft dimension'
          call handle_ncerr( nf_inq_dimlen( inp_ncid, dimid, npfts ), trim(err_msg) )
      
!------------------------------------------------------------------------------
!       ... allocate input arrays
!------------------------------------------------------------------------------
          allocate( lats(nlats),stat=astat )
          if( astat /= nf_noerr) then 
            write(*,*) 'make_wes: failed to allocate lats ; error = ',astat
            stop
          end if
          allocate( usgs_seasonal_pft(nlats,npfts,ntimes),stat=astat )
          if( astat /= nf_noerr) then 
            write(*,*) 'make_wes: failed to allocate usgs_seasonal_pft; error = ',astat
            stop
          end if

!------------------------------------------------------------------------------
!       ... read input variables
!------------------------------------------------------------------------------
          err_msg = 'make_wes: failed to get lat variable id'
          call handle_ncerr( nf_inq_varid( inp_ncid, 'lat', varid ), trim(err_msg) )
          err_msg = 'make_wes: failed to read lat variable'
          call handle_ncerr( nf_get_var_real( inp_ncid, varid, lats ), trim(err_msg) )
          err_msg = 'make_wes: failed to get usgs seasonal pft variable id'
          call handle_ncerr( nf_inq_varid( inp_ncid, 'season_wes', varid ), trim(err_msg) )
          err_msg = 'make_wes: failed to read usgs seasonal_wes variable'
          call handle_ncerr( nf_get_var_real( inp_ncid, varid, usgs_seasonal_pft ), trim(err_msg) )
!-----------------------------------------------------------------------
!     	... close input file
!-----------------------------------------------------------------------
          call handle_ncerr( nf_close( inp_ncid ), 'make_wes: Failed to close netcdf file ' )
!-----------------------------------------------------------------------
!     	... min,max input dataset diagnostic
!-----------------------------------------------------------------------
          write(*,*) ' '
          write(*,*) 'make_wes: min,max usgs seasonal wes = ',minval(usgs_seasonal_pft),maxval(usgs_seasonal_pft)
          write(*,*) ' '
        endif initial_domain
!------------------------------------------------------------------------------
!       ... allocate wrf data arrays
!------------------------------------------------------------------------------
        allocate( seasonal_pft(west_east,south_north,npfts,ntimes),stat=astat )
        if( astat /= nf_noerr) then 
          write(*,*) 'make_wes: failed to allocate seasonal_pft; error = ',astat
          stop
        end if
!------------------------------------------------------------------------------
!       ... horizontal interpolation from input to wrf grid
!------------------------------------------------------------------------------
        do j = 1,south_north
          do i = 1,west_east
            do n = 1,nlats
              if( wrf_lats(i,j,1) < lats(n) ) then
                exit
              end if
            end do
            jy(0) = min( nlats-1,max( n-1,1 ) )
            jy(1) = jy(0) + 1
            by(0) = (wrf_lats(i,j,1) - lats(jy(0)))/(lats(jy(1)) - lats(jy(0)))
            by(1) = 1. - by(0)
            seasonal_pft(i,j,:,:) = nint( usgs_seasonal_pft(jy(0),:,:)*by(1) &
                                          + usgs_seasonal_pft(jy(1),:,:)*by(0) )
          end do
        end do

!-----------------------------------------------------------------------
!     	... min,max wrf dataset diagnostic
!-----------------------------------------------------------------------
        write(*,*) ' '
        write(*,*) 'make_wes: min,max seasonal pft = ',minval(seasonal_pft),maxval(seasonal_pft)
        write(*,*) ' '
!------------------------------------------------------------------------------
!       ... create wrf data file
!------------------------------------------------------------------------------
        filename = 'wrf_season_wes_usgs_d' // numc(2:3) // '.nc'
        err_msg = 'make_wes: failed to put create file ' // trim(filename)
        call handle_ncerr( nf_create( trim( filename ), nf_clobber, wrf_ncid ), trim(err_msg) )
!-----------------------------------------------------------------------
!     	... define dimensions
!-----------------------------------------------------------------------
        err_msg = 'make_wes: failed to time dimension'
        call handle_ncerr( nf_def_dim( wrf_ncid, 'months', ntimes, times_id ), trim(err_msg) )
        err_msg = 'make_wes: failed to create lat dimension'
        call handle_ncerr( nf_def_dim( wrf_ncid, 'south_north', south_north, south_north_id ), trim(err_msg) )
        err_msg = 'make_wes: failed to create lon dimension'
        call handle_ncerr( nf_def_dim( wrf_ncid, 'west_east', west_east, west_east_id ), trim(err_msg) )
        err_msg = 'make_wes: failed to create pft dimension'
        call handle_ncerr( nf_def_dim( wrf_ncid, 'npft', npfts, dimid ), trim(err_msg) )
!-----------------------------------------------------------------------
!     	... define the variables
!-----------------------------------------------------------------------
        err_msg = 'make_wes: failed to define lat variable'
        call handle_ncerr( nf_def_var( wrf_ncid, 'XLAT', nf_float, 2, (/ west_east_id, south_north_id /), varid ), trim(err_msg) )
        err_msg = 'make_wes: failed to define lon variable'
        call handle_ncerr( nf_def_var( wrf_ncid, 'XLONG', nf_float, 2, (/ west_east_id, south_north_id /), varid ), trim(err_msg) )
        err_msg = 'make_wes: failed to define seasonal_wes variable'
        call handle_ncerr( nf_def_var( wrf_ncid, 'seasonal_wes', nf_int, 4, &
                                     (/ west_east_id, south_north_id, dimid, times_id /), varid ), trim(err_msg) )
!------------------------------------------------------------------------------
!       ... take wrf file out of define mode
!------------------------------------------------------------------------------
        err_msg = 'make_wes: failed to take wrf file out of define mode'
        call handle_ncerr( nf_enddef( wrf_ncid ), trim(err_msg) )

!-----------------------------------------------------------------------
!     	... write new wrf variables
!-----------------------------------------------------------------------
        err_msg = 'make_wes: failed to get lat id'
        call handle_ncerr( nf_inq_varid( wrf_ncid, 'XLAT', varid ), trim(err_msg) )
        err_msg = 'make_wes: failed to write lat variable'
        call handle_ncerr( nf_put_var_real( wrf_ncid, varid, wrf_lats(:,:,1) ), trim(err_msg) )
        err_msg = 'make_wes: failed to get lon id'
        call handle_ncerr( nf_inq_varid( wrf_ncid, 'XLONG', varid ), trim(err_msg) )
        err_msg = 'make_wes: failed to write lon variable'
        call handle_ncerr( nf_put_var_real( wrf_ncid, varid, wrf_lons(:,:,1) ), trim(err_msg) )
        err_msg = 'make_wes: failed to get seasonal id'
        call handle_ncerr( nf_inq_varid( wrf_ncid, 'seasonal_wes', varid ), trim(err_msg) )
        err_msg = 'make_wes: failed to write seasonal levels'
        call handle_ncerr( nf_put_var_int( wrf_ncid, varid, seasonal_pft ), trim(err_msg) )
!-----------------------------------------------------------------------
!     	... close wrf pft file
!-----------------------------------------------------------------------
        call handle_ncerr( nf_close( wrf_ncid ), 'make_wes: Failed to close netcdf file ' )
!-----------------------------------------------------------------------
!     	... deallocate work arrays
!-----------------------------------------------------------------------
        deallocate( wrf_lats, wrf_lons, seasonal_pft )
        write(*,*) 'make_wes: created file ',trim(filename)
      end do domain_loop

      deallocate( lats, usgs_seasonal_pft )

      write(*,*) ' '
      write(*,*) '============================================='
      write(*,*) 'make_wes: completed successfully'
      write(*,*) '============================================='

      end program make_wes

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
