
      program make_exo_coldens
!------------------------------------------------------------------------------------------
!	... read mozart netcdf exo column density file and interpolate
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
      integer :: ierr
      integer :: west_east, west_east_id
      integer :: south_north, south_north_id
      integer :: wrf_ntimes
      integer :: moz_nlats
      integer :: moz_nlevs
      integer :: moz_nmonths
      integer :: ndays_id
      integer :: lev_id
      integer :: i, j, n
      integer :: moz_ncid
      integer :: wrf_ncid
      integer :: dimid
      integer :: varid
      integer :: dims(4)
      integer :: jy(0:1)
      integer, allocatable :: dates(:)
      
      real    :: by(0:1)
      real, allocatable :: wrf_lats(:,:,:)
      real, allocatable :: wrf_lons(:,:,:)
      real, allocatable :: wrf_o3_coldens(:,:,:,:)
      real, allocatable :: wrf_o2_coldens(:,:,:,:)

      real, allocatable :: moz_lats(:)
      real, allocatable :: moz_levs(:)
      real, allocatable :: moz_o3_coldens(:,:,:)
      real, allocatable :: moz_o2_coldens(:,:,:)
      real, allocatable :: day_of_year(:)

      character(len=132) :: cmd
      character(len=132) :: filespec
      character(len=80) :: exo_dir = '.'
      character(len=80) :: exo_flnm = 'exo_coldens.nc'
      character(len=80) :: wrf_dir = '.'
      character(len=80) :: filename
      character(len=80) :: err_msg
      character(len=80) :: attribute
      character(len=3)  :: numc
      logical :: there


      namelist /control/ exo_dir, wrf_dir, exo_flnm, domains

!-----------------------------------------------------------------
!     	... read control variables
!-----------------------------------------------------------------
      read(*,nml=control,iostat=istat)
      if( istat /= 0 ) then
        write(*,*) 'make_exo_coldens: failed to read namelist; error = ',istat
        stop
      end if
!-----------------------------------------------------------------
!     	... check namelist inputs
!-----------------------------------------------------------------
      if( domains < 1 ) then
        write(*,*) 'make_exo_coldens: domains must be >= 1'
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
        err_msg = 'make_exo_coldens: Failed to open ' // trim(filespec)
        call handle_ncerr( nf_open( trim(filespec), nf_nowrite, wrf_ncid ), trim(err_msg) )
!------------------------------------------------------------------------------
!       ... get wrf dimensions
!------------------------------------------------------------------------------
        err_msg = 'make_exo_coldens: failed to get wrf lat id'
        call handle_ncerr( nf_inq_dimid( wrf_ncid, 'south_north', south_north_id ), trim(err_msg) )
        err_msg = 'make_exo_coldens: failed to get wrf lat dimension'
        call handle_ncerr( nf_inq_dimlen( wrf_ncid, south_north_id, south_north ), trim(err_msg) )
        err_msg = 'make_exo_coldens: failed to get wrf lon id'
        call handle_ncerr( nf_inq_dimid( wrf_ncid, 'west_east', west_east_id ), trim(err_msg) )
        err_msg = 'make_exo_coldens: failed to get wrf lon dimension'
        call handle_ncerr( nf_inq_dimlen( wrf_ncid, west_east_id, west_east ), trim(err_msg) )
        err_msg = 'make_exo_coldens: failed to get wrf Time id'
        call handle_ncerr( nf_inq_dimid( wrf_ncid, 'Time', dimid ), trim(err_msg) )
        err_msg = 'make_exo_coldens: failed to get wrf lon dimension'
        call handle_ncerr( nf_inq_dimlen( wrf_ncid, dimid, wrf_ntimes ), trim(err_msg) )
!------------------------------------------------------------------------------
!       ... allocate wrf arrays
!------------------------------------------------------------------------------
        allocate( wrf_lats(west_east,south_north,wrf_ntimes),stat=astat )
        if( astat /= nf_noerr) then 
           write(*,*) 'make_exo_coldens: failed to allocate wrf_lats ; error = ',astat
           stop
        end if
        allocate( wrf_lons(west_east,south_north,wrf_ntimes),stat=astat )
        if( astat /= nf_noerr) then 
           write(*,*) 'make_exo_coldens: failed to allocate wrf_lons ; error = ',astat
           stop
        end if
!------------------------------------------------------------------------------
!       ... read wrf latitudes
!------------------------------------------------------------------------------
        err_msg = 'make_exo_coldens: failed to get wrf lat variable id'
        call handle_ncerr( nf_inq_varid( wrf_ncid, 'XLAT', varid ), trim(err_msg) )
        err_msg = 'make_exo_coldens: failed to read wrf lat variable'
        call handle_ncerr( nf_get_var_real( wrf_ncid, varid, wrf_lats ), trim(err_msg) )
!------------------------------------------------------------------------------
!       ... read wrf longitudes
!------------------------------------------------------------------------------
        err_msg = 'make_exo_coldens: failed to get wrf lon variable id'
        call handle_ncerr( nf_inq_varid( wrf_ncid, 'XLONG', varid ), trim(err_msg) )
        err_msg = 'make_exo_coldens: failed to read wrf lon variable'
        call handle_ncerr( nf_get_var_real( wrf_ncid, varid, wrf_lons ), trim(err_msg) )
!-----------------------------------------------------------------------
!     	... close wrfinput file
!-----------------------------------------------------------------------
        call handle_ncerr( nf_close( wrf_ncid ), 'make_exo_coldens: Failed to close netcdf file ' )

initial_domain : &
        if( domain == 1 ) then
!-----------------------------------------------------------------------
!  	... open and read the mozart column density file
!-----------------------------------------------------------------------
          filespec = trim(exo_dir) // '/' // trim(exo_flnm)
          err_msg = 'make_exo_coldens: Failed to open ' // trim(filespec)
          call handle_ncerr( nf_open( trim(filespec), nf_noclobber, moz_ncid ), trim(err_msg) )
!------------------------------------------------------------------------------
!       ... get mozart dimensions
!------------------------------------------------------------------------------
          err_msg = 'make_exo_coldens: failed to get moz lat id'
          call handle_ncerr( nf_inq_dimid( moz_ncid, 'lat', dimid ), trim(err_msg) )
          err_msg = 'make_exo_coldens: failed to get moz lat dimension'
          call handle_ncerr( nf_inq_dimlen( moz_ncid, dimid, moz_nlats ), trim(err_msg) )
          err_msg = 'make_exo_coldens: failed to get moz lev id'
          call handle_ncerr( nf_inq_dimid( moz_ncid, 'lev', dimid ), trim(err_msg) )
          err_msg = 'make_exo_coldens: failed to get moz lev dimension'
          call handle_ncerr( nf_inq_dimlen( moz_ncid, dimid, moz_nlevs ), trim(err_msg) )
          err_msg = 'make_exo_coldens: failed to get moz month id'
          call handle_ncerr( nf_inq_dimid( moz_ncid, 'month', dimid ), trim(err_msg) )
          err_msg = 'make_exo_coldens: failed to get moz month dimension'
          call handle_ncerr( nf_inq_dimlen( moz_ncid, dimid, moz_nmonths ), trim(err_msg) )
          if( moz_nmonths /= 12 ) then
            write(*,*) 'make_exo_coldens: exo coldens is not annual period'
            stop
          end if
!------------------------------------------------------------------------------
!       ... allocate mozart arrays
!------------------------------------------------------------------------------
          allocate( moz_lats(moz_nlats),stat=astat )
          if( astat /= nf_noerr) then 
            write(*,*) 'make_exo_coldens: failed to allocate moz_lats ; error = ',astat
            stop
          end if
          allocate( moz_levs(moz_nlevs),stat=astat )
          if( astat /= nf_noerr) then 
            write(*,*) 'make_exo_coldens: failed to allocate moz_levs ; error = ',astat
            stop
          end if
          allocate( moz_o3_coldens(moz_nlats,moz_nlevs,moz_nmonths),stat=astat )
          if( astat /= nf_noerr) then 
            write(*,*) 'make_exo_coldens: failed to allocate moz o3 col; error = ',astat
            stop
          end if
          allocate( moz_o2_coldens(moz_nlats,moz_nlevs,moz_nmonths),stat=astat )
          if( astat /= nf_noerr) then 
            write(*,*) 'make_exo_coldens: failed to allocate moz o2 col; error = ',astat
            stop
          end if
          allocate( dates(moz_nmonths),day_of_year(moz_nmonths),stat=astat )
          if( astat /= nf_noerr) then 
            write(*,*) 'make_exo_coldens: failed to allocate dates,day_of_year; error = ',astat
            stop
          end if
!------------------------------------------------------------------------------
!       ... read mozart variables
!------------------------------------------------------------------------------
          err_msg = 'make_exo_coldens: failed to get lat variable id'
          call handle_ncerr( nf_inq_varid( moz_ncid, 'lat', varid ), trim(err_msg) )
          err_msg = 'make_exo_coldens: failed to read lat variable'
          call handle_ncerr( nf_get_var_real( moz_ncid, varid, moz_lats ), trim(err_msg) )
          err_msg = 'make_exo_coldens: failed to get level variable id'
          call handle_ncerr( nf_inq_varid( moz_ncid, 'lev', varid ), trim(err_msg) )
          err_msg = 'make_exo_coldens: failed to read cross_sections variable'
          call handle_ncerr( nf_get_var_real( moz_ncid, varid, moz_levs ), trim(err_msg) )
          err_msg = 'make_exo_coldens: failed to get moz o3 col variable id'
          call handle_ncerr( nf_inq_varid( moz_ncid, 'O3_column_density', varid ), trim(err_msg) )
          err_msg = 'make_exo_coldens: failed to read moz o3 col variable'
          call handle_ncerr( nf_get_var_real( moz_ncid, varid, moz_o3_coldens ), trim(err_msg) )
          err_msg = 'make_exo_coldens: failed to get moz o2 col variable id'
          call handle_ncerr( nf_inq_varid( moz_ncid, 'O2_column_density', varid ), trim(err_msg) )
          err_msg = 'make_exo_coldens: failed to read moz o2 col variable'
          call handle_ncerr( nf_get_var_real( moz_ncid, varid, moz_o2_coldens ), trim(err_msg) )
!-----------------------------------------------------------------------
!	... initialize the monthly day of year times
!-----------------------------------------------------------------------
          dates(:) = (/ 116, 214, 316, 415,  516,  615, &
                        716, 816, 915, 1016, 1115, 1216 /)
          do n = 1,moz_nmonths
            day_of_year(n) = caldayr( dates(n), 0 )
          end do
          deallocate( dates )
!-----------------------------------------------------------------------
!     	... close mozart file
!-----------------------------------------------------------------------
          call handle_ncerr( nf_close( moz_ncid ), 'make_exo_coldens: Failed to close netcdf file ' )
        endif initial_domain

!------------------------------------------------------------------------------
!       ... allocate wrf data arrays
!------------------------------------------------------------------------------
        allocate( wrf_o3_coldens(west_east,south_north,moz_nlevs,moz_nmonths),stat=astat )
        if( astat /= nf_noerr) then 
          write(*,*) 'make_exo_coldens: failed to allocate wrf_o3_coldens ; error = ',astat
          stop
        end if
        allocate( wrf_o2_coldens(west_east,south_north,moz_nlevs,moz_nmonths),stat=astat )
        if( astat /= nf_noerr) then 
          write(*,*) 'make_exo_coldens: failed to allocate wrf_o2_coldens ; error = ',astat
          stop
        end if
!------------------------------------------------------------------------------
!       ... horizontal interpolation from mozart to wrf grid
!------------------------------------------------------------------------------
        do j = 1,south_north
          do i = 1,west_east
            do n = 1,moz_nlats
              if( wrf_lats(i,j,1) < moz_lats(n) ) then
                exit
              end if
            end do
            jy(0) = min( moz_nlats-1,max( n-1,1 ) )
            jy(1) = jy(0) + 1
            by(0) = (wrf_lats(i,j,1) - moz_lats(jy(0)))/(moz_lats(jy(1)) - moz_lats(jy(0)))
            by(1) = 1. - by(0)
            wrf_o3_coldens(i,j,:,:) = moz_o3_coldens(jy(0),:,:)*by(1) &
                                    + moz_o3_coldens(jy(1),:,:)*by(0)
            wrf_o2_coldens(i,j,:,:) = moz_o2_coldens(jy(0),:,:)*by(1) &
                                    + moz_o2_coldens(jy(1),:,:)*by(0)
          end do
        end do

!-----------------------------------------------------------------------
!  	... create wrf exo_coldens file
!-----------------------------------------------------------------------
        filename   = 'exo_coldens_d' // numc(2:3)
        err_msg = 'make_exo_coldens: failed to put create file ' // trim(filename)
        call handle_ncerr( nf_create( trim( filename ), nf_clobber, wrf_ncid ), trim(err_msg) )
!-----------------------------------------------------------------------
!     	... define dimensions
!-----------------------------------------------------------------------
        err_msg = 'make_exo_coldens: failed to create lat dimension'
        call handle_ncerr( nf_def_dim( wrf_ncid, 'south_north', south_north, south_north_id ), trim(err_msg) )
        err_msg = 'make_exo_coldens: failed to create lon dimension'
        call handle_ncerr( nf_def_dim( wrf_ncid, 'west_east', west_east, west_east_id ), trim(err_msg) )
        err_msg = 'make_exo_coldens: failed to create ndays_of_year dimension'
        call handle_ncerr( nf_def_dim( wrf_ncid, 'ndays_of_year', moz_nmonths, ndays_id ), trim(err_msg) )
        err_msg = 'make_exo_coldens: failed to create column denisty level dimension'
        call handle_ncerr( nf_def_dim( wrf_ncid, 'coldens_levs', moz_nlevs, lev_id ), trim(err_msg) )
!-----------------------------------------------------------------------
!     	... define the variables
!-----------------------------------------------------------------------
        err_msg = 'make_exo_coldens: failed to define lat variable'
        call handle_ncerr( nf_def_var( wrf_ncid, 'XLAT', nf_float, 2, (/ west_east_id, south_north_id /), varid ), trim(err_msg) )
        err_msg = 'make_exo_coldens: failed to define lon variable'
        call handle_ncerr( nf_def_var( wrf_ncid, 'XLONG', nf_float, 2, (/ west_east_id, south_north_id /), varid ), trim(err_msg) )
        err_msg = 'make_exo_coldens: failed to define column density levels variable'
        call handle_ncerr( nf_def_var( wrf_ncid, 'coldens_levs', nf_float, 1, (/ lev_id /), varid ), trim(err_msg) )
        err_msg = 'make_exo_coldens: failed to define days of year variable'
        call handle_ncerr( nf_def_var( wrf_ncid, 'days_of_year', nf_float, 1, (/ ndays_id /), varid ), trim(err_msg) )
        dims(:) = (/ west_east_id, south_north_id, lev_id, ndays_id /)
        err_msg = 'make_exo_coldens: failed to define o3 column density variable'
        call handle_ncerr( nf_def_var( wrf_ncid, 'o3_column_density', nf_float, 4, dims, varid ), trim(err_msg) )
        err_msg = 'make_exo_coldens: failed to define o2 column density variable'
        call handle_ncerr( nf_def_var( wrf_ncid, 'o2_column_density', nf_float, 4, dims, varid ), trim(err_msg) )
!------------------------------------------------------------------------------
!       ... take wrf exo_coldens file out of define mode
!------------------------------------------------------------------------------
        err_msg = 'make_exo_coldens: failed to take wrf file out of define mode'
        call handle_ncerr( nf_enddef( wrf_ncid ), trim(err_msg) )
!-----------------------------------------------------------------------
!     	... write wrf exo_coldens file variables
!-----------------------------------------------------------------------
        err_msg = 'make_exo_coldens: failed to get lat id'
        call handle_ncerr( nf_inq_varid( wrf_ncid, 'XLAT', varid ), trim(err_msg) )
        err_msg = 'make_exo_coldens: failed to write lat variable'
        call handle_ncerr( nf_put_var_real( wrf_ncid, varid, wrf_lats(:,:,1) ), trim(err_msg) )
        err_msg = 'make_wes: failed to get lon id'
        call handle_ncerr( nf_inq_varid( wrf_ncid, 'XLONG', varid ), trim(err_msg) )
        err_msg = 'make_wes: failed to write lon variable'
        call handle_ncerr( nf_put_var_real( wrf_ncid, varid, wrf_lons(:,:,1) ), trim(err_msg) )
        err_msg = 'make_exo_coldens: failed to get days of year id'
        call handle_ncerr( nf_inq_varid( wrf_ncid, 'days_of_year', varid ), trim(err_msg) )
        err_msg = 'make_exo_coldens: failed to write days of year variable'
        call handle_ncerr( nf_put_var_real( wrf_ncid, varid, day_of_year ), trim(err_msg) )
        err_msg = 'make_exo_coldens: failed to get column density levels id'
        call handle_ncerr( nf_inq_varid( wrf_ncid, 'coldens_levs', varid ), trim(err_msg) )
        err_msg = 'make_exo_coldens: failed to write column density levels'
        call handle_ncerr( nf_put_var_real( wrf_ncid, varid, moz_levs ), trim(err_msg) )
        err_msg = 'make_exo_coldens: failed to get o3 column density id'
        call handle_ncerr( nf_inq_varid( wrf_ncid, 'o3_column_density', varid ), trim(err_msg) )
        err_msg = 'make_exo_coldens: failed to write o3 column density'
        call handle_ncerr( nf_put_var_real( wrf_ncid, varid, wrf_o3_coldens ), trim(err_msg) )
        err_msg = 'make_exo_coldens: failed to get o2 column density id'
        call handle_ncerr( nf_inq_varid( wrf_ncid, 'o2_column_density', varid ), trim(err_msg) )
        err_msg = 'make_exo_coldens: failed to write o2 column density'
        call handle_ncerr( nf_put_var_real( wrf_ncid, varid, wrf_o2_coldens ), trim(err_msg) )
!-----------------------------------------------------------------------
!     	... close wrf exo_coldens file
!-----------------------------------------------------------------------
        call handle_ncerr( nf_close( wrf_ncid ), 'make_exo_coldens: Failed to close netcdf file ' )
        deallocate( wrf_lons, wrf_lats, wrf_o3_coldens, wrf_o2_coldens )
        write(*,*) ' '
        write(*,*) 'make_exo_coldens: created file ',trim(filename)
      end do domain_loop

      deallocate( moz_lats, moz_levs, moz_o3_coldens, moz_o2_coldens, day_of_year )

      write(*,*) ' '
      write(*,*) '============================================='
      write(*,*) 'make_exo_coldens: completed successfully'
      write(*,*) '============================================='

      end program make_exo_coldens

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
