
   module data_file_utils

   use anthro_types, only : linsize

   implicit none

!---------------------------------------------------------------------
!  include files
!---------------------------------------------------------------------
   include 'netcdf.inc'

   private

   public :: data_file_init
   public :: get_src_time_ndx
   public :: read_src_data
   public :: tinterp_src_data
   public :: anthro_dir
   public :: data_yrs_offset

   integer                :: data_yrs_offset     = 0         ! years( + or - )
   character(len=linsize) :: anthro_dir

   contains

   subroutine data_file_init( data_file, start_output, cat_var_prefix, cat_var_suffix, &
                              domain, dx, ide, jde )
!---------------------------------------------------------------------
!  initialize data file type
!---------------------------------------------------------------------

   use constants_module, only : rad_per_deg, earth_radius_m
   use mapper_types
   use area_mapper,      only : area_interp_init
   use area_mapper,      only : xlong => lon, xlat => lat
   use anthro_types,     only : data_file_type, dates
   use utils,            only : diag_level

!---------------------------------------------------------------------
!  dummy arguments
!---------------------------------------------------------------------
   integer, intent(in) :: ide
   integer, intent(in) :: jde
   integer, intent(in) :: domain
   real, intent(in)    :: dx
   character(len=*), intent(in)        :: cat_var_prefix
   character(len=*), intent(in)        :: cat_var_suffix
   type(data_file_type), intent(inout) :: data_file
   type(dates), intent(in) :: start_output

!---------------------------------------------------------------------
!  local variables
!---------------------------------------------------------------------
    integer, parameter :: lower = 0
    integer, parameter :: upper = 1

    integer :: i, j
    integer :: ids, jds
    integer :: il, iu, jl, ju
    integer :: n, varid, dimid
    integer :: nlon_src, nlat_src
    integer :: ncid
    integer :: astat, ierr, status
    integer :: xndx_src(2)
    integer :: yndx_src(2)
    integer :: maxind(1)
    real    :: d1, d2, ds1
    real    :: data_dx
    real    :: wrf_lon, wrf_lat
    real    :: wrf_lon_min, wrf_lat_min
    real    :: wrf_lon_max, wrf_lat_max
    real, allocatable :: src_lons(:)
    real, allocatable :: src_lats(:)
    real, allocatable :: wrk_emis(:,:)
    character(len=132) :: varname
    character(len=132) :: message
    logical :: new_grid
    logical :: has_area_map
    logical :: has_lon_shift
    logical :: reorder_lons
    logical :: reorder_lats
    type(grid_type), pointer :: grid

    write(*,*) ' '
    write(*,*) 'data_file_init: Initializing type for src emission file ' // trim(data_file%filename)
!---------------------------------------------------------------------
!   open src emission dataset file
!---------------------------------------------------------------------
    message = 'data_file_init: Failed to open ' // trim(data_file%filespec)
    call handle_ncerr( nf_open( trim(data_file%filespec), nf_noclobber, ncid ), message )       
!---------------------------------------------------------------------
!   get src dataset dimesions
!---------------------------------------------------------------------
    message = 'data_file_init: Failed to get lon dimension id'
    call handle_ncerr( nf_inq_dimid( ncid, 'lon', dimid ), message )
    message = 'data_file_init: Failed to get lon dimension'
    call handle_ncerr( nf_inq_dimlen( ncid, dimid, nlon_src ), message )
    message = 'data_file_init: Failed to get lat dimension id'
    call handle_ncerr( nf_inq_dimid( ncid, 'lat', dimid ), message )
    message = 'data_file_init: Failed to get lat dimension'
    call handle_ncerr( nf_inq_dimlen( ncid, dimid, nlat_src ), message )
    write(*,*) 'data_file_init:  nlon_src, nlat_src = ',nlon_src,nlat_src

!---------------------------------------------------------------------
!   allocate working variable
!---------------------------------------------------------------------
    if( allocated( src_lons ) ) then
      deallocate( src_lons )
    endif
    allocate( src_lons(nlon_src),stat=astat )
    if( astat /= 0 ) then
      write(*,*) 'data_file_init: Failed to allocate src_lons; error = ',astat
      stop 'allocate failed'
    endif
!---------------------------------------------------------------------
!   read src longitude variable
!---------------------------------------------------------------------
    message = 'data_file_init: Failed to get lon variable id'
    call handle_ncerr( nf_inq_varid( ncid, 'lon', varid ), message )
    message = 'data_file_init: Failed to read lon variable'
    call handle_ncerr( nf_get_var_real( ncid, varid, src_lons ), message )
!---------------------------------------------------------------------
!   check that lons are monotinicity increasing and reorder if not
!---------------------------------------------------------------------
    reorder_lons = .false.
    if( all( src_lons(2:nlon_src) < src_lons(1:nlon_src-1) ) ) then
      reorder_lons = .true.
      do i = 1,nlon_src/2
        d1 = src_lons(i)
        iu = nlon_src - i + 1
        src_lons(i)  = src_lons(iu)
        src_lons(iu) = d1
      end do
    endif
!---------------------------------------------------------------------
!   make sure src_lons are in [-180,180]
!---------------------------------------------------------------------
   has_lon_shift = .false.
   if( all( src_lons(:) >= 0. ) ) then
     if( count( src_lons(:) >= 180. ) > 0 ) then
       d2 = 360. + src_lons(1) - src_lons(nlon_src)
       d1 = src_lons(2) - src_lons(1)
       if( abs( d2 - d1 ) > 1.e-3*d1 ) then
         write(*,*) 'data_file_init: can not mapped data longitudes to WRF grid'
         stop
       endif
       src_lons(:) = cshift( src_lons(:),nlon_src/2 ) 
       src_lons(1:nlon_src/2) = src_lons(1:nlon_src/2) - 360.
       has_lon_shift = .true.
     endif
   endif

   if( allocated( src_lats ) ) then
     deallocate( src_lats )
   endif
   allocate( src_lats(nlat_src),stat=ierr )
   if( ierr /= 0 ) then
     write(*,*) 'data_file_init: Failed to allocate src_lats; error = ',ierr
     stop 'allocate failed'
   endif
!---------------------------------------------------------------------
!   read src latitude variable
!---------------------------------------------------------------------
   message = 'data_file_init: Failed to get lat variable id'
   call handle_ncerr( nf_inq_varid( ncid, 'lat', varid ), message )
   message = 'data_file_init: Failed to read lat variable'
   call handle_ncerr( nf_get_var_real( ncid, varid, src_lats ), message )

!---------------------------------------------------------------------
!   check that lats are monotinicity increasing and reorder if not
!---------------------------------------------------------------------
    reorder_lats = .false.
    if( all( src_lats(2:nlat_src) < src_lats(1:nlat_src-1) ) ) then
      reorder_lats = .true.
      do j = 1,nlat_src/2
        d1 = src_lats(j)
        ju = nlat_src - j + 1
        src_lats(j)  = src_lats(ju)
        src_lats(ju) = d1
      end do
    endif

!---------------------------------------------------------------------
!   determine interpolation type; bilinear or area conserving
!---------------------------------------------------------------------
    data_dx = earth_radius_m * (src_lats(2) - src_lats(1)) * rad_per_deg
    has_area_map = data_dx < dx
    write(*,*) 'data_file_init: data_dx,dx,has_area_map = ',data_dx,dx,has_area_map

!-------------------------------------------------------------
!   check for match against prior datasets
!-------------------------------------------------------------
   if( grid_cnt >= grid_max ) then
     write(*,*) 'data_file_init: reached grid cache max: ',grid_max
     stop
   endif
   grid_ndx = 0
   new_grid = .true.
   do n = 1,grid_cnt
     if( grid_specs(n)%nlons /= nlon_src .or. grid_specs(n)%nlats /= nlat_src ) then
       cycle
     endif
     if( any( grid_specs(n)%lon(:) /= src_lons(:) ) ) then
       cycle
     endif
     if( any( grid_specs(n)%lat(:) /= src_lats(:) ) ) then
       cycle
     endif
     grid_ndx = n
     new_grid = .false.
     exit
   end do
!-------------------------------------------------------------
!   new data grid to cache
!-------------------------------------------------------------
has_new_grid : &
   if( new_grid ) then
     grid_cnt = grid_cnt + 1
     grid => grid_specs(grid_cnt)
     grid%nlons = nlon_src
     grid%nlats = nlat_src
     grid%has_area_map  = has_area_map
     grid%has_lon_shift = has_lon_shift
     grid%reorder_lons  = reorder_lons
     grid%reorder_lats  = reorder_lats
     grid_ndx = grid_cnt
     allocate( grid%lon(nlon_src), &
               grid%xedge(nlon_src+1),stat=ierr )
     if( ierr /= 0 ) then
       write(*,*) 'data_file_init: Failed to allocate src_lats; error = ',ierr
       stop 'allocate failed'
     endif
     allocate( grid%lat(nlat_src), &
               grid%yedge(nlat_src+1),stat=ierr )
     if( ierr /= 0 ) then
       write(*,*) 'data_file_init: Failed to allocate src_lats; error = ',ierr
       stop 'allocate failed'
     endif
     grid%lon(:) = src_lons(:)
     grid%lat(:) = src_lats(:)
     write(*,*) 'data_file_init: file ' // trim(data_file%filename),' is a new grid'

is_area_map : &
     if( has_area_map ) then
!---------------------------------------------------------------------
!   form src longitude edges
!---------------------------------------------------------------------
       grid%xedge(2:nlon_src) = .5_8*(real(src_lons(1:nlon_src-1),kind=8) + real(src_lons(2:nlon_src),kind=8))
       grid%xedge(1)          = real(src_lons(1),kind=8) - .5_8*(real(src_lons(2),kind=8) - real(src_lons(1),kind=8))
       grid%xedge(nlon_src+1) = real(src_lons(nlon_src),kind=8) + .5_8*(real(src_lons(nlon_src),kind=8) - real(src_lons(nlon_src-1),kind=8))
       write(*,'(''data_file_init: xcen_src(1,2)  = '',1p,2g22.15)') src_lons(1:2)
       write(*,'(''data_file_init: xedge(1,2) = '',1p,2g22.15)') grid_specs(grid_cnt)%xedge(1:2)
       write(*,'(''data_file_init: dx = '',1pg22.15)') int( 1./(src_lons(2) - src_lons(1)) )
!---------------------------------------------------------------------
!   form src latitude edges
!---------------------------------------------------------------------
       grid%yedge(2:nlat_src) = .5_8*(src_lats(1:nlat_src-1) + src_lats(2:nlat_src))
       grid%yedge(1)          = src_lats(1) - .5_8*(src_lats(2) - src_lats(1))
       grid%yedge(nlat_src+1) = src_lats(nlat_src) + .5_8*(src_lats(nlat_src) - src_lats(nlat_src-1))

       write(*,'(''data_file_init: nlon_src,nlat_src = '',i6,1x,i6)') nlon_src,nlat_src
       write(*,'(''data_file_init: ycen_src  = '',1p,2g22.15)') src_lats(nlat_src-1:nlat_src)
       write(*,'(''data_file_init: yedge_src = '',1p,2g22.15)') grid_specs(grid_cnt)%yedge(nlat_src:nlat_src+1)

       allocate( grid%model_area_type(ide,jde),stat=astat )
       if( astat /= 0 ) then
         write(*,*) 'proj_init; failed to allocate model_area_type: error = ',astat
         stop 'Alloc error'
       endif
       grid%model_area_type(:,:)%has_data = .false.
       grid%model_area_type(:,:)%active_dcell_cnt = 0
       grid%model_area_type(:,:)%total_dcell_cnt  = 0
       grid%model_area_type(:,:)%interior_dcell_cnt = 0
       grid%model_area_type(:,:)%partial_dcell_cnt  = 0
!---------------------------------------------------------------------
!   area conserving interpolation
!---------------------------------------------------------------------
       if( .not. allocated( wrk_emis ) ) then
         allocate( wrk_emis(ide,jde),stat=status )
         if( status /= 0 ) then
           write(*,*) 'data_file_init: allocate for wrk_emis failed; error = ',ierr
           stop 'Alloc error'
         endif
       endif
       call area_interp_init( grid, grid%model_area_type, diag_level )
     else is_area_map
!---------------------------------------------------------------------
!   form src coordinate limits
!---------------------------------------------------------------------
       ids = 1 ; jds = 1
       wrf_lon_min = minval( xlong(ids:ide,jds:jde) )
       wrf_lon_max = maxval( xlong(ids:ide,jds:jde) )
       wrf_lat_min = minval( xlat(ids:ide,jds:jde) )
       wrf_lat_max = maxval( xlat(ids:ide,jds:jde) )
       if( diag_level > 100 ) then
        write(*,*) ' '
        write(*,'('' data_file_init: model lon limits = '',1p,2g25.16)') wrf_lon_min,wrf_lon_max
        write(*,'('' data_file_init: model lat limits = '',1p,2g25.16)') wrf_lat_min,wrf_lat_max
        write(*,*) ' '
        write(*,'('' data_file_init: src lon limits = '',1p,2g25.16)') src_lons(1),src_lons(nlon_src)
        write(*,'('' data_file_init: src lat limits = '',1p,2g25.16)') src_lats(1),src_lats(nlat_src)
        write(*,*) ' '
       endif
!---------------------------------------------------------------
!     allocate memory space to store interpolation coef.
!---------------------------------------------------------------
      ierr = 0
      allocate( grid%ax(ids:ide,jds:jde,0:1), &
                grid%by(ids:ide,jds:jde,0:1), stat=status )
      ierr = ierr + status
      allocate( grid%ix(ids:ide,jds:jde,0:1), &
                grid%jy(ids:ide,jds:jde,0:1), stat=status )
      ierr = ierr + status
      if( ierr /= 0 ) then
        write(*,*) 'allocate for ax ... jy failed; error = ',ierr
      endif
!---------------------------------------------------------------------
!   set bilinear interp variables
!---------------------------------------------------------------------
lat_loop : &
      do j = jds,jde
        do i = ids,ide
!---------------------------------------------------------------------
!   longitudes
!---------------------------------------------------------------------
          wrf_lon = xlong(i,j)
          if( wrf_lon >= src_lons(nlon_src) .or. &
              wrf_lon < src_lons(1) ) then
            grid%ix(i,j,lower) = nlon_src
          else
            do n = 2,nlon_src
              if( wrf_lon < src_lons(n) ) then
                grid%ix(i,j,lower) = min( nlon_src-1, max(n-1,1) )
                exit
              endif
            end do
          endif
          grid%ix(i,j,upper) = mod( grid%ix(i,j,lower),nlon_src ) + 1
          data_dx = src_lons(grid%ix(i,j,upper)) - src_lons(grid%ix(i,j,lower))
          if( data_dx < 0. ) then
            data_dx = 360. + data_dx
          endif
          ds1 = wrf_lon - src_lons(grid%ix(i,j,lower))
          if( ds1 < 0. ) then
            ds1 = 360. + ds1
          endif
          grid%ax(i,j,lower) = ds1/data_dx
          grid%ax(i,j,upper) = 1.0 - grid%ax(i,j,lower)
!---------------------------------------------------------------------
!   latitudes
!---------------------------------------------------------------------
          wrf_lat = xlat(i,j)
          if( wrf_lat < src_lats(1) ) then
            grid%jy(i,j,0:1) = -1
            grid%by(i,j,0:1) = 0.
          elseif( wrf_lat > src_lats(nlat_src) ) then
            grid%jy(i,j,0:1) = -2
            grid%by(i,j,0:1) = 0.
          else
            do n = 1,nlat_src
              if( wrf_lat < src_lats(n) ) then
                exit
              endif
            end do
            grid%jy(i,j,lower) = min( nlat_src-1, max(n-1,1) )
            grid%jy(i,j,upper) = grid%jy(i,j,lower) + 1
            grid%by(i,j,lower) = (wrf_lat - src_lats(grid%jy(i,j,lower))) &
                             /(src_lats(grid%jy(i,j,upper)) - src_lats(grid%jy(i,j,lower)))
            grid%by(i,j,upper) = 1.0 - grid%by(i,j,lower)
          endif
        end do
      end do lat_loop
!---------------------------------------------------------------------
!   form dataset index limits
!---------------------------------------------------------------------
      write(*,*) 'data_file_init: count of points <,> data min,max lat = ',count(grid%ix(:,:,0) == nlon_src )
      xndx_src(1) = minval( grid%ix(:,:,lower) )
      xndx_src(2) = maxval( grid%ix(:,:,upper) )
      write(*,*) 'xndx_src = ',xndx_src(:)
      write(*,*) 'data_file_init: count of points < data min lat = ',count(grid%jy(:,:,0) == -1)
      write(*,*) 'data_file_init: count of points > data max lat = ',count(grid%jy(:,:,0) == -2)
      yndx_src(1) = minval( grid%jy(:,:,lower),mask=grid%jy(:,:,lower)>0 )
      yndx_src(2) = maxval( grid%jy(:,:,upper),mask=grid%jy(:,:,upper)>0 )
      write(*,*) 'yndx_src = ',yndx_src(:)
     endif is_area_map
   endif has_new_grid

!---------------------------------------------------------------------
!   if not set, get src molecular weight
!---------------------------------------------------------------------
   if( data_file%molecw == 0. ) then
     call get_molecwght( ncid, cat_var_prefix, cat_var_suffix, data_file%molecw, data_file%missing_value )
   endif

!---------------------------------------------------------------------
!   set category units
!---------------------------------------------------------------------
   if( domain == 1 ) then
     call get_units( ncid, cat_var_prefix, cat_var_suffix, data_file%molecw, data_file%con_fac )
   endif

   data_file%grid_ndx = grid_ndx

!---------------------------------------------------------------------
!   initialize timing for data file
!---------------------------------------------------------------------
   data_file%lo_tndx = 0
   data_file%hi_tndx = 0
   data_file%lo_buf_ndx = 1
   data_file%hi_buf_ndx = 2
   data_file%ncid_lo    = 0
   data_file%ncid_hi    = 0
   data_file%in_gap     = .false.
   data_file%t_interp   = .false.

!---------------------------------------------------------------------
!   initialize data file timing
!---------------------------------------------------------------------
   call data_file_timing_init( ncid, start_output, data_file ) 

   if( data_file%lo_tndx > 0 ) then
     data_file%lo_tndx = 0
     data_file%ncid_lo = ncid
     data_file%ncid_hi = ncid
   endif

   n = count( data_file%cat_active(:) )
   allocate( data_file%emis(nlon_src,nlat_src,2,n), &
             data_file%src_data(nlon_src,nlat_src),stat=astat )
   if( astat /= 0 ) then
     write(*,*) 'data_file_init: allocate for emis,src_data failed; error = ',astat
     stop 'Alloc error'
   endif

   end subroutine data_file_init

   subroutine data_file_timing_init( ncid, interp_time, data_file )
!---------------------------------------------------------------
!     initialize data file timing info
!---------------------------------------------------------------

   use anthro_types, only : data_file_type, dates
   use mo_calendar, only  : diffdat

!---------------------------------------------------------------
!     dummy arguments
!---------------------------------------------------------------
   integer, intent(in)          :: ncid
   type(data_file_type), intent(inout) :: data_file
   type(dates),          intent(in)    :: interp_time

!---------------------------------------------------------------
!     local variables
!---------------------------------------------------------------
   integer :: status, tstat
   integer :: dimid, varid
   integer :: i, j, k, n 

!---------------------------------------------------------------
!     read times
!---------------------------------------------------------------
   call read_src_times( ncid, data_file%date, data_file%secs, data_file%ntimes )

   data_file%lo_tndx = lotim( interp_time%date, interp_time%secs, data_file )
   if( data_file%lo_tndx <= 0 ) then
!---------------------------------------------------------------
!     check either prior or next data file
!---------------------------------------------------------------
     status = nf_close( ncid )
     if( data_file%lo_tndx == 0 ) then
       call next_flnm( data_file%filename, .false. )
     else
       call next_flnm( data_file%filename, .true. )
     endif
!---------------------------------------------------------------
!     open the prior or next mozart netCDF file
!---------------------------------------------------------------
     data_file%filespec = trim(anthro_dir) // '/' // trim(data_file%filename)
     status = nf_open( trim(data_file%filespec), nf_nowrite, ncid )
     if( status /= nf_noerr ) then
       write(*,*) 'data_file_timing_init : Failed to open ',trim(data_file%filespec)
       call handle_error( status )
     end if
     write(*,*) 'data_file_timing_init: opened ',trim(data_file%filespec)
     call read_src_times( ncid, data_file%date, data_file%secs, data_file%ntimes )
     data_file%lo_tndx = lotim( interp_time%date, interp_time%secs, data_file )
     if( data_file%lo_tndx == 0 ) then
       write(*,*) 'data_file_timing_init: time ',interp_time%date,' ',interp_time%secs
       write(*,*) '                 is before '
       write(*,*) '                 ',data_file%date(1),' ',data_file%secs(1)
       write(*,*) '                 the first time in the data file ',trim(data_file%filename)
       status = nf_close( ncid )
       stop
     end if
   end if

   end subroutine data_file_timing_init

   subroutine read_src_times( ncid, date, secs, ntimes )
!---------------------------------------------------------------
!     read times from current src input file
!---------------------------------------------------------------

      use mo_calendar, only : newdate

!---------------------------------------------------------------
!     dummy arguments
!---------------------------------------------------------------
      integer, intent(in) :: ncid
      integer, intent(out) :: ntimes
      integer, allocatable:: date(:)
      integer, allocatable:: secs(:)

!---------------------------------------------------------------
!     local variables
!---------------------------------------------------------------
      integer :: m
      integer :: dimid, varid
      integer :: istat, status
      integer :: sl, su
      integer :: yr, mn, dy
      real, allocatable :: times(:)
      character(len=132) :: units_text
      character(len=4)   :: time_name
      logical :: got_it

!---------------------------------------------------------------
!     get time dimension and allocate arrays
!---------------------------------------------------------------
      time_name = 'time'
      status = nf_inq_dimid( ncid, trim(time_name), dimid )
      if( status /= nf_noerr )  then
        time_name = 'date'
        status = nf_inq_dimid( ncid, trim(time_name), dimid )
        if( status /= nf_noerr )  then
          write(*,*) 'read_src_times: can not find either time or date dimension'
          call handle_error( status )
        endif
      endif

      status = nf_inq_dimlen( ncid, dimid, ntimes )
      if( status /= nf_noerr ) call handle_error( status )

      if( allocated( date ) ) then
        deallocate( date, secs )
      end if
      allocate( date(ntimes), secs(ntimes), stat=status )
      if( status /= 0 ) then
        write(*,*) 'failed to allocate date, secs; error = ',status
        stop 'Alloc err'
      end if

!---------------------------------------------------------------
!     check for "standard" date, datesec variables
!---------------------------------------------------------------
      got_it = .false.
      status = nf_inq_varid( ncid, 'date', varid )
      if( status == nf_noerr ) then
        status = nf_inq_varid( ncid, 'datesec', m )
        if( status == nf_noerr ) then
          status = nf_get_var_int( ncid, varid, date )
          if( status /= nf_noerr )  call handle_error( status )
          status = nf_get_var_int( ncid, m, secs )
          if( status /= nf_noerr )  call handle_error( status )
          got_it = .true.
        endif
      endif
!---------------------------------------------------------------
!     try to get date, datesec from time variable
!---------------------------------------------------------------
      if( .not. got_it ) then
        status = nf_inq_varid( ncid, trim(time_name), varid )
        if( status /= nf_noerr )  call handle_error( status )
        units_text = ' '
        status = nf_get_att_text( ncid, varid, 'units', units_text )
        if( status /= nf_noerr )  call handle_error( status )
        sl = scan( units_text, '0123456789' )
        if( sl > 0 ) then
          su = index( units_text, '-', back=.true. )
          if( units_text(sl+4:sl+4) == '-' .and. su == sl+7 ) then
            read(units_text(sl:sl+3),*,iostat=istat) yr
            if( istat == 0 ) then
              read(units_text(sl+5:sl+6),*,iostat=istat) mn
              if( istat == 0 ) then
                read(units_text(sl+8:sl+9),*,iostat=istat) dy
              endif
            endif
            if( istat == 0 ) then
              allocate( times(ntimes),stat=istat )
              if( istat /= 0 ) then
                write(*,*) 'read_src_times: failed to allocate times; error = ',istat
                stop 'Alloc err'
              endif
              status = nf_get_var_real( ncid, varid, times )
              if( status /= nf_noerr )  call handle_error( status )
              date(:) = dy + 100*(mn + 100*yr)
              do m = 1,ntimes
                date(m) = newdate( date(m), int( times(m) ) )
                secs(m) = 86400*(times(m) - aint(times(m)))
              end do
              got_it = .true.
            endif
          endif
        endif
      endif

      if( .not. got_it ) then
        write(*,*) 'read_src_times: failed to read date,secs'
        stop
      endif

      if( data_yrs_offset /= 0 ) then
        do m = 1,ntimes
          date(m) = date(m) + 10000*data_yrs_offset
        end do
      endif

      end subroutine read_src_times

      subroutine get_src_time_ndx( data_file, loop_time )
!---------------------------------------------------------------
!     get src time index for loop_time
!---------------------------------------------------------------

      use anthro_types, only : data_file_type, dates
      use mo_calendar, only : diffdat

!---------------------------------------------------------------
!     dummy arguments
!---------------------------------------------------------------
      type(data_file_type), intent(inout)  :: data_file
      type(dates), intent(in)              :: loop_time

!---------------------------------------------------------------
!     local variables
!---------------------------------------------------------------
      integer :: i, n
      integer :: status
      integer :: ncid
      character(len=128) :: filenm
      logical :: found

      write(*,*) ' '
      write(*,*) 'get_src_time_ndx; src_dir,src_fn = ',trim(data_file%filespec)
      write(*,*) 'get_src_time_ndx; interp_date,datesec,ntimes = ',loop_time%date,loop_time%secs,data_file%ntimes

      n = lotim( loop_time%date, loop_time%secs, data_file )
      if( n > 0 ) then
!---------------------------------------------------------------
!     loop time in present src dataset
!---------------------------------------------------------------
        write(*,*) 'get_src_time_ndx; tndx = ',n
        if( .not. data_file%in_gap ) then
          data_file%read_lo_tndx = n /= data_file%lo_tndx
          if( data_file%read_lo_tndx ) then
            data_file%lo_tndx = n
          endif
          data_file%in_gap  = .false.
          if( data_file%t_interp ) then
            data_file%read_hi_tndx = (n + 1) /= data_file%hi_tndx
            if( data_file%read_hi_tndx ) then
              data_file%hi_tndx = n + 1
            endif
          endif
        else
!---------------------------------------------------------------
!     special handling if prior time was in "gap"
!---------------------------------------------------------------
          data_file%in_gap  = .false.
          data_file%lo_tndx = n
          data_file%read_lo_tndx = .true.
          if( data_file%t_interp ) then
            data_file%hi_tndx = n + 1
            data_file%read_hi_tndx = .true.
          endif
        endif
        if( data_file%ncid_hi /= data_file%ncid_lo ) then
          status = nf_close( data_file%ncid_lo )
          if( status /= 0 ) then
            write(*,*) 'get_src_time_ndx: failed to close ',trim(data_file%filespec),' ; error = ',status
            stop
          end if
          data_file%ncid_lo = data_file%ncid_hi
        endif
        data_file%ncid_hi = data_file%ncid_lo
      else if( n < 0 ) then
!---------------------------------------------------------------
!     loop time after present src dataset
!---------------------------------------------------------------
         data_file%lo_tndx = data_file%ntimes
         call next_flnm( data_file%filename, .true. )
!---------------------------------------------------------------
!     open the input netCDF file
!---------------------------------------------------------------
         data_file%filespec = trim(anthro_dir) // '/' // trim(data_file%filename)
         status = nf_open( trim(data_file%filespec), nf_nowrite, ncid )
         if( status /= nf_noerr ) call handle_error( status )
         write(*,*) 'get_src_time_ndx: opened ',trim(data_file%filespec)
         data_file%gap_date = data_file%date(data_file%ntimes)
         data_file%gap_secs = data_file%secs(data_file%ntimes)
         call read_src_times( ncid, data_file%date, data_file%secs, data_file%ntimes )
         n = lotim( loop_time%date, loop_time%secs, data_file )
         if( n > 0 ) then
           write(*,*) 'get_src_time_ndx; src_tndx = ',n
           status = nf_close( data_file%ncid_lo )
           if( status /= 0 ) then
             write(*,*) 'get_src_time_ndx: failed to close ',trim(data_file%filespec),' ; error = ',status
             stop
           end if
           data_file%in_gap  = .false.
           data_file%ncid_lo = ncid
           data_file%ncid_hi = ncid
           data_file%lo_tndx = n
           data_file%read_lo_tndx = .true.
           if( data_file%t_interp ) then
             data_file%hi_tndx = n + 1
             data_file%read_hi_tndx = .true.
           endif
         else if( n == 0 ) then
           data_file%in_gap  = .true.
           data_file%hi_tndx = 1
           data_file%ncid_hi = ncid
           data_file%dels = data_file%dels/diffdat( data_file%gap_date, data_file%gap_secs, data_file%date(1), data_file%secs(1) )
           data_file%t_interp = .true.
           data_file%read_lo_tndx = .true.
           data_file%read_hi_tndx = .true.
         else
           write(*,*) 'get_src_time_ndx: failed to find ',loop_time%date,' : ',loop_time%secs
           write(*,*) '                  in file ',trim(data_file%filename)
           stop
         end if
      else if( data_file%in_gap ) then
        data_file%dels = data_file%dels/diffdat( data_file%gap_date, data_file%gap_secs, data_file%date(1), data_file%secs(1) )
        data_file%read_lo_tndx = .false.
        data_file%read_hi_tndx = .false.
      end if

      end subroutine get_src_time_ndx

      subroutine read_src_data( data_file, varname, bndx, tndx, ncid, cat_ndx )
!---------------------------------------------------------------
!     read source emission category
!---------------------------------------------------------------

      use anthro_types, only : data_file_type
      use mapper_types, only : grid_type, grid_specs

!---------------------------------------------------------------
!     dummy arguments
!---------------------------------------------------------------
      integer, intent(in) :: bndx              ! buffer index
      integer, intent(in) :: tndx              ! time index
      integer, intent(in) :: ncid              ! netcdf file index
      integer, intent(in) :: cat_ndx           ! category buffer index
      type(data_file_type), intent(inout) :: data_file
      character(len=*), intent(in) :: varname

!---------------------------------------------------------------
!     local variables
!---------------------------------------------------------------
      integer :: gndx
      integer :: n
      integer :: status
      integer :: varid
      integer :: nstt(3), ncnt(3)

      gndx = data_file%grid_ndx
      nstt(:) = (/ 1, 1, tndx /)
      ncnt(:) = (/ grid_specs(gndx)%nlons, grid_specs(gndx)%nlats, 1 /)
      data_file%emis(:,:,bndx,cat_ndx) = 0.

      status = nf_inq_varid( ncid, trim(varname), varid )
      if( status /= nf_noerr )  then
        write(*,*) 'read_src_data: failed to get id of ',trim(varname)
        call handle_error( status )
      end if
      status = nf_get_vara_real( ncid, varid, nstt, ncnt, data_file%emis(1,1,bndx,cat_ndx) )
      if( status /= nf_noerr ) then
        write(*,*) 'read_src_data: failed to read ',trim(varname)
        call handle_error( status )
      end if

!---------------------------------------------------------------
!     replace missing data with zero
!---------------------------------------------------------------
      where( data_file%emis(:,:,bndx,cat_ndx) == data_file%missing_value )
        data_file%emis(:,:,bndx,cat_ndx) = 0.
      endwhere

      end subroutine read_src_data

      subroutine tinterp_src_data( data_file, cat_ndx )
!---------------------------------------------------------------
!     time interpolation of source emission
!---------------------------------------------------------------

      use anthro_types, only : data_file_type
      use mapper_types, only : grid_type, grid_specs

!---------------------------------------------------------------
!     dummy arguments
!---------------------------------------------------------------
      integer, intent(in) :: cat_ndx
      type(data_file_type), intent(inout) :: data_file

!---------------------------------------------------------------
!     local variables
!---------------------------------------------------------------
      integer :: i, j, iu, ju
      integer :: lo_ndx, hi_ndx
      integer :: astat
      integer :: gndx
      real    :: d1
      real    :: dels                ! linear interp factor
      real    :: delsm1              ! linear interp factor
      real, allocatable  :: data_col(:)

      gndx = data_file%grid_ndx
      dels = data_file%dels
      lo_ndx = data_file%lo_buf_ndx
      if( data_file%t_interp ) then
        hi_ndx = data_file%hi_buf_ndx
        delsm1 = 1. - dels
        do j = 1,grid_specs(gndx)%nlats
          data_file%src_data(:,j) = data_file%emis(:,j,lo_ndx,cat_ndx) * dels &
                                  + data_file%emis(:,j,hi_ndx,cat_ndx) * delsm1
        end do
      else
        do j = 1,grid_specs(gndx)%nlats
          data_file%src_data(:,j) = data_file%emis(:,j,lo_ndx,cat_ndx)
        end do
      endif

      if( grid_specs(gndx)%reorder_lons ) then
!-------------------------------------------------------------
!  longitude reorder the src data
!-------------------------------------------------------------
        do j = 1,grid_specs(gndx)%nlats
          do i = 1,grid_specs(gndx)%nlons/2
            d1 = data_file%src_data(i,j)
            iu = grid_specs(gndx)%nlons - i + 1
            data_file%src_data(i,j)  = data_file%src_data(iu,j)
            data_file%src_data(iu,j) = d1
          end do
        end do
      endif

      if( grid_specs(gndx)%has_lon_shift ) then
!-------------------------------------------------------------
!  longitude shift the src data
!-------------------------------------------------------------
        do j = 1,grid_specs(gndx)%nlats
          data_file%src_data(:,j) = cshift( data_file%src_data(:,j),grid_specs(gndx)%nlons/2 )
        end do
      endif

      if( grid_specs(gndx)%reorder_lats ) then
!-------------------------------------------------------------
!  latitude reorder the src data
!-------------------------------------------------------------
        allocate( data_col(grid_specs(gndx)%nlons),stat=astat )
        if( astat /= 0 ) then
          write(*,*) 'tinerpt_src_data: failed to allocate data_col: error = ',astat
          stop 'Alloc error'
        endif
        do j = 1,grid_specs(gndx)%nlats/2
          ju = grid_specs(gndx)%nlats - j + 1
          data_col(:) = data_file%src_data(:,j)
          data_file%src_data(:,j)  = data_file%src_data(:,ju)
          data_file%src_data(:,ju) = data_col(:)
        end do
        deallocate( data_col )
      endif

      end subroutine tinterp_src_data

      subroutine next_flnm( filenm, incr )
!---------------------------------------------------------------
!     increment mozart filename
!---------------------------------------------------------------

!---------------------------------------------------------------
!     dummy arguments
!---------------------------------------------------------------
      character(len=*), intent(inout) :: filenm
      logical, intent(in)             :: incr

!---------------------------------------------------------------
!     local variables
!---------------------------------------------------------------
      integer :: il, iu, ios
      integer :: file_number
      character(len=5) :: numa
      logical :: found

      found = .false.
      iu    = scan( trim(filenm), '0123456789', back=.true. )
      do il = iu-1,1,-1
         if( index( '0123456789', filenm(il:il) ) == 0 ) then
            found = .true.
            exit
         end if
      end do
      if( .not. found ) then
         write(*,*) 'next_filenm: file ',trim(filenm),' is not proper format'
         stop
      end if

      il = il + 1
      write(*,*) ' '
      write(*,*) 'next_flnm; trying to increment file ',trim(filenm)
      write(*,*) 'next_flnm; il, iu = ',il,iu
      read(filenm(il:iu),*,iostat=ios) file_number
      if( ios /= 0 ) then
         write(*,*) 'next_filenm: failed to read ',filenm(il:iu),' ; error = ',ios
         stop
      end if
      write(*,*) 'next_flnm; file_number = ',file_number
      if( incr ) then
        file_number = file_number + 1
      else
        file_number = file_number - 1
      endif
      write(numa,'(i5)') file_number+10000
      filenm(il:iu) = numa(2:5)

      write(*,*) 'next_flnm; new file = ',trim(filenm)

      end subroutine next_flnm

      integer function lotim( cdate, csec, data_file )
!-----------------------------------------------------------------------
! 	... return the index of the time sample that is the lower
!           bound of the interval that contains the input date.  if
!           (cdate,csec) is earlier than the first time sample then 0 is
!           returned.  if (cdate,csec) is later than the last time sample then
!           -index is returned.  if (cdate,csec) is equal to the date of a
!           dynamics time sample then that index is returned.
!-----------------------------------------------------------------------

      use mo_calendar,  only : diffdat
      use anthro_types, only : data_file_type

      implicit none

!-----------------------------------------------------------------------
! 	... dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in) :: cdate    ! date in yyyymmdd
      integer, intent(in) :: csec     ! seconds relative to date
      type(data_file_type), intent(inout) :: data_file

!-----------------------------------------------------------------------
! 	... local variables
!-----------------------------------------------------------------------
      integer :: n
      real    :: dtime

!-----------------------------------------------------------------------
!     	... find latest date that is earlier than or equal to (date,sec)
!-----------------------------------------------------------------------
      do n = 1,data_file%ntimes
        dtime = diffdat( cdate, csec, data_file%date(n), data_file%secs(n) )
        if( dtime > 0. ) then
          lotim = n - 1
          if( lotim > 0 ) then
            data_file%dels = dtime/diffdat( data_file%date(lotim), data_file%secs(lotim), data_file%date(n), data_file%secs(n) )
            data_file%t_interp = data_file%dels /= 1.
          else
            data_file%dels = dtime
          endif
          exit
        endif
      end do

      if( n > data_file%ntimes ) then
        if( dtime == 0. ) then
          lotim = data_file%ntimes
          data_file%dels = 1.
          data_file%t_interp = .false.
        else
          lotim = -data_file%ntimes
          data_file%dels = dtime
        endif
      endif

      end function lotim

   subroutine get_molecwght( ncid, cat_var_prefix, cat_var_suffix, molecw, missing_value )
!---------------------------------------------------------------------
!	... try to get the source species molecular wght from src file
!---------------------------------------------------------------------

   use utils, only : n_sub_cats, sub_cats

!---------------------------------------------------------------------
!	... dummy args
!---------------------------------------------------------------------
   integer, intent(in) :: ncid
   real, intent(in)    :: missing_value
   real, intent(inout) :: molecw
   character(len=*), intent(in) :: cat_var_prefix
   character(len=*), intent(in) :: cat_var_suffix

!---------------------------------------------------------------------
!	... local variables
!---------------------------------------------------------------------
   integer :: m
   integer :: ncstat
   integer :: varid
   integer :: molecw_int
   integer(2) :: molecw_int2
   integer :: att_type, att_len
   real(8) :: molecw_8
   character(len=64)  :: var_nam
   character(len=132) :: message
   logical :: got_it

   got_it = .false.
!---------------------------------------------------------------------
!	... first try to get the variable molecular_weight
!---------------------------------------------------------------------
   ncstat = nf_inq_varid( ncid, 'molecular_weight', varid )
   if( ncstat == nf_noerr ) then
     ncstat = nf_get_var_real( ncid, varid, molecw )
     if( ncstat /= nf_noerr ) then
       write(*,*) 'get_molecwght: failed to read molecular weight: error = ',ncstat
       stop 'Netcdf error'
     endif
     got_it = molecw /= missing_value .and. molecw > 0. .and. molecw < 300.
   endif

   if( .not. got_it ) then
!---------------------------------------------------------------------
!	... next try to get global attribute molecular_weight
!---------------------------------------------------------------------
     ncstat = nf_inq_att( ncid, nf_global, 'molecular_weight', att_type, att_len )
     if( ncstat == nf_noerr ) then
       if( att_len == 1 ) then
         select case( att_type )
           case( nf_short )
             ncstat = nf_get_att_int2( ncid, nf_global, 'molecular_weight', molecw_int2 )
           case( nf_int )
             ncstat = nf_get_att_int( ncid, nf_global, 'molecular_weight', molecw_int )
           case( nf_float )
             ncstat = nf_get_att_real( ncid, nf_global, 'molecular_weight', molecw )
           case( nf_double )
             ncstat = nf_get_att_double( ncid, nf_global, 'molecular_weight', molecw_8 )
           case default
             ncstat = nf_noerr + 1
         end select
         if( ncstat /= nf_noerr ) then
           write(*,*) 'get_molecwght: failed to read molecular weight: error = ',ncstat
           stop 'Netcdf error'
         endif
         select case( att_type )
           case( nf_short )
             molecw = real(molecw_int2)
           case( nf_int )
             molecw = real(molecw_int)
           case( nf_double )
             molecw = real(molecw_8)
         end select
         got_it = .true.
       endif
     endif
   endif

   if( .not. got_it ) then
!---------------------------------------------------------------------
!	... finally try category attribute
!---------------------------------------------------------------------
sub_cat_loop : &
     do m = 1,n_sub_cats
       var_nam = trim( cat_var_prefix) // trim( sub_cats(m) ) // trim( cat_var_suffix )
       ncstat = nf_inq_varid( ncid, trim(var_nam), varid )
       if( ncstat == nf_noerr ) then
         ncstat = nf_inq_att( ncid, varid, 'molecular_weight', att_type, att_len )
         if( ncstat == nf_noerr ) then
           if( att_len == 1 ) then
             select case( att_type )
               case( nf_short )
                 ncstat = nf_get_att_int2( ncid, varid, 'molecular_weight', molecw_int2 )
               case( nf_int )
                 ncstat = nf_get_att_int( ncid, varid, 'molecular_weight', molecw_int )
               case( nf_float )
                 ncstat = nf_get_att_real( ncid, varid, 'molecular_weight', molecw )
               case( nf_double )
                 ncstat = nf_get_att_double( ncid, varid, 'molecular_weight', molecw_8 )
               case default
                 ncstat = nf_noerr + 1
             end select
             if( ncstat == nf_noerr ) then
               select case( att_type )
                 case( nf_short )
                   molecw = real(molecw_int2)
                 case( nf_int )
                   molecw = real(molecw_int)
                 case( nf_double )
                   molecw = real(molecw_8)
               end select
               got_it = .true.
             endif
           endif
         endif
       endif
       if( got_it ) then
         exit sub_cat_loop
       endif
     end do sub_cat_loop
   endif

   got_it = molecw /= missing_value .and. molecw > 0. .and. molecw < 300.
   if( .not. got_it ) then
     write(*,*) 'get_molecwght: failed to read molecular weight'
     call handle_error( ncstat )
     stop
   endif

   end subroutine get_molecwght

   subroutine get_units( ncid, cat_var_prefix, cat_var_suffix, molecw, con_fac )
!---------------------------------------------------------------------
!	... try to get the source species molecular wght from src file
!---------------------------------------------------------------------

   use utils, only : n_sub_cats, sub_cats
   use utils, only : upcase

!---------------------------------------------------------------------
!	... dummy args
!---------------------------------------------------------------------
   integer, intent(in) :: ncid
   real, intent(in)    :: molecw
   real, intent(inout) :: con_fac(2)
   character(len=*), intent(in) :: cat_var_prefix
   character(len=*), intent(in) :: cat_var_suffix

!---------------------------------------------------------------------
!	... local variables
!---------------------------------------------------------------------
   real, parameter :: avogadro = 6.02214129e23

   integer :: m
   integer :: ncstat
   integer :: varid
   integer :: att_type, att_len
   character(len=64)  :: var_nam
   character(len=132) :: units, wrk_string

sub_cat_loop : &
   do m = 1,n_sub_cats
     var_nam = trim( cat_var_prefix) // trim( sub_cats(m) ) // trim( cat_var_suffix )
     ncstat = nf_inq_varid( ncid, trim(var_nam), varid )
     if( ncstat == nf_noerr ) then
       ncstat = nf_inq_att( ncid, varid, 'units', att_type, att_len )
       if( ncstat == nf_noerr ) then
         if( att_type == nf_char .and. att_len /= 0 ) then
           units = ' '
           ncstat = nf_get_att_text( ncid, varid, 'units', units )
           if( ncstat == nf_noerr ) then
             call upcase( units(:len_trim(units)), wrk_string )
             select case( trim(wrk_string) )
               case( 'KG M-2 S-1','KG/M2/S','KG/M^2/S' )
                 con_fac(1) = 3.6e12 ; con_fac(2) = 1.e9
               case( 'MOLECULES/CM2/S','MOLECULES/CM^2/S','MOLECULES CM-2 S-1' )
                 con_fac(1) = 3.6e13/avogadro
                 con_fac(2) = 1.e10*molecw/avogadro
               case default
                 con_fac(1) = 3.6e12 ; con_fac(2) = 1.e9
             end select
             exit sub_cat_loop
           endif
         endif
       endif
     endif
   end do sub_cat_loop

   write(*,*) '===================================================='
   write(*,*) 'get_units: con_fac(1,2) = ',con_fac(:)
   write(*,*) '===================================================='

   end subroutine get_units

   subroutine handle_ncerr( ret, mes )
!---------------------------------------------------------------------
!	... netcdf error handling routine
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!	... dummy arguments
!---------------------------------------------------------------------
   integer, intent(in) :: ret
   character(len=*), intent(in) :: mes

   if( ret /= nf_noerr ) then
      write(*,*) nf_strerror( ret )
      stop 'netcdf error'
   endif

   end subroutine handle_ncerr

      subroutine handle_error( status )
!---------------------------------------------------------------
!     handle errors produced by calling netCDF functions
!---------------------------------------------------------------

!---------------------------------------------------------------
!     dummy arguments :
!---------------------------------------------------------------
      integer, intent(in) :: status

!---------------------------------------------------------------
!     print the error information from processing NETcdf file
!---------------------------------------------------------------
      write(*,*) nf_strerror( status )
      stop 'Netcdf error'

      end subroutine handle_error

   end module data_file_utils
