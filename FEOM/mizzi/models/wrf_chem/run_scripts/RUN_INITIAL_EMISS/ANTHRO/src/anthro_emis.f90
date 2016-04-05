
   program anthro_emis

   use anthro_types
   use utils
   use area_mapper, only : xlong => lon, xlat => lat
   use mapper_types
   use misc_definitions_module, only : PROJ_PS, PROJ_LATLON, PROJ_CASSINI
   use mo_calendar, only : diffdat, addsec2dat
   use data_file_utils, only : anthro_dir, data_yrs_offset
   use data_file_utils, only : data_file_init, get_src_time_ndx
   use data_file_utils, only : read_src_data, tinterp_src_data

   implicit none

!-----------------------------------------------------------------
!     control variables
!-----------------------------------------------------------------
   integer                :: nemis
   integer                :: emissions_zdim_stag = 10
   integer                :: domains = 1
   integer                :: output_interval = 3600      ! seconds
   character(len=linsize) :: wrf_dir
   character(len=linsize) :: src_file_prefix
   character(len=linsize) :: src_file_suffix
   character(len=linsize) :: cat_var_prefix
   character(len=linsize) :: cat_var_suffix
   character(len=3)       :: numa
   character(len=linsize) :: emis_map(linemax)
   character(len=namsize) :: sub_categories(maxsrc)
   character (LEN=19)     :: start_output_time = ' '
   character (LEN=19)     :: stop_output_time  = ' '
   character (LEN=19)     :: start_data_time = ' '
   character (LEN=19)     :: stop_data_time  = ' '
   character (LEN=16)     :: date_frmt  = ' '
   logical                :: serial_output = .false.

   namelist /control/ anthro_dir, wrf_dir, emis_map, domains, src_names, &
                      src_file_prefix, src_file_suffix, cat_var_prefix, &
                      cat_var_suffix, emissions_zdim_stag, sub_categories, &
                      date_frmt, start_output_time, stop_output_time, output_interval, &
                      start_data_time, stop_data_time, serial_output, data_yrs_offset, diag_level

!-----------------------------------------------------------------
!     local variables
!-----------------------------------------------------------------
   integer, parameter :: lower = 0
   integer, parameter :: upper = 1
   real, parameter    :: confac_gas = 3.6e12             !(kg/m^2/s -> mole/km^2/hr)
   real, parameter    :: confac_aer = 1.e9               !(kg/m^2/s -> ug/m^2/s)

   integer :: nc, ns, ns1, it
   integer :: src
   INTEGER :: icnt
   INTEGER :: ids, ide, jds, jde
   INTEGER :: i, j, n
   INTEGER :: ii, jj
   integer :: ncid
   integer :: ngatts
   integer :: cat_ndx
   integer :: map_proj
   integer :: ierr, astat, istat
   integer :: dimid, varid
   integer :: nlon_src, nlat_src
   integer :: domain
   integer :: xndx_src(2)
   integer :: yndx_src(2)
   integer, allocatable :: ix(:,:,:)                        ! index used by interpolation
   integer, allocatable :: jy(:,:,:)                        ! index used by interpolation

   real    :: scale_factor
   real    :: wrk_sum
   real    :: ds1, ds2
   real    :: xl, xu
   real    :: yl, yu, dy
   real    :: wrf_lon_min
   real    :: wrf_lon_max
   real    :: wrf_lat_min
   real    :: wrf_lat_max
   real    :: cen_lon
   real    :: cen_lat
   real    :: stand_lon
   real    :: truelat1
   real    :: truelat2
   real    :: loninc
   real    :: latinc
   real    :: knowni
   real    :: knownj
   real    :: dx
   real(8), allocatable :: xedge_src(:)
   real(8), allocatable :: yedge_src(:)
   real, allocatable :: wrk_emis(:,:)
   real, allocatable :: src_lons(:)
   real, allocatable :: src_lats(:)
   real, allocatable :: ax(:,:,:)                        ! weight coef. all domain
   real, allocatable :: by(:,:,:)                        ! weight coef. all domain

   character(len=19)   :: proj_name(0:3) = (/ 'LATLON             ', 'LAMBERT            ', &
                                              'POLAR STEREOGRAPHIC', 'MERCATOR           ' /)
   CHARACTER (LEN=192) :: filespec
   CHARACTER (LEN=132) :: varname
   CHARACTER (LEN=132) :: inpname
   CHARACTER (LEN=80)  :: message
   CHARACTER (LEN=80)  :: attribute
   CHARACTER (LEN=80)  :: units_attribute
   CHARACTER (LEN=80)  :: description_attribute
   CHARACTER (LEN=80)  :: stagger_attribute
   CHARACTER (LEN=80)  :: coor_attribute
   CHARACTER (LEN=80)  :: memord_attribute
   CHARACTER (LEN=80)  :: outpname
   CHARACTER (LEN=19)  :: Times(1)
   CHARACTER (LEN=3)   :: num

   logical :: has_area_map
   logical :: new_grid
   logical :: lexist
   logical, allocatable :: cat_active(:)

   type(dates)                :: start_output
   type(dates)                :: stop_output
   type(dates)                :: loop_date
   type(glb_att), allocatable :: attrs(:)
   type(data_file_type), allocatable :: data_file(:)

!---------------------------------------------------------------------
!  include files
!---------------------------------------------------------------------
   include 'netcdf.inc'

!---------------------------------------------------------------------
!  set namelist variables default values
!---------------------------------------------------------------------
   src_names(:) = ' '
   wrf_dir    = '.'
   anthro_dir = '.'
   src_file_prefix = ' '
   src_file_suffix = ' '
   cat_var_prefix  = ' '
   cat_var_suffix  = ' '
   sub_categories(:) = ' ' 
   emis_map(:)       = ' ' 
!-----------------------------------------------------------------
!     read control variables
!-----------------------------------------------------------------
   read(*,nml=control,iostat=istat)
   if( istat /= 0 ) then
     write(*,*) 'main_bc_wrfchem: failed to read namelist; error = ',istat
     stop
   end if
!-----------------------------------------------------------------
!     check namelist inputs
!-----------------------------------------------------------------
   if( domains < 1 ) then
     write(*,*) 'anthro_emis: domains must be >= 1'
     stop 'Namelist err'
   endif
   do n = 1,domains
     filespec = trim(wrf_dir) // '/wrfinput_d'
     write(filespec(len_trim(filespec)+1:),'(i2.2)') n
     inquire( file=trim(filespec),exist=lexist )
     if( .not. lexist ) then
       write(*,*) 'anthro_emis: wrf input file'
       write(*,*) trim(filespec)
       write(*,*) 'anthro_emis: does not exist'
       stop 'File err'
     endif
   end do
   if( emissions_zdim_stag < 1 ) then
     write(*,*) 'anthro_emis: emissions zdim must >= 1'
     stop 'Namelist err'
   endif
   if( serial_output .and. output_interval < 0 ) then
     write(*,*) 'anthro_emis: output_interval must be >= 0'
     stop 'Input parm error'
   endif

!-----------------------------------------------------------------
!     set maps
!-----------------------------------------------------------------
   call mapper( nemis, emis_map, sub_categories )

   do n = 1,maxsrc
     if( src_names(n) /= ' ' ) then
       filespec = trim(anthro_dir) // '/' // trim(src_file_prefix) &
                                   // trim(src_names(n)) // trim(src_file_suffix)
       inquire( file=trim(filespec),exist=lexist )
       if( .not. lexist ) then
         write(*,*) 'anthro_emis: anthro source file'
         write(*,*) trim(filespec)
         write(*,*) 'anthro_emis: does not exist'
         stop 'File err'
       endif
     else
       exit
     endif
   end do

   write(*,*) 'main: nemis = ',nemis
   do n = 1,nemis
     write(*,*) '================================================='
     write(*,*) 'anthro_map(',n,'):'
     write(*,*) 'src species count = ',anthro_map(n)%src_cnt
     write(*,*) 'emis species name = ',anthro_map(n)%emis_name
     write(*,*) anthro_map(n)%src_wght(:anthro_map(n)%src_cnt)
     write(*,*) anthro_map(n)%src_var(:anthro_map(n)%src_cnt)
     write(*,*) 'cat wghts'
     do ns = 1,anthro_map(n)%src_cnt
       write(*,*) anthro_map(n)%cat_wght(:,ns)
     end do
   end do
   write(*,*) ' '
   write(*,*) 'src active'
   write(*,*) src_active(:nsrc_names)
   write(*,*) 'active src = ',count(src_active(:))
   write(*,*) '================================================='

   allocate( cat_active(n_sub_cats),stat=astat )
   if( astat /= 0 ) then
     write(*,*) 'anthro_emis: failed to allocate cat_active array'
     stop 'Alloc err'
   endif
   allocate( data_file(nsrc_names),stat=astat )
   if( astat /= 0 ) then
     write(*,*) 'anthro_emis: failed to allocate data_file type'
     stop 'Alloc err'
   endif

domain_loop : &
   do domain = 1,domains
     call wrf_file
!-----------------------------------------------------------------
!     setup data, output times
!-----------------------------------------------------------------
first_domain : &
     if( domain == 1 ) then
       if( start_output_time == ' ' ) then
         start_output_time = Times(1)
       endif
       call wrf2mz_time( start_output_time, start_output%date, start_output%secs )
       if( .not. serial_output ) then
         stop_output_time = start_output_time
       elseif( stop_output_time == ' ' ) then
         stop_output_time = start_output_time
       endif
       call wrf2mz_time( stop_output_time, stop_output%date, stop_output%secs )
       if( diffdat( start_output%date, start_output%secs, stop_output%date, stop_output%secs ) < 0. ) then
         write(*,*) 'anthro_emis: start output time > stop output time'
         stop 'Input parameter error'
       endif
       do src = 1,nsrc_names
         cat_active(:) = .false.
         if( src_active(src) ) then
           do ns = 1,nemis
             do ns1 = 1,anthro_map(ns)%src_cnt
               if( trim(anthro_map(ns)%src_var(ns1)) ==  trim(src_names(src)) ) then
                 cat_active(:) = cat_active(:) .or. anthro_map(ns)%cat_wght(:,ns1) /= 0.
               endif
             end do
           end do
           allocate( data_file(src)%cat_active(n_sub_cats),stat=astat )
           if( astat /= 0 ) then
             write(*,*) 'anthro_emis: allocate for data_file%cat_active failed; error = ', astat
            stop 'Alloc err'
           endif
           data_file(src)%cat_active(:) = cat_active(:)
           write(*,*) 'will use sub cats: ',cat_active(:)
         endif
       end do
     endif first_domain

     if( .not. allocated( wrk_emis ) ) then
       allocate( wrk_emis(ide,jde),stat=astat )
       if( astat /= 0 ) then
         write(*,*) 'anthro_emis: allocate for wrk_emis failed; error = ', astat
        stop 'Alloc err'
       endif
     endif
!-----------------------------------------------------------------
!     initialize data file type
!-----------------------------------------------------------------
     do src = 1,nsrc_names
       if( src_active(src) ) then
         data_file(src)%filename = trim(src_file_prefix) // trim(src_names(src)) // trim(src_file_suffix)
         data_file(src)%filespec = trim(anthro_dir) // '/' // trim(data_file(src)%filename)
         if( domain == 1 ) then
           data_file(src)%molecw   = molecw(src)
         endif
         call data_file_init( data_file(src), start_output, cat_var_prefix, cat_var_suffix, &
                              domain, dx, ide, jde )
       endif
     end do
!-----------------------------------------------------------------
!     allocate emission array
!-----------------------------------------------------------------
     do ns = 1,nemis
       if( allocated( anthro_map(ns)%emission ) ) then
         deallocate( anthro_map(ns)%emission )
       endif
       allocate( anthro_map(ns)%emission(ide,jde),stat=astat )
       if( astat /= 0 ) then
         write(*,*) 'anthro_emis: failed to allocate emission array'
         stop 'Alloc err'
       endif
     end do

     loop_date = start_output
time_loop : &
     do
       do ns = 1,nemis
         anthro_map(ns)%emission(:,:) = 0.
       end do
src_loop : &
       do src = 1,nsrc_names
use_src : &
         if( src_active(src) ) then
           write(*,*) ' '
           write(*,*) 'will use source file for ',trim(src_names(src))

           call get_src_time_ndx( data_file(src), loop_date )

           cat_ndx = 0
cat_loop : &
           do nc = 1,n_sub_cats
use_cat :    if( data_file(src)%cat_active(nc) ) then
               varname = trim(cat_var_prefix) // trim(sub_cats(nc)) // trim(cat_var_suffix)
               cat_ndx = cat_ndx + 1
               if( data_file(src)%read_lo_tndx ) then
                 call read_src_data( data_file(src), varname, data_file(src)%lo_buf_ndx, &
                                     data_file(src)%lo_tndx, data_file(src)%ncid_lo, cat_ndx )
               endif
               if( data_file(src)%t_interp ) then
                 if( data_file(src)%read_hi_tndx ) then
                   call read_src_data( data_file(src), varname, data_file(src)%hi_buf_ndx, &
                                       data_file(src)%hi_tndx, data_file(src)%ncid_hi, cat_ndx )
                 endif
               endif
               call tinterp_src_data( data_file(src), cat_ndx )
               call map_src_emissions( data_file(src), grid_specs(data_file(src)%grid_ndx) )
               do ns = 1,nemis
                 do ns1 = 1,anthro_map(ns)%src_cnt
                   if( trim(anthro_map(ns)%src_var(ns1)) ==  trim(src_names(src)) ) then
                     if( anthro_map(ns)%cat_wght(nc,ns1) /= 0.) then
                       if( anthro_map(ns)%is_gas ) then
                         scale_factor = anthro_map(ns)%src_wght(ns1) &
                                        * anthro_map(ns)%cat_wght(nc,ns1) * data_file(src)%con_fac(1) &
                                        / data_file(src)%molecw
                       else
                         scale_factor = anthro_map(ns)%src_wght(ns1) &
                                        * anthro_map(ns)%cat_wght(nc,ns1) * data_file(src)%con_fac(2)
                       endif
                       anthro_map(ns)%emission(:,:) = anthro_map(ns)%emission(:,:) &
                                                    + scale_factor*wrk_emis(:,:)
                     endif
                   endif
                 end do
               end do
             endif use_cat
           end do cat_loop

           if( diag_level > 300 ) then
             do ns = 1,nemis
             do ns1 = 1,anthro_map(ns)%src_cnt
               if( trim(anthro_map(ns)%src_var(ns1)) ==  trim(src_names(src)) ) then
                 do nc = 1,n_sub_cats
                   if( cat_active(nc) .and. anthro_map(ns)%cat_wght(nc,ns1) /= 0.) then
                     write(*,*) 'will use sub cats ',trim(sub_cats(nc)),' with wght = ',anthro_map(ns)%cat_wght(nc,ns1)
                   endif
                 end do
               endif
             end do
             end do
           endif
         endif use_src
       end do src_loop
!---------------------------------------------------------------------
!   write emission file
!---------------------------------------------------------------------
       if( serial_output ) then
         call write_emis( 0 )
         call addsec2dat( output_interval, loop_date%date, loop_date%secs )
         if( diffdat( loop_date%date, loop_date%secs, stop_output%date, stop_output%secs ) < 0. ) then
           exit time_loop
         endif
       else
         call write_emis( 1 )
         call write_emis( 2 )
         exit time_loop
       endif
     end do time_loop
     
!---------------------------------------------------------------------
!   cleanup for next domain
!---------------------------------------------------------------------
     do src = 1,nsrc_names
       if( data_file(src)%ncid_lo /= 0 ) then
         ierr = nf_close( data_file(src)%ncid_lo )
       endif
       if( data_file(src)%ncid_hi /= 0 .and. &
           data_file(src)%ncid_hi /= data_file(src)%ncid_lo ) then
         ierr = nf_close( data_file(src)%ncid_hi )
       endif
       if( allocated( data_file(src)%emis ) ) then
         deallocate( data_file(src)%emis )
       endif
       if( allocated( data_file(src)%src_data ) ) then
         deallocate( data_file(src)%src_data )
       endif
       call cleanup_grid( grid_specs(data_file(src)%grid_ndx), ide, jde )
     end do
     grid_cnt = 0
     do ns = 1,nemis
       deallocate( anthro_map(ns)%emission )
     end do
     deallocate( wrk_emis )
   end do domain_loop

   do ns = 1,nemis
     if( allocated( anthro_map(ns)%cat_wght ) ) then
       deallocate( anthro_map(ns)%cat_wght )
     endif
   end do
   if( allocated( anthro_map ) ) then
     deallocate( anthro_map )
   endif
   if( allocated( data_file ) ) then
     deallocate( data_file )
   endif

   write(*,*) ' '
   write(*,*) '----------------------------------'
   write(*,*) 'anthro_emis completed successfully'
   write(*,*) '----------------------------------'

   CONTAINS

   subroutine wrf_file
!---------------------------------------------------------------------
!   read wrf file
!---------------------------------------------------------------------

   use area_mapper, only : proj_init
   use constants_module, only : rad_per_deg, earth_radius_m

!---------------------------------------------------------------------
!   local variables
!---------------------------------------------------------------------
   real, allocatable :: wrk(:,:)

   inpname = 'wrfinput_d'
   write(inpname(len_trim(inpname)+1:),'(i2.2)') domain
   filespec = trim( wrf_dir ) // '/' // trim( inpname )
!---------------------------------------------------------------------
!   open wrf input file
!---------------------------------------------------------------------
   message = 'wrf_file: Failed to open ' // trim(inpname)
   call handle_ncerr( nf_open( trim(filespec), nf_noclobber, ncid ), message )       
!---------------------------------------------------------------------
!   get wrf dimesions
!---------------------------------------------------------------------
   message = 'Failed to get lon dimension id'
   call handle_ncerr( nf_inq_dimid( ncid, 'west_east', dimid ), message )
   message = 'Failed to get lon dimension'
   call handle_ncerr( nf_inq_dimlen( ncid, dimid, ide ), message )
   message = 'Failed to get lat dimension id'
   call handle_ncerr( nf_inq_dimid( ncid, 'south_north', dimid ), message )
   message = 'Failed to get lat dimension'
   call handle_ncerr( nf_inq_dimlen( ncid, dimid, jde ), message )
!---------------------------------------------------------------------
!   get wrf map projection variables
!---------------------------------------------------------------------
   message = 'Failed to get map_proj'
   call handle_ncerr( nf_get_att_int( ncid, nf_global, 'MAP_PROJ', map_proj ), message )
   write(*,*) ' '
   write(*,*) 'wrf_file: MAP_PROJ is ',trim(proj_name(map_proj))
   write(*,*) ' '
   if( map_proj == PROJ_LATLON .or. map_proj == PROJ_CASSINI ) then
     message = 'Failed to get XLONG variable id'
     call handle_ncerr( nf_inq_varid( ncid, 'XLONG', varid ), message )
     allocate( wrk(ide,jde),stat=astat )
     if( astat /= 0 ) then
       write(*,*) 'wrf_file: failed to allocate wrk variable; error = ',astat
       stop 'Alloc error'
     endif
     message = 'Failed to read XLONG variable'
     call handle_ncerr( nf_get_var_real( ncid, varid, wrk ), message )
     cen_lon = wrk(1,1)
     loninc = wrk(2,1) - cen_lon
     message = 'Failed to get XLAT variable id'
     call handle_ncerr( nf_inq_varid( ncid, 'XLAT', varid ), message )
     message = 'Failed to read XLAT variable'
     call handle_ncerr( nf_get_var_real( ncid, varid, wrk ), message )
     cen_lat = wrk(1,1)
     latinc = wrk(1,2) - cen_lat
     knowni = 1.
     knownj = 1.
     deallocate( wrk )
     map_proj = PROJ_LATLON
     dx = earth_radius_m * latinc * rad_per_deg
   else
     message = 'wrf_file: Failed to get cen_lon'
     call handle_ncerr( nf_get_att_real( ncid, nf_global, 'CEN_LON', cen_lon ), message )
     write(*,*) 'wrf_file: CEN_LON = ',cen_lon
     message = 'wrf_file: Failed to get cen_lat'
     call handle_ncerr( nf_get_att_real( ncid, nf_global, 'CEN_LAT', cen_lat ), message )
     write(*,*) 'wrf_file: CEN_LAT = ',cen_lat
     message = 'wrf_file: Failed to get stand_lon'
     call handle_ncerr( nf_get_att_real( ncid, nf_global, 'STAND_LON', stand_lon ), message )
     write(*,*) 'wrf_file: STAND_LON = ',stand_lon
     message = 'wrf_file: Failed to get truelat1'
     call handle_ncerr( nf_get_att_real( ncid, nf_global, 'TRUELAT1', truelat1 ), message )
     write(*,*) 'wrf_file: TRUELAT1 = ',truelat1
     message = 'wrf_file: Failed to get truelat2'
     call handle_ncerr( nf_get_att_real( ncid, nf_global, 'TRUELAT2', truelat2 ), message )
     write(*,*) 'wrf_file: TRUELAT2 = ',truelat2
     message = 'wrf_file: Failed to get dx'
     call handle_ncerr( nf_get_att_real( ncid, nf_global, 'DX', dx ), message )
     write(*,*) 'wrf_file: DX = ',dx
   endif

!---------------------------------------------------------------------
!   initialize map projection
!---------------------------------------------------------------------
   call proj_init( map_proj, cen_lon, cen_lat, truelat1, truelat2, &
                   stand_lon, loninc, latinc, knowni, knownj, &
                   dx, ide, jde )

   ids = 1
   jds = 1

   message = 'Failed to get Times id'
   call handle_ncerr( nf_inq_varid( ncid, 'Times', varid ), message )
   message = 'Failed to read Times'
   call handle_ncerr( nf_get_var_text( ncid, varid, Times ), message )

   write(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
   write(*,*) 'wrf_file: time = ',trim(Times(1))
   write(*,*) 'wrf_file: grid dimensions'
   write(*,*) 'wrf_file: ids,ide,jds,jde'
   write(*,'(4i6)') ids,ide,jds,jde
   write(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
!---------------------------------------------------------------------
!   get wrfinput_<domain> global attributes
!---------------------------------------------------------------------
   call get_glb_atts
!---------------------------------------------------------------------
!   close wrf file
!---------------------------------------------------------------------
   message = 'wrf_file: Failed to close ' // trim(inpname)
   call handle_ncerr( nf_close( ncid ), message )       

   end subroutine wrf_file

   subroutine map_src_emissions( data_file, grid )
!---------------------------------------------------------------------
!   map src dataset to wrf grid
!---------------------------------------------------------------------

   use area_mapper, only : area_interp
   use constants_module, only : rad_per_deg, earth_radius_m
   use mapper_types

!---------------------------------------------------------------------
!   dummy arguments
!---------------------------------------------------------------------
    type(data_file_type), intent(in) :: data_file
    type(grid_type), intent(in)      :: grid

!---------------------------------------------------------------------
!   local variables
!---------------------------------------------------------------------
    integer :: il, iu, jl, ju
    integer :: n, varid
    integer :: status
    real    :: data_dx
    real    :: wrf_lon, wrf_lat
    real, allocatable :: raw_data(:,:)
    logical :: debug = .false.

is_area_map : &
    if( grid%has_area_map ) then
!---------------------------------------------------------------------
!   area conserving interpolation
!---------------------------------------------------------------------
      call area_interp( data_file, grid, grid%model_area_type, wrk_emis, diag_level )
    else is_area_map
!---------------------------------------------------------------------
!   form src coordinate limits
!---------------------------------------------------------------------
      wrf_lon_min = minval( xlong(ids:ide,jds:jde) )
      wrf_lon_max = maxval( xlong(ids:ide,jds:jde) )
      wrf_lat_min = minval( xlat(ids:ide,jds:jde) )
      wrf_lat_max = maxval( xlat(ids:ide,jds:jde) )
      if( diag_level > 100 ) then
        write(*,*) ' '
        write(*,'('' map_src_emissions: model lon limits = '',1p,2g25.16)') wrf_lon_min,wrf_lon_max
        write(*,'('' map_src_emissions: model lat limits = '',1p,2g25.16)') wrf_lat_min,wrf_lat_max
        write(*,*) ' '
        write(*,'('' map_src_emissions: src lon limits = '',1p,2g25.16)') grid%lon(1),grid%lon(grid%nlons)
        write(*,'('' map_src_emissions: src lat limits = '',1p,2g25.16)') grid%lat(1),grid%lat(grid%nlats)
        write(*,*) ' '
      endif
!---------------------------------------------------------------------
!   form dataset index limits
!---------------------------------------------------------------------
      write(*,*) 'map_src_emissions: count of points <,> data min,max lat = ',count(grid%ix(:,:,0) == grid%nlons )
      xndx_src(1) = minval( grid%ix(:,:,lower) )
      xndx_src(2) = maxval( grid%ix(:,:,upper) )
      write(*,*) 'xndx_src = ',xndx_src(:)
      write(*,*) 'map_src_emissions: count of points < data min lat = ',count(grid%jy(:,:,0) == -1)
      write(*,*) 'map_src_emissions: count of points > data max lat = ',count(grid%jy(:,:,0) == -2)
      yndx_src(1) = minval( grid%jy(:,:,lower),mask=grid%jy(:,:,lower)>0 )
      yndx_src(2) = maxval( grid%jy(:,:,upper),mask=grid%jy(:,:,upper)>0 )
      write(*,*) 'yndx_src = ',yndx_src(:)

      if( debug ) then
      write(*,*) ' '
      write(*,*) 'map_src_emissions: bilinear interp diagnostics'
      write(*,*) 'map_src_emissions: ix'
      write(*,*) grid%ix(ids,jds,:)
      write(*,*) 'map_src_emissions: ax'
      write(*,*) grid%ax(ids,jds,:)
      write(*,*) 'map_src_emissions: src lons'
      write(*,*) grid%lon(grid%ix(ids,jds,0)),grid%lon(grid%ix(ids,jds,1))
      write(*,*) 'map_src_emissions: wrf lon = ',xlong(ids,jds)
      write(*,*) 'map_src_emissions: jy'
      write(*,*) grid%jy(ids,jds,:)
      write(*,*) 'map_src_emissions: by'
      write(*,*) grid%by(ids,jds,:)
      write(*,*) 'map_src_emissions: src lats'
      write(*,*) grid%lat(grid%jy(ids,jds,0)),grid%lat(grid%jy(ids,jds,1))
      write(*,*) 'map_src_emissions: wrf lat = ',xlat(ids,jds)
      write(*,*) ' '
      do j = jds,jde
        do i = ids,ide
          if( grid%ix(i,j,lower) == grid%nlons ) then
      write(*,*) 'map_src_emissions: bilinear interp diagnostics'
      write(*,*) 'map_src_emissions: ix'
      write(*,*) grid%ix(i,j,:)
      write(*,*) 'map_src_emissions: ax'
      write(*,*) grid%ax(i,j,:)
      write(*,*) 'map_src_emissions: src lons'
      write(*,*) grid%lon(grid%ix(i,j,0)),grid%lon(grid%ix(i,j,1))
      write(*,*) 'map_src_emissions: wrf lon = ',xlong(i,j)
      write(*,*) 'map_src_emissions: jy'
      write(*,*) grid%jy(i,j,:)
      write(*,*) 'map_src_emissions: by'
      write(*,*) grid%by(i,j,:)
      write(*,*) 'map_src_emissions: src lats'
      write(*,*) grid%lat(grid%jy(i,j,0)),grid%lat(grid%jy(i,j,1))
      write(*,*) 'map_src_emissions: wrf lat = ',xlat(i,j)
            stop 'diagnostics'
          endif
        end do
      end do
      stop 'diagnostics'
      endif

!---------------------------------------------------------------------
!   allocate and read dataset variable
!---------------------------------------------------------------------
       if( .not. allocated( raw_data ) ) then
         allocate( raw_data(xndx_src(1):xndx_src(2),yndx_src(1):yndx_src(2)),stat=status )
         if( status /= 0 ) then
           write(*,*) 'map_src_emissions: allocate for raw_data failed; error = ', status
          stop 'Alloc err'
         endif
       endif
!-------------------------------------------------------------
!  transfer shifted raw data to final working array
!-------------------------------------------------------------
       do j = yndx_src(1),yndx_src(2)
         raw_data(:,j) = data_file%src_data(xndx_src(1):xndx_src(2),j)
       end do

       if( diag_level > 100 ) then
         write(*,*)  'dataset size = ',size(data_file%src_data)
         write(*,*)  'dataset min,max values = ',minval(data_file%src_data(:,:)),maxval(data_file%src_data(:,:))
       endif

!---------------------------------------------------------------------
!   set wrf anthro emission
!---------------------------------------------------------------------
       do j = jds,jde
         do i = ids,ide
           jl = grid%jy(i,j,0)
           if( jl > 0 ) then
             il = grid%ix(i,j,0)
             iu = grid%ix(i,j,1)
             ju = grid%jy(i,j,1)
             wrk_sum  = raw_data(il,jl)*grid%ax(i,j,upper)*grid%by(i,j,upper) &
                      + raw_data(il,ju)*grid%ax(i,j,upper)*grid%by(i,j,lower) &
                      + raw_data(iu,jl)*grid%ax(i,j,lower)*grid%by(i,j,upper) &
                      + raw_data(iu,ju)*grid%ax(i,j,lower)*grid%by(i,j,lower)
           else
             wrk_sum  = 0.
           endif
           wrk_emis(i,j) = wrk_sum
         end do
       end do
       if( allocated( raw_data ) ) then
         deallocate( raw_data )
       endif
    endif is_area_map

    end subroutine map_src_emissions

    subroutine write_emis( nfile )
!---------------------------------------------------------------------
!	... write the netcdf anthro emission file
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!	... dummy arguments
!---------------------------------------------------------------------
      integer, intent(in) :: nfile

!---------------------------------------------------------------------
!	... local variables
!---------------------------------------------------------------------
      integer :: lon_id
      integer :: lat_id
      integer :: time_id
      integer :: zdim_id
      integer :: string_id
      integer :: dims(4)
      integer :: start_ndx(4)
      integer :: length(4)
      integer :: astat
      integer :: m, nt, nt1
      real, allocatable :: wrk_emis(:,:,:)
      character(len=132) :: message, text
      character(len=10)  :: ctime
      character(len=8)   :: cdate
      character(len=10)  :: t_string(2)
      character(len=19)  :: wrf_time

      write(num(2:3),'(i2.2)') domain
      if( serial_output ) then
        call mz2wrf_time( wrf_time, loop_date%date, loop_date%secs )
        outpname = 'wrfchemi_d' // num(2:3) // '_' // wrf_time
      else
        if( nfile == 1 ) then
          outpname = 'wrfchemi_00z_d' // num(2:3)
        else
          outpname = 'wrfchemi_12z_d' // num(2:3)
        endif
      endif
!-----------------------------------------------------------------------
!     	... create netcdf anthro emission file and enter define mode
!-----------------------------------------------------------------------
      message = 'write_emis: Failed to create ' // trim( outpname )
      call handle_ncerr( nf_create( trim( outpname ), nf_clobber, ncid ), message )

!-----------------------------------------------------------------------
!     	... define the dimensions
!-----------------------------------------------------------------------
      call handle_ncerr( nf_def_dim( ncid, 'west_east', ide, lon_id ), &
                         'write_emis: Failed to define longitude dimension' )
      call handle_ncerr( nf_def_dim( ncid, 'south_north', jde, lat_id ), &
                         'write_emis: Failed to define latitude dimension' )
      call handle_ncerr( nf_def_dim( ncid, 'emissions_zdim_stag', emissions_zdim_stag, zdim_id ), &
                         'write_emis: Failed to define emissions_zdim_stag dimension' )
      call handle_ncerr( nf_def_dim( ncid, 'DateStrLen', 19, string_id ), &
                         'write_emis: Failed to define DateStrLen dimension' )
      call handle_ncerr( nf_def_dim( ncid, 'Time', nf_unlimited, time_id ), &
                         'write_emis: Failed to create Time dimension' )

!-----------------------------------------------------------------------
!     	... define the variables
!-----------------------------------------------------------------------
      dims(1:2) = (/ string_id, time_id /)
      call handle_ncerr( nf_def_var( ncid, 'Times', nf_char, 2, dims(1:2), varid ), &
                         'write_emis: Failed to define Times variable' )
      dims(1:2) = (/ lon_id, lat_id /)
      call handle_ncerr( nf_def_var( ncid, 'XLONG', nf_float, 2, dims(1:2), varid ), &
                         'write_emis: Failed to define XLONG variable' )
      call handle_ncerr( nf_def_var( ncid, 'XLAT', nf_float, 2, dims(1:2), varid ), &
                         'write_emis: Failed to define XLAT variable' )

      dims(:) = (/ lon_id, lat_id, zdim_id, time_id /)
      do m = 1,nemis
        message = 'write_emis: Failed to define ' // trim(anthro_map(m)%emis_name)
        call handle_ncerr( nf_def_var( ncid, 'E_'//trim(anthro_map(m)%emis_name), nf_float, 4, dims, varid ), &
                           trim(message) )
      end do

!-----------------------------------------------------------------------
!   ... define variable attributes
!-----------------------------------------------------------------------
      varname = 'XLONG'
      units_attribute       = 'degree east'
      description_attribute = 'LONGITUDE, WEST IS NEGATIVE'
      stagger_attribute     = ''
      memord_attribute      = 'XY '
      coor_attribute        = ' '
      call write_attributes

      varname = 'XLAT'
      units_attribute       = 'degree north'
      description_attribute = 'LATITUDE, SOUTH IS NEGATIVE'
      stagger_attribute     = ''
      memord_attribute      = 'XY '
      coor_attribute        = ' '
      call write_attributes

      do m = 1,nemis
        varname = 'E_' // trim(anthro_map(m)%emis_name)
        if( anthro_map(m)%is_gas ) then
          units_attribute       = 'mol km^-2 hr^-1'
        else
          units_attribute       = 'ug m^-2 s^-1'
        endif
        description_attribute = 'EMISSIONS'
        stagger_attribute     = 'Z'
        memord_attribute      = 'XYZ'
        call write_attributes
      end do

!-----------------------------------------------------------------------
!   ... define global attributes
!-----------------------------------------------------------------------
      message = 'global_attributes: Failed to write title'
      text    = 'Anthropogenic emissions'
      call handle_ncerr( nf_put_att_text( ncid, nf_global, 'Title', len_trim(text), trim(text) ), message )
      message = 'global_attributes: Failed to write History'
      call date_and_time( cdate, ctime )
      t_string(1) = cdate(1:4) // '-' // cdate(5:6) // '-' // cdate(7:8)
      t_string(2) = ctime(1:2) // ':' // ctime(3:4)
      text    = 'Created on ' // trim(t_string(1)) // ' at ' // trim(t_string(2))
      call handle_ncerr( nf_put_att_text( ncid, nf_global, 'History', len_trim(text), trim(text) ), message )
      message = 'global_attributes: Author to write Files'
      text    = 'anthro_emis'
      call handle_ncerr( nf_put_att_text( ncid, nf_global, 'Author', len_trim(text), trim(text) ), message )

      if( ngatts > 0 ) then
        call set_glb_atts
      endif

!-----------------------------------------------------------------------
!     	... leave define mode
!-----------------------------------------------------------------------
      call handle_ncerr( nf_enddef( ncid ), 'write_emis: Failed to leave define mode' )

      allocate( wrk_emis(ide,jde,emissions_zdim_stag),stat=astat )
      if( astat /= 0 ) then
        write(*,*) 'write_emis: failed to allocate wrk_emis; error = ',astat
        stop 'Alloc err'
      endif

!-----------------------------------------------------------------------
!     	... write the variables
!-----------------------------------------------------------------------
      start_ndx(1:2) = (/ 1,1 /)
      length(1:2)    = (/ 19, 1 /)
      call handle_ncerr( nf_inq_varid( ncid, 'Times', varid ), &
                         'write_emis: Failed to get Times variable id' )
      if( serial_output ) then
        call handle_ncerr( nf_put_vara_text( ncid, varid, start_ndx(:2), length(:2), wrf_time ), &
                           'write_emis: Failed to write Times variable' )
      else
        wrf_time = Times(1)
        do nt = 0,11
          start_ndx(2) = nt+1
          write(wrf_time(12:13),'(i2.2)') (nfile-1)*12 + nt
          call handle_ncerr( nf_put_vara_text( ncid, varid, start_ndx(:2), length(:2), wrf_time ), &
                             'write_emis: Failed to write Times variable' )
        end do
      endif

      call handle_ncerr( nf_inq_varid( ncid, 'XLONG', varid ), &
                         'write_emis: Failed to get xlong variable id' )
      call handle_ncerr( nf_put_var_real( ncid, varid, xlong ), &
                         'write_emis: Failed to write xlong variable' )
      call handle_ncerr( nf_inq_varid( ncid, 'XLAT', varid ), &
                         'write_emis: Failed to get xlat variable id' )
      call handle_ncerr( nf_put_var_real( ncid, varid, xlat ), &
                         'write_emis: Failed to write xlat variable' )
      start_ndx(:) = 1
      length(:)    = (/ ide, jde, emissions_zdim_stag, 1 /)
      do m = 1,nemis
        message = 'write_emis: Failed to write ' // trim(anthro_map(m)%emis_name)
        call handle_ncerr( nf_inq_varid( ncid, 'E_'//trim(anthro_map(m)%emis_name), varid ), &
                           trim(message) )
        wrk_emis(:,:,2:emissions_zdim_stag) = 0.
        wrk_emis(:,:,1)    = anthro_map(m)%emission(:,:)
        if( serial_output ) then
          start_ndx(4) = 1
          call handle_ncerr( nf_put_vara_real( ncid, varid, start_ndx, length, wrk_emis ), &
                             trim(message) )
        else
          do nt = 1,12
            start_ndx(4) = nt
            call handle_ncerr( nf_put_vara_real( ncid, varid, start_ndx, length, wrk_emis ), &
                               trim(message) )
          end do
        endif
      end do
!---------------------------------------------------------------------
!   close wrf file
!---------------------------------------------------------------------
      message = 'Failed to close ' // trim(outpname)
      call handle_ncerr( nf_close( ncid ), message )       

      deallocate( wrk_emis )

   end subroutine write_emis

   subroutine write_attributes
!---------------------------------------------------------------------
!   write common variable attributes
!---------------------------------------------------------------------

      message = 'write_attributes: Failed to get ' // trim(varname) // ' variable id'
      call handle_ncerr( nf_inq_varid( ncid, trim(varname), varid ), message )
      message = 'write_attributes: Failed to create ' // trim(varname) // ' attribute'
      call handle_ncerr( nf_put_att_text( ncid, varid, 'MemoryOrder', 3, memord_attribute ), message )
      call handle_ncerr( nf_put_att_text( ncid, varid, 'description', &
                                          len_trim(description_attribute), trim(description_attribute) ), message )
      call handle_ncerr( nf_put_att_text( ncid, varid, 'units', &
                                          len_trim(units_attribute), trim(units_attribute) ), message )
      call handle_ncerr( nf_put_att_text( ncid, varid, 'stagger', &
                                          len_trim(stagger_attribute), trim(stagger_attribute) ), message )
      if( coor_attribute /= ' ' ) then
         call handle_ncerr( nf_put_att_text( ncid, varid, 'coordinates', &
                                             len_trim(coor_attribute), trim(coor_attribute) ), message )
      endif
      ii = 104
      call handle_ncerr( nf_put_att_int( ncid, varid, 'FieldType', nf_int, 1, ii ), message )

   end subroutine write_attributes

   subroutine get_glb_atts
!---------------------------------------------------------------------
!   read the global attributes
!---------------------------------------------------------------------

   implicit none

!---------------------------------------------------------------------
!   local variables
!---------------------------------------------------------------------
   integer              :: m
   integer              :: astat
   integer              :: attr_len
   character(len=132)   :: message
   character(len=132)   :: attr_name

!---------------------------------------------------------------------
!   get global attr count
!---------------------------------------------------------------------
   message = 'glb_attr: Failed to get glb attr count'
   call handle_ncerr( nf_inq_natts( ncid, ngatts ), message )       
!---------------------------------------------------------------------
!   allocate variables
!---------------------------------------------------------------------
   if( ngatts > 0 ) then
     if( allocated( attrs ) ) then
       deallocate( attrs )
     endif
     allocate( attrs(ngatts),stat=astat )
     if( astat /= 0 ) then
       write(*,*) 'glb_attr: failed to allocate type glb_att'
       stop 'Alloc err'
     endif
     attrs(:)%name = ' '
   endif
!---------------------------------------------------------------------
!   loop over glb attributes
!---------------------------------------------------------------------
glb_attr_loop : &
   do m = 1,ngatts
     write(message,*) 'glb_attr: Failed to get glb attr # ',m,' name'
     call handle_ncerr( nf_inq_attname( ncid, nf_global, m, attr_name ), message )       
     attrs(m)%name = attr_name
     write(message,*) 'glb_attr: Failed to get glb attr # ',m,' type,len'
     call handle_ncerr( nf_inq_att( ncid, nf_global, trim(attr_name), attrs(m)%type, attr_len ), message )       
     attrs(m)%len = attr_len
     message = 'glb_attr: Failed to get ' // trim(attr_name)
     select case( attrs(m)%type )
       case( nf_byte )
         allocate( attrs(m)%attr_byte(attr_len),stat=astat )
         if( astat /= 0 ) then
           write(*,*) 'glb_attr: failed to allocate attr_byte'
           stop 'Alloc err'
         endif
         call handle_ncerr( nf_get_att_int1( ncid, nf_global, trim(attr_name), attrs(m)%attr_byte ), message )       
       case( nf_char )
         attrs(m)%attr_char = ' '
         call handle_ncerr( nf_get_att_text( ncid, nf_global, trim(attr_name), attrs(m)%attr_char ), message )       
       case( nf_short )
         allocate( attrs(m)%attr_short(attr_len),stat=astat )
         if( astat /= 0 ) then
           write(*,*) 'glb_attr: failed to allocate attr_short'
           stop 'Alloc err'
         endif
         call handle_ncerr( nf_get_att_int2( ncid, nf_global, trim(attr_name), attrs(m)%attr_short ), message )       
       case( nf_int )
         allocate( attrs(m)%attr_int(attr_len),stat=astat )
         if( astat /= 0 ) then
           write(*,*) 'glb_attr: failed to allocate attr_int'
           stop 'Alloc err'
         endif
         call handle_ncerr( nf_get_att_int( ncid, nf_global, trim(attr_name), attrs(m)%attr_int ), message )       
       case( nf_float )
         allocate( attrs(m)%attr_real(attr_len),stat=astat )
         if( astat /= 0 ) then
           write(*,*) 'glb_attr: failed to allocate attr_real'
           stop 'Alloc err'
         endif
         call handle_ncerr( nf_get_att_real( ncid, nf_global, trim(attr_name), attrs(m)%attr_real ), message )       
       case( nf_double )
         allocate( attrs(m)%attr_dbl(attr_len),stat=astat )
         if( astat /= 0 ) then
           write(*,*) 'glb_attr: failed to allocate attr_dbl'
           stop 'Alloc err'
         endif
         call handle_ncerr( nf_get_att_double( ncid, nf_global, trim(attr_name), attrs(m)%attr_dbl ), message )       
     end select
   end do glb_attr_loop

   end subroutine get_glb_atts

   subroutine set_glb_atts
!---------------------------------------------------------------------
!   set the global attributes
!---------------------------------------------------------------------

   implicit none

!---------------------------------------------------------------------
!   local variables
!---------------------------------------------------------------------
   integer              :: m
   integer              :: attr_len
   integer              :: attr_xtype
   integer              :: slen
   character(len=132)   :: message
   character(len=132)   :: attr_name

!---------------------------------------------------------------------
!   loop over glb attributes
!---------------------------------------------------------------------
glb_attr_loop : &
   do m = 1,ngatts
     attr_name = trim(attrs(m)%name)
     if( trim(attr_name) == 'TITLE' .or. trim(attr_name) == 'START_DATE' .or. &
         trim(attr_name) == 'SIMULATION_START_DATE' ) then
       cycle
     endif
     slen      = len_trim(attr_name)
     write(message,*) 'set_glb_att: Failed to define glb att ',trim(attr_name)
     attr_len   = attrs(m)%len
     attr_xtype = attrs(m)%type
     select case( attrs(m)%type )
       case( nf_byte )
         call handle_ncerr( nf_put_att_int1( ncid, nf_global, attr_name(:slen), attr_xtype, attr_len, attrs(m)%attr_byte ), message )       
       case( nf_char )
         call handle_ncerr( nf_put_att_text( ncid, nf_global, attr_name(:slen), attr_len, attrs(m)%attr_char ), message )       
       case( nf_short )
         call handle_ncerr( nf_put_att_int2( ncid, nf_global, attr_name(:slen), attr_xtype, attr_len, attrs(m)%attr_short ), message )       
       case( nf_int )
         call handle_ncerr( nf_put_att_int( ncid, nf_global, attr_name(:slen), attr_xtype, attr_len, attrs(m)%attr_int ), message )       
       case( nf_float )
         call handle_ncerr( nf_put_att_real( ncid, nf_global, attr_name(:slen), attr_xtype, attr_len, attrs(m)%attr_real ), message )       
       case( nf_double )
         call handle_ncerr( nf_put_att_double( ncid, nf_global, attr_name(:slen), attr_xtype, attr_len, attrs(m)%attr_dbl ), message )       
     end select
   end do glb_attr_loop

   end subroutine set_glb_atts

   subroutine dealloc_glb_atts
!---------------------------------------------------------------------
!   deallocate variables
!---------------------------------------------------------------------

   integer :: m

   do m = 1,ngatts
     select case( attrs(m)%type )
       case( nf_byte )
         deallocate( attrs(m)%attr_byte )
       case( nf_short )
         deallocate( attrs(m)%attr_short )
       case( nf_int )
         deallocate( attrs(m)%attr_int )
       case( nf_float )
         deallocate( attrs(m)%attr_real )
       case( nf_double )
         deallocate( attrs(m)%attr_dbl )
     end select
   end do

   end subroutine dealloc_glb_atts

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

   end program anthro_emis
