
   module srf_types

   implicit none

   private
   public :: ntypes
   public :: cover
   public :: plant_cover
   public :: fire_srf_init
   public :: fire_srf_types
   public :: fire_srf_final

   integer, parameter   :: ntypes = 4

   character(len=32)    :: plant_name(ntypes)
   character(len=8)     :: plant_type(ntypes)

   type cover
     real, pointer :: type_frac(:,:,:)
   end type cover

   type(cover), allocatable :: plant_cover(:)

   include 'netcdf.inc'

   contains

   subroutine fire_srf_init( domains )
!---------------------------------------------------------------------
!   initialize module
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!   dummy arguments
!---------------------------------------------------------------------
   integer, intent(in) :: domains

!---------------------------------------------------------------------
!   local variables
!---------------------------------------------------------------------
   integer :: astat

   allocate( plant_cover(domains),stat=astat )
   if( astat /= 0 ) then
     write(*,*) 'fire_srf_init: failed to allocate plant_cover type; error = ',astat
     stop 'Alloc error'
   endif

   end subroutine fire_srf_init

   subroutine fire_srf_types( fire_directory, domain, proj )
!---------------------------------------------------------------------
!   read srf type file
!---------------------------------------------------------------------

   use wrf_utils
   use utils, only : diag_level

!---------------------------------------------------------------------
!   dummy arguments
!---------------------------------------------------------------------
   integer, intent(in)          :: domain
   character(len=*), intent(in) :: fire_directory
   type(proj_info)              :: proj

!---------------------------------------------------------------------
!   local variables
!---------------------------------------------------------------------
   integer            :: file
   integer            :: i, j, n
   integer            :: astat
   integer            :: ncid
   integer            :: dimid, varid
   integer            :: nlon, nlat
   integer            :: lon_ndx, lat_ndx
   integer            :: lon_beg, lon_end
   integer            :: lat_beg, lat_end
   integer            :: zero_cnt
   integer(1), allocatable :: data(:,:)
   integer(1), allocatable :: xdata(:,:,:)
   integer, allocatable    :: popcnt(:,:,:)
   integer, allocatable    :: cell_cnt(:,:)
   real               :: wrf_lon, wrf_lat
   real               :: x, y
   real               :: min_x, max_x
   real               :: min_y, max_y
   real               :: minmax_wrf_lon(2)
   real               :: minmax_wrf_lat(2)
   real, allocatable  :: type_sum(:,:)
   real, allocatable  :: lon(:)
   real, allocatable  :: lat(:)
   character(len=128) :: message
   character(len=256) :: filespec
   logical            :: found

     plant_type(:) = (/ 'tempfor ', 'tropfor ', 'shrub   ', 'grass   ' /)
     plant_name(:) = (/ 'tempfor ', 'tropfor ', 'shrub   ', 'grass   ' /)
!---------------------------------------------------------------------
!   get min,max lon,lat of wrf grid
!---------------------------------------------------------------------
     minmax_wrf_lon(:) = (/ 1000.,-1000. /)
     minmax_wrf_lat(:) = (/ 1000.,-1000. /)

     do j = 1,proj%jde+1
       y = real(j) - .5
       do i = 1,proj%ide+1
         call ijll( real(i)-.5, y, proj, wrf_lat, wrf_lon )
         minmax_wrf_lon(1) = min( minmax_wrf_lon(1),wrf_lon )
         minmax_wrf_lon(2) = max( minmax_wrf_lon(2),wrf_lon )
         minmax_wrf_lat(1) = min( minmax_wrf_lat(1),wrf_lat )
         minmax_wrf_lat(2) = max( minmax_wrf_lat(2),wrf_lat )
       end do
     end do
!---------------------------------------------------------------------
!   allocate module variable
!---------------------------------------------------------------------
!  if( allocated(type_frac) ) then
!    deallocate( type_frac,stat=astat )
!    if( astat /= 0 ) then
!      write(*,*) 'fire_srf_types: failed to deallocate type_frac; error = ',astat
!      stop 'Dealloc error'
!    endif
!  endif
   allocate( plant_cover(domain)%type_frac(proj%ide,proj%jde,ntypes),type_sum(proj%ide,proj%jde),stat=astat )
   if( astat /= 0 ) then
     write(*,*) 'fire_srf_types: failed to allocate type_frac,type_sum; error = ',astat
     stop 'Alloc error'
   endif
   allocate( popcnt(proj%ide,proj%jde,ntypes),cell_cnt(proj%ide,proj%jde),stat=astat )
   if( astat /= 0 ) then
     write(*,*) 'fire_srf_types: failed to allocate popcnt,cell_cnt; error = ',astat
     stop 'Alloc error'
   endif


   cell_cnt(:,:) = 0
   min_x = .5
   max_x = real(proj%ide) + .5
   min_y = .5
   max_y = real(proj%jde) + .5

file_loop : &
     do file = 1,ntypes
       filespec = trim(fire_directory) // trim(plant_type(file)) // '_from_img.nc'
!---------------------------------------------------------------------
!   open input file
!---------------------------------------------------------------------
       message = 'fire_srf_types: Failed to open ' // trim(filespec)
       call handle_ncerr( nf_open( trim(filespec), nf_noclobber, ncid ), message )       
!---------------------------------------------------------------------
!   get wrf dimesions
!---------------------------------------------------------------------
       message = 'fire_srf_types: Failed to get lon dimension id'
       call handle_ncerr( nf_inq_dimid( ncid, 'lon', dimid ), message )
       message = 'fire_srf_types: Failed to get lon dimension'
       call handle_ncerr( nf_inq_dimlen( ncid, dimid, nlon ), message )
       message = 'fire_srf_types: Failed to get lat dimension id'
       call handle_ncerr( nf_inq_dimid( ncid, 'lat', dimid ), message )
       message = 'fire_srf_types: Failed to get lat dimension'
       call handle_ncerr( nf_inq_dimlen( ncid, dimid, nlat ), message )

!---------------------------------------------------------------------
!   allocate arrays
!---------------------------------------------------------------------
       allocate( lon(nlon),lat(nlat),stat=astat )
       if( astat /= 0 ) then
         write(*,*) 'failed to allocate lon,lat; error = ',astat
         stop 'Alloc error'
       endif

!---------------------------------------------------------------------
!   read coordinate arrays
!---------------------------------------------------------------------
       message = 'fire_srf_types: Failed to get lon varid'
       call handle_ncerr( nf_inq_varid( ncid, 'lon', varid ), message )
       message = 'fire_srf_types: Failed to read lon var'
       call handle_ncerr( nf_get_var_real( ncid, varid, lon ), message )
       message = 'fire_srf_types: Failed to get lat varid'
       call handle_ncerr( nf_inq_varid( ncid, 'lat', varid ), message )
       message = 'fire_srf_types: Failed to read lat var'
       call handle_ncerr( nf_get_var_real( ncid, varid, lat ), message )

!---------------------------------------------------------------------
!   establish lon limits for srf cover type file
!---------------------------------------------------------------------
       found = .false.
       do lon_beg = 1,nlon
         if( lon(lon_beg) >= minmax_wrf_lon(1) ) then
           found = .true.
           exit
         endif
       end do

       if( .not. found ) then
         write(*,*) 'fire_srf_types: lon ',minmax_wrf_lon(1),' not in [',lon(1),',',lon(2),']'
         stop 'Range error'
       endif
       lon_beg = max( 1,lon_beg-1 )

       found = .false.
       do lon_end = nlon,lon_beg,-1
         if( lon(lon_end) <= minmax_wrf_lon(2) ) then
           found = .true.
           exit
         endif
       end do

       if( .not. found ) then
         write(*,*) 'fire_srf_types: lon ',minmax_wrf_lon(2),' not in [',lon(1),',',lon(2),']'
         stop 'Range error'
       endif
       lon_end = min( lon_end+1,nlon )

!---------------------------------------------------------------------
!   establish lat limits for srf cover type file
!---------------------------------------------------------------------
       found = .false.
       do lat_end = nlat,1,-1
         if( lat(lat_end) > minmax_wrf_lat(1) ) then
           found = .true.
           exit
         endif
       end do

       if( .not. found ) then
         write(*,*) 'fire_srf_types: lat ',minmax_wrf_lat(1),' not in [',lat(1),',',lat(2),']'
         stop 'Range error'
       endif
       lat_end = min(lat_end + 1,nlat)

       found = .false.
       do lat_beg = 1,lat_end
         if( lat(lat_beg) < minmax_wrf_lat(2) ) then
           found = .true.
           exit
         endif
       end do

       if( .not. found ) then
         write(*,*) 'fire_srf_types: lat ',minmax_wrf_lat(2),' not in [',lat(1),',',lat(2),']'
         stop 'Range error'
       endif
       lat_beg = max(lat_beg - 1,1)

!---------------------------------------------------------------------
!   allocate and read data array
!---------------------------------------------------------------------
       allocate( data(lon_beg:lon_end,lat_beg:lat_end),stat=astat )
       if( astat /= 0 ) then
         write(*,*) 'failed to allocate data; error = ',astat
         stop 'Alloc error'
       endif
       if( diag_level >= 300 .and. file == 1 ) then
         allocate( xdata(lon_beg:lon_end,lat_beg:lat_end,ntypes),stat=astat )
         if( astat /= 0 ) then
           write(*,*) 'fire_srf_types: failed to allocate xdata; error = ',astat
           stop 'Alloc error'
         endif
       endif
       message = 'fire_srf_types: Failed to get ' // trim(plant_name(file)) // ' varid'
       call handle_ncerr( nf_inq_varid( ncid, trim(plant_name(file)), varid ), message )
       message = 'fire_srf_types: Failed to read ' // trim(plant_name(file)) // ' var'
       call handle_ncerr( nf_get_vara_int1( ncid, varid, (/ lon_beg,lat_beg /), &
                                          (/ lon_end-lon_beg+1,lat_end-lat_beg+1 /), data ), message )
       if( diag_level >= 300 ) then
         xdata(:,:,file) = data(:,:)
         i = count( data(:,:) == 1 )
         write(*,*) ' '
         write(*,*) 'fire_srf_types: count of non-zero ',trim(plant_type(file)),' type = ',i
         write(*,*) 'fire_srf_types: size of data = ',size( data )
         write(*,*) 'fire_srf_types: % non-zero ',trim(plant_type(file)),' type = ',100.*real(i)/real(size(data))
       endif

!---------------------------------------------------------------------
!   set pop count
!---------------------------------------------------------------------
       popcnt(:,:,file) = 0
       do j = lat_beg,lat_end
         do i = lon_beg,lon_end
           call llij( lat(j), lon(i), proj, x, y )
           if( min_x <= x .and. x < max_x .and. &
               min_y <= y .and. y < max_y )  then
             lon_ndx = nint(x)
             lat_ndx = nint(y)
             if( file == 1 ) then
               cell_cnt(lon_ndx,lat_ndx) = cell_cnt(lon_ndx,lat_ndx) + 1
             endif
             if( data(i,j) == 1 ) then
               popcnt(lon_ndx,lat_ndx,file) = popcnt(lon_ndx,lat_ndx,file) + 1
               if( diag_level >= 300 ) then
                 if( file > 1 .and. any( xdata(i,j,:file-1) == 1 ) ) then
                   write(*,*) 'fire_srf_types: multiple cover types @ i,j = ',i,j
                   write(*,*) 'fire_srf_types: xdata(i,j) = ',xdata(i,j,:file)
                   stop 'Data problem'
                 endif
               endif
             endif
           endif
         end do
       end do

       if( diag_level >= 300 ) then
         write(*,*) ' '
         write(*,*) 'fire_srf_types: max wrf cell ',trim(plant_type(file)),' cnt = ',maxval(popcnt(:,:,file))
         i = count(popcnt(:,:,file) /= 0)
         write(*,*) 'fire_srf_types: non-zero wrf cell ',trim(plant_type(file)),' cnt = ',i
         write(*,*) 'fire_srf_types: % non-zero wrf cell ',trim(plant_type(file)),' cnt = ',100.*real(i)/real(proj%ide*proj%jde)
       endif

!---------------------------------------------------------------------
!   close wrf file
!---------------------------------------------------------------------
       message = 'fire_srf_types: Failed to close ' // trim(filespec)
       call handle_ncerr( nf_close( ncid ), message )       
!---------------------------------------------------------------------
!   deallocate arrays
!---------------------------------------------------------------------
       deallocate( lon, lat, data, stat=astat )
       if( astat /= 0 ) then
         write(*,*) 'failed to deallocate lon,lat,data arrays; error = ',astat
         stop 'Alloc error'
       endif
     end do file_loop

     if( diag_level >= 300 ) then
       zero_cnt = 0
       do j = lat_beg,lat_end
         do i = lon_beg,lon_end
           if( all(xdata(i,j,:) == 0 ) ) then
             zero_cnt = zero_cnt + 1
           endif
         end do
       end do
       write(*,*) ' '
       write(*,*) 'fire_srf_types: zero cover type count = ',zero_cnt
     endif

!-----------------------------------------------------------------------
!  form the srf cover type variables
!-----------------------------------------------------------------------
     do j = 1,proj%jde
       do i = 1,proj%ide
!        type_sum(i,j) = real(cell_cnt(i,j))
         type_sum(i,j) = real(sum(popcnt(i,j,:)))
       end do
     end do
     do n = 1,ntypes
       where( popcnt(:,:,n) /= 0 )
         plant_cover(domain)%type_frac(:,:,n) = real(popcnt(:,:,n))/type_sum(:,:)
       elsewhere
         plant_cover(domain)%type_frac(:,:,n) = 0.
       endwhere
     end do

     if( diag_level >= 300 ) then
       write(*,*) ' '
       write(*,*) 'fire_srf_types: max cover type fraction'
       write(*,*) 'fire_srf_types: ',maxval(plant_cover(domain)%type_frac(:,:,1)), &
                                     maxval(plant_cover(domain)%type_frac(:,:,2)), &
                                     maxval(plant_cover(domain)%type_frac(:,:,3)), &
                                     maxval(plant_cover(domain)%type_frac(:,:,4))
       write(*,*) 'fire_srf_types: maxloc'
       do n = 1,ntypes
         write(*,*) 'fire_srf_types: ',trim(plant_type(n)),' max @ ',maxloc( plant_cover(domain)%type_frac(:,:,n) )
       end do
     endif

     do j = 1,proj%jde
       do i = 1,proj%ide
        if( sum(popcnt(i,j,:) ) > cell_cnt(i,j) ) then
          write(*,*) 'popcnt > cell_cnt @ i,j = ',i,j
          stop 'Data problem'
        endif
       end do
     end do

!---------------------------------------------------------------------
!   deallocate arrays
!---------------------------------------------------------------------
   if( allocated(xdata) ) then
     deallocate( xdata,stat=astat )
     if( astat /= 0 ) then
       write(*,*) 'fire_srf_types: failed to deallocate xdata; error = ',astat
       stop 'Dealloc error'
     endif
   endif
   deallocate( type_sum, popcnt, cell_cnt, stat=astat )
   if( astat /= 0 ) then
     write(*,*) 'failed to deallocate type_sum,popcnt,cell_cnt arrays; error = ',astat
     stop 'Alloc error'
   endif

   end subroutine fire_srf_types

   subroutine fire_srf_final( domains )
!---------------------------------------------------------------------
!   finalize module
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!   dummy arguments
!---------------------------------------------------------------------
   integer, intent(in) :: domains

!---------------------------------------------------------------------
!   local variables
!---------------------------------------------------------------------
   integer :: domain

   do domain = 1,domains
     deallocate( plant_cover(domain)%type_frac )
   end do
   deallocate( plant_cover )

   end subroutine fire_srf_final

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
      write(*,*) 'handle_ncerr: ',trim(mes)
      write(*,*) nf_strerror( ret )
      stop 'netcdf error'
   endif

   end subroutine handle_ncerr

   end module srf_types
