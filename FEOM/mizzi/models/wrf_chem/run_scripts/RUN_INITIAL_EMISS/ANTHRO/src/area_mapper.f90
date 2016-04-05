
   module area_mapper

   use misc_definitions_module
   use constants_module
   use mapper_types

   implicit none

   private
   public :: lon, lat
   public :: proj_init
   public :: area_interp_init
   public :: area_interp

!-----------------------------------------------------------
!	module variables
!-----------------------------------------------------------
   integer :: i, j, l
   integer :: i1, j1, j2
   integer :: dlong_ndx
   integer :: lon_s, lon_e
   integer :: lat_s, lat_e
   integer :: mlat_ndx, mlon_ndx
   integer :: m_lat_s, m_lon_s
   integer :: m_lat_e, m_lon_e
   integer :: min_lon_ndx(2)
   integer :: max_lon_ndx(2)
   integer :: lon_ndx(2,2)
   integer :: lat_ndx(2,2)
   integer :: min_lat_ndx(2)
   integer :: max_lat_ndx(2)
   integer :: dvtx_shadow_ndx(4)
   integer :: cnt_dvtx_n_mcell
   real    :: tot_area
   real    :: dlon
   real    :: mcell_area
   real    :: pole_lat
   real    :: delta_lat, delta_lon
   real    :: xc, yc
   real    :: target_lon
   real    :: target_lat
   real    :: minmax_lon(2)
   real    :: minmax_lat(2)
   real(8) :: line_parms(2)
   real    :: minmax_x(2)
   real    :: minmax_y(2)
   real    :: model_x(4)
   real    :: model_y(4)
   real    :: model_lon(2,2)
   real    :: model_lat(2,2)
   real    :: data_lon(4)
   real    :: data_lat(4)
   real, allocatable :: lon(:,:)
   real, allocatable :: lat(:,:)
   real    :: x(4), y(4)
   logical :: debugging = .false.
   logical :: has_lon_cross
   logical :: lon_cross(2,2)
   logical :: lat_mask(2,2)
   logical :: x_n_mcell(4)
   logical :: y_n_mcell(4)
   logical :: dvtx_n_mcell(4)

   type line_parameters
     real(8) :: intercept
     real(8) :: slope
     logical :: vert_line
   end type

   TYPE(line_parameters) :: lin_parms

   TYPE(proj_info) :: proj

!---------------------------------------------------------------------
!	... include files
!---------------------------------------------------------------------
   include 'netcdf.inc'

   CONTAINS

   subroutine proj_init( map_proj, lon1, lat1, truelat1, truelat2, &
                         stdlon, loninc, latinc, knowni, knownj, &
                         dx, ide, jde )
!-------------------------------------------------------------
!  ... intialize wrf map projection
!-------------------------------------------------------------

!-------------------------------------------------------------
!  ... dummy arguments
!-------------------------------------------------------------
   integer, intent(in) :: map_proj
   integer, intent(in) :: ide
   integer, intent(in) :: jde
   real, intent(in)    :: lon1
   real, intent(in)    :: lat1
   real, intent(in)    :: truelat1
   real, intent(in)    :: truelat2
   real, intent(in)    :: stdlon
   real, intent(in)    :: loninc
   real, intent(in)    :: latinc
   real, intent(in)    :: knowni
   real, intent(in)    :: knownj
   real, intent(in)    :: dx

!-------------------------------------------------------------
!  ... local variables
!-------------------------------------------------------------
   integer :: astat
   integer :: grid

   if( map_proj < PROJ_LATLON .or. map_proj > PROJ_MERC ) then
     write(*,'('' proj_init: input projection, '',i2,'', is out of bounds [1-3]'')') map_proj
     stop
   else
     write(*,*) ' '
     write(*,'('' proj_init: projection = '',i2)') map_proj
   endif

   proj%code     = map_proj
   proj%lat1     = lat1
   proj%lon1     = lon1
   proj%ixdim    = ide+1
   proj%jydim    = jde+1
   proj%dx       = dx
   if( map_proj > PROJ_LATLON ) then
     proj%truelat1 = truelat1
     proj%truelat2 = truelat2
     proj%stdlon   = stdlon
     proj%knowni   = proj%ixdim/2.
     proj%knownj   = proj%jydim/2.
     proj%re_m     = EARTH_RADIUS_M
     if( proj%truelat1 < 0. ) then
       proj%hemi = -1.0 
     else
       proj%hemi = 1.0
     endif
     proj%rebydx = proj%re_m / proj%dx
     pole_lat    = proj%hemi*90.
   else
     proj%loninc   = loninc
     proj%latinc   = latinc
     proj%knowni   = knowni
     proj%knownj   = knownj
   endif

   if( proj%code == PROJ_LC ) then
     if( abs(proj%truelat2) > 90. ) then
       proj%truelat2 = proj%truelat1
     end if
     call set_lc( proj )
   elseif( proj%code == PROJ_PS ) then
     call set_ps( proj )
   elseif( proj%code == PROJ_MERC ) then
     call set_merc( proj )
   endif

   proj%init     = .true.

!-------------------------------------------------------------
!  ... a few projection variable diagnostics
!-------------------------------------------------------------
   if( map_proj > PROJ_LATLON ) then
     write(*,*) 'proj_init: proj%hemi    = ',proj%hemi
     write(*,*) 'proj_init: proj%rebydx  = ',proj%rebydx
     write(*,*) 'proj_init: proj%polei,j = ',proj%polei,proj%polej
   endif
   write(*,'('' proj_init: west-east,south-north = '',2i5)') proj%ixdim-1,proj%jydim-1

   if( allocated( lon ) ) then
     deallocate( lon )
   endif
   allocate( lon(ide,jde),stat=astat )
   if( astat /= 0 ) then
     write(*,*) 'proj_init; failed to allocate lon: error = ',astat
     stop
   endif
   if( allocated( lat ) ) then
     deallocate( lat )
   endif
   allocate( lat(ide,jde),stat=astat )
   if( astat /= 0 ) then
     write(*,*) 'proj_init; failed to allocate lat: error = ',astat
     stop
   endif

!-------------------------------------------------------------
!  ... map wrf grid centers to lon,lat
!-------------------------------------------------------------
   do j = 1,proj%jydim-1
     do i = 1,proj%ixdim-1
       call ijll( real(i), real(j), proj, lat(i,j), lon(i,j) )
     end do
   end do

   write(*,*) ' '
   write(*,*) 'wrf domain corners'
   write(*,*) '--- ------ -------'
   write(*,'('' sw corner @ ('',1p,g15.8,'','',g15.8,'')'')') lon(1,1),lat(1,1)
   write(*,'('' se corner @ ('',1p,g15.8,'','',g15.8,'')'')') lon(ide,1),lat(ide,1)
   write(*,'('' ne corner @ ('',1p,g15.8,'','',g15.8,'')'')') lon(ide,jde),lat(ide,jde)
   write(*,'('' nw corner @ ('',1p,g15.8,'','',g15.8,'')'')') lon(1,jde),lat(1,jde)
   write(*,*) ' '

   end subroutine proj_init

   subroutine area_interp_init( grid, model_area_type, diag_level )
!-------------------------------------------------------------
!  ... area conserving interpolation from data to wrf grid
!-------------------------------------------------------------

!-------------------------------------------------------------
!  ... dummy arguments
!-------------------------------------------------------------
   integer, intent(in) :: diag_level
   type(grid_type), intent(inout) :: grid
   type(area_type), intent(inout) :: model_area_type(:,:)

!-------------------------------------------------------------
!  ... local variables
!-------------------------------------------------------------
   integer :: astat
   integer :: ierr
   integer :: dcell_cnt
   integer :: dcell_ndx
   integer :: j
   integer :: nlon, nlat
   integer, allocatable :: dcell_partial_lon_ndx(:)
   integer, allocatable :: dcell_partial_lat_ndx(:)
   real    :: dcell_area
   real    :: area_min, area_max
   real, allocatable    :: dcell_wght(:)
   logical :: dcell_outside_mcell

   nlon = grid%nlons
   nlat = grid%nlats
   m_lat_s = 1
   m_lat_e = proj%jydim-1
   m_lon_s = 1
   m_lon_e = proj%ixdim-1
   mcell_area = 0.
   area_min = 1. ; area_max = 1.
model_lat_loop : &
   do mlat_ndx = m_lat_s,m_lat_e
     model_y(1:2) = real(mlat_ndx) - .5
     model_y(3:4) = real(mlat_ndx) + .5
model_lon_loop : &
     do mlon_ndx = m_lon_s,m_lon_e
       model_x(1:4:3) = real(mlon_ndx) -.5
       model_x(2:3)   = real(mlon_ndx) +.5
!-------------------------------------------------------------
!  map model cell vertices to lat,lon
!-------------------------------------------------------------
       do j = 1,2
         do i = 1,2
           call ijll( real(i+mlon_ndx-1)-.5, real(j+mlat_ndx-1)-.5, proj, model_lat(i,j), model_lon(i,j) )
         end do
       end do

       minmax_lat(1) = minval( model_lat(:,:) )
       minmax_lat(2) = maxval( model_lat(:,:) )
       minmax_lon(1) = minval( model_lon(:,:) )
       minmax_lon(2) = maxval( model_lon(:,:) )
!-------------------------------------------------------------
!  check for model cell out of range of data grid
!-------------------------------------------------------------
       if( minmax_lat(1) >= grid%yedge(nlat+1) .or.  minmax_lat(2) <= grid%yedge(1) ) then
         cycle model_lon_loop
       endif
       
         model_area_type(mlon_ndx,mlat_ndx)%has_data = .true.
!-------------------------------------------------------------
!  check for any model grid cell edge "crossing"
!  from positive to negative longitude
!-------------------------------------------------------------
         if( proj%code == PROJ_PS ) then
           lon_cross(:,1)  = (/ (model_lon(1,1)*model_lon(2,1) < 0. .and. model_lon(1,1) > 0.), &
                                (model_lon(2,1)*model_lon(2,2) < 0. .and. model_lon(2,1) > 0.) /)
           lon_cross(:,2)  = (/ (model_lon(1,2)*model_lon(1,1) < 0. .and. model_lon(1,2) > 0.), &
                                (model_lon(2,2)*model_lon(1,2) < 0. .and. model_lon(2,2) > 0.) /)
           lat_mask(:,1)  = (/ model_lat(1,1) /= pole_lat .and. model_lat(2,1) /= pole_lat, &
                               model_lat(2,1) /= pole_lat .and. model_lat(2,2) /= pole_lat /)
           lat_mask(:,2)  = (/ model_lat(1,2) /= pole_lat .and. model_lat(1,1) /= pole_lat, &
                               model_lat(2,2) /= pole_lat .and. model_lat(1,2) /= pole_lat /)
           lon_cross(:,:) = lon_cross(:,:) .and. lat_mask(:,:)
           has_lon_cross = any( lon_cross(:,:) )
         else
           has_lon_cross = .false.
         endif

lon_cross_a : &
         if( has_lon_cross ) then
cross_loop : &
           do j = 1,2
             do i = 1,2
               if( lon_cross(i,j) ) then
                 exit cross_loop
               endif
             end do
           end do cross_loop

           delta_lat = 0.
           delta_lon = 0.
           if( i == j ) then
             if( j == 1 ) then
               delta_lon = .01
             else
               delta_lon = -.01
             endif
           else
             if( j == 1 ) then
               delta_lat = .01
             else
               delta_lat = -.01
             endif
           endif
           xc = real(i+mlon_ndx-1) - .5 + delta_lon
           yc = real(j+mlat_ndx-1) - .5 + delta_lat
           call ijll_ps( xc, yc, proj, delta_lat, delta_lon )
!-------------------------------------------------------------
!  check for crossing of date line
!-------------------------------------------------------------
           has_lon_cross = delta_lon > model_lon(i,j)
           if( has_lon_cross ) then
             minmax_lon(1) = maxval( model_lon(:,:),mask=model_lon(:,:)<0. .and. model_lat(:,:)/=pole_lat )
             minmax_lon(2) = minval( model_lon(:,:),mask=model_lon(:,:)>0. .and. model_lat(:,:)/=pole_lat )
           endif
         endif lon_cross_a

!-------------------------------------------------------------
!  find lons,lats of data cells enclosing model cell
!-------------------------------------------------------------
         target_lon = minmax_lon(1)
         do lon_s = 1,nlon+1
           if( grid%xedge(lon_s) > target_lon ) then 
             exit
           endif
         enddo
         lon_s = min( max( 1,lon_s-1 ),nlon+1 )
         target_lon = minmax_lon(2)
         do lon_e = 1,nlon+1
           if( grid%xedge(lon_e) >= target_lon ) then 
             exit
           endif
         enddo
         lon_e = min( lon_e,nlon+1)

         target_lat = minmax_lat(1)
         do lat_s = 1,nlat+1
           if( grid%yedge(lat_s) > target_lat ) then 
             exit
           endif
         enddo
         lat_s = min( max( 1,lat_s-1 ),nlat+1 )
         target_lat = minmax_lat(2)
         if( target_lat >= grid%yedge(nlat+1) ) then
           lat_e = nlat
         else
           do lat_e = lat_s,nlat+1
             if( grid%yedge(lat_e) >= target_lat ) then 
               exit
             endif
           end do
         endif
         lat_e = min( lat_e,nlat+1)
!-------------------------------------------------------------
!  check data cell long range for cross over longitude endpoint
!-------------------------------------------------------------
         if( has_lon_cross ) then
           i = lon_s
           lon_s = lon_e
           lon_e = i + nlon
         endif

         model_area_type(mlon_ndx,mlat_ndx)%lon_s = lon_s
         model_area_type(mlon_ndx,mlat_ndx)%lon_e = lon_e
         model_area_type(mlon_ndx,mlat_ndx)%lat_s = lat_s
         model_area_type(mlon_ndx,mlat_ndx)%lat_e = lat_e
         dcell_cnt = (lon_e - lon_s + 1)*(lat_e - lat_s + 1)
         mcell_area = mcell_area + real(dcell_cnt)

         allocate( dcell_partial_lon_ndx(dcell_cnt), dcell_partial_lat_ndx(dcell_cnt), &
                   dcell_wght(dcell_cnt),stat=astat )
         if( astat /= 0 ) then
           write(*,*) 'area_interp; failed to allocate dcell_partial_lon_ndx ... dcell_wght: error = ',astat
           stop
         endif

       dcell_ndx = 0
!-------------------------------------------------------------
!  loop over data cells
!-------------------------------------------------------------
data_lat_loop : &
       do j = lat_s,lat_e
         data_lat(1) = real(grid%yedge(j),kind=4)
         data_lat(2) = data_lat(1)
         data_lat(3) = real(grid%yedge(j+1),kind=4)
         data_lat(4) = data_lat(3)
data_lon_loop : &
         do dlong_ndx = lon_s,lon_e
           i = mod( (dlong_ndx-1),nlon ) + 1
           data_lon(1) = real(grid%xedge(i),kind=4)
           data_lon(4) = data_lon(1)
           data_lon(2) = real(grid%xedge(i+1),kind=4)
           data_lon(3) = data_lon(2)
           call llij( real(grid%yedge(j),kind=4), real(grid%xedge(i),kind=4), proj, x(1), y(1) )
           call llij( real(grid%yedge(j),kind=4), real(grid%xedge(i+1),kind=4), proj, x(2), y(2) )
           call llij( real(grid%yedge(j+1),kind=4), real(grid%xedge(i+1),kind=4), proj, x(3), y(3) )
           call llij( real(grid%yedge(j+1),kind=4), real(grid%xedge(i),kind=4), proj, x(4), y(4) )

           minmax_x(1) = minval( x(:) )
           minmax_x(2) = maxval( x(:) )
           minmax_y(1) = minval( y(:) )
           minmax_y(2) = maxval( y(:) )

           model_area_type(mlon_ndx,mlat_ndx)%total_dcell_cnt = &
               model_area_type(mlon_ndx,mlat_ndx)%total_dcell_cnt + 1
           dcell_outside_mcell = minmax_x(1) >= model_x(2) .or. minmax_x(2) <= model_x(1) .or. &
                                 minmax_y(1) >= model_y(3) .or. minmax_y(2) <= model_y(1)
           if( debugging ) then
             write(*,'(''area_interp: dcell_outside_mcell = '',l)') dcell_outside_mcell
           endif

inside_mcell : &
           if( .not. dcell_outside_mcell ) then
             model_area_type(mlon_ndx,mlat_ndx)%active_dcell_cnt = &
                 model_area_type(mlon_ndx,mlat_ndx)%active_dcell_cnt + 1

             do l = 1,4
               x_n_mcell(l) = model_x(1) <= x(l) .and. x(l) <= model_x(2)
               y_n_mcell(l) = model_y(1) <= y(l) .and. y(l) <= model_y(4)
             end do
             cnt_dvtx_n_mcell = count( x_n_mcell(:) .and. y_n_mcell(:) )

             dcell_ndx = dcell_ndx + 1
             dcell_partial_lon_ndx(dcell_ndx) = dlong_ndx
             dcell_partial_lat_ndx(dcell_ndx) = j
             if( cnt_dvtx_n_mcell == 4 ) then
               dcell_wght(dcell_ndx) = poly_area( 4, x, y )
             else
               dvtx_shadow_ndx(:) = 0
               do l = 1,4
                 if( x_n_mcell(l) .neqv. y_n_mcell(l) ) then
                   dvtx_shadow_ndx(l) = shadow_map( x(l),y(l), n_shadow_zone=.true. )
                 elseif( .not. (x_n_mcell(l) .or. y_n_mcell(l)) ) then
                   dvtx_shadow_ndx(l) = shadow_map( x(l),y(l), n_shadow_zone=.false. )
                 endif
               end do
               dcell_wght(dcell_ndx) = area_map( i,j,mlon_ndx,mlat_ndx )
             endif
             if( debugging ) then
               write(*,'(''area_interp: dcell_area('',i4,'','',i4,'') = '',1pg15.7)') i,j,dcell_wght(dcell_ndx)
               write(*,'(''area_interp: dcell x = '',1p4g15.7)') x(:)
               write(*,'(''area_interp: dcell y = '',1p4g15.7)') y(:)
             endif
           endif inside_mcell
         end do data_lon_loop
       end do data_lat_loop

       if( dcell_ndx > 0 ) then
         model_area_type(mlon_ndx,mlat_ndx)%partial_dcell_cnt = dcell_ndx
         allocate( model_area_type(mlon_ndx,mlat_ndx)%dcell_lon_ndx(dcell_ndx), &
                   model_area_type(mlon_ndx,mlat_ndx)%dcell_lat_ndx(dcell_ndx), &
                   model_area_type(mlon_ndx,mlat_ndx)%wght(dcell_ndx),stat=astat )
         if( astat /= 0 ) then
           write(*,*) 'area_interp; failed to allocate model type dcell_partial_lon_ndx ... dcell_wght: error = ',astat
           stop
         endif
         model_area_type(mlon_ndx,mlat_ndx)%dcell_lon_ndx(:dcell_ndx) = dcell_partial_lon_ndx(:dcell_ndx)
         model_area_type(mlon_ndx,mlat_ndx)%dcell_lat_ndx(:dcell_ndx) = dcell_partial_lat_ndx(:dcell_ndx)
         model_area_type(mlon_ndx,mlat_ndx)%wght(:dcell_ndx)          = dcell_wght(:dcell_ndx)
         area_min = min( area_min, sum(dcell_wght(:dcell_ndx)) )
         area_max = max( area_max, sum(dcell_wght(:dcell_ndx)) )
       endif
       deallocate( dcell_partial_lon_ndx, dcell_partial_lat_ndx, dcell_wght )
     end do model_lon_loop
   end do model_lat_loop

   if( diag_level > 100 ) then
     write(*,*) ' '
     write(*,*) 'area_interp: diagnostics'
     write(*,'(''area_interp: min model lat ndx = '',i6)') minval( model_area_type(:,:)%lat_s, &
                                                                      mask=model_area_type(:,:)%has_data )
     write(*,'(''area_interp: max model lat ndx = '',i6)') maxval( model_area_type(:,:)%lat_s, &
                                                                      mask=model_area_type(:,:)%has_data )
     write(*,'(''area_interp: mcells outside data domain = '',1pg22.15)') count( .not. model_area_type(:,:)%has_data )
     write(*,'(''area_interp: wght size = '',1pg22.15)') mcell_area
     write(*,'(''area_interp: total dcell cnt            = '',i9)') sum( model_area_type(:,:)%total_dcell_cnt )
     write(*,'(''area_interp: active dcell cnt           = '',i9)') sum( model_area_type(:,:)%active_dcell_cnt )
     write(*,'(''area_interp: interior dcell cnt         = '',i9)') sum( model_area_type(:,:)%interior_dcell_cnt )
     write(*,'(''area_interp: partial dcell cnt          = '',i9)') sum( model_area_type(:,:)%partial_dcell_cnt )
   endif

   end subroutine area_interp_init

   subroutine area_interp( data_file, grid, model_area_type, wrk_emis, diag_level )
!-------------------------------------------------------------
!  ... area conserving interpolation from data to wrf grid
!-------------------------------------------------------------

   use anthro_types, only : data_file_type

!-------------------------------------------------------------
!  ... dummy arguments
!-------------------------------------------------------------
   integer, intent(in) :: diag_level
   real, intent(inout) :: wrk_emis(:,:)
   type(grid_type), intent(in)      :: grid
   type(data_file_type), intent(in) :: data_file
   type(area_type), intent(in)      :: model_area_type(:,:)


!-------------------------------------------------------------
!  ... local variables
!-------------------------------------------------------------
   integer :: nlon
   integer :: nlat
   integer :: astat
   integer :: ierr
   integer :: dcell_cnt
   integer :: dcell_ndx
   integer :: j
   integer, allocatable :: dcell_partial_lon_ndx(:)
   integer, allocatable :: dcell_partial_lat_ndx(:)
   real    :: missing_value
   real    :: dcell_area
   real    :: wrk_sum
   real(8) :: lon_edge(grid%nlons+1)
   real(8) :: lat_edge(grid%nlats+1)
   logical :: dcell_outside_mcell

   nlon = grid%nlons
   nlat = grid%nlats
   missing_value = data_file%missing_value
   lon_edge(:) = grid%xedge(:)
   lat_edge(:) = grid%yedge(:)

   write(*,*) 'aera_interp: raw dataset max value   = ',maxval(data_file%src_data(:,:))
   write(*,*) 'aera_interp: raw dataset max indices = ',maxloc(data_file%src_data(:,:))

   wrk_emis(:,:) = 0.
   m_lat_s = 1
   m_lat_e = proj%jydim-1
   m_lon_s = 1
   m_lon_e = proj%ixdim-1
   mcell_area = 0.
model_lat_loop : &
   do mlat_ndx = m_lat_s,m_lat_e
     model_y(1:2) = real(mlat_ndx) - .5
     model_y(3:4) = real(mlat_ndx) + .5
model_lon_loop : &
     do mlon_ndx = m_lon_s,m_lon_e
       model_x(1:4:3) = real(mlon_ndx) -.5
       model_x(2:3)   = real(mlon_ndx) +.5
!-------------------------------------------------------------
!  map model cell vertices to lat,lon
!-------------------------------------------------------------
       do j = 1,2
         do i = 1,2
           call ijll( real(i+mlon_ndx-1)-.5, real(j+mlat_ndx-1)-.5, proj, model_lat(i,j), model_lon(i,j) )
         end do
       end do

       minmax_lat(1) = minval( model_lat(:,:) )
       minmax_lat(2) = maxval( model_lat(:,:) )
       minmax_lon(1) = minval( model_lon(:,:) )
       minmax_lon(2) = maxval( model_lon(:,:) )
!-------------------------------------------------------------
!  check for model cell out of range of data grid
!-------------------------------------------------------------
       if( minmax_lat(1) >= lat_edge(nlat+1) .or.  minmax_lat(2) <= lat_edge(1) ) then
         cycle model_lon_loop
       endif
       
       lon_s = model_area_type(mlon_ndx,mlat_ndx)%lon_s
       lon_e = model_area_type(mlon_ndx,mlat_ndx)%lon_e
       lat_s = model_area_type(mlon_ndx,mlat_ndx)%lat_s
       lat_e = model_area_type(mlon_ndx,mlat_ndx)%lat_e

       dcell_ndx = 0
       wrk_sum   = 0.
!-------------------------------------------------------------
!  loop over data cells
!-------------------------------------------------------------
data_lat_loop : &
       do j = lat_s,lat_e
         data_lat(1) = real(lat_edge(j),kind=4)
         data_lat(2) = data_lat(1)
         data_lat(3) = real(lat_edge(j+1),kind=4)
         data_lat(4) = data_lat(3)
data_lon_loop : &
         do dlong_ndx = lon_s,lon_e
           i = mod( (dlong_ndx-1),nlon ) + 1
           data_lon(1) = real(lon_edge(i),kind=4)
           data_lon(4) = data_lon(1)
           data_lon(2) = real(lon_edge(i+1),kind=4)
           data_lon(3) = data_lon(2)
           call llij( real(lat_edge(j),kind=4), real(lon_edge(i),kind=4), proj, x(1), y(1) )
           call llij( real(lat_edge(j),kind=4), real(lon_edge(i+1),kind=4), proj, x(2), y(2) )
           call llij( real(lat_edge(j+1),kind=4), real(lon_edge(i+1),kind=4), proj, x(3), y(3) )
           call llij( real(lat_edge(j+1),kind=4), real(lon_edge(i),kind=4), proj, x(4), y(4) )

           minmax_x(1) = minval( x(:) )
           minmax_x(2) = maxval( x(:) )
           minmax_y(1) = minval( y(:) )
           minmax_y(2) = maxval( y(:) )

           dcell_outside_mcell = minmax_x(1) >= model_x(2) .or. minmax_x(2) <= model_x(1) .or. &
                                 minmax_y(1) >= model_y(3) .or. minmax_y(2) <= model_y(1)
           if( debugging ) then
             write(*,'(''area_interp: dcell_outside_mcell = '',l)') dcell_outside_mcell
           endif
inside_mcell : &
           if( .not. dcell_outside_mcell ) then
             dcell_ndx  = dcell_ndx + 1
             dcell_area = model_area_type(mlon_ndx,mlat_ndx)%wght(dcell_ndx)
             wrk_sum    = wrk_sum + dcell_area*data_file%src_data(i,j)
             if( debugging ) then
               write(*,'(''area_interp: dcell_area('',i4,'','',i4,'') = '',1pg15.7)') i,j,dcell_area
               write(*,'(''area_interp: dcell x = '',1p4g15.7)') x(:)
               write(*,'(''area_interp: dcell y = '',1p4g15.7)') y(:)
             endif
           endif inside_mcell
         end do data_lon_loop
       end do data_lat_loop

       wrk_emis(mlon_ndx,mlat_ndx) = wrk_sum

     end do model_lon_loop
   end do model_lat_loop

   if( diag_level > 100 ) then
     write(*,*) ' '
     write(*,*)  'raw dataset size = ',size(data_file%src_data)
     write(*,*)  'raw dataset min,max values = ',minval(data_file%src_data(:,:)),maxval(data_file%src_data(:,:))
     write(*,*)  'raw dataset max    indices = ',maxloc(data_file%src_data(:,:))
     write(*,*)  'raw dataset missing value count = ',count(data_file%src_data(:,:) == missing_value )
     write(*,*)  '% valid data = ',100.* real( count(data_file%src_data(:,:) /= missing_value ) ) /real(size(data_file%src_data))
   endif

   if( diag_level > 100 ) then
     write(*,*) ' '
     write(*,*)  'emis dataset size = ',size(wrk_emis)
     write(*,*)  'emis dataset min,max values = ',minval(wrk_emis(:,:)),maxval(wrk_emis(:,:))
     write(*,*)  'emis dataset missing value count = ',count(wrk_emis(:,:) == missing_value )
     write(*,*)  '% valid data = ',100.* real( count(wrk_emis(:,:) /= missing_value ) ) /real(size(wrk_emis))
     write(*,*)  'emis dataset surf integral       = ',sum(wrk_emis(:,:))
     write(*,*)  'emis dataset min value           = ',minval(wrk_emis(:,:))
     write(*,*)  'emis dataset min indices         = ',minloc(wrk_emis(:,:))
     write(*,*)  'emis dataset max value           = ',maxval(wrk_emis(:,:))
     write(*,*)  'emis dataset max indices         = ',maxloc(wrk_emis(:,:))

     write(*,*) ' '
     write(*,*) 'area_interp: diagnostics'
     write(*,'(''area_interp: min model lat ndx = '',i6)') minval( model_area_type(:,:)%lat_s, &
                                                                      mask=model_area_type(:,:)%has_data )
     write(*,'(''area_interp: max model lat ndx = '',i6)') maxval( model_area_type(:,:)%lat_s, &
                                                                      mask=model_area_type(:,:)%has_data )
     write(*,'(''area_interp: mcells outside data domain = '',1pg22.15)') count( .not. model_area_type(:,:)%has_data )
     write(*,'(''area_interp: wght size = '',1pg22.15)') mcell_area
     write(*,'(''area_interp: total dcell cnt            = '',i9)') sum( model_area_type(:,:)%total_dcell_cnt )
     write(*,'(''area_interp: active dcell cnt           = '',i9)') sum( model_area_type(:,:)%active_dcell_cnt )
     write(*,'(''area_interp: interior dcell cnt         = '',i9)') sum( model_area_type(:,:)%interior_dcell_cnt )
     write(*,'(''area_interp: partial dcell cnt          = '',i9)') sum( model_area_type(:,:)%partial_dcell_cnt )
   endif

!  stop 'diagnostics'

   end subroutine area_interp

   subroutine llij( lat, lon, proj, i, j)
!-----------------------------------------------------------
! Subroutine to convert geographical latitude,longitude values to cartesian x/y
!-----------------------------------------------------------

!-----------------------------------------------------------
! dummy arguments
!-----------------------------------------------------------
   real, intent(in)              :: lat      ! Latitude (-90->90 deg N)
   real, intent(in)              :: lon      ! Longitude (-180->180 E)
   real, intent(out)             :: i        ! Cartesian X coordinate
   real, intent(out)             :: j        ! Cartesian Y coordinate
   type(proj_info)   :: proj

   if( proj%code == PROJ_LATLON ) then
     call llij_latlon( lat, lon, proj, i, j)
   elseif( proj%code == PROJ_LC ) then
     call llij_lc( lat, lon, proj, i, j)
   elseif( proj%code == PROJ_PS ) then
     call llij_ps( lat, lon, proj, i, j)
   elseif( proj%code == PROJ_MERC ) then
     call llij_merc( lat, lon, proj, i, j)
   endif

   end subroutine llij

   subroutine ijll( i, j, proj, lat, lon )
!-----------------------------------------------------------
! Subroutine to convert cartesian x/y values to longitude,latitude values
!-----------------------------------------------------------

!-----------------------------------------------------------
! dummy arguments
!-----------------------------------------------------------
   real, intent(in)              :: i        ! Cartesian X coordinate
   real, intent(in)              :: j        ! Cartesian Y coordinate
   real, intent(out)             :: lat      ! Latitude (-90->90 deg N)
   real, intent(out)             :: lon      ! Longitude (-180->180 E)
   type(proj_info)   :: proj

   if( proj%code == PROJ_LATLON ) then
     call ijll_latlon( i, j, proj, lat, lon )
   elseif( proj%code == PROJ_LC ) then
     call ijll_lc( i, j, proj, lat, lon )
   elseif( proj%code == PROJ_PS ) then
     call ijll_ps( i, j, proj, lat, lon )
   elseif( proj%code == PROJ_MERC ) then
     call ijll_merc( i, j, proj, lat, lon )
   endif

   end subroutine ijll

   SUBROUTINE set_lc(proj)
      ! Initialize the remaining items in the proj structure for a
      ! lambert conformal grid.

      IMPLICIT NONE
      
      TYPE(proj_info), INTENT(INOUT)     :: proj
  
      REAL                               :: arg
      REAL                               :: deltalon1
      REAL                               :: tl1r
      REAL                               :: ctl1r
  
      ! Compute cone factor
      CALL lc_cone( proj%truelat1, proj%truelat2, proj%cone )
  
      ! Compute longitude differences and ensure we stay out of the
      ! forbidden "cut zone"
      deltalon1 = proj%lon1 - proj%stdlon
      IF (deltalon1 > 180.) then
        deltalon1 = deltalon1 - 360.
      elseIF (deltalon1 < -180.) then
        deltalon1 = deltalon1 + 360.
      endif
  
      ! Convert truelat1 to radian and compute COS for later use
      tl1r = proj%truelat1 * rad_per_deg
      ctl1r = COS(tl1r)
  
      ! Compute the radius to our known lower-left (SW) corner
      proj%rsw = proj%rebydx * ctl1r/proj%cone * &
             (TAN((90.*proj%hemi-proj%lat1)*rad_per_deg/2.) / &
              TAN((90.*proj%hemi-proj%truelat1)*rad_per_deg/2.))**proj%cone
  
      ! Find pole point
      arg = proj%cone*(deltalon1*rad_per_deg)
      proj%polei = proj%hemi*proj%knowni - proj%hemi * proj%rsw * SIN(arg)
      proj%polej = proj%hemi*proj%knownj + proj%rsw * COS(arg)  
  
   END SUBROUTINE set_lc                             

   SUBROUTINE lc_cone(truelat1, truelat2, cone)
 
   ! Subroutine to compute the cone factor of a Lambert Conformal projection

      IMPLICIT NONE
      
      ! Input Args
      REAL, INTENT(IN)             :: truelat1  ! (-90 -> 90 degrees N)
      REAL, INTENT(IN)             :: truelat2  !   "   "  "   "     "
  
      ! Output Args
      REAL, INTENT(OUT)            :: cone
  
      ! Locals
  
      ! BEGIN CODE
  
      ! First, see if this is a secant or tangent projection.  For tangent
      ! projections, truelat1 = truelat2 and the cone is tangent to the 
      ! Earth's surface at this latitude.  For secant projections, the cone
      ! intersects the Earth's surface at each of the distinctly different
      ! latitudes
      IF (ABS(truelat1-truelat2) > 0.1) THEN
         cone = ALOG10(COS(truelat1*rad_per_deg)) - &
                ALOG10(COS(truelat2*rad_per_deg))
         cone = cone /(ALOG10(TAN((45.0 - ABS(truelat1)/2.0) * rad_per_deg)) - &
                ALOG10(TAN((45.0 - ABS(truelat2)/2.0) * rad_per_deg)))        
      ELSE
         cone = SIN(ABS(truelat1)*rad_per_deg )  
      ENDIF

      RETURN

   END SUBROUTINE lc_cone

   SUBROUTINE ijll_lc( i, j, proj, lat, lon)
 
   ! Subroutine to convert from the (i,j) cartesian coordinate to the 
   ! geographical latitude and longitude for a Lambert Conformal projection.
 
   ! History:
   ! 25 Jul 01: Corrected by B. Shaw, NOAA/FSL
   ! 
      IMPLICIT NONE
  
      ! Input Args
      REAL, INTENT(IN)              :: i        ! Cartesian X coordinate
      REAL, INTENT(IN)              :: j        ! Cartesian Y coordinate
      TYPE(proj_info),INTENT(IN)    :: proj     ! Projection info structure
  
      ! Output Args                 
      REAL, INTENT(OUT)             :: lat      ! Latitude (-90->90 deg N)
      REAL, INTENT(OUT)             :: lon      ! Longitude (-180->180 E)
  
      ! Locals 
      REAL                          :: inew
      REAL                          :: jnew
      REAL                          :: r
      REAL                          :: chi,chi1,chi2
      REAL                          :: r2
      REAL                          :: xx
      REAL                          :: yy
  
      ! BEGIN CODE
  
      chi1 = (90. - proj%hemi*proj%truelat1)*rad_per_deg
      chi2 = (90. - proj%hemi*proj%truelat2)*rad_per_deg
  
      ! See if we are in the southern hemispere and flip the indices
      ! if we are. 
      inew = proj%hemi * i
      jnew = proj%hemi * j
  
      ! Compute radius**2 to i/j location
      xx = inew - proj%polei
      yy = proj%polej - jnew
      r2 = (xx*xx + yy*yy)
      r = SQRT(r2)/proj%rebydx
     
      ! Convert to lat/lon
      IF (r2 == 0.) THEN
         lat = proj%hemi * 90.
         lon = proj%stdlon
      ELSE
         
         ! Longitude
         lon = proj%stdlon + deg_per_rad * ATAN2(proj%hemi*xx,yy)/proj%cone
         lon = AMOD(lon+360., 360.)
   
         ! Latitude.  Latitude determined by solving an equation adapted 
         ! from:
         !  Maling, D.H., 1973: Coordinate Systems and Map Projections
         ! Equations #20 in Appendix I.  
           
         IF (chi1 .EQ. chi2) THEN
            chi = 2.0*ATAN( ( r/TAN(chi1) )**(1./proj%cone) * TAN(chi1*0.5) )
         ELSE
            chi = 2.0*ATAN( (r*proj%cone/SIN(chi1))**(1./proj%cone) * TAN(chi1*0.5)) 
         ENDIF
         lat = (90.0-chi*deg_per_rad)*proj%hemi
  
      ENDIF
  
      IF (lon .GT. +180.) lon = lon - 360.
      IF (lon .LT. -180.) lon = lon + 360.
 
   END SUBROUTINE ijll_lc

   SUBROUTINE llij_lc( lat, lon, proj, i, j)
 
   ! Subroutine to compute the geographical latitude and longitude values
   ! to the cartesian x/y on a Lambert Conformal projection.
     
      IMPLICIT NONE
  
      ! Input Args
      REAL, INTENT(IN)              :: lat      ! Latitude (-90->90 deg N)
      REAL, INTENT(IN)              :: lon      ! Longitude (-180->180 E)
      TYPE(proj_info),INTENT(IN)      :: proj     ! Projection info structure
  
      ! Output Args                 
      REAL, INTENT(OUT)             :: i        ! Cartesian X coordinate
      REAL, INTENT(OUT)             :: j        ! Cartesian Y coordinate
  
      ! Locals 
      REAL                          :: arg
      REAL                          :: deltalon
      REAL                          :: tl1r
      REAL                          :: rm
      REAL                          :: ctl1r
      
  
      ! BEGIN CODE
      
      ! Compute deltalon between known longitude and standard lon and ensure
      ! it is not in the cut zone
      deltalon = lon - proj%stdlon
      IF (deltalon .GT. +180.) deltalon = deltalon - 360.
      IF (deltalon .LT. -180.) deltalon = deltalon + 360.
      
      ! Convert truelat1 to radian and compute COS for later use
      tl1r = proj%truelat1 * rad_per_deg
      ctl1r = COS(tl1r)     
     
      ! Radius to desired point
      rm = proj%rebydx * ctl1r/proj%cone * &
           (TAN((90.*proj%hemi-lat)*rad_per_deg/2.) / &
            TAN((90.*proj%hemi-proj%truelat1)*rad_per_deg/2.))**proj%cone
  
      arg = proj%cone*(deltalon*rad_per_deg)
      i = proj%polei + proj%hemi * rm * SIN(arg)
      j = proj%polej - rm * COS(arg)
  
      ! Finally, if we are in the southern hemisphere, flip the i/j
      ! values to a coordinate system where (1,1) is the SW corner
      ! (what we assume) which is different than the original NCEP
      ! algorithms which used the NE corner as the origin in the 
      ! southern hemisphere (left-hand vs. right-hand coordinate?)
      i = proj%hemi * i  
      j = proj%hemi * j

   END SUBROUTINE llij_lc

   SUBROUTINE set_merc(proj)
   
      ! Sets up the remaining basic elements for the mercator projection
  
      IMPLICIT NONE
      TYPE(proj_info), INTENT(INOUT)       :: proj
      REAL                                 :: clain
  
  
      !  Preliminary variables
  
      clain = COS(rad_per_deg*proj%truelat1)
      proj%dlon = proj%dx / (proj%re_m * clain)
  
      ! Compute distance from equator to origin, and store in the 
      ! proj%rsw tag.
  
      proj%rsw = 0.
      IF (proj%lat1 .NE. 0.) THEN
         proj%rsw = (ALOG(TAN(0.5*((proj%lat1+90.)*rad_per_deg))))/proj%dlon
      ENDIF

   END SUBROUTINE set_merc

   SUBROUTINE llij_merc(lat, lon, proj, i, j)
 
      ! Compute i/j coordinate from lat lon for mercator projection
    
      IMPLICIT NONE
      REAL, INTENT(IN)              :: lat
      REAL, INTENT(IN)              :: lon
      TYPE(proj_info),INTENT(IN)    :: proj
      REAL,INTENT(OUT)              :: i
      REAL,INTENT(OUT)              :: j
      REAL                          :: deltalon
  
      deltalon = lon - proj%lon1
      IF (deltalon > -180.) deltalon = deltalon + 360.
      IF (deltalon > 180.) deltalon = deltalon - 360.
      i = proj%knowni + (deltalon/(proj%dlon*deg_per_rad))
      j = proj%knownj + (ALOG(TAN(0.5*((lat + 90.) * rad_per_deg)))) / &
             proj%dlon - proj%rsw
  
   END SUBROUTINE llij_merc

   SUBROUTINE ijll_merc(i, j, proj, lat, lon)
 
      ! Compute the lat/lon from i/j for mercator projection
  
      IMPLICIT NONE
      REAL,INTENT(IN)               :: i
      REAL,INTENT(IN)               :: j    
      TYPE(proj_info),INTENT(IN)    :: proj
      REAL, INTENT(OUT)             :: lat
      REAL, INTENT(OUT)             :: lon 
  
  
      lat = 2.0*ATAN(EXP(proj%dlon*(proj%rsw + j-proj%knownj)))*deg_per_rad - 90.
      lon = (i-proj%knowni)*proj%dlon*deg_per_rad + proj%lon1
      IF (lon > 180.) lon = lon - 360.
      IF (lon < -180.) lon = lon + 360.

   END SUBROUTINE ijll_merc

   SUBROUTINE set_ps(proj)
      ! Initializes a polar-stereographic map projection from the partially
      ! filled proj structure. This routine computes the radius to the
      ! southwest corner and computes the i/j location of the pole for use
      ! in llij_ps and ijll_ps.
      IMPLICIT NONE
   
      ! Declare args
      TYPE(proj_info), INTENT(INOUT)    :: proj
  
      ! Local vars
      REAL                              :: ala1
      REAL                              :: alo1
      REAL                              :: reflon
      REAL                              :: scale_top
  
      ! Executable code
      reflon = proj%stdlon + 90.
  
      ! Compute numerator term of map scale factor
      scale_top = 1. + proj%hemi * SIN(proj%truelat1 * rad_per_deg)
  
      ! Compute radius to lower-left (SW) corner
      ala1 = proj%lat1 * rad_per_deg
      proj%rsw = proj%rebydx*COS(ala1)*scale_top/(1.+proj%hemi*SIN(ala1))
  
      ! Find the pole point
      alo1 = (proj%lon1 - reflon) * rad_per_deg
      proj%polei = proj%knowni - proj%rsw * COS(alo1)
      proj%polej = proj%knownj - proj%hemi * proj%rsw * SIN(alo1)

   END SUBROUTINE set_ps

   SUBROUTINE llij_ps(lat,lon,proj,i,j,debug)
      ! Given latitude (-90 to 90), longitude (-180 to 180), and the
      ! standard polar-stereographic projection information via the 
      ! public proj structure, this routine returns the i/j indices which
      ! if within the domain range from 1->nx and 1->ny, respectively.
  
      IMPLICIT NONE
  
      ! Delcare input arguments
      REAL, INTENT(IN)               :: lat
      REAL, INTENT(IN)               :: lon
      logical, optional, INTENT(IN)  :: debug
      TYPE(proj_info),INTENT(IN)     :: proj
  
      ! Declare output arguments     
      REAL, INTENT(OUT)              :: i !(x-index)
      REAL, INTENT(OUT)              :: j !(y-index)
  
      ! Declare local variables
      
      REAL                           :: reflon
      REAL                           :: scale_top
      REAL                           :: ala
      REAL                           :: alo
      REAL                           :: rm
  
      ! BEGIN CODE
    
      reflon = proj%stdlon + 90.
     
      ! Compute numerator term of map scale factor
  
      scale_top = 1. + proj%hemi * SIN(proj%truelat1 * rad_per_deg)
  
      ! Find radius to desired point
      ala = lat * rad_per_deg
      rm = proj%rebydx * COS(ala) * scale_top/(1. + proj%hemi *SIN(ala))
      alo = (lon - reflon) * rad_per_deg
      i = proj%polei + rm * COS(alo)
      j = proj%polej + proj%hemi * rm * SIN(alo)
      if( present( debug ) ) then
        if( debug ) then
          write(*,*) 'llij_ps: lat,lon   = ',lat,lon
          write(*,*) 'llij_ps: ala, alo  = ',ala,alo
          write(*,*) 'llij_ps: scale_top,hemi = ',scale_top,proj%hemi
          write(*,*) 'llij_ps: reflon,rm = ',reflon,rm
        endif
      endif
   
   END SUBROUTINE llij_ps

   SUBROUTINE ijll_ps(i, j, proj, lat, lon)
 
      ! This is the inverse subroutine of llij_ps.  It returns the 
      ! latitude and longitude of an i/j point given the projection info 
      ! structure.  
  
      IMPLICIT NONE
  
      ! Declare input arguments
      REAL, INTENT(IN)                    :: i    ! Column
      REAL, INTENT(IN)                    :: j    ! Row
      TYPE (proj_info), INTENT(IN)        :: proj
      
      ! Declare output arguments
      REAL, INTENT(OUT)                   :: lat     ! -90 -> 90 north
      REAL, INTENT(OUT)                   :: lon     ! -180 -> 180 East
  
      ! Local variables
      REAL                                :: reflon
      REAL                                :: scale_top
      REAL                                :: xx,yy
      REAL                                :: gi2, r2
      REAL                                :: arccos
  
      ! Begin Code
  
      ! Compute the reference longitude by rotating 90 degrees to the east
      ! to find the longitude line parallel to the positive x-axis.
      reflon = proj%stdlon + 90.
     
      ! Compute numerator term of map scale factor
      scale_top = 1. + proj%hemi * SIN(proj%truelat1 * rad_per_deg)
  
      ! Compute radius to point of interest
      xx = i - proj%polei
      yy = (j - proj%polej) * proj%hemi
      r2 = xx**2 + yy**2
  
      ! Now the magic code
      IF (r2 .EQ. 0.) THEN 
         lat = proj%hemi * 90.
         lon = reflon
      ELSE
         gi2 = (proj%rebydx * scale_top)**2.
         lat = deg_per_rad * proj%hemi * ASIN((gi2-r2)/(gi2+r2))
         arccos = ACOS(xx/SQRT(r2))
         IF (yy .GT. 0) THEN
            lon = reflon + deg_per_rad * arccos
         ELSE
            lon = reflon - deg_per_rad * arccos
         ENDIF
      ENDIF
    
      ! Convert to a -180 -> 180 East convention
      IF (lon > 180.) then
        lon = lon - 360.
      ELSEIF (lon < -180.) then
        lon = lon + 360.
      ENDIF

   END SUBROUTINE ijll_ps

   SUBROUTINE llij_latlon( lat, lon, proj, i, j )
  
      ! Compute the i/j location of a lat/lon on a LATLON grid.
      IMPLICIT NONE
      REAL, INTENT(IN)             :: lat
      REAL, INTENT(IN)             :: lon
      TYPE(proj_info), INTENT(IN)  :: proj
      REAL, INTENT(OUT)            :: i
      REAL, INTENT(OUT)            :: j
  
      REAL                         :: deltalat
      REAL                         :: deltalon
  
      ! Compute deltalat and deltalon as the difference between the input 
      ! lat/lon and the origin lat/lon
      deltalat = lat - proj%lat1
      deltalon = lon - proj%lon1      
      
      ! Compute i/j
      i = deltalon/proj%loninc
      j = deltalat/proj%latinc

      i = i + proj%knowni
      j = j + proj%knownj
  
   END SUBROUTINE llij_latlon

   SUBROUTINE ijll_latlon( i, j, proj, lat, lon )
  
      ! Compute the lat/lon location of an i/j on a LATLON grid.
      IMPLICIT NONE
      REAL, INTENT(IN)             :: i
      REAL, INTENT(IN)             :: j
      TYPE(proj_info), INTENT(IN)  :: proj
      REAL, INTENT(OUT)            :: lat
      REAL, INTENT(OUT)            :: lon
  
      REAL                         :: i_work, j_work
      REAL                         :: deltalat
      REAL                         :: deltalon
  
      i_work = i - proj%knowni
      j_work = j - proj%knownj

      ! Compute deltalat and deltalon 
      deltalat = j_work*proj%latinc
      deltalon = i_work*proj%loninc
  
      lat = proj%lat1 + deltalat
      lon = proj%lon1 + deltalon
  
   END SUBROUTINE ijll_latlon

   integer FUNCTION shadow_map( x, y, n_shadow_zone )
!---------------------------------------------------------------
!  map data vertex to shadow "zone" index
!---------------------------------------------------------------

!---------------------------------------------------------------
!  dummy arguments
!---------------------------------------------------------------
   real, intent(in) :: x
   real, intent(in) :: y
   logical, intent(in) :: n_shadow_zone  ! .t. => + shd_ndx, .f. => - shd_ndx

   if( n_shadow_zone ) then
     if( x < model_x(1) ) then
       shadow_map = 4
     elseif( x > model_x(2) ) then
       shadow_map = 2
     elseif( y < model_y(1) ) then
       shadow_map = 1
     elseif( y > model_y(4) ) then
       shadow_map = 3
     endif
   else
     if( x < model_x(1) ) then
       if( y < model_y(1) ) then
         shadow_map = 1
       else
         shadow_map = 4
       endif
     elseif( x > model_x(2) ) then
       if( y < model_y(1) ) then
         shadow_map = 2
       else
         shadow_map = 3
       endif
     endif
     shadow_map = -shadow_map
   endif

   end FUNCTION shadow_map

   real FUNCTION area_map( i,j,mlon_ndx,mlat_ndx )
!---------------------------------------------------------------
!  calculate the intersection area of data cell and model grid cell
!---------------------------------------------------------------

   integer, intent(in) :: i,j,mlon_ndx,mlat_ndx
!---------------------------------------------------------------
!  local variables
!---------------------------------------------------------------
   integer :: k, l, lp1, linc
   integer :: m, mp1, n
   integer :: ml, mu, minc
   integer :: nv
   integer :: shd_ndx
   integer :: vtx_cnt
   integer :: ndx(1)
   integer :: wrk_cnt(4)
   integer :: wrk1_cnt(4)
   real    :: dcell_area
   real    :: mvtx_lon(4)
   real    :: mvtx_lat(4)
   real    :: vtx_x(8)
   real    :: vtx_y(8)
   logical :: match_pnt
   logical :: mvtx_n_dcell(4)

!---------------------------------------------------------------
!  data vertex in model cell
!---------------------------------------------------------------
   dvtx_n_mcell(:) = x_n_mcell(:) .and. y_n_mcell(:)
!---------------------------------------------------------------
!  determine if model cell vertices are in data cell
!---------------------------------------------------------------
   mvtx_lon(:) = (/ model_lon(1,1), model_lon(2,1), model_lon(2,2), model_lon(1,2) /)
   mvtx_lat(:) = (/ model_lat(1,1), model_lat(2,1), model_lat(2,2), model_lat(1,2) /)
   do l = 1,4
     mvtx_n_dcell(l) = pnt_in_quad( (/model_x(l), model_y(l) /), x, y )
   end do

   if( all( mvtx_n_dcell(:) ) ) then
     area_map = poly_area( 4, model_x, model_y )
     return
   endif

   dcell_area = poly_area( 4, x, y )
valid_data_cell : &
   if( dcell_area > 0. ) then
     nv = 0
vertex_loop : &
     do l = 1,4
       lp1 = mod( l,4 ) + 1
       if( l /= 4 ) then
         linc = 1
       else
         linc = -3
       endif
       vtx_cnt = count( x_n_mcell(l:lp1:linc) ) + count( y_n_mcell(l:lp1:linc) )
       wrk_cnt(l) = vtx_cnt
       if( maxval( abs(dvtx_shadow_ndx(l:lp1:linc)) ) == 0 ) then
!---------------------------------------------------------------
!  both dline endpoints are in mdl cell
!---------------------------------------------------------------
         nv = nv + 1
         vtx_x(nv) = x(l)
         vtx_y(nv) = y(l)
         wrk1_cnt(l) = nv
       elseif( minval( abs(dvtx_shadow_ndx(l:lp1:linc)) ) == 0 ) then
!---------------------------------------------------------------
!  one dline endpoint is in mdl cell
!  exactly one xsect with mdl cell boundary
!---------------------------------------------------------------
         CALL set_vertices( l, lp1, nv, vtx_x, vtx_y )
       elseif( abs(dvtx_shadow_ndx(l)) /= abs(dvtx_shadow_ndx(lp1)) ) then
!---------------------------------------------------------------
!  neither dline endpoint is in mdl cell
!  either 0,1, or 2 xsects with mdl cell boundary
!---------------------------------------------------------------
         if( minval( dvtx_shadow_ndx(l:lp1:linc) ) > 0 ) then
!---------------------------------------------------------------
!  both dcell edge endpnts in + shd zone
!---------------------------------------------------------------
           CALL set_vertices( l, lp1, nv, vtx_x, vtx_y )
         elseif( maxval( dvtx_shadow_ndx(l:lp1:linc) ) < 0 ) then
           if( mod( sum( abs(dvtx_shadow_ndx(l:lp1:linc)) ),2 ) == 0 ) then
!---------------------------------------------------------------
!  both dcell edge endpnts in - shd zone
!---------------------------------------------------------------
             CALL set_vertices( l, lp1, nv, vtx_x, vtx_y )
           endif
         else
!---------------------------------------------------------------
!  dcell edge endpnts in +,- shd zones
!---------------------------------------------------------------
           if( dvtx_shadow_ndx(l) > 0 ) then
             do n = 2,3
               shd_ndx = dvtx_shadow_ndx(l) + n
               if( shd_ndx > 4 ) shd_ndx = mod( shd_ndx,4 )
               if( abs(dvtx_shadow_ndx(lp1)) == shd_ndx ) then
                 CALL set_vertices( l, lp1, nv, vtx_x, vtx_y )
                 exit
               endif
             end do
           else
             do n = 1,2
               shd_ndx = abs(dvtx_shadow_ndx(l)) + n
               if( shd_ndx > 4 ) shd_ndx = mod( shd_ndx,4 )
               if( dvtx_shadow_ndx(lp1) == shd_ndx ) then
                 CALL set_vertices( l, lp1, nv, vtx_x, vtx_y )
                 exit
               endif
             end do
           endif
         endif
       endif
     wrk1_cnt(l) = nv
   end do vertex_loop

   vtx_cnt = count( mvtx_n_dcell(:) )
!---------------------------------------------------------------
!  check for model cell vertices in data cell
!---------------------------------------------------------------
   if( vtx_cnt > 0 ) then
     ml = 1; mu = 4; minc = 1  
     if( nv == 2 .and. (mvtx_n_dcell(1) .and. mvtx_n_dcell(4)) ) then
       ml = 4; mu = 1; minc = -1  
     endif
     do k = ml,mu,minc
       if( mvtx_n_dcell(k) ) then
         do l = 1,nv
           match_pnt = vtx_x(l) == model_x(k) .and. vtx_y(l) == model_y(k)
           if( match_pnt ) then
             exit
           endif
         end do
         if( .not. match_pnt ) then
           CALL insert_mvtx( k, nv, vtx_x, vtx_y )
         endif
       endif
     end do
   endif

   if( nv > 2 ) then
     area_map = poly_area( nv, vtx_x, vtx_y )
     if( area_map < 0. ) then
       if( abs( area_map ) > .01 * abs( dcell_area ) ) then
         write(*,*) ' '
         write(*,'(''area_map: xsecting data  cell ('',i5,'','',i5,'') area = '',1pg15.8)') i,j,area_map
         write(*,'(''area_map: data  cell area = '',1pg15.8)') dcell_area
         write(*,'(''area_map: model cell ('',i5,'','',i5,'')'')') mlon_ndx,mlat_ndx
         write(*,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         write(*,*) 'area_map: model cell vertices(x,y)'
         do m = 1,4
           write(*,'('' ('',g15.8,'','',g15.8,'')'')') model_x(m),model_y(m)
         end do
         write(*,*) ' '
         write(*,*) 'area_map: model cell vertices(lon,lat)'
         do m = 1,4
           write(*,'('' ('',g15.8,'','',g15.8,'')'')') mvtx_lon(m),mvtx_lat(m)
         end do
         write(*,*) ' '
         write(*,*) 'area_map: data cell vertices(x,y)'
         do m = 1,4
           write(*,'('' ('',g15.8,'','',g15.8,'')'')') x(m),y(m)
         end do
         write(*,*) ' '
         write(*,*) 'area_map: data cell vertices(lon,lat)'
         do m = 1,4
           write(*,'('' ('',g15.8,'','',g15.8,'')'')') data_lon(m),data_lat(m)
         end do
         write(*,*) ' '
         write(*,'(''area_map: dvtx_n_mcell = ('',3(l1,'',''),l1,'')'')') dvtx_n_mcell(:)
         write(*,'(''area_map: mvtx_n_dcell = ('',3(l1,'',''),l1,'')'')') mvtx_n_dcell(:)
         write(*,'(''area_map: shadow_cnt   = ('',3(i1,'',''),i1,'')'')') wrk_cnt(:)
         write(*,'(''area_map: shadow_ndx   = ('',3(i2,'',''),i2,'')'')') dvtx_shadow_ndx(:)
         write(*,'(''area_map: vtx cnt      = ('',3(i1,'',''),i1,'')'')') wrk1_cnt(:)
         write(*,*) ' '
         do m = 1,nv
           write(*,'('' ('',g15.8,'','',g15.8,'')'')') vtx_x(m),vtx_y(m)
         end do
         write(*,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
       endif
     endif
     area_map = max( area_map, 0. )
   else
     area_map = 0.
   endif
   else valid_data_cell
     area_map = 0.
   endif valid_data_cell

!  END FUNCTION area_map

   CONTAINS

   SUBROUTINE set_vertices( v, vp1, nv, vtx_x, vtx_y )
!---------------------------------------------------------------
!  calculate intersection of line from data vertices
!  (x(l),y(l)) -> (x(lp1),y(lp1)) and the appropriate
!  model cell edge
!---------------------------------------------------------------

!---------------------------------------------------------------
!  dummy arguments
!---------------------------------------------------------------
  integer, intent(in)    :: v
  integer, intent(in)    :: vp1
  integer, intent(inout) :: nv
  real, intent(inout)    :: vtx_x(:)
  real, intent(inout)    :: vtx_y(:)

!---------------------------------------------------------------
!  local variables
!---------------------------------------------------------------
   integer :: n1, n2
   integer :: ndx
   integer :: ndxs, ndxe
   integer :: vndx
   integer :: shadow_ndx
   integer :: shd_ndx_sum
   real(8) :: slope
   real(8) :: intercept
   real    :: xs_pnt(2)
   logical :: found
   logical :: has_xsect

!---------------------------------------------------------------
!  linear eqn between data vertices
!---------------------------------------------------------------
   CALL lin_eqn( x(v), x(vp1), y(v), y(vp1), slope, intercept )

vtx_in_model_cell : &
   if( dvtx_n_mcell(v) .or. dvtx_n_mcell(vp1) ) then
!---------------------------------------------------------------
!  one dvtx in mcell the other not
!  set ndx to index of pnt out of mcell
!---------------------------------------------------------------
     if( dvtx_n_mcell(v) ) then
       ndx = vp1
     else
       ndx = v
     endif
     shadow_ndx = dvtx_shadow_ndx(ndx)
     ndx = abs( shadow_ndx )
in_shadow_zone : &
     if( shadow_ndx > 0 ) then
!---------------------------------------------------------------
!  only one xsect possible
!---------------------------------------------------------------
       call set_xsect( ndx, xs_pnt, slope, intercept, has_xsect )
       if( .not. has_xsect ) then
         write(*,'(''set_vertices: Something very wrong'')')
         write(*,'(''set_vertices: i,j,mlon,mlat_ndx = '',4i6)') i,j,mlon_ndx,mlat_ndx
       endif
       if( dvtx_n_mcell(v) ) then
         nv = nv + 1
         vtx_x(nv) = x(v)
         vtx_y(nv) = y(v)
!        if( xs_pnt(1) /= x(v) .or. xs_pnt(2) /= y(v) ) then
         if( .not. pnts_match( xs_pnt(1), xs_pnt(2), x(v), y(v) ) ) then
           nv = nv + 1
           vtx_x(nv) = xs_pnt(1)
           vtx_y(nv) = xs_pnt(2)
         endif
       elseif( dvtx_n_mcell(vp1) ) then
!        if( xs_pnt(1) /= x(vp1) .or. xs_pnt(2) /= y(vp1) ) then
!        if( .not. pnt_on_bndy( x(vp1), y(vp1), model_x, model_y ) ) then
         if( vp1 /= 4 .or. .not. pnts_match( xs_pnt(1), xs_pnt(2), x(vp1), y(vp1) ) ) then
           nv = nv + 1
           vtx_x(nv) = xs_pnt(1)
           vtx_y(nv) = xs_pnt(2)
         endif
       endif
     else in_shadow_zone
!---------------------------------------------------------------
!  two possible xsects, only one valid
!---------------------------------------------------------------
         n1 = ndx - 1
         if( n1 < 1 ) n1 = n1 + 4
         call set_xsect( n1, xs_pnt, slope, intercept, has_xsect )
         if( .not. has_xsect ) then
           n1 = ndx
           call set_xsect( n1, xs_pnt, slope, intercept, has_xsect )
         endif
         if( dvtx_n_mcell(v) ) then
           nv = nv + 1
           vtx_x(nv) = x(v)
           vtx_y(nv) = y(v)
           if( .not. pnts_match( xs_pnt(1), xs_pnt(2), x(v), y(v) ) ) then
             nv = nv + 1
             vtx_x(nv) = xs_pnt(1)
             vtx_y(nv) = xs_pnt(2)
           endif
         elseif( dvtx_n_mcell(vp1) ) then
           if( vp1 /= 4 .or. .not. pnts_match( xs_pnt(1), xs_pnt(2), x(vp1), y(vp1) ) ) then
             nv = nv + 1
             vtx_x(nv) = xs_pnt(1)
             vtx_y(nv) = xs_pnt(2)
           endif
         endif
     endif in_shadow_zone
   else vtx_in_model_cell
!---------------------------------------------------------------
!  both dvtx in shadow zone
!  either 0, 1, or 2 itersections
!---------------------------------------------------------------
     ndxs = dvtx_shadow_ndx(v)
     ndxe = dvtx_shadow_ndx(vp1)

     if( abs(ndxs) /= abs(ndxe) ) then
       shd_ndx_sum = abs(ndxs) + abs(ndxe)
       if( min( ndxs,ndxe) > 0 ) then
!---------------------------------------------------------------
!  both shd zones > 0
!---------------------------------------------------------------
         if( mod( shd_ndx_sum,2 ) == 0 ) then
!---------------------------------------------------------------
!  exactly two xsects
!---------------------------------------------------------------
           call set_xsect( ndxs, xs_pnt, slope, intercept, has_xsect )
           nv = nv + 1
           vtx_x(nv) = xs_pnt(1)
           vtx_y(nv) = xs_pnt(2)
           call set_xsect( ndxe, xs_pnt, slope, intercept, has_xsect )
           nv = nv + 1
           vtx_x(nv) = xs_pnt(1)
           vtx_y(nv) = xs_pnt(2)
         else
!---------------------------------------------------------------
!  either 0, 1, or 2 xsects
!---------------------------------------------------------------
           call set_xsect( ndxs, xs_pnt, slope, intercept, has_xsect )
           if( has_xsect ) then
             nv = nv + 1
             vtx_x(nv) = xs_pnt(1)
             vtx_y(nv) = xs_pnt(2)
             call set_xsect( ndxe, xs_pnt, slope, intercept, has_xsect )
             if( has_xsect ) then
               nv = nv + 1
               vtx_x(nv) = xs_pnt(1)
               vtx_y(nv) = xs_pnt(2)
             endif
           endif
         endif
       elseif( max( ndxs,ndxe ) < 0 .and.  mod( shd_ndx_sum,2 ) == 0 ) then
         found = .false.
         do n2 = abs(ndxs)-1,abs(ndxs)
           n1 = n2
           if( n1 == 0 ) n1 = 4
           call set_xsect( n1, xs_pnt, slope, intercept, has_xsect )
           if( has_xsect ) then
             nv = nv + 1
             vtx_x(nv) = xs_pnt(1)
             vtx_y(nv) = xs_pnt(2)
             found = .true.
             exit
           endif
         end do
         if( found ) then
           do n2 = abs(ndxe)-1,abs(ndxe)
             n1 = n2
             if( n1 == 0 ) n1 = 4
             call set_xsect( n1, xs_pnt, slope, intercept, has_xsect )
             if( has_xsect ) then
               nv = nv + 1
               vtx_x(nv) = xs_pnt(1)
               vtx_y(nv) = xs_pnt(2)
               exit
             endif
           end do
         endif
!---------------------------------------------------------------
!  both shd zones < 0
!---------------------------------------------------------------
       else
!---------------------------------------------------------------
!  one shd zone > 0, the other < 0
!---------------------------------------------------------------
         if( ndxs > 0 ) then
           call set_xsect( ndxs, xs_pnt, slope, intercept, has_xsect )
           if( has_xsect ) then
             nv = nv + 1
             vtx_x(nv) = xs_pnt(1)
             vtx_y(nv) = xs_pnt(2)
             do n2 = abs(ndxe)-1,abs(ndxe)
               n1 = n2
               if( n1 == 0 ) n1 = 4
               call set_xsect( n1, xs_pnt, slope, intercept, has_xsect )
               if( has_xsect ) then
                 nv = nv + 1
                 vtx_x(nv) = xs_pnt(1)
                 vtx_y(nv) = xs_pnt(2)
                 exit
               endif
             end do
           endif
         else
           found = .false.
           do n2 = abs(ndxs)-1,abs(ndxs)
             n1 = n2
             if( n1 == 0 ) n1 = 4
             call set_xsect( n1, xs_pnt, slope, intercept, has_xsect )
             if( has_xsect ) then
               nv = nv + 1
               vtx_x(nv) = xs_pnt(1)
               vtx_y(nv) = xs_pnt(2)
               found = .true.
               exit
             endif
           end do
           if( found ) then
             do n2 = abs(ndxe),abs(ndxe)+1
               n1 = n2
               if( n1 > 4 ) n1 = mod( n1,4 )
               call set_xsect( n1, xs_pnt, slope, intercept, has_xsect )
               if( has_xsect ) then
                 nv = nv + 1
                 vtx_x(nv) = xs_pnt(1)
                 vtx_y(nv) = xs_pnt(2)
                 exit
               endif
             end do
           endif
         endif
       endif
     endif
   endif vtx_in_model_cell

   end SUBROUTINE set_vertices

   END FUNCTION area_map

   SUBROUTINE set_xsect( shd_ndx, xs_pnt, slope, intercept, has_xsect )
!---------------------------------------------------------------
!  set xsect coordinates
!---------------------------------------------------------------
!---------------------------------------------------------------
!  dummy arguments
!---------------------------------------------------------------
   integer, intent(in) :: shd_ndx                ! shadow zone index( > 0 )
   real(8), intent(in) :: slope
   real(8), intent(in) :: intercept
   real, intent(inout) :: xs_pnt(2)
   logical, intent(out) :: has_xsect

!---------------------------------------------------------------
!  local variables
!---------------------------------------------------------------
   real    :: xl, xu
   real    :: yl, yu

   if( mod( shd_ndx,2 ) /= 0 ) then
     xs_pnt(:) = (/ 0., model_y(shd_ndx) /)
   else
     xs_pnt(:) = (/ model_x(shd_ndx), 0. /)
   endif

   CALL xs_coord( xs_pnt, slope, intercept )

   if( mod( shd_ndx,2 ) /= 0 ) then
     xl = model_x(1)
     xu = model_x(2)
     has_xsect = xs_pnt(1) >= xl .and.  xs_pnt(1) <= xu
   else
     yl = model_y(1)
     yu = model_y(4)
     has_xsect = xs_pnt(2) >= yl .and. xs_pnt(2) <= yu
   endif

   end SUBROUTINE set_xsect

   SUBROUTINE xs_coord( xs, slope, intercept )
!---------------------------------------------------------------
!  compute intersection coordinate
!---------------------------------------------------------------
!---------------------------------------------------------------
!  dummy arguments
!---------------------------------------------------------------
   real(8), intent(in) :: slope
   real(8), intent(in) :: intercept
   real, intent(inout) :: xs(2)

   if( xs(1) == 0. ) then
     if( slope /= 0._8 ) then
       xs(1) = real((real(xs(2),kind=8) - intercept)/slope,kind=4)
     else
       xs(1) = real(intercept,kind=4)
     endif
   else
     xs(2) = real(slope*real(xs(1),kind=8) + intercept,kind=4)
   endif

   end SUBROUTINE xs_coord

   real FUNCTION poly_area( nv, x, y )
!---------------------------------------------------------------
!  calculate the area of polynomial with nv vertices
!---------------------------------------------------------------

!---------------------------------------------------------------
!  dummy arguments
!---------------------------------------------------------------
   integer, intent(in) :: nv
   real, intent(in)    :: x(nv)
   real, intent(in)    :: y(nv)

!---------------------------------------------------------------
!  local variables
!---------------------------------------------------------------
   integer :: i, im1, ip1
   real    :: wrk(nv)

   do i = 1,nv
     ip1 = mod( i,nv ) + 1
     im1 = i - 1
     if( im1 == 0 ) im1 = nv
     wrk(i) = (x(ip1) - x(im1))*y(i)
   end do

   poly_area = -.5*sum( wrk(:) )


   END FUNCTION poly_area

   SUBROUTINE insert_mvtx( l, nv, vtx_x, vtx_y )
!---------------------------------------------------------------
!  insert model vertex in data polygon
!---------------------------------------------------------------
!---------------------------------------------------------------
!  dummy arguments
!---------------------------------------------------------------
   integer, intent(in)    :: l
   integer, intent(inout) :: nv
   real, intent(inout)    :: vtx_x(:)
   real, intent(inout)    :: vtx_y(:)

!---------------------------------------------------------------
!  local variables
!---------------------------------------------------------------
   integer :: k
   real    :: wrk_x(nv)
   real    :: wrk_y(nv)

   select case( l )
     case( 1,3 )
       do k = 1,nv
         if( model_x(l) == vtx_x(k) ) then
           exit
         endif
       end do
     case( 2,4 )
       do k = 1,nv
         if( model_y(l) == vtx_y(k) ) then
           exit
         endif
       end do
   end select

   if( k <= nv ) then
     k = k + 1
   
     if( k <= nv ) then
       wrk_x(k:nv) = vtx_x(k:nv)
       wrk_y(k:nv) = vtx_y(k:nv)
     endif
     vtx_x(k)    = model_x(l)
     vtx_y(k)    = model_y(l)
     if( k <= nv ) then
       vtx_x(k+1:nv+1) = wrk_x(k:nv)
       vtx_y(k+1:nv+1) = wrk_y(k:nv)
     endif
     nv = nv + 1
   endif

   end SUBROUTINE insert_mvtx

   logical FUNCTION pnts_match( pnt_x, pnt_y, x, y )
!---------------------------------------------------------------
!  do pnts have same coordinates?
!---------------------------------------------------------------

!---------------------------------------------------------------
!  dummy arguments
!---------------------------------------------------------------
   real, intent(in) :: pnt_x
   real, intent(in) :: pnt_y
   real, intent(in) :: x
   real, intent(in) :: y

!---------------------------------------------------------------
!  local variables
!---------------------------------------------------------------
   real, parameter :: delta = 1.e-7

   if( abs(x - pnt_x) <= delta*abs(pnt_x) .and. &
       abs(y - pnt_y) <= delta*abs(pnt_y) ) then
     pnts_match = .true.
   else
     pnts_match = .false.
   endif

   end FUNCTION pnts_match

   logical FUNCTION pnt_in_rect( pnt_x, pnt_y, x, y )
!---------------------------------------------------------------
!  determine whether point is in a rectangle
!---------------------------------------------------------------

!---------------------------------------------------------------
!  dummy arguments
!---------------------------------------------------------------
   real, intent(in) :: pnt_x
   real, intent(in) :: pnt_y
   real, intent(in) :: x(:)
   real, intent(in) :: y(:)

!---------------------------------------------------------------
!  local variables
!---------------------------------------------------------------
   integer :: l
   real    :: pnt(2)
   logical :: out_of_rect

   out_of_rect = pnt_x <= minval( x(:) ) .or. pnt_x >= maxval( x(:) ) .or. &
                 pnt_y <= minval( y(:) ) .or. pnt_y >= maxval( y(:) )

   pnt_in_rect = .not. out_of_rect

   end FUNCTION pnt_in_rect

   logical FUNCTION pnt_in_quad( pnt, x, y )
!---------------------------------------------------------------
!  determine whether input point is in quadrilateral
!---------------------------------------------------------------

!---------------------------------------------------------------
!  dummy arguments
!---------------------------------------------------------------
   real, intent(in) :: pnt(:)
   real, intent(in) :: x(:)
   real, intent(in) :: y(:)

!---------------------------------------------------------------
!  local variables
!---------------------------------------------------------------
   real :: tri_x(4)
   real :: tri_y(4)
 
   tri_x(1:3) = x(1:3)
   tri_y(1:3) = y(1:3)
   tri_x(4)   = pnt(1)
   tri_y(4)   = pnt(2)
   pnt_in_quad = pnt_in_triangle( tri_x, tri_y )
   if( .not. pnt_in_quad ) then
     tri_x(1:3) = x((/1,3,4/))
     tri_y(1:3) = y((/1,3,4/))
     pnt_in_quad = pnt_in_triangle( tri_x, tri_y )
   endif

   end FUNCTION pnt_in_quad

   logical FUNCTION pnt_in_triangle( x, y )
!---------------------------------------------------------------
!  determine whether input point is in triangle
!---------------------------------------------------------------

!---------------------------------------------------------------
!  dummy arguments
!---------------------------------------------------------------
   real, intent(in) :: x(:)
   real, intent(in) :: y(:)

!---------------------------------------------------------------
!  local variables
!---------------------------------------------------------------
   real    :: a, b, c

   a = (x(1) - x(4))*(y(2) - y(4)) - (x(2) - x(4))*(y(1) - y(4))
   b = (x(2) - x(4))*(y(3) - y(4)) - (x(3) - x(4))*(y(2) - y(4))
   c = (x(3) - x(4))*(y(1) - y(4)) - (x(1) - x(4))*(y(3) - y(4))

   pnt_in_triangle = (sign( 1.,a) == sign( 1.,b )) .and. (sign( 1.,b ) == sign( 1.,c ))

   end FUNCTION pnt_in_triangle

   SUBROUTINE lin_eqn( xs, xe, ys, ye, slope, intercept )
!---------------------------------------------------------------
!  calculate slope, intercept for linear equation between
!  (xs,ys) and (xe,ye)
!---------------------------------------------------------------
   
!---------------------------------------------------------------
!  dummy arguments
!---------------------------------------------------------------
   real, intent(in)  :: xs, xe 
   real, intent(in)  :: ys, ye 
   real(8), intent(out) :: slope
   real(8), intent(out) :: intercept

   if( xs /= xe ) then
     slope = real((ye - ys),kind=8)/real((xe - xs),kind=8)
     intercept = real(ys,kind=8) - slope*real(xs,kind=8)
     lin_parms%slope = slope
     lin_parms%intercept = intercept
     lin_parms%vert_line = .false.
   else
     slope = 0._8
     intercept = real(xs,kind=8)
     lin_parms%slope = slope
     lin_parms%intercept = intercept
     lin_parms%vert_line = .true.
   endif

   END SUBROUTINE lin_eqn

   logical function pnt_on_bndy( x, y, model_x, model_y )
!---------------------------------------------------------------------
!  is the point (x,y) on the boundary of a model cell
!---------------------------------------------------------------------

   real, intent(in) :: x
   real, intent(in) :: y
   real, intent(in) :: model_x(:)
   real, intent(in) :: model_y(:)

   pnt_on_bndy = (x == model_x(1) .or. x == model_x(2)) .or. &
                 (y == model_y(1) .or. y == model_y(3))

   end function pnt_on_bndy

   logical function pnt_is_vtx( nv, x, y, vtx_x, vtx_y )
!---------------------------------------------------------------------
!  is the point (x,y) a vertex
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  dummy arguments
!---------------------------------------------------------------------
   integer, intent(in) :: nv
   real, intent(in) :: x
   real, intent(in) :: y
   real, intent(in) :: vtx_x(:)
   real, intent(in) :: vtx_y(:)

!---------------------------------------------------------------------
!  local variables
!---------------------------------------------------------------------
   integer :: n

   pnt_is_vtx = .false.
   do n = 1,nv
     if( x == vtx_x(n) .and. y == vtx_y(n) ) then
       pnt_is_vtx = .true.
       exit
     endif
   end do

   end function pnt_is_vtx

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

   end module area_mapper
